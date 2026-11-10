#!/usr/bin/env python3
"""Capture every conda environment the project could have run in.

    python scripts/45_capture_envs.py                  # report only
    python scripts/45_capture_envs.py --into LATEST    # write into the
                                                       # newest freeze
    python scripts/45_capture_envs.py --into /path/to/98_freeze/<tag>

WHY THIS IS NOT ONE FILE

  44_freeze.py exported the active environment and called it
  ENVIRONMENT.txt. That is wrong for this project, in a way that only
  shows up when someone tries to rebuild it.

  The pipeline does not run in one environment. DIAMOND is in the BUSCO
  environment and 24_call_orfs.py goes looking for it there; the
  embedding jobs needed torch and transformers on the GPU nodes; the
  differential expression needed pydeseq2. Whichever environment happened
  to be active when the freeze ran, the export describes that one and
  silently omits the rest. A reproducibility section built on it would
  send someone to install the wrong set of packages.

  So: enumerate every environment, export each one three ways, and record
  which of them actually contains each tool the pipeline calls. That last
  table is the useful artefact — it answers "where does jackhmmer live",
  which is the question a person rebuilding this will have.

WHAT IS WRITTEN

  ENVIRONMENTS/INDEX.md              which env has which tool
  ENVIRONMENTS/<env>.export.txt      conda list --export  (name=ver=build)
  ENVIRONMENTS/<env>.explicit.txt    conda list --explicit (exact URLs)
  ENVIRONMENTS/<env>.pip.txt         pip freeze from that env's python

  The original ENVIRONMENT.txt is left exactly as it was. It is a true
  record of what the freeze captured, and overwriting it would hide that
  the first attempt was incomplete.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Tools the pipeline shells out to, and the import names it needs.
TOOLS = ["diamond", "jackhmmer", "hmmsearch", "hmmbuild", "hisat2",
         "stringtie", "samtools", "foldseek", "mmseqs", "prodigal",
         "TransDecoder.LongOrfs", "busco", "fasterq-dump", "zstd"]

# Groups of synonyms that must be present SOMEWHERE across the
# environments. Each group counts as satisfied if any member is found.
REQUIRED = [("pytorch", "torch"), ("transformers",), ("pydeseq2",),
            ("biopython",), ("numpy",), ("scipy",)]

PACKAGES = ["pytorch", "torch", "transformers", "sentencepiece", "pydeseq2",
            "biopython", "numpy", "scipy", "pandas", "matplotlib",
            "scikit-learn", "statsmodels", "pyhmmer", "python-docx",
            "openpyxl", "pyyaml", "goatools", "anndata"]


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def run(cmd: list[str], timeout: int = 300) -> tuple[int, str, str]:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout)
        return r.returncode, r.stdout, r.stderr
    except (OSError, subprocess.SubprocessError) as e:
        return 1, "", f"{e.__class__.__name__}: {e}"


# ---------------------------------------------------------------------------

def discover() -> list[tuple[str, Path]]:
    """Every conda environment, by name and prefix.

    Tried in order of reliability: the JSON API, the text listing, then
    the file conda keeps of known prefixes. The last one survives a conda
    that will not run at all, which is exactly the situation that
    produced the bad export in the first place.
    """
    envs: dict[str, Path] = {}

    if shutil.which("conda"):
        rc, out, _ = run(["conda", "env", "list", "--json"])
        if rc == 0:
            try:
                for p in json.loads(out).get("envs", []):
                    envs[Path(p).name] = Path(p)
            except ValueError:
                pass
        if not envs:
            rc, out, _ = run(["conda", "env", "list"])
            if rc == 0:
                for line in out.splitlines():
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.replace("*", " ").split()
                    if parts and Path(parts[-1]).is_dir():
                        envs[parts[0] if len(parts) > 1
                             else Path(parts[-1]).name] = Path(parts[-1])

    known = Path.home() / ".conda" / "environments.txt"
    if known.is_file():
        for line in known.read_text().splitlines():
            p = Path(line.strip())
            if p.is_dir():
                envs.setdefault(p.name, p)

    # The active environment, even if nothing above listed it. This is the
    # case that bit us: `conda activate pmbio` reported
    # EnvironmentNameNotFound while the prompt showed (pmbio), so the env
    # was live at a prefix conda would not resolve by name.
    pre = os.environ.get("CONDA_PREFIX")
    if pre and Path(pre).is_dir():
        envs.setdefault(Path(pre).name, Path(pre))

    return sorted(envs.items())


def env_python(prefix: Path) -> Path | None:
    for c in (prefix / "bin" / "python", prefix / "python.exe"):
        if c.is_file():
            return c
    return None


def capture(name: str, prefix: Path, out: Path) -> dict:
    """Export one environment three ways. Partial failure is recorded,
    never fatal — a missing pip is not a reason to lose the conda list."""
    info = {"name": name, "prefix": str(prefix), "conda": 0, "pip": 0,
            "notes": []}

    for flag, suffix in (("--export", "export"), ("--explicit", "explicit")):
        rc, txt, err = run(["conda", "list", "-p", str(prefix), flag])
        if rc == 0 and txt.strip():
            (out / f"{name}.{suffix}.txt").write_text(txt)
            if suffix == "export":
                info["conda"] = sum(
                    1 for l in txt.splitlines()
                    if l.strip() and not l.startswith("#"))
        else:
            info["notes"].append(f"conda list {flag} failed")

    py = env_python(prefix)
    if py:
        rc, txt, err = run([str(py), "-m", "pip", "freeze"])
        if rc == 0 and txt.strip():
            (out / f"{name}.pip.txt").write_text(txt)
            info["pip"] = len([l for l in txt.splitlines() if l.strip()])
        else:
            info["notes"].append("pip freeze failed")
    else:
        info["notes"].append("no python in this prefix")

    # which tools are on this env's PATH
    binp = prefix / "bin"
    info["tools"] = sorted(t for t in TOOLS if (binp / t).exists())

    # which key packages this env has, from whichever export worked
    have = set()
    for suffix in ("export", "pip"):
        f = out / f"{name}.{suffix}.txt"
        if not f.is_file():
            continue
        low = f.read_text().lower()
        for pkg in PACKAGES:
            if re.search(rf"(^|[\s=/]){re.escape(pkg.lower())}[\s=<>-]", low,
                         re.M):
                have.add(pkg)
    info["packages"] = sorted(have)
    return info


def write_index(path: Path, infos: list[dict], active: str | None):
    L = [f"# Conda environments — captured "
         f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}", "",
         f"Active at capture time: `{active or 'none'}`", "",
         "This supersedes `ENVIRONMENT.txt` in the same directory, which "
         "recorded only the environment that happened to be active. That "
         "file is left in place rather than overwritten, because it is a "
         "true record of what the freeze captured.", "",
         "## Environments", "",
         "| env | conda pkgs | pip pkgs | prefix |", "|---|---:|---:|---|"]
    for i in infos:
        L.append(f"| `{i['name']}` | {i['conda']} | {i['pip']} | "
                 f"`{i['prefix']}` |")

    L += ["", "## Where each tool lives", "",
          "The question someone rebuilding this will actually ask.", "",
          "| tool | environment(s) |", "|---|---|"]
    for t in TOOLS:
        where = [i["name"] for i in infos if t in i["tools"]]
        L.append(f"| `{t}` | {', '.join(f'`{w}`' for w in where) or '—'} |")

    L += ["", "## Where each package lives", "",
          "| package | environment(s) |", "|---|---|"]
    for p in PACKAGES:
        where = [i["name"] for i in infos if p in i["packages"]]
        L.append(f"| `{p}` | {', '.join(f'`{w}`' for w in where) or '—'} |")

    # Synonym groups: conda ships it as `pytorch`, pip as `torch`, and
    # reporting one of them missing because the other was installed is a
    # false alarm that trains people to ignore this section.
    missing = [" / ".join(grp) for grp in REQUIRED
               if not any(p in i["packages"] for p in grp for i in infos)]
    if missing:
        L += ["", "### Not found in any environment", "",
              "These were used by the pipeline but are absent from every "
              "environment captured here. Either the environment was "
              "deleted, or it lives somewhere conda does not list — a "
              "container, a module, or a prefix outside the search paths. "
              "This needs resolving before the methods section claims a "
              "reproducible environment.", ""]
        for p in missing:
            L.append(f"- `{p}`")

    notes = [(i["name"], n) for i in infos for n in i["notes"]]
    if notes:
        L += ["", "## Capture problems", ""]
        for n, msg in notes:
            L.append(f"- `{n}`: {msg}")

    L += ["", "## Rebuilding", "", "```bash",
          "conda create --name <env> --file <env>.explicit.txt",
          "# then, if that env has a pip file:",
          "pip install -r <env>.pip.txt", "```", ""]
    path.write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------------

def unlock(d: Path) -> list[tuple[Path, int]]:
    """Make a frozen directory writable, remembering the old modes."""
    saved = []
    for p in [d] + sorted(d.rglob("*")):
        try:
            m = p.stat().st_mode
            saved.append((p, m))
            p.chmod(m | stat.S_IWUSR)
        except OSError:
            pass
    return saved


def relock(saved):
    for p, m in reversed(saved):
        try:
            p.chmod(m)
        except OSError:
            pass


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--into", default=None,
                    help="freeze directory to write into, or LATEST for "
                         "the newest under the project's 98_freeze")
    args = ap.parse_args()

    envs = discover()
    if not envs:
        print("no conda environments found — is conda on PATH?",
              file=sys.stderr)
        return 1
    active = os.environ.get("CONDA_DEFAULT_ENV") or (
        Path(os.environ["CONDA_PREFIX"]).name
        if os.environ.get("CONDA_PREFIX") else None)
    log(f"{len(envs)} environment(s); active: {active or 'none'}")
    for n, p in envs:
        log(f"  {n:<22}{p}")

    target = None
    if args.into:
        if args.into.upper() == "LATEST":
            sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
            from pmdark import config as cfgmod
            fr = cfgmod.load().p("share") / "98_freeze"
            cand = sorted((d for d in fr.glob("*") if d.is_dir()),
                          key=lambda d: d.name)
            if not cand:
                print(f"no freeze directories under {fr}", file=sys.stderr)
                return 1
            target = cand[-1]
        else:
            target = Path(args.into).expanduser().resolve()
        if not target.is_dir():
            print(f"not a directory: {target}", file=sys.stderr)
            return 1

    if target is None:
        print("\nREPORT ONLY — pass --into LATEST to write the capture.")
        out = Path(os.environ.get("TMPDIR", "/tmp")) / "pmdark_envs"
    else:
        out = target / "ENVIRONMENTS"
    saved = unlock(target) if target else []
    out.mkdir(parents=True, exist_ok=True)

    infos = []
    for name, prefix in envs:
        log(f"capturing {name}")
        infos.append(capture(name, prefix, out))

    write_index(out / "INDEX.md", infos, active)
    log(f"wrote {out}/INDEX.md")

    if target:
        # the new subdirectory is frozen too, then the parent is relocked
        for p in sorted(out.rglob("*")) + [out]:
            try:
                p.chmod(p.stat().st_mode & ~0o222)
            except OSError:
                pass
        relock(saved)
        log("freeze directory re-locked")

    print("\n" + "=" * 78)
    print(f"{'env':<22}{'conda':>7}{'pip':>6}   tools")
    print("=" * 78)
    for i in infos:
        print(f"{i['name']:<22}{i['conda']:>7}{i['pip']:>6}   "
              f"{', '.join(i['tools']) or '—'}")
    print()
    for grp in REQUIRED:
        where = sorted({i["name"] for p in grp for i in infos
                        if p in i["packages"]})
        flag = "" if where else "   <-- NOT FOUND IN ANY ENVIRONMENT"
        print(f"  {' / '.join(grp):<22}{', '.join(where) or '—'}{flag}")
    print(f"\n  written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
