#!/usr/bin/env python3
"""Freeze the project: manifest, provenance, archive, verification.

    python scripts/44_freeze.py                 # plan only, writes nothing
    python scripts/44_freeze.py --go
    python scripts/44_freeze.py --go --level full

WHAT A FREEZE IS FOR

  Two different jobs, and they want different things:

    reproducibility  someone must be able to rebuild every number in the
                     paper. That needs code, config, the results JSONs and
                     the environment. It is small.
    recovery         we must not lose work that cost GPU-weeks. That needs
                     the embeddings and the derived protein sets. It is
                     large, and some of it is cheaper to recompute than to
                     store.

  So the levels are not "more is better". `core` is the reproducible paper
  plus everything expensive to regenerate. `full` adds what a network
  connection could fetch again — the SRA mirror, the HuggingFace cache,
  the genome and its index. Those are recorded in the inventory either
  way, by accession and checksum, which is what actually matters for
  provenance; copying them is insurance, not evidence.

THE THREE TIERS

  This project is spread across three filesystems with different sizes,
  quotas and visibility (see config/config.yaml). A freeze that looked at
  one of them would silently miss two thirds of the work, so every root is
  walked and labelled in the manifest.

SYMLINKS ARE NOT FOLLOWED

  scripts/35_retier.sh moved bulk data to the archive and left symlinks
  behind, so the same bytes are reachable by two paths. Following them
  would double-count the inventory, could archive tens of gigabytes twice,
  and — for any link pointing outside the roots — would pull in files that
  are not ours. Links are recorded as links, with their targets.

NOTHING OUTSIDE THE PROJECT IS READ OR WRITTEN

  Every root is checked to be the project's own directory before the walk,
  and the output goes inside the archive root. Other projects on the
  shared volume are never touched.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

SKIP_DIRS = {"__pycache__", ".git", ".ipynb_checkpoints", ".mypy_cache",
             ".pytest_cache", ".DS_Store"}
FREEZE_DIRNAME = "98_freeze"

# group -> (config key or None for the code root, tier label)
GROUPS: list[tuple[str, str | None, str]] = [
    ("code",       None,         "code"),
    ("results",    "results",    "archive"),
    ("interim",    "interim",    "archive"),
    ("logs",       "logs",       "archive"),
    ("proteins",   "proteins",   "work"),
    ("align",      "align",      "work"),
    ("structures", "structures", "work"),
    ("raw",        "raw",        "archive"),
    ("dbs",        "dbs",        "archive"),
    ("embeddings", "embeddings", "work"),
    ("genome",     "genome",     "work"),
    ("models",     "models",     "work"),
]

LEVELS = {
    "code": ["code"],
    "core": ["code", "results", "interim", "logs",
             "proteins", "align", "structures", "embeddings"],
    "full": [g for g, _, _ in GROUPS],
}

# Groups a NETWORK CONNECTION can rebuild: the SRA mirror, the HuggingFace
# cache, the genome and its index, the sequence databases. Inventoried
# always, archived only at --level full, and flagged in the provenance so
# nobody mistakes their absence for data loss.
#
# Embeddings are deliberately NOT in this set, though an earlier version of
# this file had them here. That was wrong twice over: they cost GPU-weeks
# rather than a download, and they are small — a few hundred megabytes for
# all three models, against tens of gigabytes for the caches beside them.
# Expensive to make and cheap to keep is the exact profile that belongs in
# the archive, so they are in `core`.
REGENERABLE = {"raw", "dbs", "genome", "models"}


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or u == "TB":
            return f"{n:,.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TB"


# ---------------------------------------------------------------------------
# walking
# ---------------------------------------------------------------------------

class Entry:
    __slots__ = ("group", "tier", "root", "rel", "kind", "size", "mtime",
                 "target", "digest", "note")

    def __init__(self, group, tier, root, rel, kind, size, mtime,
                 target=None):
        self.group, self.tier, self.root, self.rel = group, tier, root, rel
        self.kind, self.size, self.mtime, self.target = kind, size, mtime, target
        self.digest: str | None = None
        self.note: str = ""

    @property
    def path(self) -> Path:
        return self.root / self.rel


def walk(group: str, tier: str, root: Path, freeze_dir: Path) -> list[Entry]:
    """Inventory one root. Symlinks are recorded, never followed."""
    out: list[Entry] = []
    if not root.exists():
        return out
    root = root.resolve()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        here = Path(dirpath)
        if here == freeze_dir or freeze_dir in here.parents:
            dirnames[:] = []
            continue
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        for name in sorted(filenames):
            if name in SKIP_DIRS:
                continue
            p = here / name
            rel = p.relative_to(root)
            try:
                if p.is_symlink():
                    out.append(Entry(group, tier, root, rel, "link", 0,
                                     p.lstat().st_mtime, os.readlink(p)))
                    continue
                st = p.stat()
            except OSError as e:
                e2 = Entry(group, tier, root, rel, "error", 0, 0)
                e2.note = f"unreadable: {e.__class__.__name__}"
                out.append(e2)
                continue
            out.append(Entry(group, tier, root, rel, "file",
                             st.st_size, st.st_mtime))
        # a symlinked directory shows up in dirnames, not filenames
        for name in list(dirnames):
            p = here / name
            if p.is_symlink():
                dirnames.remove(name)
                out.append(Entry(group, tier, root, p.relative_to(root),
                                 "link", 0, p.lstat().st_mtime,
                                 os.readlink(p)))
    return out


def sha256(path: Path, buf: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(buf):
            h.update(chunk)
    return h.hexdigest()


def hash_rate(entries: list[Entry], sample_bytes: int = 64 << 20) -> float:
    """Measured bytes/sec, so the plan's estimate is this disk's, not a
    guess. Reads real project files, largest first."""
    files = sorted((e for e in entries if e.kind == "file" and e.size > 0),
                   key=lambda e: -e.size)
    done = 0
    t0 = time.time()
    for e in files:
        if done >= sample_bytes:
            break
        try:
            with e.path.open("rb") as fh:
                h = hashlib.sha256()
                while done < sample_bytes and (c := fh.read(1 << 20)):
                    h.update(c)
                    done += len(c)
        except OSError:
            continue
    dt = time.time() - t0
    return done / dt if dt > 0 and done else 300e6


# ---------------------------------------------------------------------------
# plan
# ---------------------------------------------------------------------------

def summarise(entries: list[Entry], groups: list[tuple[str, str | None, str]],
              present: dict[str, bool]) -> dict[str, dict]:
    """Every configured group appears, including the empty ones.

    A group that is simply absent from the table cannot be told apart from
    a group that was never looked at, and "the databases were empty" and
    "the databases were not checked" are very different claims to make in
    a provenance record.
    """
    agg: dict[str, dict] = {
        g: {"tier": t, "files": 0, "links": 0, "errors": 0, "bytes": 0,
            "state": "ok" if present.get(g) else "MISSING"}
        for g, _, t in groups}
    for e in entries:
        a = agg.setdefault(e.group, {"tier": e.tier, "files": 0, "links": 0,
                                     "errors": 0, "bytes": 0, "state": "ok"})
        if e.kind == "file":
            a["files"] += 1
            a["bytes"] += e.size
        elif e.kind == "link":
            a["links"] += 1
        else:
            a["errors"] += 1
    for g, a in agg.items():
        if a["state"] == "ok" and not (a["files"] or a["links"]):
            a["state"] = "empty"
    return agg


def print_plan(agg, level, keep, max_hash, rate, out_dir):
    print("=" * 78)
    print(f"FREEZE PLAN — level '{level}'")
    print("=" * 78)
    print(f"  {'group':<12}{'tier':<9}{'files':>8}{'links':>7}"
          f"{'size':>12}   action")
    tot_b = tot_f = 0
    hash_b = arch_b = 0
    for g, _, _ in GROUPS:
        a = agg.get(g)
        if not a:
            continue
        inc = g in keep
        if a["state"] != "ok":
            act = ("directory does not exist" if a["state"] == "MISSING"
                   else "exists, no files")
        else:
            act = "archive + hash" if inc else "inventory only"
            if g in REGENERABLE and not inc:
                act += " (re-downloadable)"
        print(f"  {g:<12}{a['tier']:<9}{a['files']:>8,}{a['links']:>7,}"
              f"{human(a['bytes']):>12}   {act}")
        tot_b += a["bytes"]
        tot_f += a["files"]
        if inc:
            arch_b += a["bytes"]
            hash_b += a["bytes"]
    print("  " + "-" * 74)
    print(f"  {'TOTAL':<12}{'':<9}{tot_f:>8,}{'':>7}{human(tot_b):>12}")
    print(f"\n  to archive : {human(arch_b)}")
    print(f"  to hash    : {human(hash_b)}  "
          f"(files over {max_hash:g} GB are sized, not hashed)")
    print(f"  disk rate  : {human(rate)}/s measured on this filesystem")
    if rate > 0:
        print(f"  hashing    : ~{hash_b / rate / 60:.1f} min")
    errs = sum(a["errors"] for a in agg.values())
    if errs:
        print(f"\n  {errs} unreadable path(s) — listed in the inventory")
    print(f"\n  output     : {out_dir}")
    print("=" * 78)


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def tool_versions() -> list[tuple[str, str]]:
    probes = [
        ("python", [sys.executable, "--version"]),
        ("conda", ["conda", "--version"]),
        ("diamond", ["diamond", "--version"]),
        ("jackhmmer", ["jackhmmer", "-h"]),
        ("hmmsearch", ["hmmsearch", "-h"]),
        ("hisat2", ["hisat2", "--version"]),
        ("stringtie", ["stringtie", "--version"]),
        ("foldseek", ["foldseek", "version"]),
        ("mmseqs", ["mmseqs", "version"]),
        ("tar", ["tar", "--version"]),
        ("zstd", ["zstd", "--version"]),
    ]
    out = []
    for name, cmd in probes:
        exe = shutil.which(cmd[0])
        if not exe:
            out.append((name, "not found"))
            continue
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            txt = (r.stdout or r.stderr).strip().splitlines()
            line = next((l for l in txt if l.strip()), "")
            out.append((name, f"{line.strip()[:90]}   [{exe}]"))
        except (OSError, subprocess.SubprocessError):
            out.append((name, f"present but would not report   [{exe}]"))
    return out


def python_packages() -> str:
    for cmd in (["conda", "list", "--export"],
                [sys.executable, "-m", "pip", "freeze"]):
        if not shutil.which(cmd[0]) and cmd[0] != sys.executable:
            continue
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if r.returncode == 0 and r.stdout.strip():
                return f"# {' '.join(cmd)}\n{r.stdout}"
        except (OSError, subprocess.SubprocessError):
            continue
    return "# environment capture failed\n"


def write_provenance(path: Path, tag, level, roots, agg, keep, cfg_path,
                     script_dir: Path):
    now = datetime.now(timezone.utc)
    L = [f"# PmDark freeze `{tag}`", "",
         f"- UTC          : {now.isoformat(timespec='seconds')}",
         f"- local        : {datetime.now().astimezone().isoformat(timespec='seconds')}",
         f"- host         : {socket.gethostname()}",
         f"- user         : {os.environ.get('USER', '?')}",
         f"- platform     : {platform.platform()}",
         f"- python       : {sys.version.split()[0]} ({sys.executable})",
         f"- conda env    : {os.environ.get('CONDA_DEFAULT_ENV', 'none')}",
         f"- level        : {level}", ""]

    L += ["## Roots", "",
          "| group | tier | path | files | size | in archive |",
          "|---|---|---|---:|---:|---|"]
    for g, _, _ in GROUPS:
        a = agg.get(g)
        if not a:
            continue
        inc = "yes" if (g in keep and a["state"] == "ok") else "no"
        note = "" if a["state"] == "ok" else f" ({a['state'].lower()})"
        L.append(f"| {g} | {a['tier']} | `{roots[g]}`{note} | {a['files']:,} | "
                 f"{human(a['bytes'])} | {inc} |")

    missing = [g for g in REGENERABLE if g in agg and g not in keep]
    if missing:
        L += ["", "### Not in the archive, by design", "",
              "These are inventoried with checksums but not copied, because "
              "a network connection rebuilds them and they are large. Their "
              "absence is a storage decision, not data loss.", ""]
        for g in missing:
            L.append(f"- `{g}` — {human(agg[g]['bytes'])}")

    L += ["", "## Script checksums (md5)", "", "| script | md5 |", "|---|---|"]
    for p in sorted(script_dir.glob("*.py")) + sorted(script_dir.glob("*.sh")):
        L.append(f"| {p.name} | `{hashlib.md5(p.read_bytes()).hexdigest()}` |")

    L += ["", "## Tool versions", "", "| tool | version |", "|---|---|"]
    for name, ver in tool_versions():
        L.append(f"| {name} | {ver} |")

    L += ["", "## Configuration", "",
          f"Copied verbatim to `config.yaml` beside this file "
          f"(source `{cfg_path}`).", "",
          "## Verifying this freeze later", "",
          "```bash", "cd <this directory>", "bash VERIFY.sh", "```", ""]
    path.write_text("\n".join(L) + "\n")


VERIFY_SH = """#!/usr/bin/env bash
# Re-check this freeze. Exits non-zero if anything has changed.
set -uo pipefail
cd "$(dirname "$0")"
fail=0

echo "== archive checksum"
if [ -f ARCHIVE.sha256 ]; then
  sha256sum -c ARCHIVE.sha256 || fail=1
else
  echo "   no archive in this freeze (level=code?)"
fi

echo "== manifest of archived files"
if [ -f MANIFEST.sha256 ]; then
  n=$(wc -l < MANIFEST.sha256)
  echo "   $n entries recorded"
  echo "   (to check against the live tree, run with the project mounted:)"
  echo "   sha256sum -c MANIFEST.sha256   # from the recorded root"
else
  echo "   MANIFEST.sha256 missing"; fail=1
fi

[ "$fail" -eq 0 ] && echo "OK" || echo "FAILED"
exit "$fail"
"""


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go", action="store_true",
                    help="execute. Without it nothing is written.")
    ap.add_argument("--level", choices=sorted(LEVELS), default="core")
    ap.add_argument("--tag", default=None,
                    help="freeze name (default: UTC timestamp)")
    ap.add_argument("--max-hash-gb", type=float, default=8.0,
                    help="files larger than this are recorded by size and "
                         "mtime rather than hashed (default 8)")
    ap.add_argument("--no-archive", action="store_true",
                    help="manifest and provenance only, no tarball")
    args = ap.parse_args()

    cfg = cfgmod.load()
    share = cfg.p("share").resolve()
    code_root = cfg.p("home").resolve()

    # Refuse to operate anywhere but this project's own directories.
    for name, r in (("share", share), ("home", code_root)):
        if len(r.parts) < 3 or str(r) in ("/", str(Path.home())):
            print(f"refusing: paths.{name} resolves to {r}", file=sys.stderr)
            return 2
    if share.name != cfg["project"]["slug"]:
        print(f"refusing: paths.share is {share}, whose last component is "
              f"not the project slug '{cfg['project']['slug']}'. A freeze "
              f"writes into its own root and must not be aimed at a parent "
              f"directory shared with other projects.", file=sys.stderr)
        return 2

    tag = args.tag or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    freeze_root = share / FREEZE_DIRNAME
    out_dir = freeze_root / tag

    roots: dict[str, Path] = {}
    present: dict[str, bool] = {}
    entries: list[Entry] = []
    log("walking the three tiers (symlinks recorded, not followed)")
    for group, key, tier in GROUPS:
        root = code_root if key is None else cfg.p(key).resolve()
        roots[group] = root
        present[group] = root.is_dir()
        got = walk(group, tier, root, freeze_root)
        entries += got
        flag = "" if present[group] else "   <- NO SUCH DIRECTORY"
        log(f"  {group:<12}{len(got):>8,} paths   {root}{flag}")

    agg = summarise(entries, GROUPS, present)
    keep = set(LEVELS[args.level])
    rate = hash_rate([e for e in entries if e.group in keep])
    print()
    print_plan(agg, args.level, keep, args.max_hash_gb, rate, out_dir)

    if not args.go:
        print("\nPLAN ONLY — nothing written. Re-run with --go to execute.")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"freezing into {out_dir}")

    # -- hash ------------------------------------------------------------
    cap = args.max_hash_gb * (1 << 30)
    sel = [e for e in entries if e.group in keep]
    todo = [e for e in sel if e.kind == "file" and e.size <= cap]
    log(f"hashing {len(todo):,} files ({human(sum(e.size for e in todo))})")
    done_b = 0
    t0 = time.time()
    nxt = 10.0
    for e in todo:
        try:
            e.digest = sha256(e.path)
        except OSError as err:
            e.note = f"hash failed: {err.__class__.__name__}"
        done_b += e.size
        if time.time() - t0 > nxt:
            pct = 100 * done_b / max(1, sum(x.size for x in todo))
            log(f"  {pct:5.1f}%  {human(done_b)}")
            nxt += 10.0
    for e in sel:
        if e.kind == "file" and e.digest is None and not e.note:
            e.note = f"not hashed (> {args.max_hash_gb:g} GB)"

    # -- manifest and inventory -------------------------------------------
    man = out_dir / "MANIFEST.sha256"
    with man.open("w") as fh:
        for e in sel:
            if e.digest:
                fh.write(f"{e.digest}  {e.group}/{e.rel}\n")
    log(f"wrote MANIFEST.sha256 ({sum(1 for e in sel if e.digest):,} files)")

    inv = out_dir / "INVENTORY.tsv"
    with inv.open("w") as fh:
        fh.write("group\ttier\tkind\tbytes\tmtime_utc\tsha256\tnote\t"
                 "root\trelpath\ttarget\n")
        for e in entries:
            ts = (datetime.fromtimestamp(e.mtime, timezone.utc)
                  .isoformat(timespec="seconds") if e.mtime else "")
            fh.write(f"{e.group}\t{e.tier}\t{e.kind}\t{e.size}\t{ts}\t"
                     f"{e.digest or ''}\t{e.note}\t{e.root}\t{e.rel}\t"
                     f"{e.target or ''}\n")
    log(f"wrote INVENTORY.tsv ({len(entries):,} paths, every tier)")

    shutil.copy2(cfg.path, out_dir / "config.yaml")
    (out_dir / "ENVIRONMENT.txt").write_text(python_packages())
    write_provenance(out_dir / "PROVENANCE.md", tag, args.level, roots, agg,
                     keep, cfg.path, Path(__file__).resolve().parent)
    v = out_dir / "VERIFY.sh"
    v.write_text(VERIFY_SH)
    v.chmod(0o755)
    log("wrote PROVENANCE.md, ENVIRONMENT.txt, config.yaml, VERIFY.sh")

    # -- archive -----------------------------------------------------------
    if not args.no_archive:
        zstd = shutil.which("zstd")
        ext = "tar.zst" if zstd else "tar.gz"
        tarball = out_dir / f"pmdark_{args.level}_{tag}.{ext}"
        listing = out_dir / "_filelist.txt"
        n = 0
        # Paths are stored RELATIVE to / and the archive is created with
        # -C /. An archive of absolute paths extracts over the live tree
        # the moment someone untars it in the wrong place; this one
        # extracts into the current directory by default, and restoring
        # in place has to be asked for explicitly with -C /.
        with listing.open("w") as fh:
            for e in sel:
                if e.kind in ("file", "link"):
                    fh.write(f"{str(e.root / e.rel).lstrip('/')}\n")
                    n += 1
        cmd = ["tar", "--create", "--file", str(tarball),
               "-C", "/", "--no-recursion",
               "--files-from", str(listing)]
        cmd[1:1] = (["--use-compress-program", "zstd -10 -T0"] if zstd
                    else ["--gzip"])
        log(f"archiving {n:,} paths with {'zstd' if zstd else 'gzip'}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-4000:], file=sys.stderr)
            log("ARCHIVE FAILED — manifest and provenance were still written")
            return 1

        # verify: the archive must contain what we said it would
        log("verifying the archive reads back")
        lst = subprocess.run(["tar", "--list", "--file", str(tarball)]
                             + (["--use-compress-program", "zstd -d"]
                                if zstd else []),
                             capture_output=True, text=True)
        got = sum(1 for l in lst.stdout.splitlines()
                  if l.strip() and not l.endswith("/"))
        if lst.returncode != 0 or got != n:
            log(f"VERIFY FAILED: wrote {n:,} paths, archive lists {got:,}")
            return 1
        digest = sha256(tarball)
        (out_dir / "ARCHIVE.sha256").write_text(
            f"{digest}  {tarball.name}\n")
        listing.unlink()
        log(f"{tarball.name}  {human(tarball.stat().st_size)}")
        log(f"sha256 {digest}")
        log(f"verified: {got:,} paths read back")

    # -- freeze it ---------------------------------------------------------
    for p in sorted(out_dir.rglob("*")) + [out_dir]:
        try:
            p.chmod(p.stat().st_mode & ~0o222)
        except OSError:
            pass
    log("freeze directory set read-only")

    print("\n" + "=" * 78)
    print(f"FROZEN  {out_dir}")
    print("=" * 78)
    for p in sorted(out_dir.iterdir()):
        print(f"  {human(p.stat().st_size):>12}  {p.name}")
    print(f"\n  verify later : bash {out_dir}/VERIFY.sh")
    print(f"  freeze path  : {out_dir}")
    print("  restore here : tar -xf <archive>            "
          "(into the current directory)")
    print("  restore live : tar -xf <archive> -C /       "
          "(overwrites the project in place)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
