#!/usr/bin/env python3
"""Build the GitHub repository and the Zenodo deposit.

    python scripts/46_release.py --init          # write release.yaml stub
    python scripts/46_release.py                 # plan (writes nothing)
    python scripts/46_release.py --go

WHY THE METADATA IS A FILE AND NOT ARGUMENTS

  A Zenodo DOI is permanent. A record published with a misspelt author or
  a missing ORCID cannot be quietly corrected: the DOI is minted, indexed
  and cited. So the author list, the affiliations and the repository name
  come from `config/release.yaml`, which a person fills in and reads back,
  rather than from flags typed once at a prompt.

  Every field starts as REPLACE_ME and this script REFUSES to run while
  any remain, listing each one. That is deliberate. Placeholder text that
  reaches a published artefact is the single most repeated failure in
  this project's history; a stub that silently became "Author, A." in a
  citation record would be the worst instance of it.

WHAT IS BUILT

  97_release/<tag>/repo/      exactly what goes in the GitHub repository
  97_release/<tag>/zenodo/    the deposit archive and its metadata JSON
  97_release/<tag>/NEXT_STEPS.md

  The repository gets code, documentation, the conda environment exports,
  the figures and the small result files. Any single file above
  --max-repo-mb is left out and recorded instead in the data-availability
  section as living in the Zenodo record, because GitHub rejects files
  over 100 MB and a repository nobody can clone is not a release.

NO NUMBER IS WRITTEN INTO THE README

  The same rule as 42_consolidate.py and 43_figures.py. Every figure in
  the README is read from a results JSON at build time, so the repository
  front page cannot drift from the paper. A number whose source file is
  absent is omitted rather than guessed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import date, datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

PLACEHOLDER = "REPLACE_ME"

STUB = """# Release metadata for GitHub and Zenodo.
#
# Every REPLACE_ME must be filled in. scripts/46_release.py refuses to run
# while any remain — a Zenodo DOI is permanent and placeholder text in a
# citation record cannot be taken back.

github:
  owner: REPLACE_ME            # your GitHub username or organisation
  repo: pmdark                 # repository name
  private: true                # flip to false when the paper is accepted
  default_branch: main

zenodo:
  # "reserve" mints a DOI you can cite in the submitted manuscript while
  # the record stays unpublished. "publish" makes it live immediately.
  mode: reserve
  sandbox_first: true          # rehearse on sandbox.zenodo.org first
  upload_type: software
  access_right: restricted     # becomes "open" when you publish
  keywords:
    - protein language model
    - function annotation
    - dark proteome
    - Penaeus monodon
    - whitening
    - Gene Ontology

licence:
  code: MIT
  data: CC-BY-4.0

authors:
  # One block per author, in the order they appear on the paper.
  # orcid may be left as an empty string if an author has none; it must
  # not be left as REPLACE_ME.
  - name: REPLACE_ME           # "Family, Given" — Zenodo's expected form
    affiliation: REPLACE_ME
    orcid: REPLACE_ME          # "0000-0000-0000-0000", or "" if none

corresponding: REPLACE_ME      # email for the repository README
"""


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:,.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} GB"


import re                                                      # noqa: E402

ORCID_RE = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$")

# Text that looks like it came from an example rather than from a person.
# REPLACE_ME alone is not enough: the usual way placeholder text reaches a
# permanent record is that someone pastes the documented example command,
# which is a perfectly sensible thing to do. The example must therefore be
# rejected as firmly as the stub.
LOOKS_LIKE_EXAMPLE = re.compile(
    r"^\s*<.*>\s*$"                      # <anything in angle brackets>
    r"|^(yourname|your[-_ ]?name|username|your[-_ ]?github[-_ ]?user"
    r"|youruser|myname|name|surname|author|firstname|lastname|test|foo)\s*$"
    r"|surname\s*,\s*(given|a\.?\s*b\.?|initials?)"
    r"|^author\s*,\s*a\.?"
    r"|\byour\s+(university|institution|affiliation|organisation"
    r"|organization|lab)\b"
    r"|^you@|@example\.(com|org|net)|@institution\b"
    r"|replace[-_ ]?me",
    re.I)


def orcid_valid(s: str) -> bool:
    """ORCID checksum, ISO 7064 MOD 11-2.

    This is what catches 0000-0000-0000-0000 — the number used in every
    piece of ORCID documentation, and therefore the one most likely to be
    pasted by accident. Its check digit should be 1, so the all-zero form
    is not merely unassigned, it is arithmetically impossible.
    """
    if not ORCID_RE.match(s):
        return False
    digits = s.replace("-", "")
    total = 0
    for ch in digits[:15]:
        total = (total + int(ch)) * 2
    expected = (12 - total % 11) % 11
    return ("X" if expected == 10 else str(expected)) == digits[15].upper()

# Which flag sets which field, so the error tells you what to type rather
# than which line of a file to go and edit.
FLAG_FOR = {
    "github.owner": "--github-user YOURNAME",
    "github.repo": "--repo NAME",
    "corresponding": "--email you@example.org",
    "authors[].name": '--author "Surname, Given"',
    "authors[].affiliation": '--affiliation "Your University"',
    "authors[].orcid": "--orcid 0000-0000-0000-0000   (or --orcid none)",
}


def apply_flags(meta: dict, args) -> list[str]:
    """Fold command-line flags into the metadata. Returns what changed."""
    changed = []
    if args.github_user:
        meta.setdefault("github", {})["owner"] = args.github_user
        changed.append("github.owner")
    if args.repo:
        meta.setdefault("github", {})["repo"] = args.repo
        changed.append("github.repo")
    if args.email:
        meta["corresponding"] = args.email
        changed.append("corresponding")
    if args.public is not None:
        meta.setdefault("github", {})["private"] = not args.public
        meta.setdefault("zenodo", {})["access_right"] = (
            "open" if args.public else "restricted")
        changed.append("visibility")
    if getattr(args, "sandbox", None) is not None:
        meta.setdefault("zenodo", {})["sandbox_first"] = args.sandbox
        changed.append("zenodo.sandbox_first")

    if args.author:
        affs, orcs = args.affiliation, args.orcid
        # One affiliation for several authors is the common case and
        # saves retyping the same institution.
        if len(affs) == 1 and len(args.author) > 1:
            affs = affs * len(args.author)
        if len(orcs) == 1 and len(args.author) > 1 and orcs[0].lower() in (
                "none", ""):
            orcs = orcs * len(args.author)
        if affs and len(affs) != len(args.author):
            raise SystemExit(
                f"{len(args.author)} --author but {len(affs)} "
                f"--affiliation. Give one each, in the same order, or a "
                f"single --affiliation for all of them.")
        if orcs and len(orcs) != len(args.author):
            raise SystemExit(
                f"{len(args.author)} --author but {len(orcs)} --orcid. "
                f"Give one each, in the same order, using --orcid none "
                f"for an author without one.")
        people = []
        for i, name in enumerate(args.author):
            oc = orcs[i].strip() if i < len(orcs) else PLACEHOLDER
            # Validate before anything is written. Saving first and
            # rejecting after leaves a bad ORCID sitting in the file for
            # the next run to trip over.
            if oc and oc.lower() not in ("none", "-") \
                    and oc != PLACEHOLDER and not orcid_valid(oc):
                why = ("is not in the form 0000-0000-0000-0000"
                       if not ORCID_RE.match(oc) else
                       "fails the ORCID checksum, so it is not a real "
                       "ORCID (0000-0000-0000-0000 is the documentation "
                       "example and fails this check by design)")
                raise SystemExit(
                    f"--orcid '{oc}' for {name} {why}.\n"
                    f"A wrong ORCID in a minted DOI record points at "
                    f"someone else, or nobody. Nothing was saved.")
            people.append({
                "name": name,
                "affiliation": affs[i] if i < len(affs) else PLACEHOLDER,
                "orcid": "" if oc.lower() in ("none", "-") else oc,
            })
        meta["authors"] = people
        changed.append(f"{len(people)} author(s)")
    return changed


def save_meta(path: Path, meta: dict) -> None:
    import yaml
    body = yaml.safe_dump(meta, sort_keys=False, allow_unicode=True,
                          default_flow_style=False)
    path.write_text(
        "# Release metadata for GitHub and Zenodo.\n"
        "# Written by scripts/46_release.py from its command-line flags —\n"
        "# pass them again to change anything. No need to edit this by hand.\n"
        "\n" + body)


def find_placeholders(node, trail=""):
    """Every REPLACE_ME still in the metadata, by dotted path."""
    out = []
    if isinstance(node, dict):
        for k, v in node.items():
            out += find_placeholders(v, f"{trail}.{k}" if trail else str(k))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            out += find_placeholders(v, f"{trail}[{i}]")
    elif isinstance(node, str) and PLACEHOLDER in node:
        out.append(trail)
    return out


# ---------------------------------------------------------------------------
# numbers for the README, read from the results at build time
# ---------------------------------------------------------------------------

def readme_numbers(res: Path) -> list[tuple[str, str]]:
    """(label, value) pairs. A row whose source is missing is dropped,
    never estimated."""
    def J(name):
        p = res / name
        try:
            return json.loads(p.read_text()) if p.exists() else None
        except ValueError:
            return None

    rows: list[tuple[str, str]] = []
    boot = J("bootstrap_ci.json") or {}
    cells = boot.get("cells", {})

    def cell(key):
        c = cells.get(key)
        if not c or c.get("fmax") is None:
            return None
        ci = c.get("fmax_ci") or [None, None]
        s = f"{c['fmax']:.3f}"
        if ci[0] is not None:
            s += f" [{ci[0]:.3f}, {ci[1]:.3f}]"
        return s, c.get("n")

    t0 = cell("esm2/whiten/T0_full")
    if t0:
        rows.append(("Fmax, all benchmark queries (ESM-2, whitened)", t0[0]))
        if t0[1]:
            rows.append(("Benchmark queries", f"{t0[1]:,}"))
    orp = cell("esm2/whiten/ORPHAN_nohit")
    if orp:
        rows.append(("Fmax on proteins with no BLAST-class hit", orp[0]))
        if orp[1]:
            rows.append(("Proteins with no BLAST-class hit", f"{orp[1]:,}"))

    h = J("orphan_hardening.json") or {}
    hard = (h.get("strata") or {}).get("HARDENED") or {}
    if hard.get("fmax") is not None:
        ci = hard.get("ci") or [None, None]
        s = f"{hard['fmax']:.3f}"
        if ci[0] is not None:
            s += f" [{ci[0]:.3f}, {ci[1]:.3f}]"
        rows.append(("Fmax after iterative-profile hardening", s))
        if hard.get("n"):
            rows.append(("Proteins with no hit from either search",
                         f"{hard['n']:,}"))

    an = J("annotation_summary.json") or {}
    if an.get("n_targets"):
        rows.append(("STAT-dependent uncharacterised targets annotated",
                     f"{an['n_targets']:,}"))
    if an.get("n_orphan") is not None:
        rows.append(("…of which have no reachable homologue",
                     f"{an['n_orphan']:,}"))

    mt = J("matched_enrichment.json") or {}
    tr = mt.get("trend") or {}
    if tr.get("spearman_rho") is not None:
        pv = tr.get("spearman_p", tr.get("p"))
        s = f"{tr['spearman_rho']:+.2f}"
        if pv is not None:
            s += f" (p = {pv:.3f})" if pv >= 1e-3 else " (p < 0.001)"
        rows.append(("Spearman |δ| vs expression decile", s))
    return rows


# ---------------------------------------------------------------------------
# generated files
# ---------------------------------------------------------------------------

MIT = """MIT License

Copyright (c) {year} {holders}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

GITIGNORE = """__pycache__/
*.py[cod]
.ipynb_checkpoints/
.DS_Store
.env
.venv/
*.egg-info/

# Large or machine-specific artefacts — these live in the Zenodo record
*.npy
*.npz
*.bam
*.sam
*.fastq
*.fastq.gz
*.tar.zst
*.tar.gz
"""


def cff(meta, title, version, doi) -> str:
    L = ["cff-version: 1.2.0",
         "message: If you use this software, please cite it as below.",
         f'title: "{title}"',
         "type: software",
         f"version: {version}",
         f"date-released: {date.today().isoformat()}",
         f"license: {meta['licence']['code']}",
         "authors:"]
    for a in meta["authors"]:
        fam, _, giv = a["name"].partition(",")
        L.append(f'  - family-names: "{fam.strip()}"')
        if giv.strip():
            L.append(f'    given-names: "{giv.strip()}"')
        L.append(f'    affiliation: "{a["affiliation"]}"')
        if a.get("orcid"):
            L.append(f'    orcid: "https://orcid.org/{a["orcid"]}"')
    if doi:
        L += ["identifiers:", "  - type: doi", f"    value: {doi}",
              "    description: Archived release"]
    return "\n".join(L) + "\n"


def zenodo_json(meta, title, description) -> str:
    creators = []
    for a in meta["authors"]:
        c = {"name": a["name"], "affiliation": a["affiliation"]}
        if a.get("orcid"):
            c["orcid"] = a["orcid"]
        creators.append(c)
    z = {
        "title": title,
        "upload_type": meta["zenodo"]["upload_type"],
        "description": description,
        "creators": creators,
        "keywords": meta["zenodo"]["keywords"],
        "access_right": meta["zenodo"]["access_right"],
        "license": meta["licence"]["code"].lower().replace(" ", "-"),
        "language": "eng",
    }
    return json.dumps(z, indent=2, ensure_ascii=False) + "\n"


def readme(meta, title, rows, figs, excluded, envs, tag) -> str:
    gh = meta["github"]
    L = [f"# {title}", "",
         "Calibrated protein-language-model function transfer, evaluated "
         "as homologous evidence is withdrawn, and applied to the "
         "STAT-dependent uncharacterised fraction of the *Penaeus "
         "monodon* hemocyte proteome.", "",
         "> **Status.** Manuscript in preparation. This repository is the "
         "code and the derived results; the archived release carries the "
         "DOI cited in the paper.", ""]

    if rows:
        L += ["## Headline results", "",
              "Every figure below is read from the result files in "
              "`results/` when this README is generated, so it cannot "
              "drift from the paper.", "",
              "| quantity | value |", "|---|---|"]
        # A pipe inside a cell ends the cell. "Spearman |δ|" would split
        # the row into four columns and render as garbage.
        esc = lambda s: str(s).replace("|", "\\|")
        L += [f"| {esc(k)} | {esc(v)} |" for k, v in rows]
        L.append("")

    L += ["## Layout", "", "```",
          "scripts/     numbered pipeline, run in order",
          "src/pmdark/  importable package: config, embedding, transfer, "
          "GO evaluation",
          "config/      config.yaml — every path and parameter",
          "slurm/       batch scripts for the GPU stages",
          "envs/        conda exports, one per environment (see "
          "envs/INDEX.md)",
          "figures/     the paper's figures, in PNG, PDF, SVG and TIFF",
          "results/     result files; every number in the paper comes "
          "from these",
          "docs/        design limits and assembly strategy",
          "```", ""]

    L += ["## Reproducing", "",
          "The pipeline runs across more than one conda environment — "
          "`envs/INDEX.md` records which environment holds which tool, "
          "because DIAMOND, the PLM stack and the differential-expression "
          "stack are not installed together.", "", "```bash"]
    for e in envs[:3]:
        L.append(f"conda create --name {e} --file envs/{e}.explicit.txt")
    L += ["```", "",
          "Then edit `config/config.yaml` to point at your own storage and "
          "run the numbered scripts in order. Scripts 31 onward reproduce "
          "every figure and table from the files already in `results/`, "
          "without re-running the GPU stages.", ""]

    if figs:
        L += ["## Figures", ""]
        for f in figs:
            L.append(f"- `figures/{f}`")
        L.append("")

    L += ["## Data availability", "",
          f"Sequencing data are from SRA study "
          f"`{meta.get('_sra', 'see config/config.yaml')}`, generated by a "
          f"previous study and not by this work.", "",
          "Derived results, figures and tables are archived on Zenodo "
          "(DOI below). This repository carries the same files, except "
          "those too large for Git.", ""]
    if excluded:
        L += ["Available only from the Zenodo record:", ""]
        L += [f"- `{p}` ({human(s)})" for p, s in excluded]
        L.append("")

    L += ["## Citation", "",
          "See `CITATION.cff`. The archived release DOI is minted on "
          "Zenodo; until the paper is accepted the record is "
          f"`{meta['zenodo']['access_right']}`.", "",
          "## Licence", "",
          f"Code is {meta['licence']['code']} (`LICENSE`). Results, "
          f"figures and tables are {meta['licence']['data']}.", "",
          "## Contact", "", f"{meta['corresponding']}", "",
          "---", "",
          f"Built from freeze `{tag}` on "
          f"{datetime.now(timezone.utc).strftime('%Y-%m-%d')}.", ""]
    return "\n".join(L)


# ---------------------------------------------------------------------------

def next_steps(meta, tag, out: Path, zip_path: Path, repo: Path) -> str:
    gh = meta["github"]
    slug = f"{gh['owner']}/{gh['repo']}"
    vis = "--private" if gh.get("private", True) else "--public"
    sand = meta["zenodo"].get("sandbox_first", True)
    host = "sandbox.zenodo.org" if sand else "zenodo.org"
    L = [f"# Release steps — freeze `{tag}`", "",
         "## 1. GitHub", "",
         "From a machine with your GitHub credentials, with the `repo` "
         "directory copied across:", "", "```bash",
         f"cd {repo.name}", "git init -b " + gh.get("default_branch", "main"),
         "git add -A",
         'git commit -m "Initial release: calibrated PLM function transfer"',
         f"gh repo create {slug} {vis} --source=. --push",
         "```", "",
         "Without the `gh` CLI, create the repository in the web interface "
         "and then:", "", "```bash",
         f"git remote add origin git@github.com:{slug}.git",
         f"git push -u origin {gh.get('default_branch', 'main')}",
         "```", "",
         "## 2. Zenodo", ""]
    if sand:
        L += ["Rehearse on the sandbox first. It behaves identically and "
              "the DOIs are throwaway, so a mistake costs nothing — which "
              "is the whole point, because a real DOI cannot be "
              "withdrawn.", ""]
    L += [f"Get a token at https://{host}/account/settings/applications/ "
          "with the `deposit:write` and `deposit:actions` scopes, then:",
          "",
          "The token is read by a prompt rather than typed into the "
          "command. Pasting a secret onto a command line leaves it in "
          "the shell history, and a command with a gap in it is a "
          "command that gets run with the gap still in it.",
          "", "```bash",
          "read -r -s -p 'Zenodo token: ' ZENODO_TOKEN; echo; "
          "export ZENODO_TOKEN",
          f"export ZENODO_HOST={host}", "",
          "# create an empty deposition and keep its id",
          'DEP=$(curl -s -X POST "https://$ZENODO_HOST/api/deposit/depositions" '
          '-H "Authorization: Bearer $ZENODO_TOKEN" '
          '-H "Content-Type: application/json" -d "{}" '
          "| python -c 'import json,sys; print(json.load(sys.stdin)[\"id\"])')",
          'echo "deposition $DEP"', "",
          "# attach the metadata",
          'curl -s -X PUT "https://$ZENODO_HOST/api/deposit/depositions/$DEP" '
          '-H "Authorization: Bearer $ZENODO_TOKEN" '
          '-H "Content-Type: application/json" '
          f'-d @zenodo/zenodo_metadata.json', "",
          "# upload the archive",
          'BUCKET=$(curl -s "https://$ZENODO_HOST/api/deposit/depositions/$DEP" '
          '-H "Authorization: Bearer $ZENODO_TOKEN" '
          "| python -c 'import json,sys; print(json.load(sys.stdin)"
          '["links"]["bucket"])\')',
          f'curl -s --upload-file zenodo/{zip_path.name} '
          f'"$BUCKET/{zip_path.name}?access_token=$ZENODO_TOKEN" > /dev/null',
          "", "# reserve the DOI without publishing",
          'curl -s "https://$ZENODO_HOST/api/deposit/depositions/$DEP" '
          '-H "Authorization: Bearer $ZENODO_TOKEN" '
          "| python -c 'import json,sys; "
          'print(json.load(sys.stdin)["metadata"]["prereserve_doi"]["doi"])\'',
          "```", "",
          "That last command prints the DOI to put in the manuscript's "
          "data-availability statement. The record stays unpublished until "
          "you press publish, and the DOI resolves once you do.", "",
          "## 3. On acceptance", "", "```bash",
          f"gh repo edit {slug} --visibility public --accept-visibility-change-consequences",
          "```", "",
          "and in the Zenodo web interface set access to Open and press "
          "Publish. Then re-run this script so `CITATION.cff` carries the "
          "final DOI, and commit that one file.", ""]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--author", action="append", default=[], metavar="NAME",
                    help='author as "Surname, Given". Repeat once per '
                         'author, in paper order.')
    ap.add_argument("--affiliation", action="append", default=[],
                    metavar="TEXT",
                    help="affiliation for the author at the same position. "
                         "Give one per --author, or one in total to apply "
                         "to all of them.")
    ap.add_argument("--orcid", action="append", default=[], metavar="ID",
                    help='ORCID for the author at the same position, as '
                         '0000-0000-0000-0000. Use "none" for an author '
                         'who has none.')
    ap.add_argument("--github-user", default=None, metavar="NAME")
    ap.add_argument("--repo", default=None, metavar="NAME",
                    help="repository name (default pmdark)")
    ap.add_argument("--email", default=None, metavar="ADDR",
                    help="corresponding address for the README")
    ap.add_argument("--public", dest="public", action="store_true",
                    default=None,
                    help="repository public and Zenodo open. Bioinformatics "
                         "requires this at submission.")
    ap.add_argument("--private", dest="public", action="store_false",
                    help="keep both closed for now")
    ap.add_argument("--envs", default=None, metavar="A,B,C",
                    help="which conda environments to ship. Default: all "
                         "captured except the conda base install. Pass a "
                         "comma-separated list to narrow it.")
    ap.add_argument("--sandbox", dest="sandbox", action="store_true",
                    default=None,
                    help="write the Zenodo steps against "
                         "sandbox.zenodo.org, to rehearse. A sandbox DOI "
                         "is throwaway; a real one cannot be withdrawn.")
    ap.add_argument("--no-sandbox", dest="sandbox", action="store_false",
                    help="write them against the real zenodo.org")
    ap.add_argument("--init", action="store_true",
                    help="write config/release.yaml and stop. Not needed "
                         "if you pass the flags above.")
    ap.add_argument("--go", action="store_true")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--max-repo-mb", type=float, default=25.0,
                    help="single files above this stay out of the "
                         "repository (default 25)")
    ap.add_argument("--doi", default=None,
                    help="final DOI, once Zenodo has minted it")
    args = ap.parse_args()

    import yaml
    cfg = cfgmod.load()
    share = cfg.p("share").resolve()
    code = cfg.p("home").resolve()
    res = cfg.p("results").resolve()
    meta_path = code / "config" / "release.yaml"

    if args.init:
        if meta_path.exists():
            print(f"{meta_path} already exists — not overwriting.")
            return 1
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta_path.write_text(STUB)
        print(f"wrote {meta_path}\n\nFill in every REPLACE_ME, then re-run.")
        return 0

    # Metadata comes from flags; release.yaml is where it is kept so the
    # flags only have to be typed once. No editor is involved.
    meta = (yaml.safe_load(meta_path.read_text()) if meta_path.exists()
            else yaml.safe_load(STUB))
    changed = apply_flags(meta, args)
    if changed:
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        save_meta(meta_path, meta)
        log(f"saved {len(changed)} setting(s) to {meta_path}")

    gaps = find_placeholders(meta)
    if gaps:
        print(f"\nrefusing: {len(gaps)} thing(s) not set yet.\n",
              file=sys.stderr)
        for g in gaps:
            key = re.sub(r"\[\d+\]", "[]", g)
            print(f"  {g:<26}{FLAG_FOR.get(key, '')}", file=sys.stderr)
        print("\nPass them as flags, for example:\n", file=sys.stderr)
        print('  python scripts/46_release.py \\\n'
              '      --github-user YOURNAME \\\n'
              '      --author "Surname, Given" '
              '--affiliation "Your University" \\\n'
              '      --orcid 0000-0000-0000-0000 \\\n'
              '      --email you@kmutt.ac.th --public\n', file=sys.stderr)
        print("Repeat --author/--affiliation/--orcid once per author, in "
              "paper order.\nThey are saved, so you only type them once.\n"
              "A Zenodo DOI is permanent, which is why nothing is guessed.",
              file=sys.stderr)
        return 2

    # Text that came from the example command rather than from a person.
    suspect: list[tuple[str, str, str]] = []
    for label, val, flag in (
            ("github.owner", meta["github"].get("owner", ""),
             "--github-user"),
            ("github.repo", meta["github"].get("repo", ""), "--repo"),
            ("corresponding", meta.get("corresponding", ""), "--email")):
        if isinstance(val, str) and LOOKS_LIKE_EXAMPLE.search(val):
            suspect.append((label, val, flag))
    for i, a in enumerate(meta["authors"]):
        for key, flag in (("name", "--author"),
                          ("affiliation", "--affiliation")):
            v = a.get(key) or ""
            if LOOKS_LIKE_EXAMPLE.search(v):
                suspect.append((f"authors[{i}].{key}", v, flag))
        oc = (a.get("orcid") or "").strip()
        if oc and not orcid_valid(oc):
            suspect.append((f"authors[{i}].orcid", oc, "--orcid"))
    if suspect:
        print(f"\nrefusing: {len(suspect)} value(s) look like the example "
              f"text, not real details.\n", file=sys.stderr)
        for label, val, flag in suspect:
            print(f"  {label:<24}{val!r:<42}set with {flag}",
                  file=sys.stderr)
        print("\nThese came from the example command in the instructions. "
              "Re-run with\nyour own details — the same flags, real "
              "values. Nothing was built.\n"
              "A Zenodo DOI is permanent; this is the last point at which "
              "it is free to fix.", file=sys.stderr)
        return 2

    title = cfg.raw.get("paper", {}).get("title") or (
        "A structurally coherent, STAT-regulated dark proteome in "
        "Penaeus monodon hemocytes")
    meta["_sra"] = cfg.raw.get("data", {}).get("sra_study", "")
    tag = args.tag or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = share / "97_release" / tag
    repo = out / "repo"
    zdir = out / "zenodo"

    # ---- choose what goes where ----------------------------------------
    cap = args.max_repo_mb * (1 << 20)
    plan: list[tuple[Path, str]] = []        # (source, repo-relative)
    for sub in ("scripts", "src", "config", "slurm", "docs"):
        d = code / sub
        if not d.is_dir():
            continue
        for p in sorted(d.rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts:
                plan.append((p, str(Path(sub) / p.relative_to(d))))

    # Walk the results tree in full. Globbing only the top level and
    # figures/ silently dropped every other subdirectory — the freeze
    # counted 156 files where the release staged 76, and a deposit that
    # omits the files behind the numbers cannot be checked against the
    # paper, which is the only reason the deposit exists.
    excluded: list[tuple[str, int]] = []
    deposit_only: list[tuple[Path, str]] = []
    figs: list[str] = []
    n_seen = 0
    # Directories whose names begin with "_" are working scratch — the
    # ones this project has held superseded results, e.g. the copy taken
    # before the I20 tier was rebuilt. A DOI is permanent and citable,
    # so shipping a stale results file beside the real one invites
    # somebody to quote a number the paper withdrew. Skipped, counted,
    # and named in the summary rather than dropped quietly.
    skipped_dirs: dict[str, int] = {}
    for p in sorted(res.rglob("*")):
        if not p.is_file() or p.is_symlink():
            continue
        if "__pycache__" in p.parts:
            continue
        scratch = next((part for part in p.relative_to(res).parts[:-1]
                        if part.startswith("_")), None)
        if scratch:
            skipped_dirs[scratch] = skipped_dirs.get(scratch, 0) + 1
            continue
        n_seen += 1
        rel = p.relative_to(res)
        dest = ("figures/" + str(rel.relative_to("figures"))
                if rel.parts and rel.parts[0] == "figures"
                else f"results/{rel}")
        if p.stat().st_size > cap:
            # Out of the repository, but NOT out of the deposit. The
            # source path is kept so the archive can carry it: the
            # README tells readers these files are "available only from
            # the Zenodo record", and before this the archive was built
            # from the repository directory alone, so the cap dropped
            # them from both and the record quietly failed the promise
            # the README makes for it.
            excluded.append((dest, p.stat().st_size))
            deposit_only.append((p, dest))
            continue
        plan.append((p, dest))
        if rel.parts and rel.parts[0] == "figures" and len(rel.parts) == 2 \
                and p.suffix.lower() in (".png", ".pdf", ".svg", ".tif",
                                         ".docx", ".md"):
            figs.append(p.name)

    envs: list[str] = []
    fr = share / "98_freeze"
    envsrc = None
    for d in sorted((x for x in fr.glob("*") if x.is_dir()),
                    key=lambda x: x.name, reverse=True):
        if (d / "ENVIRONMENTS" / "INDEX.md").is_file():
            envsrc = d / "ENVIRONMENTS"
            break
    BASE_ENVS = {"miniconda3", "anaconda3", "miniforge3", "base",
                 "mambaforge"}
    want = ({e.strip() for e in args.envs.split(",") if e.strip()}
            if args.envs else None)
    skipped_envs: list[str] = []
    if envsrc:
        for p in sorted(envsrc.glob("*")):
            if not p.is_file():
                continue
            stem = p.name.split(".")[0]
            if p.name != "INDEX.md":
                # The conda base install is not a project environment,
                # and an unrelated project's env is noise in a citation
                # artefact.
                if want is not None and stem not in want:
                    if stem not in skipped_envs:
                        skipped_envs.append(stem)
                    continue
                if want is None and stem in BASE_ENVS:
                    if stem not in skipped_envs:
                        skipped_envs.append(stem)
                    continue
            plan.append((p, f"envs/{p.name}"))
            if p.name.endswith(".explicit.txt"):
                envs.append(p.name[:-len(".explicit.txt")])

    total = sum(p.stat().st_size for p, _ in plan)
    print("=" * 78)
    print(f"RELEASE PLAN — {tag}")
    print("=" * 78)
    by: dict[str, list[int]] = {}
    for p, rel in plan:
        k = rel.split("/")[0]
        b = by.setdefault(k, [0, 0])
        b[0] += 1
        b[1] += p.stat().st_size
    for k in sorted(by):
        print(f"  {k:<12}{by[k][0]:>6} files   {human(by[k][1]):>12}")
    print("  " + "-" * 60)
    print(f"  {'repo total':<12}{len(plan):>6} files   {human(total):>12}")
    if excluded:
        print(f"\n  too large for Git (> {args.max_repo_mb:g} MB each), "
              f"Zenodo only:")
        for name, sz in excluded:
            print(f"    {human(sz):>12}  {name}")
    print(f"\n  environments : {', '.join(envs) or 'none found'}")
    if skipped_envs:
        print(f"  not shipped  : {', '.join(skipped_envs)}"
              f"   (--envs to change)")
    print(f"  results tree : {n_seen} files seen under {res.name}")
    if skipped_dirs:
        for nm, cnt in sorted(skipped_dirs.items()):
            print(f"    skipped    : {cnt} file(s) in {nm}/ "
                  f"(scratch; a DOI is permanent and citable)")
    if not envsrc:
        print("    WARNING: no ENVIRONMENTS capture found under 98_freeze.")
        print("    Run scripts/45_capture_envs.py --into LATEST first.")
    print(f"  github       : {meta['github']['owner']}/"
          f"{meta['github']['repo']}  "
          f"({'private' if meta['github'].get('private', True) else 'public'})")
    print(f"  zenodo       : {meta['zenodo']['mode']}, "
          f"{meta['zenodo']['access_right']}"
          f"{', sandbox first' if meta['zenodo'].get('sandbox_first') else ''}")
    print(f"  authors      : {len(meta['authors'])}")
    print(f"  output       : {out}")
    print("=" * 78)

    if not args.go:
        print("\nPLAN ONLY — nothing written. Re-run with --go.")
        return 0

    repo.mkdir(parents=True, exist_ok=True)
    zdir.mkdir(parents=True, exist_ok=True)
    for src, rel in plan:
        dst = repo / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    log(f"staged {len(plan)} files into {repo}")

    rows = readme_numbers(res)
    holders = "; ".join(a["name"] for a in meta["authors"])
    (repo / "LICENSE").write_text(
        MIT.format(year=date.today().year, holders=holders)
        if meta["licence"]["code"].upper() == "MIT"
        else f"{meta['licence']['code']}\n\nCopyright (c) "
             f"{date.today().year} {holders}\n")
    (repo / ".gitignore").write_text(GITIGNORE)
    (repo / "README.md").write_text(
        readme(meta, title, rows, figs, excluded, envs, tag))
    version = f"0.1.0-{tag[:8]}"
    (repo / "CITATION.cff").write_text(cff(meta, title, version, args.doi))
    desc = ("<p>Code and derived results for a calibrated "
            "protein-language-model approach to protein function "
            "annotation, evaluated under progressive withdrawal of "
            "homologous evidence and applied to the STAT-dependent "
            "uncharacterised fraction of the <i>Penaeus monodon</i> "
            "hemocyte proteome.</p>")
    if rows:
        desc += "<ul>" + "".join(f"<li>{k}: {v}</li>" for k, v in rows) \
                + "</ul>"
    (zdir / "zenodo_metadata.json").write_text(
        json.dumps({"metadata": json.loads(
            zenodo_json(meta, title, desc))}, indent=2, ensure_ascii=False)
        + "\n")
    log("wrote README.md, LICENSE, CITATION.cff, .gitignore, "
        "zenodo_metadata.json")
    if not rows:
        log("WARNING: no result files found — README has no numbers table")

    # ---- the deposit ----------------------------------------------------
    zip_base = zdir / f"pmdark_{tag}"
    log("building the Zenodo archive")
    shutil.make_archive(str(zip_base), "zip", root_dir=repo)
    zip_path = Path(f"{zip_base}.zip")
    # The deposit is the complete record; the repository is the part of
    # it Git should carry. Everything the cap held back is appended here
    # at the same path it would have had in the repository, so the two
    # differ only in what Git would struggle with.
    if deposit_only:
        import zipfile as _zf
        with _zf.ZipFile(zip_path, "a", _zf.ZIP_DEFLATED) as z:
            for src, rel in deposit_only:
                z.write(src, rel)
        log(f"appended {len(deposit_only)} file(s) held back from the "
            f"repository into the archive")
    import hashlib
    h = hashlib.sha256()
    with zip_path.open("rb") as fh:
        while c := fh.read(1 << 20):
            h.update(c)
    (zdir / "ARCHIVE.sha256").write_text(f"{h.hexdigest()}  {zip_path.name}\n")

    import zipfile
    with zipfile.ZipFile(zip_path) as z:
        bad = z.testzip()
        n = len(z.namelist())
    if bad:
        log(f"ARCHIVE CORRUPT at {bad}")
        return 1
    log(f"{zip_path.name}  {human(zip_path.stat().st_size)}  "
        f"{n} entries, CRC verified")

    (out / "NEXT_STEPS.md").write_text(
        next_steps(meta, tag, out, zip_path, repo))
    log("wrote NEXT_STEPS.md")

    print("\n" + "=" * 78)
    print(f"RELEASE BUILT  {out}")
    print("=" * 78)
    print(f"  repo/              {len(plan) + 4} files, "
          f"{human(total)} plus generated")
    print(f"  zenodo/            {zip_path.name}  "
          f"{human(zip_path.stat().st_size)}")
    print(f"  NEXT_STEPS.md      git and Zenodo commands, filled in")
    print(f"\n  read {out}/NEXT_STEPS.md and work through it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
