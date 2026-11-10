#!/usr/bin/env python3
"""Find and free space on the quota'd work tier. Reports by default.

    python scripts/34_space.py                 # measure only, change nothing
    python scripts/34_space.py --free           # show the plan, still no action
    python scripts/34_space.py --free --yes     # actually do it

WHY

  Every ProtT5 and ProstT5 worker died on its first shard write:

      OSError: 2048000 requested and 0 written
        np.save(fh, matrix)   src/pmdark/shards.py:156

  "0 written" is the filesystem refusing, not numpy failing. $HOME on
  /storage/slurm carries a ~100 GB quota shared with five other projects,
  and `df` reports the whole 9.5 TB filesystem, so it looked like there
  were 6 TB free the entire time. ESM-2 embedded successfully before the
  tier also had to hold two more ~2.8 GB model caches, a 3.3 GB HISAT2
  index, a decompressed Swiss-Prot, and the re-alignment outputs.

WHAT THIS CLASSIFIES

  REGENERABLE  rebuilt by a documented command in seconds or minutes.
               Safe to delete; the command is printed next to each one.
  ARCHIVABLE   needed later but not by a compute node, so it belongs on
               the archive tier (login-visible, 6 TB free) rather than
               under quota. Moved, never deleted.
  KEEP         required by work in progress or expensive to rebuild.

  Nothing outside this project is touched. The archive destination is only
  ever this project's own directory — the other five projects on that
  volume are never read, moved or deleted.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

# relative to the work tier
PLAN = [
    ("align/_unmapped", "REGENERABLE",
     "unmapped reads for the deferred de novo arm; "
     "re-fetch with scripts/21_align_assemble.py --save-unmapped"),
    ("proteins/transdecoder/uniprot_sprot.fasta", "REGENERABLE",
     "decompressed Swiss-Prot; 24_call_orfs.py re-expands it from the "
     "archive .gz on demand"),
    ("proteins/transdecoder/homology_both.faa", "REGENERABLE",
     "concatenated homology DB; 24_call_orfs.py rebuilds it"),
    ("proteins/ref_proteome.faa", "REGENERABLE",
     "staged copy of the reference proteome; 27_recovery_exact.py restages"),
    ("proteins/merge_sweep", "REGENERABLE",
     "merge-sweep GTFs; 23b_merge_sweep.py regenerates (minutes)"),
    ("proteins/compare", "REGENERABLE",
     "diamond databases for 31_compare_protein_sets.py"),
    ("genome/pmon.*.ht2", "KEEP",
     "HISAT2 index, 3.3 GB, 17 minutes to rebuild — keep unless desperate"),
    ("models", "KEEP",
     "HF weight cache for all PLMs; deleting forces a re-download"),
    ("embeddings", "KEEP", "the actual results"),
    ("proteins/transdecoder", "KEEP",
     "the predicted proteome and its genome GFF3"),
]

ARCHIVE_CANDIDATES = [
    ("align", "per-library GTFs and abundance tables — only "
              "30_expressed_denominator.py reads them, on the login node"),
]


def human(n: float) -> str:
    for u in ("B", "K", "M", "G", "T"):
        if n < 1024 or u == "T":
            return f"{n:.1f}{u}"
        n /= 1024
    return f"{n:.1f}T"


def usage(p: Path) -> int:
    if not p.exists():
        return 0
    if p.is_file():
        return p.stat().st_size
    tot = 0
    for f in p.rglob("*"):
        try:
            if f.is_file():
                tot += f.stat().st_size
        except OSError:
            pass
    return tot


def quota_report() -> None:
    for cmd in (["quota", "-s"], ["lfs", "quota", "-h", str(Path.home())]):
        if shutil.which(cmd[0]):
            print(f"\n$ {' '.join(cmd)}")
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=60)
                print((r.stdout or r.stderr).rstrip() or "  (no output)")
            except (OSError, subprocess.TimeoutExpired) as e:
                print(f"  failed: {e}")
            return
    print("\nNo quota command available. `df` shows the whole filesystem,")
    print("not your allowance, which is why this was invisible.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--free", action="store_true",
                    help="print the plan; combine with --yes to execute")
    ap.add_argument("--yes", action="store_true")
    ap.add_argument("--archive-align", action="store_true",
                    help="also move align/ to the archive tier")
    args = ap.parse_args()

    cfg = cfgmod.load()
    work = Path(str(cfg["paths"]["work"]).replace("~", str(Path.home())))
    share = cfg.p("share")

    print("=" * 78)
    print(f"WORK TIER  {work}")
    print("=" * 78)
    quota_report()

    print("\nTop-level usage:")
    rows = []
    for child in sorted(work.iterdir()) if work.exists() else ():
        n = usage(child)
        rows.append((n, child.name))
    for n, name in sorted(rows, reverse=True):
        print(f"  {human(n):>9}  {name}")
    print(f"  {human(sum(n for n, _ in rows)):>9}  TOTAL")

    print("\n" + "=" * 78)
    print("CLASSIFIED")
    print("=" * 78)
    freeable = 0
    actions: list[tuple[Path, str]] = []
    for rel, kind, why in PLAN:
        if "*" in rel:
            paths = sorted(work.glob(rel))
            n = sum(usage(p) for p in paths)
        else:
            paths = [work / rel]
            n = usage(paths[0])
        if n == 0:
            continue
        print(f"  [{kind:<11}] {human(n):>9}  {rel}")
        print(f"                {why}")
        if kind == "REGENERABLE":
            freeable += n
            actions += [(p, "delete") for p in paths if p.exists()]

    if args.archive_align:
        for rel, why in ARCHIVE_CANDIDATES:
            p = work / rel
            n = usage(p)
            if n:
                print(f"  [ARCHIVABLE  ] {human(n):>9}  {rel}")
                print(f"                {why}")
                freeable += n
                actions.append((p, "archive"))

    print("-" * 78)
    print(f"  recoverable now: {human(freeable)}")
    print("\n  Embeddings themselves are small — 19,304 x 1024 floats is")
    print("  ~79 MB per model for the reference set, ~4.5 MB for the query.")
    print("  The tier does not need much; it needs some.")

    if not args.free:
        print("\n(measurement only; --free to see the plan)")
        return 0

    print("\n" + "=" * 78)
    print("PLAN" + ("" if args.yes else "  (dry run — add --yes to execute)"))
    print("=" * 78)
    dest = share / "90_offloaded"
    for p, how in actions:
        if how == "delete":
            print(f"  delete   {p}  ({human(usage(p))})")
        else:
            print(f"  move     {p}  ->  {dest / p.name}  ({human(usage(p))})")
    if not args.yes:
        return 0

    dest.mkdir(parents=True, exist_ok=True)
    done = 0
    for p, how in actions:
        n = usage(p)
        try:
            if how == "delete":
                shutil.rmtree(p) if p.is_dir() else p.unlink()
            else:
                target = dest / p.name
                if target.exists():
                    shutil.rmtree(target) if target.is_dir() else target.unlink()
                shutil.move(str(p), str(target))
            done += n
            print(f"  ok  {p.name}  (+{human(n)})")
        except OSError as e:
            print(f"  FAILED  {p}: {e}")
    print(f"\nfreed {human(done)}")
    print("\nNow re-run the embedding fan-outs. Completed shards are kept, so")
    print("nothing already computed is lost.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
