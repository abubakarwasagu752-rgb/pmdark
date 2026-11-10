#!/usr/bin/env python3
"""Progress of a fan-out embedding run across all partitions.

    python scripts/embed_status.py --model esm2
    python scripts/embed_status.py --model esm2 --watch
    python scripts/embed_status.py --model esm2 --clear-stale
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod          # noqa: E402
from pmdark.claim import ClaimRegistry        # noqa: E402
from pmdark.shards import fasta_checksum, plan_shards, read_fasta  # noqa: E402


def render(out_dir: Path, n_shards: int, stale_after: float) -> str:
    is_done = lambda i: (out_dir / f"shard_{i:05d}.done").exists()  # noqa: E731
    reg = ClaimRegistry(out_dir, stale_after=stale_after)
    st = reg.status(n_shards, is_done)

    width = 46
    filled = int(width * st["done"] / max(1, st["total"]))
    active = int(width * st["in_progress"] / max(1, st["total"]))
    bar = "#" * filled + "+" * active + "." * (width - filled - active)
    pct = 100.0 * st["done"] / max(1, st["total"])

    lines = [
        f"[{bar}] {pct:5.1f}%",
        f"  done {st['done']}/{st['total']}   running {st['in_progress']}   "
        f"stale {st['stale']}   unclaimed {st['unclaimed']}",
    ]
    if st["workers"]:
        lines.append("  active workers:")
        for key, n in sorted(st["workers"].items(), key=lambda kv: -kv[1]):
            lines.append(f"    {key:<34} {n} shard(s)")
    if st["stale"]:
        lines.append(f"  {st['stale']} stale claim(s) — workers killed by the "
                     "walltime cap; they are reclaimed automatically")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--fasta", type=Path, default=None)
    ap.add_argument("--stale-after", type=float, default=3600.0)
    ap.add_argument("--watch", action="store_true")
    ap.add_argument("--interval", type=float, default=30.0)
    ap.add_argument("--clear-stale", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    fasta = args.fasta or (cfg.p("proteins") / "benchmark" / "reference.faa")
    # The worker writes to embeddings/<model>/<fasta stem>/ so query and
    # reference do not collide. Mirror that here or the bar reads 0%.
    out_dir = cfg.p("embeddings", args.model, Path(fasta).stem)
    if not out_dir.exists():
        print(f"no shard directory yet: {out_dir}")
        return 1

    if fasta.exists():
        records = read_fasta(fasta)
        n_shards = len(plan_shards(records, cfg["compute"]["shard_size"]))
        ck = fasta_checksum(records)
    else:
        n_shards = len(list(out_dir.glob("shard_*.npy"))) or 1
        ck = "?"
        print(f"(FASTA not found at {fasta}; inferring shard count)")

    if args.clear_stale:
        n = ClaimRegistry(out_dir, stale_after=args.stale_after).clear_stale()
        print(f"cleared {n} stale claim(s)")
        return 0

    header = f"model {args.model}   shards {n_shards}   input {ck}\n{out_dir}"

    if not args.watch:
        print(header)
        print(render(out_dir, n_shards, args.stale_after))
        return 0

    try:
        while True:
            print("\033[2J\033[H", end="")
            print(time.strftime("%F %T"))
            print(header)
            print(render(out_dir, n_shards, args.stale_after))
            time.sleep(args.interval)
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
