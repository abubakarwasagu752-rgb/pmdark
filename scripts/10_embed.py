#!/usr/bin/env python3
"""Embed a protein FASTA with one PLM, resumably, under a hard walltime.

Typical use (as a Slurm array task):

    python scripts/10_embed.py \
        --fasta  $SHARE/01_interim/all_proteins.faa \
        --model  esm2 \
        --out    $SHARE/04_embeddings/esm2 \
        --walltime-min 720

Re-running the same command picks up exactly where it stopped.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                      # noqa: E402
from pmdark.embed import (                                # noqa: E402
    Deadline,
    EmbedConfig,
    Embedder,
    autoscale_batch,
    gpu_report,
    run_shards,
    slice_for_task,
)
from pmdark.shards import ShardStore, fasta_checksum, plan_shards, read_fasta  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--model", required=True,
                    help="key in config models: esm2 | esm2_3b | prott5 | prostt5 | saprot")
    ap.add_argument("--out", required=True, type=Path,
                    help="shard output directory")
    ap.add_argument("--walltime-min", type=float, default=720.0,
                    help="job walltime in minutes (default 720 = 12 h cap)")
    ap.add_argument("--margin-min", type=float, default=12.0)
    ap.add_argument("--batch", type=int, default=None, help="override batch size")
    ap.add_argument("--max-tokens", type=int, default=16384)
    ap.add_argument("--shard-size", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="plan only: print shard counts and exit")
    args = ap.parse_args()

    cfg = cfgmod.load()
    models = cfg["models"]
    if args.model not in models:
        ap.error(f"--model must be one of {sorted(models)}")

    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    n_tasks = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", 1))

    records = read_fasta(args.fasta)
    checksum = fasta_checksum(records)
    shard_size = args.shard_size or cfg["compute"]["shard_size"]
    shards = plan_shards(records, shard_size)

    print(f"fasta      : {args.fasta}")
    print(f"sequences  : {len(records):,}  (checksum {checksum})")
    print(f"shards     : {len(shards)} x {shard_size}")
    print(f"model      : {args.model} -> {models[args.model]}")
    print(f"array task : {task_id + 1}/{n_tasks}")
    print(f"device     : {gpu_report()}")

    if args.dry_run:
        lengths = [r.length for r in records]
        print(f"len min/median/max: {min(lengths)}/"
              f"{sorted(lengths)[len(lengths)//2]}/{max(lengths)}")
        return 0

    store = ShardStore(args.out, checksum=checksum, n_shards=len(shards))
    assigned = slice_for_task(len(shards), task_id, n_tasks)
    pending = store.pending(assigned)
    print(f"assigned   : {len(assigned)} shards, {len(pending)} pending")

    if not pending:
        done, total = store.progress()
        print(f"nothing to do — {done}/{total} complete overall")
        return 0

    default_batch = cfg["compute"]["batch"].get(args.model, {}).get("l40", 32)
    batch = args.batch or autoscale_batch(default_batch)
    print(f"batch size : {batch}")

    embedder = Embedder(EmbedConfig(
        model_key=args.model,
        model_id=models[args.model],
        cache_dir=cfg.p("models"),
        max_len=cfg["compute"]["max_len"],
        max_batch=batch,
        max_tokens=args.max_tokens,
        dtype=cfg["compute"]["dtype"],
    ))

    deadline = Deadline(args.walltime_min, args.margin_min)
    stats = run_shards(shards, store, embedder, pending, deadline)

    done, total = store.progress()
    print(json.dumps({
        "completed": len(stats["completed"]),
        "skipped": len(stats["skipped"]),
        "stopped_early": stats["stopped_early"],
        "global_progress": f"{done}/{total}",
    }, indent=2))

    # Exit 99 signals "resubmit me" to the wrapper script
    return 99 if stats["stopped_early"] or done < total else 0


if __name__ == "__main__":
    raise SystemExit(main())
