#!/usr/bin/env python3
"""Pre-registered hyperparameter selection for the calibration benchmark.

    python scripts/32_sweep.py --model esm2
    python scripts/32_sweep.py --model esm2 --transform whiten

WHY THIS IS A SEPARATE SCRIPT

  k, sim_floor and power were chosen on general grounds, not fitted to this
  data. The first run showed that matters: sim_floor=0.5 is sensible for
  centred embeddings but far too high for whitened ones, where the median
  top similarity at the hardest tier is 0.292 — so 1,002 of 1,087 queries
  got no prediction at all and Fmax collapsed to 0.093. That is a
  mis-specified threshold, not a property of the method.

  But tuning against the hard tiers and then REPORTING the hard tiers is
  circular, and a reviewer will say so. So:

      SELECTION  happens on T0_full and I50 only — the easy tiers.
      The chosen values are then FROZEN and written to disk.
      EVALUATION on I30 / Inone / ORPHAN uses those frozen values and is
      never allowed to influence them.

  This script does the selection half and records it. 31_calibrate.py then
  reads the frozen parameters, so the split is enforced by the code rather
  than by discipline.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                          # noqa: E402
from pmdark.goeval import GOGraph                             # noqa: E402
from pmdark.transfer import center, whiten                    # noqa: E402

# Selection runs ONLY on these. Never add a hard tier here.
SELECTION_TIERS = ("T0_full", "I50")

GRID_K = (10, 25, 50, 100)
GRID_FLOOR = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
GRID_POWER = (1.0, 2.0, 3.0, 5.0)


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="esm2")
    ap.add_argument("--transform", default="center",
                    choices=["none", "center", "whiten"])
    ap.add_argument("--min-coverage", type=float, default=0.90,
                    help="reject settings that refuse to predict for more "
                         "than this fraction of queries; a high Fmax on 8%% "
                         "of proteins is not a usable annotator")
    ap.add_argument("--block", type=int, default=512)
    args = ap.parse_args()

    # Import here so the module docstring above stays the first thing read.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "cal", Path(__file__).resolve().parent / "31_calibrate.py")
    cal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cal)

    cfg = cfgmod.load()
    bench = cfg.p("proteins", "benchmark")

    with gzip.open(bench / "benchmark.json.gz", "rt") as fh:
        meta = json.load(fh)
    ident_path = bench / "identity_tiers.json.gz"
    with gzip.open(ident_path, "rt") as fh:
        ident_excl = json.load(fh)["exclusions"]

    graph = GOGraph.from_obo(cfg.p("raw", "benchmark") / "go-basic.obo")
    q_ids, q_emb = cal.load_embeddings(cfg, args.model, bench / "query.faa")
    r_ids, r_emb = cal.load_embeddings(cfg, args.model, bench / "reference.faa")

    if args.transform == "center":
        r_emb, mu = center(r_emb.astype(np.float32))
        q_emb, _ = center(q_emb.astype(np.float32), mu)
    elif args.transform == "whiten":
        r_emb, st = whiten(r_emb.astype(np.float32))
        q_emb, _ = whiten(q_emb.astype(np.float32), st)

    q_meta, r_meta = meta["query"], meta["reference"]
    q_ids = [a for a in q_ids if a in q_meta]
    r_ids = [a for a in r_ids if a in r_meta]

    log(f"model {args.model}  transform {args.transform}")
    log(f"selection tiers: {', '.join(SELECTION_TIERS)}  "
        f"(hard tiers deliberately excluded)")
    log(f"grid: {len(GRID_K)}x{len(GRID_FLOOR)}x{len(GRID_POWER)} = "
        f"{len(GRID_K)*len(GRID_FLOOR)*len(GRID_POWER)} settings "
        f"x {len(SELECTION_TIERS)} tiers")

    rows = []
    total = len(GRID_K) * len(GRID_FLOOR) * len(GRID_POWER)
    n = 0
    for k in GRID_K:
        for floor in GRID_FLOOR:
            for power in GRID_POWER:
                n += 1
                scores, covs = [], []
                for tier in SELECTION_TIERS:
                    res, diag = cal.run_tier(
                        tier, q_ids, q_emb, q_meta, r_ids, r_emb, r_meta,
                        graph, k, floor, power, args.block, ident_excl)
                    scores.append(res["overall"]["fmax"])
                    covs.append(res["overall"]["coverage"])
                row = {
                    "k": k, "sim_floor": floor, "power": power,
                    "mean_fmax": float(np.mean(scores)),
                    "min_coverage": float(min(covs)),
                    "per_tier": dict(zip(SELECTION_TIERS, scores)),
                }
                rows.append(row)
                if n % 10 == 0 or n == total:
                    log(f"  {n}/{total}  best so far "
                        f"{max(r['mean_fmax'] for r in rows if r['min_coverage'] >= args.min_coverage) if any(r['min_coverage'] >= args.min_coverage for r in rows) else float('nan'):.4f}")

    usable = [r for r in rows if r["min_coverage"] >= args.min_coverage]
    if not usable:
        log(f"WARNING: no setting reached coverage >= {args.min_coverage}; "
            "relaxing the constraint")
        usable = rows
    best = max(usable, key=lambda r: r["mean_fmax"])

    print("\n" + "=" * 74)
    print(f"SELECTION — {args.model} / {args.transform}")
    print("=" * 74)
    print(f"{'k':>5}{'floor':>8}{'power':>7}{'meanF':>9}{'minCov':>9}   tiers")
    print("-" * 74)
    for r in sorted(usable, key=lambda r: -r["mean_fmax"])[:12]:
        mark = " <-- selected" if r is best else ""
        tiers = "  ".join(f"{t}={v:.3f}" for t, v in r["per_tier"].items())
        print(f"{r['k']:>5}{r['sim_floor']:>8.2f}{r['power']:>7.1f}"
              f"{r['mean_fmax']:>9.4f}{r['min_coverage']:>9.2f}   {tiers}{mark}")
    print("-" * 74)

    frozen = cfg.p("results") / f"params_{args.model}_{args.transform}.json"
    frozen.parent.mkdir(parents=True, exist_ok=True)
    frozen.write_text(json.dumps({
        "model": args.model,
        "transform": args.transform,
        "selected_on": list(SELECTION_TIERS),
        "min_coverage": args.min_coverage,
        "k": best["k"], "sim_floor": best["sim_floor"], "power": best["power"],
        "selection_mean_fmax": best["mean_fmax"],
        "grid": {"k": list(GRID_K), "sim_floor": list(GRID_FLOOR),
                 "power": list(GRID_POWER)},
        "all_rows": rows,
    }, indent=2))

    print(f"\nFROZEN: k={best['k']}  sim_floor={best['sim_floor']}  "
          f"power={best['power']}")
    print(f"wrote {frozen}")
    print("\nThese were selected on T0_full and I50 ONLY. Evaluate with:")
    print(f"  python scripts/31_calibrate.py --model {args.model} "
          f"--transform {args.transform} \\")
    print(f"      --k {best['k']} --sim-floor {best['sim_floor']} "
          f"--power {best['power']} \\")
    print("      --tiers T0_full,I50,I30,Inone")
    print("\nReport the frozen values in the methods and state that the hard")
    print("tiers played no part in choosing them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
