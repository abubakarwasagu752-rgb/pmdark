#!/usr/bin/env python3
"""Fmax on the structurally hardened orphan subset.

    python scripts/50_fmax_structural.py
    python scripts/50_fmax_structural.py --model prott5

WHY THIS EXISTS

  The manuscript reports Fmax 0.625 [0.542, 0.706] on the 70 proteins
  with no hit from DIAMOND or jackhmmer. 47_structure.py then found that
  18 of those 70 have a structural relative in the same lookup set, so
  the set with no relative by ANY of the three searches is 52, not 70,
  and 0.625 is the accuracy of a set we no longer claim.

  This recomputes it on the 52 and, in the same run and the same table,
  on the 70 and the 94, so the three are directly comparable and the
  attrition is visible rather than asserted.

  Expect the number to FALL. The 18 removed are the ones with detectable
  relatives, which are the easier ones; dropping them leaves a harder
  set. A figure that went up would mean something was wrong, and the
  script says so rather than leaving the reader to notice.

HOW IT SCORES

  Through 36_bootstrap_baseline.py's own Bench.register_subset and
  predictions path — the same code that produced every other number in
  the paper. Nothing about transfer, thresholds or the bootstrap is
  reimplemented here, because a subset scored by a second implementation
  would not be comparable with the first.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

import numpy as np                                               # noqa: E402

from pmdark import config as cfgmod                              # noqa: E402


def script_roots() -> list[Path]:
    """Every directory that might hold the numbered scripts.

    The project has two: the code tree named by paths.home in
    config.yaml, and the copy under the shared project directory. New
    scripts have been landing in the second and the older ones live in
    the first, so a script that assumes its own directory holds its
    dependencies fails depending on where it was put.
    """
    roots = [HERE]
    try:
        cfg = cfgmod.load()
        for key in ("home", "share"):
            d = cfg.p(key) / "scripts"
            if d.is_dir() and d not in roots:
                roots.append(d)
    except Exception:
        pass
    return roots


def load_mod(name: str, filename: str):
    tried = []
    for root in script_roots():
        path = root / filename
        tried.append(str(path))
        if not path.is_file():
            continue
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            continue
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        if root != HERE:
            print(f"  (loaded {filename} from {root})")
        return mod
    raise SystemExit(
        f"cannot find {filename}. Looked in:\n  " + "\n  ".join(tried) +
        "\n\nThe numbered scripts are split across two directories. "
        "Copy the missing one\nalongside this script, or run this script "
        "from the directory that has it.")


def main() -> int:
    boot = load_mod("pmdark_boot", "36_bootstrap_baseline.py")
    ap = boot.build_parser(__doc__)
    ap.add_argument("--model", default="esm2")
    ap.add_argument("--transform", default="whiten")
    ap.add_argument("--ids", default=None,
                    help="identifier list of the structurally hardened "
                         "set (default: results/structure_hardened_ids.txt)")
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results")
    calib = boot.load_calib()
    bench = boot.Bench(cfg, calib)

    ids_path = Path(args.ids) if args.ids else \
        res / "structure_hardened_ids.txt"
    if not ids_path.exists():
        print(f"missing {ids_path}\nRun scripts/48_structure_stats.py "
              f"first — it writes the surviving identifiers.",
              file=sys.stderr)
        return 2
    struct_ids = [l.strip() for l in ids_path.read_text().split()
                  if l.strip()]

    orphans = bench.orphans()

    # The jackhmmer-hardened 70. If 37 did not leave an identifier list,
    # derive it exactly as 47 does — orphans with no profile hit at the
    # transfer threshold, iteration suffix stripped — rather than drop
    # the middle row and leave the attrition looking like one step.
    hard_path = res / "orphan_hardening_ids.txt"
    jack_ids = None
    if hard_path.exists():
        jack_ids = [l.strip() for l in hard_path.read_text().split()
                    if l.strip()]
    else:
        st = load_mod("pmdark_struct", "47_structure.py")
        oc = cfg.p("proteins") / "benchmark" / "orphan_check"
        tbl = st.find_tbl(oc) if oc.is_dir() else None
        if tbl:
            hit, diag = st.parse_tbl(tbl, args.transfer_evalue)
            jack_ids = [i for i in orphans if i not in hit]
            print(f"  derived the jackhmmer-hardened set from "
                  f"{tbl.name}: {len(orphans)} - "
                  f"{len(orphans) - len(jack_ids)} = {len(jack_ids)}"
                  f"   (self-hits dropped: {diag['self_rows']})")

    # Every identifier must be a real orphan. A set assembled from the
    # wrong namespace would score fine and mean nothing.
    stray = [i for i in struct_ids if i not in set(orphans)]
    if stray:
        print(f"refusing: {len(stray)} of {len(struct_ids)} identifiers "
              f"are not in the orphan set; first {stray[0]}",
              file=sys.stderr)
        return 2

    print("=" * 78)
    print(f"Fmax ACROSS THE THREE ORPHAN DEFINITIONS  "
          f"({args.model}/{args.transform})")
    print("=" * 78)
    print(f"  PARAMS k={args.k} sim_floor={args.sim_floor} "
          f"power={args.power}  (must match the frozen values)")

    # The nested definitions, then the two complements. 52 is a SUBSET
    # of 70, so if Fmax rises when 18 proteins are removed, those 18
    # must have been harder than the average. Scoring them directly
    # turns an anomaly into an explanation; leaving them out would make
    # the rise look like noise we were choosing not to investigate.
    nested: list[tuple[str, str, list[str]]] = [
        ("ORPHAN_nohit", "no BLAST-class hit", list(orphans)),
    ]
    if jack_ids:
        nested.append(("HARDENED", "no hit, either sequence search",
                       jack_ids))
    nested.append(("STRUCT_HARDENED",
                   "no hit from sequence or structure", struct_ids))

    complements: list[tuple[str, str, list[str]]] = []
    if jack_ids:
        rescued_profile = [i for i in orphans if i not in set(jack_ids)]
        if len(rescued_profile) >= 10:
            complements.append(("RESCUED_PROFILE",
                                "  ...rescued by profile search only",
                                rescued_profile))
        rescued_struct = [i for i in jack_ids if i not in set(struct_ids)]
        if len(rescued_struct) >= 10:
            complements.append(("RESCUED_STRUCTURE",
                                "  ...rescued by structure only",
                                rescued_struct))
    panels = nested + complements

    out: dict = {"model": args.model, "transform": args.transform,
                 "k": args.k, "sim_floor": args.sim_floor,
                 "power": args.power, "replicates": args.replicates,
                 "seed": args.seed, "panels": {}}

    print(f"\n  {'definition':<36}{'n':>5}{'Fmax':>9}{'95% CI':>20}")
    rows = []
    for label, pretty, ids in panels:
        if len(ids) < 10:
            print(f"  {pretty:<36}{len(ids):>5}   too few to bootstrap")
            continue
        bench.register_subset(label, "T0_full", ids)
        order, truth, prop = boot.predictions(
            bench, args.model, args.transform, label, args.k,
            args.sim_floor, args.power, args.block, not args.no_cache)
        if not order:
            print(f"  {pretty:<36}{len(ids):>5}   no predictions")
            continue
        A = boot.threshold_counts(order, truth, prop)
        r = boot.single_ci(A, np.arange(len(order)), args.replicates,
                           args.seed, args.alpha)
        lo, hi = r.get("ci", [None, None])
        print(f"  {pretty:<36}{len(order):>5}{r['fmax']:>9.3f}"
              f"   [{lo:.3f}, {hi:.3f}]")
        out["panels"][label] = {"pretty": pretty, "n": len(order),
                                "fmax": r["fmax"], "ci": [lo, hi]}
        rows.append((label, pretty, len(order), r["fmax"], lo, hi))
        if label.startswith("RESCUED"):
            rows.pop()
            out["panels"][label]["complement"] = True

    # ---- the direction check -------------------------------------------
    # An explanation for a rise, if the complements supply one.
    rp = out["panels"].get("RESCUED_STRUCTURE")
    sh = out["panels"].get("STRUCT_HARDENED")
    hd = out["panels"].get("HARDENED")
    if rp and sh and hd and rp["fmax"] < hd["fmax"]:
        out["rise_explained"] = (
            f"the {rp['n']} proteins removed by structural hardening "
            f"score {rp['fmax']:.3f}, below the {hd['n']}-protein set's "
            f"{hd['fmax']:.3f}, so removing them raises the mean. A "
            f"remote structural relative does not make a protein easier "
            f"to annotate by embedding transfer; fold-level similarity "
            f"is often functionally uninformative.")
        print(f"\n  EXPLANATION  {out['rise_explained']}")

    note = None
    if len(rows) >= 2:
        prev = None
        for label, pretty, n, f, lo, hi in rows:
            if prev is not None and f > prev[3] + 1e-9:
                note = (f"Fmax ROSE from {prev[3]:.3f} ({prev[1]}, "
                        f"n={prev[2]}) to {f:.3f} ({pretty}, n={n}). "
                        f"Removing proteins that have detectable "
                        f"relatives should leave a HARDER set, so a rise "
                        f"needs explaining before it is reported.")
            prev = (label, pretty, n, f, lo, hi)
        first, last = rows[0], rows[-1]
        out["attrition"] = {"from_n": first[2], "to_n": last[2],
                            "from_fmax": first[3], "to_fmax": last[3],
                            "delta": last[3] - first[3]}
        print(f"\n  attrition {first[2]} -> {last[2]} proteins, "
              f"Fmax {first[3]:.3f} -> {last[3]:.3f} "
              f"({last[3] - first[3]:+.3f})")
        ov = (last[4] <= first[5] and first[4] <= last[5])
        print(f"  intervals {'overlap' if ov else 'DO NOT overlap'}"
              f" — the hardened figure is "
              f"{'consistent with' if ov else 'distinguishable from'}"
              f" the original")
        out["intervals_overlap"] = bool(ov)
    if note:
        print(f"\n  WARNING  {note}")
        out["warning"] = note

    out["note"] = ("Scored through 36_bootstrap_baseline.py's own "
                   "register_subset and predictions path, so these are "
                   "comparable with every other figure in the paper.")
    p = cfg.p("results") / "fmax_structural.json"
    p.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  wrote {p}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
