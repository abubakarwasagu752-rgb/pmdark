#!/usr/bin/env python3
"""Replace the DIAMOND baseline with an iterative-profile one (jackhmmer).

    python scripts/38_profile_baseline.py --threads 16
    python scripts/38_profile_baseline.py --transfer-evalue 0.01

WHY THE DIAMOND BASELINE IS NOT ENOUGH

  scripts/37_harden_orphans.py settled it: jackhmmer finds homologs in the
  reference set for 51 of the 94 queries DIAMOND called orphans (24 of them
  at the conventional E<1e-3). DIAMOND is therefore demonstrably the weaker
  search on exactly this data, and every accuracy number measured against
  it -- +0.017 at T0, +0.021 at I50, +0.051 at I30 -- is measured against a
  baseline we have already shown to be beatable by a better search.

  A reviewer who knows this will ask why the comparison was not run against
  a profile method. This runs it.

WHAT IS HELD FIXED

  Everything except the search. The tier masks, the top-k vote, the
  ancestor propagation, the Fmax code and the paired bootstrap are the
  same objects used for the DIAMOND arm, imported from
  scripts/36_bootstrap_baseline.py rather than reimplemented, so the only
  difference between the two baselines is which tool found the neighbours.

  Hits are used for transfer at a conventional significance threshold
  (--transfer-evalue, default 1e-3) rather than at jackhmmer's permissive
  reporting threshold. Transferring GO from E=10 hits would inject noise
  and flatter our own method; the baseline gets its conventional operating
  point.

A LIMITATION THIS SURFACES

  The identity tiers are defined by DIAMOND percent-identity, so they can
  only exclude relatives DIAMOND could detect. We now know it misses some.
  A reference that DIAMOND never aligned carries no recorded identity, is
  therefore in no exclusion list, and stays permitted at every I-tier even
  if it is a genuine close relative. This leaks equally into BOTH arms --
  the PLM arm may use such a reference too -- so the comparison stays fair,
  but "no homolog above X% identity" is a weaker statement for the I-tiers
  than it sounds. The ORPHAN stratum is the one that has been hardened
  against this, which is another reason it carries the paper's claim.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "src"))


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_mod(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot import {filename}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_pairs(tbl: Path, suffix_re, max_evalue: float,
                drop_self: bool = True) -> tuple[dict, int, int]:
    """query -> {subject: bitscore}, for hits at or below `max_evalue`.

    Column indices are the ones verified against real HMMER output in
    scripts/37: target(0), query(2), full-sequence E-value(4), score(5).
    The query column may carry jackhmmer's per-iteration "-iN" suffix
    depending on the build, so it is stripped the same way.

    Keeps the best (highest) score per pair, since later iterations
    re-report the same targets.
    """
    out: dict[str, dict[str, float]] = {}
    kept = dropped = 0
    with tbl.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 6:
                continue
            target = f[0]
            query = suffix_re.sub("", f[2])
            if drop_self and target == query:
                continue
            try:
                ev, score = float(f[4]), float(f[5])
            except ValueError:
                continue
            if ev > max_evalue:
                dropped += 1
                continue
            d = out.setdefault(query, {})
            if score > d.get(target, -1.0):
                d[target] = score
            kept += 1
    return out, kept, dropped


def sensitivity_table(cfg, new_report: dict, label: str,
                      evalue: float) -> None:
    """Print jackhmmer against DIAMOND, tier by tier.

    Only valid when both arms used the SAME transfer threshold. The first
    version of this table compared jackhmmer at E<=1e-3 against DIAMOND
    unfiltered (E<=10) and so reported jackhmmer as LESS sensitive at the
    permissive tiers, which is backwards. The thresholds are checked here
    and the table refuses to print if they disagree.
    """
    # The DIAMOND baseline was written as homology_baseline.json before the
    # files were namespaced by source, so accept either name rather than
    # silently skipping the comparison on an older results directory.
    old = None
    for name in ("homology_baseline_DIAMOND.json", "homology_baseline.json"):
        path = cfg.p("results") / name
        if not path.exists():
            continue
        try:
            old = json.loads(path.read_text())
            break
        except ValueError:
            continue
    if old is None:
        return

    old_e = old.get("transfer_evalue", "unfiltered")
    if old_e != evalue:
        print("\n" + "-" * 78)
        print("SEARCH SENSITIVITY  (skipped)")
        print("-" * 78)
        print(f"  The DIAMOND baseline on disk used transfer threshold "
              f"{old_e}, this run used {evalue:g}.")
        print("  Comparing counts across different thresholds would make the")
        print("  weaker search look more sensitive. Regenerate the DIAMOND")
        print("  arm at a matched threshold first:")
        print(f"    python scripts/36_bootstrap_baseline.py --part baseline "
              f"--transfer-evalue {evalue:g}")
        return

    print("\n" + "-" * 78)
    print(f"SEARCH SENSITIVITY  (queries with any permitted hit, "
          f"both at E<={evalue:g})")
    print("-" * 78)
    print(f"  {'tier':<14} {'DIAMOND':>18} {'  ' + label:>18}   {'gain':>8}")
    for tier, row in new_report.get("tiers", {}).items():
        n = row.get("best", {})
        o = old.get("tiers", {}).get(tier, {}).get("best", {})
        if "answerable" not in n or "answerable" not in o:
            continue
        tot = n.get("n_scored") or 1
        gain = n["answerable"] - o["answerable"]
        print(f"  {tier:<14} {o['answerable']:>7}/{tot:<5} "
              f"({100 * o['answerable'] / tot:>4.1f}%) "
              f"{n['answerable']:>7}/{tot:<5} "
              f"({100 * n['answerable'] / tot:>4.1f}%)   {gain:>+8}")
    print("\n  A positive gain is homology that DIAMOND missed. Where the")
    print("  gain is large, the DIAMOND baseline was understating what")
    print("  homology transfer can do, and the profile row is the one to")
    print("  report.")


def main() -> int:
    boot = load_mod("pmdark_boot", "36_bootstrap_baseline.py")
    hard = load_mod("pmdark_hard", "37_harden_orphans.py")

    ap = boot.build_parser(__doc__)
    ap.add_argument("--iterations", type=int, default=5)
    ap.add_argument("--report-evalue", type=float, default=10.0,
                    help="jackhmmer reporting threshold; permissive so the "
                         "transfer threshold can be varied without "
                         "re-searching")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--jackhmmer", default=None)
    ap.add_argument("--force-search", action="store_true")
    args = ap.parse_args()

    cfg = boot.cfgmod.load()
    calib = boot.load_calib()
    bench = boot.Bench(cfg, calib)

    # Pin the query universe to the embedded set, as the other parts do,
    # so both arms are scored over identical proteins.
    first = args.cmp_models.split(",")[0].strip()
    bench.embeddings(first)
    q_ids, _ = bench.ids()
    log(f"query universe: {len(q_ids)} proteins")

    exe = hard.find_jackhmmer(args.jackhmmer)
    have_py = importlib.util.find_spec("pyhmmer") is not None
    if exe:
        log(f"jackhmmer: {exe}")
    elif have_py:
        log("jackhmmer CLI not found; using the bundled pyhmmer engine")
    else:
        print("\nNo HMMER available. Either will do:")
        print("    conda install -n pmbio -c bioconda hmmer")
        print("    pip install pyhmmer")
        return 2

    out_dir = bench.dir / "profile_search"
    out_dir.mkdir(parents=True, exist_ok=True)
    tbl = out_dir / f"jackhmmer_all_N{args.iterations}.tbl"
    qfa, rfa = bench.dir / "query.faa", bench.dir / "reference.faa"

    if not (tbl.exists() and tbl.stat().st_size > 0) or args.force_search:
        log(f"searching all {len(q_ids)} queries against the reference set")
        log("  this is ~11x the orphan-only search; expect a few minutes")
    if exe:
        ok = hard.run_jackhmmer(exe, qfa, rfa, tbl, args.iterations,
                                args.report_evalue, args.threads,
                                args.force_search)
    else:
        ok = hard.run_pyhmmer(qfa, rfa, tbl, args.iterations, args.threads,
                              args.force_search)
    if not ok:
        return 1

    hits, kept, dropped = parse_pairs(tbl, hard.ITER_SUFFIX,
                                      args.transfer_evalue)
    log(f"transfer set: {kept:,} alignments at E<={args.transfer_evalue:g} "
        f"for {len(hits):,} queries ({dropped:,} weaker hits excluded)")
    if not hits:
        print("\nNo hits survived the transfer threshold. Nothing to compare.")
        return 1

    print("\n" + "=" * 78)
    print(f"PROFILE BASELINE   jackhmmer -N {args.iterations}, "
          f"transfer at E<={args.transfer_evalue:g}")
    print("=" * 78)

    rep = boot.run_baseline(bench, args, hits=hits, label="jackhmmer")
    sensitivity_table(cfg, rep, "jackhmmer", args.transfer_evalue)
    boot.run_compare(bench, args, hits=hits, label="jackhmmer")

    leak = []
    for tier, row in rep.get("tiers", {}).items():
        b = row.get("best", {})
        if tier in ("I20", "Inone") and b.get("answerable_frac", 0) > 0.25:
            leak.append((tier, b))
    if leak:
        print("\n" + "-" * 78)
        print("TIER LEAKAGE — THESE TIERS ARE NOT HOMOLOG-FREE")
        print("-" * 78)
        for tier, b in leak:
            print(f"  {tier:<8} profile search still answers "
                  f"{b['answerable']}/{b['n_scored']} "
                  f"({100 * b['answerable_frac']:.1f}%)")
        print()
        print("  The identity tiers are built from DIAMOND percent-identity,")
        print("  so they can only exclude relatives DIAMOND DETECTED. A")
        print("  reference it never aligned has no recorded identity, sits")
        print("  in no exclusion list, and stays permitted even at Inone.")
        print("  So Inone removes DIAMOND-detectable homologs, not all")
        print("  homologs, and must not be described as an orphan")
        print("  condition. It leaks into both arms equally, so the")
        print("  comparison stays fair.")
        print()
        print("  The hardened ORPHAN stratum from scripts/37 is the only")
        print("  condition verified against profile search, which is why it")
        print("  should carry the paper's central claim.")

    print("\n" + "-" * 78)
    print("WHICH BASELINE TO REPORT")
    print("-" * 78)
    print("  Report the profile baseline as the primary comparison. It is")
    print("  the stronger search, and we have already demonstrated on the")
    print("  orphan stratum that DIAMOND misses homologs it finds, so a")
    print("  DIAMOND-only comparison overstates the method's advantage.")
    print("  Keep the DIAMOND row as the BLAST-class reference point that")
    print("  most readers will have in mind, clearly labelled as such.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
