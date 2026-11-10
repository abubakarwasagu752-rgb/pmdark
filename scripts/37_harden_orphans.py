#!/usr/bin/env python3
"""Test whether the orphan stratum is actually orphan, with jackhmmer.

    python scripts/37_harden_orphans.py
    python scripts/37_harden_orphans.py --threads 16
    python scripts/37_harden_orphans.py --harden-at 1e-3   # looser definition

WHY THIS IS THE CRITICAL PATH

  Every claim about the P. monodon dark proteome rests on 94 benchmark
  proteins, and "orphan" currently means one thing only:

      no DIAMOND hit at e-value 10 under --ultra-sensitive

  DIAMOND is fast but it is not the most sensitive thing available. If an
  iterative profile search finds relatives for a meaningful share of those
  94, then homology transfer WOULD have worked on them, the stratum is
  contaminated, and Fmax 0.664 is not the orphan number. That is the
  cheapest attack a reviewer can make on the paper, so it is worth making
  it ourselves first.

WHAT IS AND IS NOT BEING TESTED

  The search is against the BENCHMARK REFERENCE SET, not UniProt at large,
  and that is deliberate. The question is not "is this protein known
  anywhere" but "could homology transfer have annotated it from the
  reference set we transfer from". A relative with no GO annotation, or
  one outside the reference set, cannot help homology transfer and so
  cannot contaminate the stratum.

  jackhmmer is run with the default inclusion threshold (--incE 1e-3) so
  the profile does not drift on noise, but reports down to a very
  permissive E-value so weak hits are visible. Orphan status is then
  tabulated across a range of cuts rather than at one, because "orphan" is
  a threshold-dependent claim and the paper should state which threshold
  it used.

  The default --harden-at is 10.0, matching DIAMOND's permissive e-value.
  That is the most hostile choice to our own claim: any query with even a
  marginal profile hit is thrown out. If the number survives that, it
  survives anything.

WHAT IT DOES NOT COVER

  Structural search. Foldseek would need predicted structures for the
  benchmark queries, which do not exist yet (ESMFold is still deferred).
  A profile search is the strongest sequence-level test available without
  new compute, and it is what HHblits-class criticism usually reduces to.
  State the gap in the limitations rather than implying it was checked.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "src"))

import numpy as np                                                 # noqa: E402

from pmdark.shards import read_fasta, write_fasta                  # noqa: E402

CUTS = [1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 1.0, 10.0]


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_mod(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot import {filename}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def find_jackhmmer(explicit: str | None) -> str | None:
    """jackhmmer may be in any of the project envs, or nowhere."""
    if explicit:
        return explicit if Path(explicit).exists() else None
    env = os.environ.get("JACKHMMER")
    if env and Path(env).exists():
        return env
    found = shutil.which("jackhmmer")
    if found:
        return found
    home = Path.home()
    for base in (home / "pmdark_work" / "envs", home / "miniconda3" / "envs"):
        if not base.is_dir():
            continue
        for envdir in sorted(base.iterdir()):
            cand = envdir / "bin" / "jackhmmer"
            if cand.exists():
                return str(cand)
    return None


def run_jackhmmer(exe: str, query: Path, db: Path, tbl: Path,
                  iterations: int, evalue: float, threads: int,
                  force: bool) -> bool:
    if tbl.exists() and tbl.stat().st_size > 0 and not force:
        log(f"reusing {tbl.name}")
        return True
    cmd = [exe,
           "--cpu", str(threads),
           "-N", str(iterations),
           "-E", str(evalue),          # report weak hits
           "--domE", str(evalue),
           "--noali",
           "--tblout", str(tbl),
           str(query), str(db)]
    log("running " + " ".join(Path(c).name if "/" in c else c for c in cmd))
    nq = sum(1 for ln in query.open() if ln.startswith(">"))
    log(f"  ({nq} queries x {iterations} iterations against {db.name})")
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL)
    if r.returncode != 0:
        print(f"jackhmmer exited {r.returncode}")
        return False
    return True


ITER_SUFFIX = re.compile(r"-i\d+$")


def run_pyhmmer(query: Path, db: Path, tbl: Path, iterations: int,
                threads: int, force: bool) -> bool:
    """Fallback when the jackhmmer CLI is absent but pyhmmer is installed.

    pyhmmer bundles the same HMMER engine and writes the same tblout, so
    parse_tbl() is shared between both paths. HMMER's default reporting
    E-value is 10 and its default inclusion threshold is 1e-3, which are
    the values this script wants, so no thresholds are overridden here.
    """
    if tbl.exists() and tbl.stat().st_size > 0 and not force:
        log(f"reusing {tbl.name}")
        return True
    try:
        import pyhmmer
        from pyhmmer.easel import Alphabet, SequenceFile
    except ImportError:
        return False
    log(f"running pyhmmer jackhmmer ({iterations} iterations, "
        f"{threads} cpus)")
    abc = Alphabet.amino()
    try:
        with SequenceFile(str(query), digital=True, alphabet=abc) as fh:
            queries = fh.read_block()
        with SequenceFile(str(db), digital=True, alphabet=abc) as fh:
            targets = fh.read_block()
        first = True
        with tbl.open("wb") as out:
            for res in pyhmmer.hmmer.jackhmmer(
                    queries, targets, max_iterations=iterations,
                    cpus=threads):
                res.hits.write(out, format="targets", header=first)
                first = False
    except Exception as e:                       # API drift, bad input
        print(f"pyhmmer failed: {type(e).__name__}: {e}")
        return False
    return True


def parse_tbl(tbl: Path, drop_self: bool = True) -> dict[str, float]:
    """query -> best (smallest) full-sequence E-value over all its hits.

    HMMER tblout columns, verified against real output rather than from
    memory: target(0), target_acc(1), query(2), query_acc(3),
    full-sequence E-value(4).

    THE TRAP: jackhmmer renames the query per iteration, so the query
    column holds "P12345-i1", "P12345-i2", ... and never the bare
    accession. Looking up bare accessions against this column matches
    nothing, every query looks hit-free, and the script reports that the
    whole orphan stratum survived -- a check that always passes and
    appears to confirm the claim. The suffix must be stripped.

    Later iterations re-report the same targets, so the minimum E-value
    across the file is the best evidence found at any iteration.
    """
    best: dict[str, float] = {}
    seen_suffix = False
    with tbl.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 5:
                continue
            target, raw = f[0], f[2]
            query = ITER_SUFFIX.sub("", raw)
            if query != raw:
                seen_suffix = True
            if drop_self and target == query:
                continue
            try:
                ev = float(f[4])
            except ValueError:
                continue
            if query not in best or ev < best[query]:
                best[query] = ev
    if best and not seen_suffix:
        log("  NOTE: no -iN iteration suffixes seen in the query column; "
            "jackhmmer output format may have changed — verify the parse")
    return best


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="esm2")
    ap.add_argument("--transform", default="whiten",
                    choices=["none", "center", "whiten"])
    ap.add_argument("--iterations", type=int, default=5)
    ap.add_argument("--report-evalue", type=float, default=10.0,
                    help="jackhmmer reporting threshold; permissive on "
                         "purpose so weak hits are visible")
    ap.add_argument("--harden-at", type=float, default=10.0,
                    help="a query with a jackhmmer hit at or below this "
                         "E-value is dropped from the orphan stratum. The "
                         "default is deliberately hostile to our claim.")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--replicates", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--sim-floor", type=float, default=0.0)
    ap.add_argument("--power", type=float, default=3.0)
    ap.add_argument("--block", type=int, default=512)
    ap.add_argument("--min-stratum", type=int, default=20)
    ap.add_argument("--jackhmmer", default=None)
    ap.add_argument("--force", action="store_true",
                    help="re-run jackhmmer even if the table exists")
    ap.add_argument("--no-cache", action="store_true")
    args = ap.parse_args()

    boot = load_mod("pmdark_boot", "36_bootstrap_baseline.py")
    cfg = boot.cfgmod.load()
    calib = boot.load_calib()
    bench = boot.Bench(cfg, calib)

    # pin the query universe to the embedded set, as 36 does
    bench.embeddings(args.model)
    orphans = bench.orphans()
    log(f"orphan stratum: {len(orphans)} queries with no DIAMOND hit "
        f"(e-value {bench.ident_evalue})")

    exe = find_jackhmmer(args.jackhmmer)
    have_py = importlib.util.find_spec("pyhmmer") is not None
    if exe:
        log(f"jackhmmer: {exe}")
    elif have_py:
        log("jackhmmer CLI not found; using the bundled pyhmmer engine")
    else:
        print("\nNo HMMER available. Either will do:")
        print("    conda install -n pmbio -c bioconda hmmer")
        print("    pip install pyhmmer")
        print("  then re-run, or pass --jackhmmer /path/to/jackhmmer.")
        print("\nNothing was changed. The orphan stratum remains")
        print("DIAMOND-defined, which must be stated as a limitation")
        print("until this check runs.")
        return 2

    # ---- orphan subset FASTA -------------------------------------------
    qfa = bench.dir / "query.faa"
    records = read_fasta(qfa)
    want = set(orphans)
    sub = [r for r in records if r.pid in want]
    if len(sub) != len(orphans):
        missing = want - {r.pid for r in sub}
        raise SystemExit(
            f"{len(missing)} orphan accessions are not in {qfa.name}: "
            f"{sorted(missing)[:5]}")
    out_dir = bench.dir / "orphan_check"
    out_dir.mkdir(parents=True, exist_ok=True)
    ofa = out_dir / "orphans.faa"
    write_fasta(sub, ofa)
    log(f"wrote {ofa.name} ({len(sub)} sequences)")

    tbl = out_dir / f"jackhmmer_N{args.iterations}.tbl"
    ref = bench.dir / "reference.faa"
    if exe:
        ok = run_jackhmmer(exe, ofa, ref, tbl, args.iterations,
                           args.report_evalue, args.threads, args.force)
    else:
        ok = run_pyhmmer(ofa, ref, tbl, args.iterations, args.threads,
                         args.force)
    if not ok:
        return 1

    best = parse_tbl(tbl)
    log(f"jackhmmer found at least one hit for {len(best)} of {len(orphans)}")

    # ---- survival across thresholds ------------------------------------
    print("\n" + "=" * 78)
    print("ORPHAN SURVIVAL vs jackhmmer SENSITIVITY")
    print("=" * 78)
    print(f"  {'E-value cut':>12}  {'with hit':>9}  {'still orphan':>13}  "
          f"{'% orphan':>9}")
    survival = {}
    n = len(orphans)
    for cut in CUTS:
        hit = sum(1 for a in orphans if best.get(a, float("inf")) <= cut)
        survival[f"{cut:g}"] = {"with_hit": hit, "still_orphan": n - hit}
        print(f"  {cut:>12g}  {hit:>9}  {n - hit:>13}  "
              f"{100 * (n - hit) / n:>8.1f}%")
    print("\n  A hit here means homology transfer could have reached that")
    print("  protein from the reference set, so DIAMOND's miss was a")
    print("  sensitivity failure rather than an absence of relatives.")

    # ---- the hardened stratum ------------------------------------------
    survivors = [a for a in orphans
                 if best.get(a, float("inf")) > args.harden_at]
    dropped = [a for a in orphans if a not in set(survivors)]
    print("\n" + "=" * 78)
    print(f"HARDENED STRATUM  (dropping any hit at E <= {args.harden_at:g})")
    print("=" * 78)
    print(f"  started with {n}, dropped {len(dropped)}, "
          f"{len(survivors)} survive")

    report: dict = {
        "model": args.model, "transform": args.transform,
        "jackhmmer": {"exe": exe, "iterations": args.iterations,
                      "report_evalue": args.report_evalue},
        "diamond_evalue": bench.ident_evalue,
        "n_orphan_diamond": n,
        "survival": survival,
        "harden_at": args.harden_at,
        "n_survivors": len(survivors), "n_dropped": len(dropped),
        "dropped": sorted(dropped),
        "strata": {},
    }

    if len(survivors) < args.min_stratum:
        print(f"\n  Only {len(survivors)} survive, below --min-stratum "
              f"{args.min_stratum}.")
        print("  The orphan claim cannot be supported at this threshold.")
        print("  Report the survival table and say so plainly; do not")
        print("  quote an Fmax from a stratum this small.")
    else:
        # Score the original, hardened and dropped strata the same way, so
        # the only thing that differs between them is which queries are in.
        tag = hashlib.md5(
            ("|".join(sorted(survivors))).encode()).hexdigest()[:8]
        names = {
            "diamond_orphans": ("T0_full", orphans),
            f"hardened_{tag}": ("T0_full", survivors),
        }
        if len(dropped) >= args.min_stratum:
            names[f"dropped_{tag}"] = ("T0_full", dropped)

        print()
        for label, (mask, ids) in names.items():
            bench.register_subset(label, mask, ids)
            order, truth, prop = boot.predictions(
                bench, args.model, args.transform, label, args.k,
                args.sim_floor, args.power, args.block, not args.no_cache)
            if not order:
                continue
            A = boot.threshold_counts(order, truth, prop)
            res = boot.single_ci(A, np.arange(len(order)), args.replicates,
                                 args.seed, args.alpha)
            # The per-ontology panel and the Fmax-maximising threshold are
            # what scripts/41 needs to attach honest confidence to each
            # prediction. Without them it falls back on the un-hardened
            # 94-protein calibration, which overstates these proteins.
            from pmdark.goeval import evaluate as _evaluate
            panel = _evaluate(truth, prop, bench.graph)
            ov = panel["overall"]
            res["threshold"] = float(ov["threshold"])
            res["precision"] = float(ov["precision"])
            res["recall"] = float(ov["recall"])
            for _ns in ("BP", "MF", "CC"):
                _v = panel.get(_ns, {}).get("fmax")
                res[_ns] = float(_v) if _v is not None else None
            short = ("DIAMOND-orphan" if label.startswith("diamond")
                     else "HARDENED" if label.startswith("hardened")
                     else "dropped (had a profile hit)")
            report["strata"][short] = res
            print(f"  {short:<30} n={res['n']:>4}  Fmax {res['fmax']:.3f}  "
                  f"95% CI [{res['ci'][0]:.3f}, {res['ci'][1]:.3f}]")

        a = report["strata"].get("DIAMOND-orphan")
        b = report["strata"].get("HARDENED")
        if a and b:
            shift = b["fmax"] - a["fmax"]
            print(f"\n  Hardening moves Fmax by {shift:+.3f} "
                  f"({a['fmax']:.3f} -> {b['fmax']:.3f}).")
            # The two strata are nested, so this is not a paired test and
            # no p-value is reported for it. What matters is whether the
            # hardened point estimate sits inside the original interval.
            inside = a["ci"][0] <= b["fmax"] <= a["ci"][1]
            if inside:
                print("  The hardened estimate falls inside the original")
                print("  interval, so the orphan result is not an artifact")
                print("  of DIAMOND's sensitivity. Report the HARDENED")
                print("  number and the survival table together.")
            else:
                print("  The hardened estimate falls OUTSIDE the original")
                print("  interval. The orphan number was driven by queries")
                print("  that a profile search can reach. Report the")
                print("  HARDENED number only, and say why it changed.")
            print("\n  These strata are nested, not independent, so no")
            print("  paired p-value is given for the shift; the comparison")
            print("  is of a point estimate against an interval.")

    out = cfg.p("results") / "orphan_hardening.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    log(f"wrote {out}")

    print("\n  Still unchecked: structural homology. Foldseek needs")
    print("  predicted structures for these queries, which do not exist")
    print("  yet. State that as a limitation rather than implying a")
    print("  structure-level search was done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
