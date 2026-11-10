#!/usr/bin/env python3
"""Compare protein sets the only fair way: identical metrics, both sets.

    python scripts/31_compare_protein_sets.py --theirs upstream_proteins.faa
    python scripts/31_compare_protein_sets.py            # ours alone

WHY THE PROTEIN-COUNT COMPARISON HAD TO BE ABANDONED

  I proposed "our 18,929 proteins at >= 100 aa against their 28,423" as the
  paper's opening claim. It fails three ways, each fatal on its own:

    1. Their 28,423 sits within 0.6% of the 28,261 protein-coding genes
       NCBI annotates in this genome. A ratio against a number that
       coincides with the reference gene count cannot distinguish "they
       inflated" from "we are incomplete".

    2. Neither count is gene-level. Our 18,929 proteins collapse onto
       9,507 genes — 1.99 per gene, with 49.2% of detected genes carrying
       more than one ORF. Their 586,928 unigenes collapsed to 28,423
       proteins by a different route. It compares two differently
       inflated numbers.

    3. Our set is incomplete by construction: 36.7% of annotated genes
       detected, ~26% below 300 aa. A count drawn from an incomplete set
       measures sampling depth, not assembly quality.

  A count is the wrong instrument. What distinguishes a fragmented
  assembly from a compact one is REDUNDANCY and COMPLETENESS per gene, and
  both are measurable on either set without privileging ours.

THE METRIC: BEST RECIPROCAL HITS AGAINST THE REFERENCE PROTEOME

  Reference protein P and query protein Q are paired only if P's best hit
  is Q and Q's best hit is P. Reciprocity is what makes this paralog-safe
  without needing genome coordinates — which matters, because their
  assembly is de novo and has none. The single-direction search is exactly
  what made 28_recovery_by_family.py credit one assembled actin to 99
  annotated actin genes.

  Then, for each protein set, with the same thresholds and the same code:

    proteins                    how many it contains
    reference genes paired      how many distinct genes it actually reaches
    proteins per paired gene    REDUNDANCY — the fragmentation measure
    median query coverage       COMPLETENESS of the proteins it does have
    >= 80% coverage             how many it gets nearly full-length
    no hit to the reference     its own dark-set candidate pool

  A fragmented assembly shows high redundancy and low coverage. A compact
  one shows redundancy near 1.0 and high coverage. That comparison is
  interpretable in a way a count never was, and it does not depend on
  either set being complete.

WHAT STANDS WITHOUT THIS SCRIPT AT ALL

  586,928 unigenes against 28,261 annotated protein-coding genes in the
  same species is 20.8x, at a transcript N50 of 417 bp. That is arithmetic
  on their published figures and the public reference annotation, with no
  input from our pipeline, and it is the claim that survives regardless of
  how our own recovery turns out. This script tests whether that
  over-assembly propagated into their protein set.

GETTING THEIR PROTEIN SET

  Their repository is reachable from the cluster login node but not from
  here. On the cluster:

      cd ~ && git clone --depth 1 \\
        https://github.com/LuiguiGallardo/p_monodon_transcriptome upstream
      find ~/upstream -size +1M \\
        \\( -name '*.faa' -o -name '*.pep' -o -name '*.fasta*' \\) -ls

  Then pass whichever file is the predicted proteome with --theirs. If the
  repo holds no protein FASTA, the script still reports ours, and the
  comparison can be made against the paper's supplementary tables instead.
  No contact with the authors is needed or attempted.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402


def _load(name: str):
    sp = importlib.util.spec_from_file_location(
        "m" + name[:2], Path(__file__).resolve().parent / name)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


_rx = _load("27_recovery_exact.py")
_qc = _load("25_protein_qc.py")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def best_hits(dm: str, query: Path, db: Path, out: Path, threads: int,
              evalue: float, force: bool) -> dict[str, tuple[str, float, float]]:
    """query id -> (best subject, identity, query coverage)."""
    if force or not out.exists():
        rc = subprocess.run([
            dm, "blastp", "--query", str(query), "--db", str(db),
            "--out", str(out), "--outfmt", "6", "qseqid", "sseqid",
            "pident", "qcovhsp", "bitscore",
            "--max-target-seqs", "5", "--evalue", str(evalue),
            "--threads", str(threads), "--quiet",
        ]).returncode
        if rc != 0:
            raise SystemExit(f"diamond blastp failed ({rc}) for {query.name}")
    best: dict[str, tuple[str, float, float, float]] = {}
    with out.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                ident, qcov, bits = float(f[2]), float(f[3]), float(f[4])
            except ValueError:
                continue
            cur = best.get(f[0])
            if cur is None or bits > cur[3]:
                best[f[0]] = (f[1], ident, qcov, bits)
    return {k: (v[0], v[1], v[2]) for k, v in best.items()}


def makedb(dm: str, fa: Path, db: Path, threads: int) -> None:
    if not Path(str(db) + ".dmnd").exists():
        log(f"diamond makedb {fa.name}")
        subprocess.run([dm, "makedb", "--in", str(fa), "-d", str(db),
                        "--threads", str(threads)], capture_output=True)


def profile(name: str, fa: Path, refp: Path, prot2gene: dict,
            work: Path, dm: str, args) -> dict:
    """Best-reciprocal-hit profile of one protein set against the reference."""
    log(f"--- {name}: {fa.name}")
    lens = _qc.read_fasta_lengths(fa)
    tag = name.lower().replace(" ", "_")
    qdb = work / f"db_{tag}"
    makedb(dm, fa, qdb, args.threads)

    fwd = best_hits(dm, fa, work / "db_ref", work / f"{tag}_to_ref.tsv",
                    args.threads, args.evalue, args.force_search)
    rev = best_hits(dm, refp, qdb, work / f"ref_to_{tag}.tsv",
                    args.threads, args.evalue, args.force_search)

    # Reciprocity: our Q pairs with reference P only if each is the other's
    # best hit. This is what the single-direction search got wrong.
    pairs: list[tuple[str, str, float, float]] = []
    for q, (s, ident, qcov) in fwd.items():
        back = rev.get(s)
        if back and back[0] == q:
            pairs.append((q, s, ident, qcov))
    passing = [(q, s, i, c) for (q, s, i, c) in pairs
               if i >= args.min_identity]

    genes = defaultdict(set)
    for q, s, i, c in passing:
        g = prot2gene.get(s)
        if g:
            genes[g].add(q)
    covs = [c for (_, _, _, c) in passing]
    # Coverage here is of the QUERY (the set's own protein). For the
    # reference-side view, use the reverse search's coverage.
    rev_cov = [v[2] for k, v in rev.items()
               if v[1] >= args.min_identity and prot2gene.get(k)]

    nhit = sum(1 for q in lens if q in fwd)
    res = {
        "name": name,
        "file": str(fa),
        "proteins": len(lens),
        "median_length": int(statistics.median(lens.values())) if lens else 0,
        "with_any_ref_hit": nhit,
        "no_ref_hit": len(lens) - nhit,
        "reciprocal_pairs": len(passing),
        "genes_paired": len(genes),
        "proteins_per_paired_gene": round(
            len(passing) / max(1, len(genes)), 2),
        "genes_with_multiple_proteins": sum(1 for v in genes.values()
                                            if len(v) > 1),
        "median_qcov_of_reference": round(
            statistics.median(rev_cov), 1) if rev_cov else None,
        "ref_ge80_qcov": sum(1 for c in rev_cov if c >= 80.0),
    }
    log(f"    {res['proteins']:,} proteins -> {res['genes_paired']:,} "
        f"reference genes, {res['proteins_per_paired_gene']} per gene")
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--theirs", type=Path, default=None,
                    help="the upstream predicted proteome (see the note above)")
    ap.add_argument("--ours", type=Path, default=None)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--min-identity", type=float, default=50.0,
                    help="reciprocal pairs below this are not called the "
                         "same gene product")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--force-search", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out, res_dir = cfg.p("proteins"), cfg.p("results")
    qc = json.loads((res_dir / "protein_qc.json").read_text())
    ours = args.ours or Path(qc["pep"])
    refp = next(iter(cfg.p("genome").glob("*_protein.faa")), None)
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if refp is None or ref_gff is None:
        raise SystemExit("reference proteome / GFF not staged — run 27 first")

    dm = _rx.find_diamond()
    if dm is None:
        raise SystemExit("diamond not found")

    log("mapping reference proteins to genes")
    A = _rx.annotation(ref_gff)
    prot2gene = A["prot2gene"]
    n_coding = len(A["coding"])
    log(f"  {len(prot2gene):,} proteins over {n_coding:,} coding genes")

    work = out / "compare"
    work.mkdir(parents=True, exist_ok=True)
    makedb(dm, refp, work / "db_ref", args.threads)

    sets = [profile("ours", ours, refp, prot2gene, work, dm, args)]
    if args.theirs:
        if not args.theirs.exists():
            raise SystemExit(f"{args.theirs} not found")
        sets.append(profile("upstream", args.theirs, refp, prot2gene,
                            work, dm, args))

    print("\n" + "=" * 80)
    print("PROTEIN SETS, IDENTICAL METRICS")
    print(f"reference: {refp.name}, {n_coding:,} coding genes")
    print("=" * 80)
    rows = [
        ("proteins in the set", "proteins", "{:,}"),
        ("median protein length (aa)", "median_length", "{:,}"),
        ("reference genes paired (BRH)", "genes_paired", "{:,}"),
        ("proteins per paired gene", "proteins_per_paired_gene", "{:.2f}"),
        ("genes with >1 protein", "genes_with_multiple_proteins", "{:,}"),
        ("median coverage of ref protein", "median_qcov_of_reference", "{}%"),
        ("ref proteins at >=80% coverage", "ref_ge80_qcov", "{:,}"),
        ("no hit to the reference at all", "no_ref_hit", "{:,}"),
    ]
    w = 34
    print(f"  {'metric':<{w}}" + "".join(f"{s['name']:>18}" for s in sets))
    print("-" * 80)
    for lab, key, fmt in rows:
        cells = ""
        for s in sets:
            v = s.get(key)
            cells += f"{(fmt.format(v) if v is not None else 'n/a'):>18}"
        print(f"  {lab:<{w}}{cells}")
    print("-" * 80)
    for s in sets:
        print(f"  {s['name']}: {100*s['genes_paired']/n_coding:.1f}% of coding "
              f"genes reached")

    print("\n  HOW TO READ THIS")
    print("  Redundancy near 1.0 with high coverage = a compact assembly.")
    print("  Redundancy well above 1.0 with low coverage = a fragmented one,")
    print("  where several partial proteins stand in for one real gene.")
    print("  Neither set is complete, and this comparison does not assume")
    print("  either is: both are measured against the same reference by the")
    print("  same code at the same thresholds.")
    if len(sets) == 1:
        print("\n  Only ours was profiled. Pass --theirs to compare; the")
        print("  docstring has the two commands that fetch their set.")

    rp = res_dir / "compare_protein_sets.json"
    rp.write_text(json.dumps({"reference_coding_genes": n_coding,
                              "min_identity": args.min_identity,
                              "sets": sets}, indent=2))
    print(f"\nwrote {rp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
