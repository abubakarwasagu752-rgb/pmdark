#!/usr/bin/env python3
"""Identify what the 32% unmapped reads actually are.

    python scripts/22_diagnose_unmapped.py --n 200000
    python scripts/22_diagnose_unmapped.py --runs SRR28572183,SRR28572173

THE PROBLEM

  Genome-guided alignment to the chromosome-level NSTDA_Pmon_1 assembly
  gives a mean overall rate of 68.2% (range 56.5-78.2%), against the >=85%
  gate in docs/ASSEMBLY_STRATEGY.md. Before merging transcripts we need to
  know whether the missing third is biology, contamination, or error,
  because each implies a different fix:

    WSSV reads            -> expected, quantifiable, and scientifically
                             interesting (per-sample viral load)
    bacterial / other     -> contamination; filter and move on
    rRNA / non-coding     -> library prep; harmless to the protein set
    shrimp-like protein   -> the REFERENCE GENOME is incomplete, and the
                             de novo arm becomes essential rather than
                             optional
    nothing at all        -> low-complexity, adapter residue, or junk

  The pattern already argues against WSSV being the main cause: the WORST
  library (56.5%) is dsRNA-only with no virus, and the BEST (78.2%) is
  infected. So this looks for all five.

METHOD

  Subsample unmapped reads, DIAMOND blastx them against Swiss-Prot, and
  summarise the taxonomic origin of the hits. blastx rather than blastn
  because we have a protein database to hand and translated search is far
  more sensitive across species boundaries.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def find_diamond() -> str | None:
    import os
    for cand in (os.environ.get("DIAMOND"), shutil.which("diamond")):
        if cand and Path(cand).exists():
            return cand
    home = Path.home()
    for env in ("busco", "pmbio"):
        p = home / "pmdark_work" / "envs" / env / "bin" / "diamond"
        if p.exists():
            return str(p)
    return None


def subsample(fq_gz: Path, n: int, out: Path) -> int:
    """Take the first n reads. Not random, but RNA-seq FASTQ order is
    effectively arbitrary with respect to origin, and a head is ~100x
    faster than a reservoir sample over a 300 MB gzip."""
    written = 0
    with gzip.open(fq_gz, "rt") as fh, out.open("w") as o:
        for i, line in enumerate(fh):
            if i >= n * 4:
                break
            if i % 4 == 0:
                o.write(">" + line[1:].split()[0] + "\n")
                written += 1
            elif i % 4 == 1:
                o.write(line)
    return written


def basic_stats(fa: Path) -> dict:
    lens, gc, ns = [], 0, 0
    total = 0
    with fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            s = line.strip().upper()
            lens.append(len(s))
            gc += s.count("G") + s.count("C")
            ns += s.count("N")
            total += len(s)
    lens.sort()
    return {
        "n_reads": len(lens),
        "mean_len": round(sum(lens) / max(1, len(lens)), 1),
        "gc_percent": round(100 * gc / max(1, total), 1),
        "n_percent": round(100 * ns / max(1, total), 2),
    }


def build_db(dmnd: str, sprot: Path, dbdir: Path, threads: int) -> Path | None:
    dbdir.mkdir(parents=True, exist_ok=True)
    db = dbdir / "sprot"
    if (dbdir / "sprot.dmnd").exists():
        log(f"  reusing {db}.dmnd")
        return db
    if not sprot.exists():
        log(f"  Swiss-Prot not found at {sprot}")
        return None
    log(f"  building DIAMOND db from {sprot.name} (a few minutes)")
    rc = subprocess.run([dmnd, "makedb", "--in", str(sprot),
                         "-d", str(db), "--threads", str(threads)]).returncode
    return db if rc == 0 else None


def classify(hits: Path, sprot_headers: dict[str, str]) -> dict:
    """Summarise blastx hits by the organism in the Swiss-Prot header.

    Swiss-Prot FASTA headers carry OS=<organism>, so species attribution
    needs no taxonomy download.
    """
    best: dict[str, tuple[float, str]] = {}
    with hits.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            q, s, pid, bits = p[0], p[1], float(p[2]), float(p[-1])
            if q not in best or bits > best[q][0]:
                best[q] = (bits, s)

    org = Counter()
    for _, (_, s) in best.items():
        acc = s.split("|")[1] if "|" in s else s
        org[sprot_headers.get(acc, "unknown")] += 1
    return {"n_with_hit": len(best), "organisms": org}


def load_headers(sprot: Path, limit: int = 0) -> dict[str, str]:
    """accession -> organism, parsed from OS= in Swiss-Prot headers."""
    out: dict[str, str] = {}
    opener = gzip.open if str(sprot).endswith(".gz") else open
    with opener(sprot, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) > 2 else line[1:].split()[0]
            os_i = line.find("OS=")
            if os_i < 0:
                continue
            rest = line[os_i + 3:]
            for tag in (" OX=", " GN=", " PE=", " SV="):
                j = rest.find(tag)
                if j >= 0:
                    rest = rest[:j]
            out[acc] = rest.strip()
            if limit and len(out) >= limit:
                break
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=None,
                    help="comma list; default = 3 spanning the rate range")
    ap.add_argument("--n", type=int, default=200000,
                    help="reads to subsample per library")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--evalue", type=float, default=1e-5)
    args = ap.parse_args()

    cfg = cfgmod.load()
    unmapped = cfg.p("align") / "_unmapped"
    if not unmapped.is_dir():
        raise SystemExit(f"no unmapped reads at {unmapped}")

    dmnd = find_diamond()
    if not dmnd:
        raise SystemExit("diamond not found; see find_diamond()")
    log(f"diamond: {dmnd}")

    # Default picks span the observed alignment-rate range so we can see
    # whether the unmapped content DIFFERS between good and bad libraries.
    runs = (args.runs.split(",") if args.runs else
            ["SRR28572183", "SRR28572178", "SRR28572173"])
    log(f"libraries: {', '.join(runs)}")

    sprot = cfg.p("raw", "uniprot") / "uniprot_sprot.fasta.gz"
    work = cfg.p("align") / "_diagnosis"
    work.mkdir(parents=True, exist_ok=True)

    db = build_db(dmnd, sprot, work / "db", args.threads)
    if db is None:
        return 1

    log("parsing Swiss-Prot headers for organism attribution")
    headers = load_headers(sprot)
    log(f"  {len(headers):,} accessions")

    report: dict = {"runs": {}}
    for run in runs:
        fq = unmapped / f"{run}_unmapped_1.fq.gz"
        if not fq.exists():
            log(f"  {run}: no unmapped file, skipping")
            continue

        log(f"--- {run} ---")
        fa = work / f"{run}.sample.fa"
        n = subsample(fq, args.n, fa)
        stats = basic_stats(fa)
        log(f"  sampled {n:,} reads  mean_len {stats['mean_len']}  "
            f"GC {stats['gc_percent']}%  N {stats['n_percent']}%")

        hits = work / f"{run}.sprot.tsv"
        if not hits.exists() or hits.stat().st_size == 0:
            log("  diamond blastx vs Swiss-Prot")
            subprocess.run([
                dmnd, "blastx", "-q", str(fa), "-d", str(db),
                "-o", str(hits), "--outfmt", "6", "qseqid", "sseqid",
                "pident", "evalue", "bitscore",
                "--evalue", str(args.evalue), "--max-target-seqs", "1",
                "--threads", str(args.threads), "--quiet",
            ])

        cls = classify(hits, headers)
        frac = cls["n_with_hit"] / max(1, n)
        log(f"  {cls['n_with_hit']:,}/{n:,} ({100*frac:.1f}%) hit a protein")
        log("  top organisms:")
        for name, c in cls["organisms"].most_common(12):
            log(f"    {c:>7,}  {100*c/max(1,cls['n_with_hit']):5.1f}%  {name}")
        log(f"  {100*(1-frac):.1f}% matched NO protein "
            "(rRNA, non-coding, low-complexity or adapter)")

        report["runs"][run] = {
            "stats": stats,
            "n_with_protein_hit": cls["n_with_hit"],
            "fraction_with_hit": round(frac, 4),
            "top_organisms": cls["organisms"].most_common(25),
        }

    out = cfg.p("results") / "unmapped_diagnosis.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    print("\n" + "=" * 74)
    print("HOW TO READ THIS")
    print("=" * 74)
    print("  mostly Penaeus/crustacean proteins -> the reference genome is")
    print("     incomplete; the de novo arm is REQUIRED, not optional, and")
    print("     those loci are prime dark-proteome candidates")
    print("  mostly bacteria                    -> contamination; filter")
    print("  WSSV proteins                      -> viral load; quantify it,")
    print("     it is a real covariate for the STAT analysis")
    print("  mostly no hit at all               -> rRNA/non-coding; benign")
    print("     for a protein-level project")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
