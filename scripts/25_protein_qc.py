#!/usr/bin/env python3
"""Partition the predicted proteins against the reference, and test completeness.

    python scripts/25_protein_qc.py --threads 16
    python scripts/25_protein_qc.py --no-busco          # partition only, ~2 min

WHY THIS RUNS BEFORE THE PROTEINNESS FILTER

  ORF calling produced 18,929 proteins at >= 100 aa against the upstream
  study's 28,423 — a ratio of 0.67. I had pre-stated that this comparison
  would quantify how much of the upstream protein set was assembly
  fragments. It does not, and the reason is a coincidence worth stating
  plainly: the NCBI annotation of this same genome has 28,261
  protein-coding genes, so the upstream count of 28,423 sits almost exactly
  ON the reference gene count. Our 0.67 is therefore equally well read as

    (A) the upstream set was inflated by fragmentation, or
    (B) our set is missing a third of the proteome.

  A protein count alone cannot separate those. Two measurements can, and
  both are cheap:

    1. RECOVERY. How many DISTINCT annotated protein-coding genes do our
       proteins hit? Hemocytes are one tissue. A single tissue expressing
       60-75% of an animal's genes is unremarkable, so if recovery is in
       that range and the misses are tissue-specific genes, reading (B)
       dissolves.

    2. COMPLETENESS. BUSCO on the protein set. Conserved single-copy
       orthologues are expressed in essentially every tissue, so they are
       the part of the proteome we have no excuse for missing.

  PRE-STATED DECISION RULE, fixed before looking at the output:

    BUSCO complete >= 90%  -> the compact protein set is complete; the
                              0.67 ratio reflects upstream inflation plus
                              single-tissue sampling, and Table 1 rests on
                              the assembly metrics (586,928 unigenes and a
                              417 bp N50 versus 28,612 loci), NOT on the
                              protein count.
    BUSCO complete 80-90%  -> report the shortfall; the protein-count
                              comparison is not usable as evidence.
    BUSCO complete < 80%   -> the assembly or ORF calling is losing real
                              genes. Stop and fix it before the dark
                              proteome is defined on top of it.

WHAT ELSE THIS SETTLES

  The dark-set denominator. "Dark" is only meaningful relative to what is
  already annotated, so each protein is placed against the reference:
  sense overlap with a protein-coding gene (known), sense overlap with a
  non-coding gene, antisense-only overlap, or no overlap at all (a novel
  locus). The 70-99 aa band is reported the same way, because a short
  protein at a locus NCBI never annotated is the AMP candidate Aim 3 is
  built around.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

BIN = 100_000  # genomic index bin size


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------- sequences
def read_fasta_lengths(fa: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    cur, n = None, 0
    with fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    out[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip().replace("*", ""))
    if cur:
        out[cur] = n
    return out


# --------------------------------------------------------------- annotation
def _attr(s: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", s)
    return m.group(1) if m else None


def reference_genes(gff: Path) -> tuple[dict, dict]:
    """Binned gene index, and the id->biotype map.

    index[chrom][bin] -> list of (start, end, strand, gene_id, coding)
    """
    idx: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    biotype: dict[str, str] = {}
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            s, e = int(p[3]), int(p[4])
            gid = _attr(p[8], "ID") or _attr(p[8], "gene") or f"{p[0]}:{s}"
            bt = _attr(p[8], "gene_biotype") or "unknown"
            biotype[gid] = bt
            coding = bt == "protein_coding"
            rec = (s, e, p[6], gid, coding)
            for b in range(s // BIN, e // BIN + 1):
                idx[p[0]][b].append(rec)
    return idx, biotype


def predicted_spans(gff3: Path) -> dict[str, tuple]:
    """ORF id -> (chrom, start, end, strand) from TransDecoder's genome GFF3."""
    span: dict[str, list] = {}
    with gff3.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("mRNA", "CDS"):
                continue
            ident = _attr(p[8], "ID") if p[2] == "mRNA" else _attr(p[8], "Parent")
            if not ident:
                continue
            # IDs appear as 'MSTRG.1.1.p1' or 'GENE.MSTRG.1.1~~MSTRG.1.1.p1'
            if "~~" in ident:
                ident = ident.split("~~")[-1]
            ident = ident.split(",")[0]
            s, e = int(p[3]), int(p[4])
            if ident in span:
                r = span[ident]
                r[1] = min(r[1], s)
                r[2] = max(r[2], e)
            else:
                span[ident] = [p[0], s, e, p[6]]
    return {k: tuple(v) for k, v in span.items()}


def classify(idx: dict, chrom: str, start: int, end: int, strand: str):
    """-> (category, hit_gene_ids_sense_coding)"""
    sense_pc, sense_nc, anti = [], [], []
    bins = idx.get(chrom)
    if bins:
        seen = set()
        for b in range(start // BIN, end // BIN + 1):
            for (s, e, st, gid, coding) in bins.get(b, ()):
                if gid in seen or e < start or s > end:
                    continue
                seen.add(gid)
                if st == strand:
                    (sense_pc if coding else sense_nc).append(gid)
                else:
                    anti.append(gid)
    if sense_pc:
        return "known_pc", sense_pc
    if sense_nc:
        return "known_nc", []
    if anti:
        return "antisense", []
    return "novel", []


# -------------------------------------------------------------------- BUSCO
def find_busco() -> str | None:
    import os
    for cand in (os.environ.get("BUSCO"), shutil.which("busco")):
        if cand and Path(cand).exists():
            return cand
    p = Path.home() / "pmdark_work" / "envs" / "busco" / "bin" / "busco"
    return str(p) if p.exists() else None


def parse_busco(out_path: Path, tag: str) -> dict | None:
    for js in sorted(out_path.glob(f"**/short_summary*{tag}*.json")):
        try:
            d = json.loads(js.read_text())
            r = d.get("results", d)
            return {k: r.get(k) for k in
                    ("Complete percentage", "Single copy percentage",
                     "Multi copy percentage", "Fragmented percentage",
                     "Missing percentage", "n_markers")}
        except (OSError, ValueError):
            continue
    for txt in sorted(out_path.glob("**/short_summary*.txt")):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],"
                      r"F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", txt.read_text())
        if m:
            return {"Complete percentage": float(m.group(1)),
                    "Single copy percentage": float(m.group(2)),
                    "Multi copy percentage": float(m.group(3)),
                    "Fragmented percentage": float(m.group(4)),
                    "Missing percentage": float(m.group(5)),
                    "n_markers": int(m.group(6))}
    return None


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pep", type=Path, default=None)
    ap.add_argument("--genome-gff3", type=Path, default=None)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--no-busco", action="store_true")
    ap.add_argument("--busco-lineage", default="arthropoda_odb10")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    out = cfg.p("proteins")
    res = cfg.p("results")

    orf = res / "orf_calling.json"
    meta = json.loads(orf.read_text()) if orf.exists() else {}
    pep = args.pep or (Path(meta["pep_file"]) if meta.get("pep_file") else None)
    gff3 = args.genome_gff3 or (Path(meta["genome_gff3"])
                                if meta.get("genome_gff3") else None)
    if not pep or not pep.exists():
        raise SystemExit("protein FASTA not found — run scripts/24_call_orfs.py")
    if not gff3 or not gff3.exists():
        raise SystemExit("genome GFF3 not found — rerun 24_call_orfs.py; the "
                         "partition needs genome coordinates")

    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if ref_gff is None:
        raise SystemExit("reference GFF not found in the work tier")

    log(f"proteins : {pep.name}")
    log(f"coords   : {gff3.name}")
    log(f"reference: {ref_gff.name}")

    lens = read_fasta_lengths(pep)
    log(f"  {len(lens):,} proteins")

    log("indexing reference genes")
    idx, biotype = reference_genes(ref_gff)
    n_pc = sum(1 for b in biotype.values() if b == "protein_coding")
    log(f"  {len(biotype):,} genes, {n_pc:,} protein-coding")

    log("reading predicted ORF genome coordinates")
    spans = predicted_spans(gff3)
    log(f"  {len(spans):,} ORFs placed")
    missing = [k for k in lens if k not in spans]
    if missing:
        log(f"  WARNING {len(missing):,} proteins have no genomic span "
            f"(e.g. {missing[:3]}) — counted as unplaced")

    log("classifying against the annotation")
    cats: dict[str, int] = defaultdict(int)
    band: dict[str, int] = defaultdict(int)
    hit_pc: set[str] = set()
    for pid, L in lens.items():
        sp = spans.get(pid)
        if sp is None:
            cats["unplaced"] += 1
            if 70 <= L < 100:
                band["unplaced"] += 1
            continue
        chrom, s, e, strand = sp
        cat, genes = classify(idx, chrom, s, e, strand)
        cats[cat] += 1
        hit_pc.update(genes)
        if 70 <= L < 100:
            band[cat] += 1

    order = ["known_pc", "known_nc", "antisense", "novel", "unplaced"]
    label = {"known_pc": "annotated protein-coding (sense)",
             "known_nc": "annotated non-coding (sense)",
             "antisense": "antisense overlap only",
             "novel": "NO overlap — unannotated locus",
             "unplaced": "no genomic coordinates"}

    print("\n" + "=" * 74)
    print("PROTEIN SET vs REFERENCE ANNOTATION")
    print("=" * 74)
    tot = len(lens)
    for k in order:
        if cats.get(k):
            print(f"  {label[k]:<36} {cats[k]:>7,}  ({100*cats[k]/tot:>4.1f}%)")
    print(f"  {'total':<36} {tot:>7,}")

    print("\n  70-99 aa band (the AMP-candidate space)")
    btot = sum(band.values())
    for k in order:
        if band.get(k):
            print(f"    {label[k]:<34} {band[k]:>7,}  ({100*band[k]/max(1,btot):>4.1f}%)")
    print(f"    {'band total':<34} {btot:>7,}")

    recov = 100 * len(hit_pc) / max(1, n_pc)
    print("\n" + "=" * 74)
    print("RECOVERY OF ANNOTATED PROTEIN-CODING GENES")
    print("=" * 74)
    print(f"  distinct annotated protein-coding genes hit  {len(hit_pc):,}"
          f" / {n_pc:,}")
    print(f"  recovery                                     {recov:.1f}%")
    print("  a single tissue expressing 60-75% of an animal's genes is")
    print("  unremarkable; this is the number that tests whether 0.67x means")
    print("  'upstream was inflated' or 'we lost a third of the proteome'")

    busco = None
    if not args.no_busco:
        bb = find_busco()
        if bb is None:
            log("busco not found — skipping completeness "
                "(install or pass --no-busco to silence)")
        else:
            bdir = out / "busco"
            bdir.mkdir(parents=True, exist_ok=True)
            tag = args.busco_lineage
            busco = parse_busco(bdir, tag) if not args.force else None
            if busco is None:
                log(f"running BUSCO ({tag}, protein mode) — "
                    "this is the completeness test")
                cmd = [bb, "-i", str(pep), "-m", "proteins",
                       "-l", tag, "-o", f"pmdark_{tag}",
                       "--out_path", str(bdir), "-c", str(args.threads)]
                if args.force:
                    cmd.append("-f")
                log("$ " + " ".join(cmd))
                rc = subprocess.run(cmd).returncode
                if rc != 0:
                    log(f"busco exited {rc}; continuing without completeness")
                busco = parse_busco(bdir, tag)

    if busco:
        c = busco.get("Complete percentage")
        print("\n" + "=" * 74)
        print(f"BUSCO COMPLETENESS  ({args.busco_lineage}, protein mode)")
        print("=" * 74)
        print(f"  complete    {c}%  "
              f"(single {busco.get('Single copy percentage')}%, "
              f"duplicated {busco.get('Multi copy percentage')}%)")
        print(f"  fragmented  {busco.get('Fragmented percentage')}%")
        print(f"  missing     {busco.get('Missing percentage')}%")
        print(f"  markers     {busco.get('n_markers')}")
        print("-" * 74)
        # The rule was fixed in this file's docstring before the run.
        if c is None:
            verdict = "UNKNOWN — could not parse the summary"
        elif c >= 90:
            verdict = ("PASS — the protein set is complete. The 0.67x protein "
                       "ratio is\n  upstream inflation plus single-tissue "
                       "sampling, and Table 1 rests on\n  the assembly "
                       "metrics, not the protein count.")
        elif c >= 80:
            verdict = ("MARGINAL — report the shortfall. The protein-count\n"
                       "  comparison cannot be used as evidence.")
        else:
            verdict = ("FAIL — real genes are being lost upstream of here.\n"
                       "  Fix the assembly or ORF calling before defining the "
                       "dark proteome.")
        print(f"  {verdict}")

    report = {
        "pep": str(pep), "genome_gff3": str(gff3),
        "n_proteins": tot,
        "categories": dict(cats),
        "band_70_99": dict(band),
        "reference": {"n_genes": len(biotype), "n_protein_coding": n_pc},
        "recovery": {"genes_hit": len(hit_pc),
                     "pct_of_protein_coding": round(recov, 2)},
        "busco": busco,
        "busco_lineage": args.busco_lineage,
    }
    rp = res / "protein_qc.json"
    rp.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rp}")

    # The novel set is the dark-proteome candidate pool; give it its own file.
    novel_ids = set()
    for pid in lens:
        sp = spans.get(pid)
        if sp and classify(idx, *sp)[0] == "novel":
            novel_ids.add(pid)
    nf = out / "proteins_novel_loci.faa"
    w = False
    with pep.open() as fh, nf.open("w") as o:
        for line in fh:
            if line.startswith(">"):
                w = line[1:].split()[0] in novel_ids
            if w:
                o.write(line)
    print(f"wrote {nf.name} ({len(novel_ids):,} sequences)")

    print("\nNext: the proteinness filter on these sets "
          "(ESM-2 pseudo-log-likelihood vs matched decoys).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
