#!/usr/bin/env python3
"""Annotate the STAT-dependent uncharacterised genes with the calibrated method.

    python scripts/41_annotate_targets.py              # prepare, then transfer
    python scripts/41_annotate_targets.py --stage prepare
    python scripts/41_annotate_targets.py --hmm        # also harden with jackhmmer

WHAT THIS IS

  The paper's deliverable. scripts/40 produced 328 genes that are
  STAT-dependent under WSSV infection and carry no functional description
  in the reference annotation. This assigns each one GO terms using
  exactly the procedure that scripts/31-38 calibrated, and attaches to
  each the accuracy that was measured for its own situation rather than a
  single global number.

WHY THE TRANSFER IS NOT REIMPLEMENTED HERE

  The calibration only transfers to these proteins if the procedure is
  identical down to the tie-breaks. So this calls run_tier() from
  scripts/31_calibrate.py with the targets as queries and synthetic
  (empty) query metadata, and takes the propagated predictions it
  returns. The metrics run_tier computes alongside them are meaningless
  without ground truth and are discarded. Nothing about the neighbour
  search, the weighting, the 0.01 vote cutoff or the ancestor-max
  propagation is restated in this file.

  tier "T0_full" masks only ref_acc == query_acc, and no P. monodon
  accession appears in the benchmark reference, so nothing is excluded.

PER-PROTEIN CONFIDENCE, NOT ONE NUMBER

  Measured accuracy depends on whether a protein has a reachable homolog:

    has a homolog in the reference  -> T0_full      Fmax 0.815
    no DIAMOND and no jackhmmer hit -> orphan set   Fmax 0.625

  Each target is searched against the same reference set the GO comes
  from and assigned the stratum that applies to it. The score threshold
  used to emit terms is the one that MAXIMISED Fmax in that stratum,
  taken from the calibration file rather than chosen here, and the
  precision and recall measured at that threshold are reported next to
  every prediction.

  Terms are predictions, not findings. The honest reading of a BP term
  at the orphan threshold is "about as reliable as Fmax 0.37 on a
  held-out benchmark", which is why the per-ontology numbers are carried
  through to the output.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "src"))

from pmdark import config as cfgmod                              # noqa: E402
from pmdark.goeval import GOGraph                                # noqa: E402
from pmdark.shards import Record, read_fasta, write_fasta        # noqa: E402
from pmdark.transfer import whiten                               # noqa: E402

GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_mod(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot import {filename}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def go_names(obo: Path) -> dict[str, str]:
    """GOGraph keeps namespaces but not names; a term id alone is unreadable
    in a results table, so parse the names too."""
    names: dict[str, str] = {}
    tid = None
    with obo.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                tid = None
            elif line.startswith("id: "):
                tid = line[4:]
            elif line.startswith("name: ") and tid:
                names[tid] = line[6:]
                tid = None
    return names


def gene_to_protein(gff: Path) -> dict[str, tuple[str, int]]:
    """gene_id -> (protein_id, CDS length) for its longest coding isoform.

    Walks mRNA (ID=rna-*, Parent=gene-*) then CDS (Parent=rna-*,
    protein_id=XP_*), because protein ids live on CDS features and the
    gene link lives on mRNA features.
    """
    rna_gene: dict[str, str] = {}
    rna_prot: dict[str, str] = {}
    rna_len: dict[str, int] = defaultdict(int)
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat, start, end, attrs = f[2], f[3], f[4], f[8]
            if feat == "mRNA":
                a = dict(GFF_ATTR.findall(attrs))
                rid, par = a.get("ID"), a.get("Parent", "")
                if rid and par.startswith("gene-"):
                    rna_gene[rid] = par
            elif feat == "CDS":
                a = dict(GFF_ATTR.findall(attrs))
                par, pid = a.get("Parent", ""), a.get("protein_id", "")
                if not par:
                    continue
                if pid:
                    rna_prot.setdefault(par, pid)
                try:
                    rna_len[par] += int(end) - int(start) + 1
                except ValueError:
                    pass
    best: dict[str, tuple[str, int]] = {}
    for rid, gid in rna_gene.items():
        pid, L = rna_prot.get(rid), rna_len.get(rid, 0)
        if not pid or L <= 0:
            continue
        if L > best.get(gid, ("", 0))[1]:
            best[gid] = (pid, L)
    return best


def diamond_best(exe: str, query: Path, db_fasta: Path, out: Path,
                 threads: int, force: bool) -> dict[str, float]:
    """target -> best percent identity against the reference set."""
    db = out.parent / "ref_for_targets"
    if not db.with_suffix(".dmnd").exists():
        log("building DIAMOND db of the benchmark reference")
        if subprocess.run([exe, "makedb", "--in", str(db_fasta),
                           "-d", str(db), "--threads", str(threads)
                           ]).returncode != 0:
            raise SystemExit("diamond makedb failed")
    if force or not out.exists() or out.stat().st_size == 0:
        log("searching targets against the reference (--ultra-sensitive)")
        rc = subprocess.run([
            exe, "blastp", "-q", str(query), "-d", str(db),
            "-o", str(out), "--outfmt", "6", "qseqid", "sseqid", "pident",
            "evalue", "bitscore", "--evalue", "10", "--max-target-seqs", "50",
            "--ultra-sensitive", "--threads", str(threads)]).returncode
        if rc != 0:
            raise SystemExit(f"diamond blastp failed ({rc})")
    else:
        log(f"reusing {out.name}")
    best: dict[str, float] = {}
    with out.open() as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 3:
                continue
            try:
                pid = float(p[2])
            except ValueError:
                continue
            if pid > best.get(p[0], -1.0):
                best[p[0]] = pid
    return best


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", default=None,
                    help="gene id list (default: "
                         "results/stat_dependent_uncharacterised.txt)")
    ap.add_argument("--model", default="esm2")
    ap.add_argument("--transform", default="whiten",
                    choices=["none", "center", "whiten"])
    ap.add_argument("--stage", default="all",
                    choices=["all", "prepare", "transfer"])
    ap.add_argument("--k", type=int, default=25)
    ap.add_argument("--sim-floor", type=float, default=0.0)
    ap.add_argument("--power", type=float, default=3.0)
    ap.add_argument("--block", type=int, default=512)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--hmm", action="store_true",
                    help="also run jackhmmer, so the orphan call matches the "
                         "hardened definition from scripts/37")
    ap.add_argument("--min-score", type=float, default=None,
                    help="override the calibrated score threshold "
                         "(not recommended; it breaks the calibration)")
    ap.add_argument("--diamond", default=None,
                    help="explicit path to the diamond binary")
    ap.add_argument("--force-search", action="store_true")
    args = ap.parse_args()

    cfg = cfgmod.load()
    calib = load_mod("pmdark_calib", "31_calibrate.py")
    rec27 = load_mod("pmdark_rec", "27_recovery_exact.py")

    bench = cfg.p("proteins", "benchmark")
    work = cfg.p("proteins", "annotate")
    work.mkdir(parents=True, exist_ok=True)

    tfile = Path(args.targets) if args.targets else (
        cfg.p("results") / "stat_dependent_uncharacterised.txt")
    if not tfile.exists():
        raise SystemExit(f"missing {tfile} — run scripts/40 first")
    genes = [ln.strip() for ln in tfile.read_text().splitlines() if ln.strip()]
    log(f"{len(genes)} target genes from {tfile.name}")

    # ---- gene -> protein ----------------------------------------------
    gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    if gff is None:
        raise SystemExit("no reference GFF in the genome directory")
    log(f"mapping genes to proteins via {gff.name}")
    g2p = gene_to_protein(gff)
    pairs = [(g, g2p[g][0], g2p[g][1]) for g in genes if g in g2p]
    missing = [g for g in genes if g not in g2p]
    log(f"  {len(pairs)} genes have a coding isoform; {len(missing)} do not")
    if missing[:5]:
        log(f"  e.g. {', '.join(missing[:5])}")
    if not pairs:
        raise SystemExit("no target gene has a protein; nothing to annotate")

    # ---- target FASTA -------------------------------------------------
    asm = gff.name.replace("_genomic.gff", "")
    prot = rec27.stage_proteome(cfg, asm)
    if prot is None or not prot.exists():
        raise SystemExit(
            "could not stage the reference proteome; "
            "run scripts/27_recovery_exact.py once to fetch it")
    log(f"reading {prot.name}")
    by_pid = {r.pid.split()[0]: r for r in read_fasta(prot)}
    want = {pid: g for g, pid, _ in pairs}
    recs: list[Record] = []
    for pid, g in want.items():
        r = by_pid.get(pid)
        if r is not None:
            # name the sequence by GENE, so every downstream table keys on
            # the thing scripts/40 selected rather than an isoform id
            recs.append(Record(pid=g, seq=r.seq.rstrip("*")))
    log(f"  {len(recs)} target proteins extracted")
    fasta = work / "targets.faa"
    write_fasta(recs, fasta)
    log(f"wrote {fasta}")

    # ---- stratify: does each target have a reachable homolog? ---------
    # diamond is installed in the busco env, not pmbio, so shutil.which()
    # alone finds nothing and every target would silently be labelled as
    # having a homolog -- overstating confidence on precisely the proteins
    # this paper is about. Reuse the finder that already knows the layout.
    orf = load_mod("pmdark_orf", "24_call_orfs.py")
    exe = orf.find_diamond(args.diamond)
    ident: dict[str, float] = {}
    if exe:
        ident = diamond_best(exe, fasta, bench / "reference.faa",
                             work / "targets_vs_reference.tsv",
                             args.threads, args.force_search)
        log(f"DIAMOND: {len(ident)} of {len(recs)} targets have a hit")
    else:
        print("\n  diamond NOT FOUND. Without it no target can be shown to")
        print("  lack a homolog, so all of them would be scored against the")
        print("  T0 threshold and their reliability OVERSTATED. Refusing.")
        print("    conda install -n pmbio -c bioconda diamond")
        print("    or pass --diamond /path/to/diamond")
        return 2

    hmm_hit: set[str] = set()
    if args.hmm:
        hard = load_mod("pmdark_hard", "37_harden_orphans.py")
        jh = hard.find_jackhmmer(None)
        have_py = importlib.util.find_spec("pyhmmer") is not None
        tbl = work / "targets_jackhmmer.tbl"
        ok = False
        if jh:
            ok = hard.run_jackhmmer(jh, fasta, bench / "reference.faa", tbl,
                                    5, 10.0, args.threads, args.force_search)
        elif have_py:
            ok = hard.run_pyhmmer(fasta, bench / "reference.faa", tbl, 5,
                                  args.threads, args.force_search)
        if ok:
            best = hard.parse_tbl(tbl)
            hmm_hit = {q for q, e in best.items() if e <= 1e-3}
            log(f"jackhmmer: {len(hmm_hit)} of {len(recs)} have a hit "
                f"at E<=1e-3")
        else:
            log("jackhmmer unavailable; orphan calls rest on DIAMOND alone")

    orphan = {r.pid for r in recs
              if r.pid not in ident and r.pid not in hmm_hit}
    log(f"{len(orphan)} targets have no reachable homolog "
        f"({100 * len(orphan) / max(1, len(recs)):.1f}%)")

    # ---- embeddings ---------------------------------------------------
    try:
        q_ids, q_emb = calib.load_embeddings(cfg, args.model, fasta)
    except SystemExit as e:
        print("\n" + "=" * 78)
        print("EMBEDDINGS NEEDED")
        print("=" * 78)
        print(str(e))
        print(f"\n  bash slurm/embed_fanout.sh {args.model} {fasta}")
        print("\n  Then re-run this script; everything above is cached.")
        return 3
    r_ids, r_emb = calib.load_embeddings(cfg, args.model, bench / "reference.faa")
    log(f"embeddings: targets {q_emb.shape}, reference {r_emb.shape}")

    if args.stage == "prepare":
        print("\n(--stage prepare: stopping before transfer)")
        return 0

    # ---- transform, exactly as calibrated -----------------------------
    if args.transform == "whiten":
        r_t, st = whiten(r_emb.astype(np.float32))
        q_t, _ = whiten(q_emb.astype(np.float32), st)
    elif args.transform == "center":
        from pmdark.transfer import center
        r_t, mu = center(r_emb.astype(np.float32))
        q_t, _ = center(q_emb.astype(np.float32), mu)
    else:
        q_t, r_t = q_emb.astype(np.float32), r_emb.astype(np.float32)
    log(f"embeddings {args.transform}ed on reference statistics")

    with __import__("gzip").open(bench / "benchmark.json.gz", "rt") as fh:
        meta = json.load(fh)
    r_meta = meta["reference"]
    r_keep = [a for a in r_ids if a in r_meta]
    ridx = {a: i for i, a in enumerate(r_ids)}
    r_t = r_t[[ridx[a] for a in r_keep]]

    obo = cfg.p("raw", "benchmark") / "go-basic.obo"
    log(f"loading {obo.name}")
    graph = GOGraph.from_obo(obo)
    names = go_names(obo)

    # Synthetic query metadata: run_tier needs these keys but T0_full masks
    # on accession only, and empty "go" simply makes its metrics empty.
    q_meta = {a: {"taxid": "pmon", "order": "", "phylum": False, "go": []}
              for a in q_ids}
    log("transferring GO (run_tier, T0_full, nothing masked)")
    out = calib.run_tier("T0_full", q_ids, q_t, q_meta, r_keep, r_t, r_meta,
                         graph, args.k, args.sim_floor, args.power,
                         args.block, None, return_predictions=True)
    if len(out) != 4:
        raise SystemExit("31_calibrate.py needs the return_predictions patch")
    _m, diag, _truth, prop = out
    log(f"  median top similarity {diag['median_top_sim']:.3f}; "
        f"{diag['no_neighbour']} targets with no usable neighbour")

    # ---- calibrated thresholds ----------------------------------------
    cal_path = cfg.p("results") / \
        f"calibration_{args.model}_{args.transform}.json"
    if not cal_path.exists():
        raise SystemExit(
            f"missing {cal_path.name}; the per-protein confidence comes from "
            "it.\n  python scripts/31_calibrate.py --model "
            f"{args.model} --transform {args.transform} "
            "--tiers T0_full,I90,I50,I30,I20,Inone")
    cal = json.loads(cal_path.read_text())
    strata = {}
    for key, label in (("ORPHAN_nohit", "orphan"), ("T0_full", "homolog")):
        blk = cal["results"].get(key)
        if not blk:
            continue
        ov = blk["metrics"]["overall"]
        strata[label] = {
            "fmax": ov["fmax"], "threshold": ov["threshold"],
            "precision": ov["precision"], "recall": ov["recall"],
            "BP": blk["metrics"].get("BP", {}).get("fmax"),
            "MF": blk["metrics"].get("MF", {}).get("fmax"),
            "CC": blk["metrics"].get("CC", {}).get("fmax"),
            "n_benchmark": ov.get("n"),
        }
    # The hardened orphan number is the one that applies. Our orphan
    # targets are selected the same way scripts/37 selected its survivors
    # -- no DIAMOND hit at e-value 10 AND no jackhmmer hit at E<=1e-3 --
    # so the 94-protein ORPHAN_nohit figure in the calibration file is
    # the WRONG, optimistic one for them.
    hard_path = cfg.p("results") / "orphan_hardening.json"
    if hard_path.exists():
        try:
            hj = json.loads(hard_path.read_text())
            h = hj.get("strata", {}).get("HARDENED")
            if h and h.get("fmax") is not None:
                prev = strata.get("orphan", {}).get("fmax")
                strata["orphan"] = {
                    "fmax": h["fmax"],
                    "ci": h.get("ci"),
                    "threshold": h.get("threshold",
                                       strata.get("orphan", {}).get(
                                           "threshold", 0.3)),
                    "precision": h.get("precision",
                                       strata.get("orphan", {}).get(
                                           "precision")),
                    "recall": h.get("recall"),
                    "BP": h.get("BP"), "MF": h.get("MF"), "CC": h.get("CC"),
                    "n_benchmark": h.get("n"),
                    "source": "hardened (DIAMOND + jackhmmer)",
                    "harden_at": hj.get("harden_at"),
                }
                log(f"orphan calibration: using the HARDENED set "
                    f"(n={h.get('n')}, Fmax {h['fmax']:.3f})"
                    + (f", not the un-hardened {prev:.3f}"
                       if prev is not None else ""))
                if strata["orphan"].get("BP") is None:
                    log("  per-ontology numbers absent from "
                        "orphan_hardening.json; re-run scripts/37 to add "
                        "them, otherwise BP/MF/CC are left blank rather "
                        "than filled from the optimistic set")
        except (ValueError, KeyError, TypeError) as e:
            log(f"could not read {hard_path.name} ({e}); "
                "falling back to the un-hardened orphan calibration")
    else:
        log("WARNING: no orphan_hardening.json. The orphan calibration is "
            "the un-hardened 94-protein figure, which OVERSTATES these "
            "proteins. Run scripts/37_harden_orphans.py.")

    if "homolog" not in strata:
        raise SystemExit("calibration file has no T0_full result")
    if "orphan" not in strata:
        log("WARNING: no ORPHAN_nohit in the calibration; orphan targets "
            "will be scored against the T0 threshold, which OVERSTATES "
            "their reliability")
        strata["orphan"] = dict(strata["homolog"], overstated=True)

    print("\n" + "=" * 78)
    print(f"CALIBRATION APPLIED  ({args.model} + {args.transform})")
    print("=" * 78)
    for lab, c in strata.items():
        print(f"  {lab:<9} Fmax {c['fmax']:.3f}  threshold {c['threshold']:.2f}"
              f"  precision {c['precision']:.3f}  recall {c['recall']:.3f}")
        print(f"            BP {c['BP']:.3f}  MF {c['MF']:.3f}  "
              f"CC {c['CC']:.3f}" if c.get("BP") else "")

    # ---- is the calibration transferable to these proteins? -----------
    # The calibration was measured on benchmark proteins. If our targets
    # sit at systematically lower neighbour similarity, their neighbours
    # are worse and the measured accuracy does not carry over.
    cal_sim = None
    for key in ("ORPHAN_nohit", "T0_full"):
        blk = cal["results"].get(key)
        if blk and blk.get("diagnostics", {}).get("median_top_sim"):
            cal_sim = (key, blk["diagnostics"]["median_top_sim"])
            if key == "ORPHAN_nohit":
                break
    print("\n" + "-" * 78)
    print("IS THE CALIBRATION TRANSFERABLE?")
    print("-" * 78)
    print(f"  targets, median top similarity      : "
          f"{diag['median_top_sim']:.3f}")
    if cal_sim:
        print(f"  calibration set ({cal_sim[0]}) : {cal_sim[1]:.3f}")
        shift = diag["median_top_sim"] - cal_sim[1]
        print(f"  shift                               : {shift:+.3f}")
        if shift < -0.05:
            print("\n  The targets sit LOWER than the proteins the accuracy")
            print("  was measured on, so their nearest neighbours are worse")
            print("  and the reported Fmax is OPTIMISTIC for them. Say so,")
            print("  and treat the numbers as an upper bound.")
        elif shift > 0.05:
            print("\n  The targets sit higher than the calibration set, so")
            print("  the reported Fmax is conservative for them.")
        else:
            print("\n  Comparable. The calibration transfers.")
    else:
        print("  calibration set: no diagnostic recorded; cannot compare")

    # ---- emit ---------------------------------------------------------
    long_rows = []
    summary = []
    for g in q_ids:
        lab = "orphan" if g in orphan else "homolog"
        c = strata[lab]
        thr = args.min_score if args.min_score is not None else c["threshold"]
        terms = {t: s for t, s in prop.get(g, {}).items() if s >= thr}
        per_ns = {"BP": [], "MF": [], "CC": []}
        for t, s in sorted(terms.items(), key=lambda kv: -kv[1]):
            ns = graph.namespace.get(t, "")
            long_rows.append({
                "gene": g, "stratum": lab, "go_id": t,
                "namespace": ns, "go_name": names.get(t, ""),
                "score": round(float(s), 4),
                "score_threshold": thr,
                "stratum_fmax": c["fmax"],
                "stratum_precision": c["precision"],
                "ontology_fmax": c.get(ns),
            })
            if ns in per_ns:
                per_ns[ns].append((t, s))
        summary.append({
            "gene": g, "stratum": lab,
            "best_identity_pct": round(ident.get(g, float("nan")), 1)
            if g in ident else None,
            "jackhmmer_hit": g in hmm_hit,
            "n_terms": len(terms),
            "n_BP": len(per_ns["BP"]), "n_MF": len(per_ns["MF"]),
            "n_CC": len(per_ns["CC"]),
            "top_MF": per_ns["MF"][0][0] if per_ns["MF"] else "",
            "top_MF_name": names.get(per_ns["MF"][0][0], "")
            if per_ns["MF"] else "",
            "top_BP": per_ns["BP"][0][0] if per_ns["BP"] else "",
            "top_BP_name": names.get(per_ns["BP"][0][0], "")
            if per_ns["BP"] else "",
            "expected_fmax": c["fmax"],
        })

    res_dir = cfg.p("results")
    res_dir.mkdir(parents=True, exist_ok=True)
    try:
        import pandas as pd
        pd.DataFrame(long_rows).to_csv(
            res_dir / "annotation_targets_long.csv", index=False)
        pd.DataFrame(summary).to_csv(
            res_dir / "annotation_targets.csv", index=False)
    except ImportError:
        import csv
        for fn, rows in (("annotation_targets_long.csv", long_rows),
                         ("annotation_targets.csv", summary)):
            if not rows:
                continue
            with (res_dir / fn).open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)

    annotated = sum(1 for s in summary if s["n_terms"] > 0)
    n_orph = sum(1 for s in summary if s["stratum"] == "orphan")
    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    print(f"  {len(summary)} targets, {annotated} received at least one term "
          f"({100 * annotated / max(1, len(summary)):.1f}%)")
    print(f"  {n_orph} are orphans (expected Fmax "
          f"{strata['orphan']['fmax']:.3f}), "
          f"{len(summary) - n_orph} have a homolog "
          f"(expected Fmax {strata['homolog']['fmax']:.3f})")
    print(f"  {len(long_rows):,} gene-term predictions above the calibrated "
          f"thresholds")
    with_mf = sum(1 for s in summary if s["top_MF"])
    with_bp = sum(1 for s in summary if s["top_BP"])
    cc_only = sum(1 for s in summary
                  if s["n_terms"] > 0 and not s["top_MF"] and not s["top_BP"])
    n = max(1, len(summary))
    print("\n  USEFUL YIELD — report this, not the 'any term' figure:")
    print(f"    molecular function : {with_mf:>4} / {len(summary)} "
          f"({100 * with_mf / n:.1f}%)")
    print(f"    biological process : {with_bp:>4} / {len(summary)} "
          f"({100 * with_bp / n:.1f}%)")
    print(f"    cellular comp. ONLY: {cc_only:>4} / {len(summary)} "
          f"({100 * cc_only / n:.1f}%)")
    print("  A gene with only CC terms has been localised, not functionally")
    print("  annotated. CC is the easiest ontology (0.81 on orphans against")
    print("  0.44 for MF), so 'annotated' without this breakdown overstates")
    print("  what was achieved.")

    (res_dir / "annotation_summary.json").write_text(json.dumps({
        "model": args.model, "transform": args.transform,
        "params": {"k": args.k, "sim_floor": args.sim_floor,
                   "power": args.power},
        "n_targets_requested": len(genes),
        "n_without_protein": len(missing),
        "n_targets": len(summary), "n_annotated": annotated,
        "n_orphan": n_orph, "n_predictions": len(long_rows),
        "yield": {"molecular_function": with_mf,
                  "biological_process": with_bp,
                  "cellular_component_only": cc_only},
        "transferability": {
            "target_median_top_sim": diag["median_top_sim"],
            "calibration_median_top_sim": cal_sim[1] if cal_sim else None,
            "calibration_stratum": cal_sim[0] if cal_sim else None},
        "calibration": strata,
        "diamond_hits": len(ident), "jackhmmer_hits": len(hmm_hit),
        "median_top_similarity": diag["median_top_sim"],
        "no_neighbour": diag["no_neighbour"],
    }, indent=2))
    log(f"wrote {res_dir / 'annotation_targets.csv'}")
    log(f"wrote {res_dir / 'annotation_targets_long.csv'}")

    print("\n  These are PREDICTIONS carrying measured error rates, not")
    print("  findings. Every row reports the stratum it was scored in and")
    print("  that stratum's Fmax, including the per-ontology numbers: BP is")
    print("  by far the weakest and should never be quoted at the overall")
    print("  Fmax. Nothing here is experimental evidence of function.")
    if not args.hmm:
        print("\n  Run again with --hmm to make the orphan calls match the")
        print("  hardened definition. Without it the orphan set rests on")
        print("  DIAMOND alone, which scripts/37 showed misses homologs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
