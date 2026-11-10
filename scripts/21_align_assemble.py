#!/usr/bin/env python3
"""Per-sample align + transcript assembly, self-contained on a compute node.

ARCHITECTURE NOTE — this is why the script looks the way it does:

  /shared/a2-public is NOT mounted on compute nodes. A batch job cannot read
  the SRA mirror we downloaded. But compute nodes DO have internet
  (ebi.ac.uk -> 200), so each job fetches its own library straight from ENA
  into node-local /tmp, does all heavy I/O there on local disk, and writes
  back only small outputs to the quota'd HOME tier.

  ENA  --download-->  /tmp/pmdark_$JOBID/   (local disk, 64 GB, fast)
                          fastp
                          hisat2 | samtools sort
                          stringtie
                          unmapped read extraction
  /tmp  --copy-->     $HOME/pmdark_work/align/<run>/   (GTF, counts, stats)
  /tmp  wiped on exit

Run as a Slurm array, one task per library:

    sbatch --array=0-11 slurm/21_align.sbatch
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

ENA_FASTQ_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run(cmd: list[str] | str, shell: bool = False, check: bool = True) -> int:
    pretty = cmd if isinstance(cmd, str) else " ".join(cmd)
    log(f"$ {pretty[:300]}")
    rc = subprocess.run(cmd, shell=shell).returncode
    if check and rc != 0:
        raise RuntimeError(f"command failed (rc={rc}): {pretty[:300]}")
    return rc


def working_perl() -> str | None:
    """Find a Perl that can actually start.

    conda bakes the environment prefix into Perl's Config_heavy.pl. On this
    cluster the prefix contains an '@':

        /storage/slurm/home/68090800301@cpe.kmutt.ac.th/pmdark_work/envs/pmbio

    and Perl parses '@cpe' as an array variable inside those strings:

        Global symbol "@cpe" requires explicit package name ...
        Compilation failed in require at .../Config.pm line 105

    so conda's Perl will not run AT ALL in this env. hisat2 and TransDecoder
    are both Perl wrapper scripts, so both die. The system Perl lives at a
    path with no '@' and works, hence this probe.
    """
    for cand in ("/usr/bin/perl", "/bin/perl", shutil.which("perl")):
        if not cand or not Path(cand).exists():
            continue
        try:
            r = subprocess.run([cand, "-e", "print 1"],
                               capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and r.stdout.strip() == "1":
                return cand
        except (OSError, subprocess.TimeoutExpired):
            continue
    return None


def perl_tool(name: str, perl: str | None) -> list[str]:
    """Build the invocation for a tool that may be a Perl script."""
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"{name} not on PATH")
    try:
        first = Path(path).read_bytes()[:64].decode("utf-8", "ignore")
    except OSError:
        return [path]
    if first.startswith("#!") and "perl" in first.split("\n")[0]:
        if perl is None:
            raise RuntimeError(
                f"{name} is a Perl script but no working Perl was found. "
                "See working_perl() for why conda's Perl is broken here."
            )
        return [perl, path]
    return [path]


def read_samples(path: Path) -> list[dict]:
    rows = []
    header = None
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
    return rows


def ena_urls(acc: str) -> list[str]:
    import urllib.request
    url = ENA_FASTQ_API.format(acc=acc)
    req = urllib.request.Request(url, headers={"User-Agent": "pmdark/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        lines = r.read().decode().strip().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"ENA returned no file report for {acc}")
    cols = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    urls = [u for u in cols.get("fastq_ftp", "").split(";") if u]
    if not urls:
        raise RuntimeError(f"no fastq_ftp for {acc}")
    return [u if u.startswith("http") else "https://" + u for u in urls]


def fetch_reads(acc: str, dest: Path, threads: int, archive: Path | None) -> list[Path]:
    """Prefer the archive copy if this process can see it (login node);
    otherwise download from ENA (the compute-node path)."""
    dest.mkdir(parents=True, exist_ok=True)

    if archive and archive.is_dir():
        local = sorted(archive.glob(f"{acc}_*.fastq.gz"))
        if len(local) >= 2:
            log(f"archive copy visible — linking {len(local)} files")
            out = []
            for src in local:
                dst = dest / src.name
                if not dst.exists():
                    try:
                        os.symlink(src, dst)
                    except OSError:
                        shutil.copy2(src, dst)
                out.append(dst)
            return out

    log(f"downloading {acc} from ENA")
    out = []
    for url in ena_urls(acc):
        dst = dest / Path(url).name
        if dst.exists() and dst.stat().st_size > 0:
            out.append(dst)
            continue
        if shutil.which("aria2c"):
            run(["aria2c", "-c", "-x", "8", "-s", "8", "--console-log-level=warn",
                 "-d", str(dest), "-o", dst.name, url])
        else:
            run(["curl", "-fsSL", "-C", "-", "-o", str(dst), url])
        out.append(dst)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=Path, default=None,
                    help="HISAT2 index prefix (default: work/genome/pmon)")
    ap.add_argument("--sample", default=None, help="run accession; default from array id")
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--force", action="store_true",
                    help="redo a library even if its DONE marker exists")
    ap.add_argument("--guided", action="store_true", default=True,
                    help="StringTie -G (guide, not -e). On by default; see "
                         "the correction note at step 4.")
    ap.add_argument("--no-guided", dest="guided", action="store_false")
    ap.add_argument("--ref-quant", action="store_true", default=True,
                    help="second StringTie pass with -e -G to quantify the "
                         "reference annotation — the expression denominator")
    ap.add_argument("--no-ref-quant", dest="ref_quant", action="store_false")
    ap.add_argument("--keep-bam", action="store_true",
                    help="copy the sorted BAM back to WORK (large — off by default)")
    ap.add_argument("--save-unmapped", action="store_true", default=True,
                    help="keep unmapped reads for the de novo arm")
    args = ap.parse_args()

    cfg = cfgmod.load()
    samples = read_samples(Path(__file__).resolve().parents[1] / "config" / "samples.tsv")

    if args.sample:
        row = next((r for r in samples if r["run_accession"] == args.sample), None)
        if row is None:
            print(f"unknown sample {args.sample}")
            return 1
    else:
        idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
        if idx >= len(samples):
            print(f"array index {idx} >= {len(samples)} samples — nothing to do")
            return 0
        row = samples[idx]

    acc = row["run_accession"]
    group = row["group"]
    index = args.index or (cfg.p("genome") / "pmon")

    out_dir = cfg.p("align") / acc
    out_dir.mkdir(parents=True, exist_ok=True)
    done_marker = out_dir / "DONE"
    if done_marker.exists() and not args.force:
        # The marker from the un-guided run means "aligned and assembled",
        # not "assembled the way we now need". --force re-does the library;
        # the re-run is what adds -G and the -e -G denominator pass.
        log(f"{acc} already complete — {done_marker}")
        log("  (--force to redo with -G and the reference-quantification pass)")
        return 0

    # Scratch sizing. Empirically ~0.1 GB of compressed FASTQ per million
    # reads in this study. Peak concurrent usage with the cleanup below is
    # trimmed + BAM + sort temp, so budget ~0.35 GB per million reads with a
    # floor for small libraries.
    reads_m = int(row.get("reads", 0) or 0) / 1e6
    need_gb = max(6.0, reads_m * 0.35)

    log(f"sample   : {acc}  ({row['sample_title']}, group={group})")
    log(f"reads    : {reads_m:.1f} M  -> scratch budget {need_gb:.0f} GB")
    scratch = cfgmod.scratch(cfg, need_gb=need_gb)
    log(f"index    : {index}")
    log(f"out      : {out_dir}")
    log(f"threads  : {args.threads}")

    if not Path(str(index) + ".1.ht2").exists():
        log(f"ERROR: HISAT2 index not found at {index}.*.ht2")
        log("Build it first: sbatch slurm/20_index.sbatch")
        return 1

    free_gb = shutil.disk_usage(scratch).free / 1e9
    if free_gb < need_gb:
        log(f"WARNING: only {free_gb:.0f} GB free, budget was {need_gb:.0f} GB. "
            "Proceeding with aggressive cleanup; may still fail.")

    t_start = time.time()
    stats: dict = {"run": acc, "group": group, "title": row["sample_title"]}

    # ---- 1. reads ------------------------------------------------------
    archive = cfg.p("raw", "sra")
    reads = fetch_reads(acc, scratch / "reads", args.threads,
                        archive if archive.exists() else None)
    if len(reads) < 2:
        log(f"ERROR: expected paired FASTQ, got {len(reads)}")
        return 1
    r1, r2 = reads[0], reads[1]
    stats["input_gb"] = round(sum(p.stat().st_size for p in reads) / 1e9, 2)

    # ---- 2. QC ---------------------------------------------------------
    qc = scratch / "qc"
    qc.mkdir(exist_ok=True)
    t1, t2 = qc / f"{acc}_1.fq.gz", qc / f"{acc}_2.fq.gz"
    run(["fastp", "-i", str(r1), "-I", str(r2), "-o", str(t1), "-O", str(t2),
         "--detect_adapter_for_pe", "--thread", str(min(args.threads, 16)),
         "--json", str(out_dir / "fastp.json"),
         "--html", str(out_dir / "fastp.html"),
         "--length_required", "36"])

    # Free the raw reads immediately. On cpu01 /tmp has only ~9.7 GB, so
    # holding raw and trimmed copies at once is what tips a large library
    # over. Symlinks (archive-visible case) are unlinked, not the originals.
    for p in reads:
        if p.is_symlink():
            p.unlink()
        else:
            p.unlink(missing_ok=True)
    log(f"freed raw reads ({stats['input_gb']} GB); "
        f"{shutil.disk_usage(scratch).free / 1e9:.0f} GB now free")

    # ---- 3. align ------------------------------------------------------
    bam = scratch / f"{acc}.sorted.bam"
    un_prefix = scratch / f"{acc}_unmapped"
    hisat_log = out_dir / "hisat2.log"

    perl = working_perl()
    if perl:
        log(f"perl     : {perl}")
    hisat_cmd = [
        *perl_tool("hisat2", perl),
        "-p", str(args.threads), "-x", str(index),
        "-1", str(t1), "-2", str(t2),
        "--dta",                      # required downstream by StringTie
        "--no-unal" if not args.save_unmapped else "--un-conc-gz",
    ]
    if args.save_unmapped:
        hisat_cmd.append(f"{un_prefix}_%.fq.gz")

    sort_threads = max(1, args.threads // 2)
    cmd = (f"{' '.join(hisat_cmd)} 2> {hisat_log} | "
           f"samtools sort -@ {sort_threads} -m 1G "
           f"-T {scratch}/sorttmp -o {bam} -")
    run(cmd, shell=True)
    run(["samtools", "index", "-@", str(sort_threads), str(bam)])

    if hisat_log.exists():
        txt = hisat_log.read_text()
        stats["hisat2_log"] = txt.strip().splitlines()
        for line in txt.splitlines():
            if "overall alignment rate" in line:
                stats["alignment_rate"] = line.strip().split("%")[0]

    # Trimmed reads are no longer needed once the BAM exists.
    for p in (t1, t2):
        p.unlink(missing_ok=True)
    log(f"freed trimmed reads; "
        f"{shutil.disk_usage(scratch).free / 1e9:.0f} GB now free")

    # ---- 4. transcript assembly ----------------------------------------
    # CORRECTION. This step originally ran with no -G, on my reasoning that
    # a reference guide "biases toward known genes" and novel loci are the
    # point. That conflated two StringTie flags:
    #
    #   -e   RESTRICTS output to the reference. Novel loci are lost.
    #   -G   supplies the reference as a GUIDE. Novel transcripts and
    #        novel loci are still assembled and reported.
    #
    # -G without -e costs nothing in novel discovery and lets StringTie
    # complete transcripts at annotated loci that sparse coverage would
    # otherwise leave as sub-threshold fragments. Diagnosis in
    # 29_fragmentation.py: coverage of the annotated protein falls from
    # 100% to 56% across the protein-length bands.
    ref_gff = next(iter(cfg.p("genome").glob("*_genomic.gff")), None)
    gtf = out_dir / f"{acc}.gtf"
    st = ["stringtie", str(bam), "-p", str(args.threads), "-o", str(gtf),
          "-l", acc, "-A", str(out_dir / f"{acc}.gene_abund.tsv")]
    if args.guided and ref_gff:
        st[2:2] = ["-G", str(ref_gff)]
        log(f"guided assembly: -G {ref_gff.name} (no -e, novel loci kept)")
    elif args.guided:
        log("WARNING --guided requested but no reference GFF found")
    run(st)

    # ---- 4b. reference-gene expression, independent of our assembly -----
    # Every recovery figure so far divides by all 24,011 annotated genes,
    # including the ones a hemocyte never transcribes, because we have
    # never measured expression at the READ level. A gene with ten reads is
    # expressed but will not assemble, and our assembly cannot tell those
    # apart from a gene with zero reads.
    #
    # -e -G quantifies the reference annotation and nothing else, so the
    # gene table it writes is a true denominator. This is the one place -e
    # is correct: it is a measurement, not the assembly.
    if args.ref_quant and ref_gff:
        rq = out_dir / f"{acc}.ref_expr.gtf"
        ra = out_dir / f"{acc}.ref_gene_abund.tsv"
        log("quantifying the reference annotation (-e -G) for the denominator")
        run(["stringtie", str(bam), "-p", str(args.threads),
             "-e", "-G", str(ref_gff), "-o", str(rq), "-A", str(ra)])

    # ---- 5. unmapped reads -> archive-bound, for the de novo arm --------
    if args.save_unmapped:
        um = sorted(scratch.glob(f"{acc}_unmapped_*.fq.gz"))
        if um:
            total = sum(p.stat().st_size for p in um) / 1e9
            stats["unmapped_gb"] = round(total, 2)
            keep = cfg.p("align") / "_unmapped"
            keep.mkdir(parents=True, exist_ok=True)
            # These feed rnaSPAdes later. They are the only large thing we
            # keep on the quota'd tier, so warn loudly if they are big.
            if total > 3.0:
                log(f"WARNING: unmapped reads {total:.1f} GB — quota pressure")
            for p in um:
                shutil.copy2(p, keep / p.name)
            log(f"kept {len(um)} unmapped files ({total:.2f} GB) in {keep}")

    if args.keep_bam:
        shutil.copy2(bam, out_dir / bam.name)
        shutil.copy2(str(bam) + ".bai", out_dir / (bam.name + ".bai"))

    # ---- 6. stats -------------------------------------------------------
    stats["elapsed_min"] = round((time.time() - t_start) / 60, 1)
    stats["bam_gb"] = round(bam.stat().st_size / 1e9, 2)
    stats["gtf_transcripts"] = sum(
        1 for l in gtf.read_text().splitlines()
        if not l.startswith("#") and "\ttranscript\t" in l
    )
    (out_dir / "stats.json").write_text(json.dumps(stats, indent=2))
    done_marker.write_text(json.dumps(stats, indent=2))

    log(f"DONE {acc}: {stats['gtf_transcripts']:,} transcripts "
        f"in {stats['elapsed_min']} min")
    log(json.dumps(stats, indent=2))

    # scratch is node-local and wiped by Slurm; remove explicitly anyway so a
    # node with a sticky /tmp does not accumulate across array tasks
    shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
