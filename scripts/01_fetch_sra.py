#!/usr/bin/env python3
"""Resolve and download SRA study SRP499992, resumably.

The study holds the 12 hemocyte RNA-seq libraries (WSSV n=5,
PmSTAT dsRNA n=3, PmSTAT dsRNA + WSSV n=4). Replicates are UNBALANCED and
there is no uninfected control arm — see docs/DESIGN_LIMITS.md.

Each run is downloaded and converted independently, with a .done marker, so
the job can be killed by the 24 h cap and resumed indefinitely.

    python scripts/01_fetch_sra.py --list          # resolve runs only
    python scripts/01_fetch_sra.py --threads 8     # download + convert
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

ENA_FIELDS = (
    "run_accession,sample_accession,experiment_accession,sample_title,"
    "sample_alias,library_name,read_count,base_count,fastq_ftp,fastq_md5,"
    "fastq_bytes,scientific_name"
)


def resolve_runs(study: str, timeout: int = 60) -> list[dict]:
    """Resolve a study accession to its runs via the ENA portal API.

    ENA is used rather than NCBI eutils because it returns direct FASTQ URLs
    and MD5 checksums in one call, which makes verification trivial.
    """
    params = urllib.parse.urlencode({
        "accession": study,
        "result": "read_run",
        "fields": ENA_FIELDS,
        "format": "json",
        "limit": "0",
    })
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{params}"
    req = urllib.request.Request(url, headers={"User-Agent": "pmdark/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def have(tool: str) -> bool:
    return shutil.which(tool) is not None


def run(cmd: list[str], log=print) -> int:
    log("  $ " + " ".join(cmd))
    return subprocess.run(cmd).returncode


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    import hashlib
    h = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def download_url(url: str, dest: Path, log=print, retries: int = 3) -> bool:
    """Download with resume. Prefers aria2c, falls back to curl, then wget."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not url.startswith(("http://", "https://", "ftp://")):
        url = "https://" + url

    for attempt in range(1, retries + 1):
        if have("aria2c"):
            rc = run(["aria2c", "-c", "-x", "8", "-s", "8", "--console-log-level=warn",
                      "-d", str(dest.parent), "-o", dest.name, url], log)
        elif have("curl"):
            rc = run(["curl", "-fsSL", "-C", "-", "-o", str(dest), url], log)
        else:
            rc = run(["wget", "-c", "-q", "-O", str(dest), url], log)
        if rc == 0 and dest.exists() and dest.stat().st_size > 0:
            return True
        log(f"  attempt {attempt}/{retries} failed (rc={rc})")
        time.sleep(5 * attempt)
    return False


def fetch_via_ena(rec: dict, outdir: Path, log=print) -> bool:
    """Direct FASTQ download from ENA, with MD5 verification."""
    ftp = (rec.get("fastq_ftp") or "").split(";")
    md5s = (rec.get("fastq_md5") or "").split(";")
    ftp = [u for u in ftp if u]
    if not ftp:
        return False

    for i, url in enumerate(ftp):
        dest = outdir / Path(url).name
        expect = md5s[i] if i < len(md5s) else None

        if dest.exists() and expect:
            if md5sum(dest) == expect:
                log(f"  {dest.name} present and verified")
                continue
            log(f"  {dest.name} checksum mismatch — refetching")
            dest.unlink()

        if not download_url(url, dest, log):
            return False
        if expect and md5sum(dest) != expect:
            log(f"  {dest.name} FAILED md5 after download")
            return False
        log(f"  {dest.name} ok ({dest.stat().st_size / 1e9:.2f} GB)")
    return True


def fetch_via_sratools(acc: str, outdir: Path, threads: int, log=print) -> bool:
    """Fallback: NCBI sra-tools. Slower and needs scratch space for the .sra."""
    if not have("fasterq-dump"):
        log("  fasterq-dump not on PATH")
        return False
    tmp = outdir / f".tmp_{acc}"
    tmp.mkdir(parents=True, exist_ok=True)
    if have("prefetch"):
        run(["prefetch", "--max-size", "100G", "-O", str(tmp), acc], log)
    rc = run(["fasterq-dump", "--split-files", "--threads", str(threads),
              "--temp", str(tmp), "-O", str(outdir), acc], log)
    shutil.rmtree(tmp, ignore_errors=True)
    if rc != 0:
        return False
    for fq in sorted(outdir.glob(f"{acc}_*.fastq")):
        if have("pigz"):
            run(["pigz", "-p", str(threads), "-f", str(fq)], log)
        else:
            run(["gzip", "-f", str(fq)], log)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", default=None, help="override study accession")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--list", action="store_true", help="resolve and print runs, then exit")
    ap.add_argument("--force-sratools", action="store_true",
                    help="skip ENA, use prefetch/fasterq-dump")
    args = ap.parse_args()

    cfg = cfgmod.load()
    study = args.study or cfg["data"]["sra_study"]
    outdir = cfg.p("raw", "sra")
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"study  : {study}")
    print(f"outdir : {outdir}")

    meta_path = outdir / "runinfo.json"
    if meta_path.exists():
        runs = json.loads(meta_path.read_text())
        print(f"runs   : {len(runs)} (cached)")
    else:
        print("resolving via ENA portal API ...")
        runs = resolve_runs(study)
        if not runs:
            print(f"ERROR: no runs returned for {study}")
            return 1
        meta_path.write_text(json.dumps(runs, indent=2))
        print(f"runs   : {len(runs)} (wrote {meta_path.name})")

    total_gb = sum(
        sum(int(b) for b in (r.get("fastq_bytes") or "0").split(";") if b)
        for r in runs
    ) / 1e9
    total_reads = sum(int(r.get("read_count") or 0) for r in runs)

    print(f"\n{'run':<13}{'reads':>14}  {'GB':>7}  title")
    print("-" * 78)
    for r in runs:
        gb = sum(int(b) for b in (r.get("fastq_bytes") or "0").split(";") if b) / 1e9
        title = (r.get("sample_title") or r.get("library_name")
                 or r.get("sample_alias") or "")[:34]
        print(f"{r['run_accession']:<13}{int(r.get('read_count') or 0):>14,}  {gb:>7.2f}  {title}")
    print("-" * 78)
    print(f"{'TOTAL':<13}{total_reads:>14,}  {total_gb:>7.2f} GB")

    expected = cfg["data"]["n_libraries"]
    if len(runs) != expected:
        print(f"\nWARNING: expected {expected} libraries, resolved {len(runs)}. "
              "Check the accession and the group assignment before proceeding.")

    print("\nNOTE: sample titles above are the ONLY source for group assignment\n"
          "      (WSSV / PmSTAT_dsRNA / PmSTAT_dsRNA_WSSV). Record the mapping in\n"
          "      config/samples.tsv before running differential expression.")

    if args.list:
        return 0

    free_gb = shutil.disk_usage(outdir).free / 1e9
    if free_gb < total_gb * 1.3:
        print(f"\nERROR: need ~{total_gb * 1.3:.0f} GB, share has {free_gb:.0f} GB free")
        return 1

    print(f"\ndownloading ({free_gb:.0f} GB free) ...")
    failures: list[str] = []
    for i, rec in enumerate(runs, 1):
        acc = rec["run_accession"]
        done = outdir / f"{acc}.done"
        if done.exists():
            print(f"[{i}/{len(runs)}] {acc} already complete")
            continue

        print(f"[{i}/{len(runs)}] {acc}")
        ok = False
        if not args.force_sratools:
            ok = fetch_via_ena(rec, outdir, log=print)
        if not ok:
            print("  falling back to sra-tools")
            ok = fetch_via_sratools(acc, outdir, args.threads, log=print)

        if ok:
            done.write_text(json.dumps({"accession": acc,
                                        "reads": rec.get("read_count")}) + "\n")
        else:
            failures.append(acc)
            print(f"  FAILED: {acc}")

    if failures:
        print(f"\n{len(failures)} run(s) failed: {', '.join(failures)}")
        print("Re-run this script; completed runs are skipped.")
        return 1

    print(f"\nAll {len(runs)} runs complete in {outdir}")
    print("Next: scripts/02_fetch_refs.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
