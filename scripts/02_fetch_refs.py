#!/usr/bin/env python3
"""Acquire references via ENA + UniProt only.

NCBI is UNREACHABLE from the CPE cluster (verified 2026-11-07):
    api.ncbi.nlm.nih.gov   TIMEOUT
    ftp.ncbi.nlm.nih.gov   TIMEOUT
while these answer 200:
    www.ebi.ac.uk          (ENA portal + browser API)
    rest.uniprot.org
    huggingface.co

ENA mirrors every INSDC assembly, so GCA_* accessions are all retrievable
there. This script therefore never touches NCBI, and the `datasets` CLI is
not used.

    python scripts/02_fetch_refs.py --list
    python scripts/02_fetch_refs.py --all
    python scripts/02_fetch_refs.py --only genome,wssv
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

TAXID_PMONODON = 6687

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api"
ENA_BROWSER = "https://www.ebi.ac.uk/ena/browser/api"
UNIPROT = "https://rest.uniprot.org"

# Complete WSSV genomes, all INSDC so all mirrored by ENA.
# The shrimp in this study are Thai, so WSSV-TH is the primary reference;
# the other two are kept for ORFeome comparison and to avoid isolate bias
# in the host-virus screen.
WSSV_ISOLATES = {
    "AF440570": "WSSV-TH (Thailand, 292,967 bp)   <-- primary",
    "AF332093": "WSSV-CN (China, 305,107 bp; = RefSeq NC_003225)",
    "AF369029": "WSSV-TW (Taiwan, 307,287 bp)",
}


def fetch_bytes(url: str, timeout: int = 120) -> bytes | None:
    req = urllib.request.Request(url, headers={"User-Agent": "pmdark/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception as exc:
        print(f"    {type(exc).__name__}: {exc}")
        return None


def fetch_to(url: str, dest: Path, timeout: int = 900) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"    {dest.name} already present ({dest.stat().st_size/1e6:.1f} MB)")
        return True
    if shutil.which("curl"):
        rc = subprocess.run(
            ["curl", "-fsSL", "-C", "-", "--max-time", str(timeout), "-o", str(dest), url]
        ).returncode
    else:
        rc = subprocess.run(["wget", "-c", "-q", "-O", str(dest), url]).returncode
    ok = rc == 0 and dest.exists() and dest.stat().st_size > 0
    if ok:
        print(f"    {dest.name}  {dest.stat().st_size/1e6:.1f} MB")
    else:
        dest.unlink(missing_ok=True)
    return ok


# ----------------------------------------------------------------------
# genome, via ENA
# ----------------------------------------------------------------------
def list_assemblies() -> list[dict]:
    fields = ("accession,assembly_name,assembly_level,assembly_title,"
              "genome_representation,base_count,scientific_name,"
              "study_accession,sample_accession,version,wgs_set")
    q = urllib.parse.urlencode({
        "result": "assembly",
        "query": f"tax_eq({TAXID_PMONODON})",
        "fields": fields,
        "format": "json",
        "limit": "0",
    })
    data = fetch_bytes(f"{ENA_PORTAL}/search?{q}")
    if not data:
        return []
    try:
        rows = json.loads(data.decode())
    except ValueError:
        print("    could not parse ENA assembly response")
        return []

    def rank(r: dict) -> tuple:
        lvl = (r.get("assembly_level") or "").lower()
        order = {"chromosome": 0, "complete genome": 0, "scaffold": 1, "contig": 2}
        return (order.get(lvl, 3), -int(r.get("base_count") or 0))

    return sorted(rows, key=rank)


def download_genome(acc: str, wgs_set: str, outdir: Path) -> bool:
    """Try each ENA route in turn; report which one worked."""
    outdir.mkdir(parents=True, exist_ok=True)
    marker = outdir / f"{acc}.done"
    if marker.exists():
        print(f"    {acc} already downloaded")
        return True

    dest = outdir / f"{acc}.fasta.gz"
    routes = [
        (f"{ENA_BROWSER}/fasta/{acc}?download=true&gzip=true",
         "browser API, assembly accession"),
    ]
    if wgs_set:
        routes.append(
            (f"{ENA_BROWSER}/fasta/{wgs_set}?download=true&gzip=true",
             f"browser API, WGS set {wgs_set}")
        )
    routes.append(
        (f"{ENA_PORTAL}/files?result=assembly&accession={acc}&format=tsv",
         "portal file listing (manual follow-up)")
    )

    for url, label in routes[:-1]:
        print(f"    trying: {label}")
        if fetch_to(url, dest):
            # Guard against a 20-byte error page saved as a 'success'
            if dest.stat().st_size < 1_000_000:
                print(f"    suspiciously small ({dest.stat().st_size} B) — discarding")
                dest.unlink(missing_ok=True)
                continue
            marker.write_text(f"{acc}\nvia: {label}\n")
            return True

    print("\n    All direct routes failed. Fall back to the portal listing:")
    print(f"      curl -s '{routes[-1][0]}'")
    print("    then fetch the submitted FASTA URL it names.")
    return False


# ----------------------------------------------------------------------
# WSSV, via ENA
# ----------------------------------------------------------------------
def download_wssv(outdir: Path) -> bool:
    outdir.mkdir(parents=True, exist_ok=True)
    ok_any = False

    for acc, desc in WSSV_ISOLATES.items():
        print(f"  {acc}  {desc}")
        fa = outdir / f"{acc}.fasta"
        embl = outdir / f"{acc}.embl"
        fetch_to(f"{ENA_BROWSER}/fasta/{acc}?download=true", fa, timeout=180)
        fetch_to(f"{ENA_BROWSER}/embl/{acc}?download=true", embl, timeout=180)

        # Extract CDS translations -> the viral ORFeome for the co-folding screen
        faa = outdir / f"{acc}_proteins.faa"
        if embl.exists() and not faa.exists():
            try:
                from Bio import SeqIO
                n = 0
                with faa.open("w") as out:
                    for rec in SeqIO.parse(str(embl), "embl"):
                        for feat in rec.features:
                            if feat.type != "CDS":
                                continue
                            tr = feat.qualifiers.get("translation", [None])[0]
                            if not tr:
                                continue
                            pid = (feat.qualifiers.get("protein_id")
                                   or feat.qualifiers.get("locus_tag")
                                   or [f"{acc}_cds_{n}"])[0]
                            prod = feat.qualifiers.get("product", [""])[0]
                            out.write(f">{pid} {prod}\n{tr}\n")
                            n += 1
                print(f"    {faa.name}: {n} proteins")
                ok_any = ok_any or n > 0
            except Exception as exc:
                print(f"    protein extraction failed: {exc}")
    return ok_any


# ----------------------------------------------------------------------
# UniProt
# ----------------------------------------------------------------------
def download_swissprot(outdir: Path) -> bool:
    outdir.mkdir(parents=True, exist_ok=True)
    ok = True
    # Swiss-Prot via the UniProt REST stream (ftp.uniprot.org is a separate
    # host and was not tested; rest.uniprot.org returned 200).
    ok &= fetch_to(
        f"{UNIPROT}/uniprotkb/stream?format=fasta&compressed=true"
        "&query=%28reviewed%3Atrue%29",
        outdir / "uniprot_sprot.fasta.gz", timeout=1800,
    )
    # Experimentally-evidenced GO annotations — the ground truth for the
    # calibration benchmark in stage 04.
    ok &= fetch_to(
        f"{UNIPROT}/uniprotkb/stream?format=tsv&compressed=true"
        "&fields=accession,id,go_id,organism_id,protein_name,length"
        "&query=%28reviewed%3Atrue%29+AND+%28go_exp%3A*%29",
        outdir / "goa_experimental.tsv.gz", timeout=1800,
    )
    return ok


def download_proteomes(outdir: Path) -> bool:
    """Targeted panel for phylostratigraphy. Full reference proteomes are
    ~100 GB; this panel is a few hundred MB and spans the strata we need."""
    outdir.mkdir(parents=True, exist_ok=True)
    panel = [
        (6687, "Penaeus_monodon"),
        (6689, "Penaeus_vannamei"),
        (6693, "Penaeus_japonicus"),
        (6690, "Penaeus_chinensis"),
        (6706, "Procambarus_clarkii"),
        (6665, "Daphnia_pulex"),
        (7227, "Drosophila_melanogaster"),
        (7070, "Tribolium_castaneum"),
        (7460, "Apis_mellifera"),
        (6945, "Ixodes_scapularis"),
        (6239, "Caenorhabditis_elegans"),
        (9606, "Homo_sapiens"),
    ]
    ok = True
    for taxid, label in panel:
        dest = outdir / f"{label}_{taxid}.fasta.gz"
        url = (f"{UNIPROT}/uniprotkb/stream?format=fasta&compressed=true"
               f"&query=%28taxonomy_id%3A{taxid}%29")
        print(f"  {label} ({taxid})")
        if not fetch_to(url, dest, timeout=1200):
            print(f"    WARNING: {label} not retrieved")
            ok = False
    return ok


# ----------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--accession", default=None)
    args = ap.parse_args()

    cfg = cfgmod.load()
    raw = cfg.p("raw")

    targets = ({"genome", "wssv", "swissprot", "proteomes"} if args.all
               else {t.strip() for t in args.only.split(",") if t.strip()})

    if args.list or "genome" in targets:
        print("Penaeus monodon assemblies in ENA (NCBI is unreachable here)\n")
        rows = list_assemblies()
        if not rows:
            print("  ENA returned nothing. Check connectivity:")
            print(f"    curl -s '{ENA_PORTAL}/search?result=assembly"
                  f"&query=tax_eq({TAXID_PMONODON})&fields=accession,assembly_name"
                  "&format=tsv&limit=0'")
            return 1

        print(f"{'accession':<20}{'level':<16}{'Gb':>6}  {'wgs':<8} name")
        print("-" * 78)
        for r in rows[:12]:
            bc = int(r.get("base_count") or 0) / 1e9
            print(f"{r.get('accession',''):<20}{(r.get('assembly_level') or ''):<16}"
                  f"{bc:>6.2f}  {(r.get('wgs_set') or '-'):<8} "
                  f"{r.get('assembly_name','')}")
        print("-" * 78)

        best = rows[0]
        print(f"\nselected: {best.get('accession')}  "
              f"({best.get('assembly_name')}, {best.get('assembly_level')})")

        if args.list:
            print("\nRe-run with --only genome to download.")
            return 0

        acc = args.accession or best.get("accession")
        print(f"\ndownloading genome {acc}")
        if not download_genome(acc, best.get("wgs_set") or "", raw / "genome"):
            return 1

    if "wssv" in targets:
        print("\nWSSV genomes (ENA)")
        download_wssv(raw / "wssv")

    if "swissprot" in targets:
        print("\nUniProt Swiss-Prot + experimental GOA")
        download_swissprot(raw / "uniprot")

    if "proteomes" in targets:
        print("\nPhylostratigraphy proteome panel")
        download_proteomes(raw / "proteomes")

    if not targets and not args.list:
        ap.print_help()
        return 1

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
