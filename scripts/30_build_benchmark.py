#!/usr/bin/env python3
"""Build the calibration benchmark — the credibility core of the paper.

THE QUESTION THIS ANSWERS

  When a protein has no annotated homolog, can a protein language model
  still recover its function, and how often is it right? Published PLM
  annotation pipelines report benchmark numbers from the literature. We
  measure it ourselves, on a split designed to mimic our actual problem.

THE DESIGN

  Query set: reviewed (Swiss-Prot) ARTHROPOD proteins carrying GO terms
  with experimental evidence codes. Experimental only — electronic
  annotations would make the ground truth circular, since they are
  themselves produced by homology transfer.

  Reference set: every other reviewed protein with experimental GO.

  Four tiers of increasing difficulty, each removing more of the query's
  relatives from the reference it is allowed to learn from:

    T0 full        everything except the query itself
    T1 no_species  drop the query's own species
    T2 no_order    drop the query's taxonomic order (e.g. all Decapoda)
    T3 no_phylum   drop Arthropoda entirely        <-- the real test

  T3 is the condition our shrimp dark proteins are actually in: no
  annotated relative anywhere near them. Whatever precision we measure at
  T3 is the precision we are entitled to claim for the dark proteome. If it
  collapses, that is the paper's headline finding and we report it.

  Lineage comes from UniProt's own `lineage` field, so tier membership is a
  substring test rather than a taxonomy-tree download.

    python scripts/30_build_benchmark.py --fetch     # download (~30-60 min)
    python scripts/30_build_benchmark.py --build     # make the splits
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

UNIPROT = "https://rest.uniprot.org"
GO_OBO = "https://current.geneontology.org/ontology/go-basic.obo"

FIELDS = "accession,id,organism_id,organism_name,lineage,go_id,sequence,length"
QUERY_EXP = "%28reviewed%3Atrue%29+AND+%28go_exp%3A*%29"

QUERY_PHYLUM = "Arthropoda"
MIN_LEN, MAX_LEN = 50, 1022     # MAX matches the ESM-2 positional limit
MIN_GO_TERMS = 1


def gzip_ok(path: Path) -> bool:
    """Is this a complete gzip stream? A truncated download decompresses
    partially and then raises EOFError, which previously got mistaken for
    a valid cached file on the next run."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 22):
                pass
        return True
    except (OSError, EOFError):
        return False


def fetch(url: str, dest: Path, timeout: int = 3600, retries: int = 4) -> bool:
    """Download with validation.

    UniProt's /stream endpoint drops long HTTP/2 transfers with
        curl: (92) HTTP/2 stream 1 was not closed cleanly: INTERNAL_ERROR
    which leaves a partial .gz behind. Forcing HTTP/1.1 avoids it; the
    integrity check catches whatever still slips through.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size > 0:
        if not str(dest).endswith(".gz") or gzip_ok(dest):
            print(f"  {dest.name} present ({dest.stat().st_size/1e6:.1f} MB)")
            return True
        print(f"  {dest.name} is TRUNCATED ({dest.stat().st_size/1e6:.1f} MB) "
              "— discarding and refetching")
        dest.unlink()

    for attempt in range(1, retries + 1):
        print(f"  fetching {dest.name} (attempt {attempt}/{retries})")
        rc = subprocess.run([
            "curl", "-fsSL",
            "--http1.1",                 # the actual fix for error 92
            "--max-time", str(timeout),
            "--retry", "3", "--retry-delay", "10", "--retry-all-errors",
            "--speed-time", "120", "--speed-limit", "1024",  # abort if stalled
            "-o", str(dest), url,
        ]).returncode

        if rc != 0 or not dest.exists() or dest.stat().st_size == 0:
            print(f"    curl rc={rc}")
            dest.unlink(missing_ok=True)
            continue

        if str(dest).endswith(".gz") and not gzip_ok(dest):
            print(f"    downloaded {dest.stat().st_size/1e6:.1f} MB "
                  "but gzip stream is incomplete — retrying")
            dest.unlink()
            continue

        print(f"    ok {dest.stat().st_size/1e6:.1f} MB")
        return True

    print(f"    FAILED after {retries} attempts")
    return False


def do_fetch(cfg) -> int:
    raw = cfg.p("raw", "benchmark")
    ok = True
    # One stream carries sequence, lineage and GO together, so no separate
    # FASTA parse and no taxonomy tree download is needed.
    ok &= fetch(
        f"{UNIPROT}/uniprotkb/stream?format=tsv&compressed=true"
        f"&fields={FIELDS}&query={QUERY_EXP}",
        raw / "sprot_exp_go.tsv.gz",
    )
    ok &= fetch(GO_OBO, raw / "go-basic.obo")
    return 0 if ok else 1


def load_records(path: Path) -> list[dict]:
    out = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ["Entry", "Organism (ID)", "Taxonomic lineage",
                "Gene Ontology IDs", "Sequence", "Length"]
        missing = [c for c in need if c not in idx]
        if missing:
            print(f"  column names differ; got: {header}")
            # UniProt field labels drift; fall back to positional order,
            # which follows the `fields` parameter we requested.
            idx = {k: i for i, k in enumerate(
                ["Entry", "Entry Name", "Organism (ID)", "Organism",
                 "Taxonomic lineage", "Gene Ontology IDs", "Sequence", "Length"]
            )}

        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(idx):
                continue
            go_raw = p[idx["Gene Ontology IDs"]].strip()
            if not go_raw:
                continue
            terms = {t.strip() for t in go_raw.split(";") if t.strip().startswith("GO:")}
            if len(terms) < MIN_GO_TERMS:
                continue
            seq = p[idx["Sequence"]].strip().upper()
            if not (MIN_LEN <= len(seq) <= MAX_LEN):
                continue
            out.append({
                "acc": p[idx["Entry"]],
                "taxid": p[idx["Organism (ID)"]],
                "lineage": p[idx["Taxonomic lineage"]],
                "go": terms,
                "seq": seq,
            })
    return out


def lineage_tokens(lineage: str) -> list[str]:
    # UniProt renders lineage as "cellular organisms, Eukaryota, ..., Decapoda"
    # sometimes with rank hints in parentheses.
    return [t.strip().split(" (")[0] for t in lineage.split(",") if t.strip()]


# Ranks above order that the -poda / -ptera suffix heuristic otherwise
# swallows. "Panarthropoda" ends in -poda and is not literally equal to
# "Arthropoda", so the first version of this function assigned 267 queries
# an "order" of Panarthropoda — which would have made T2 remove every
# arthropod and silently collapse T2 into T3.
_NOT_ORDERS = {
    "Panarthropoda", "Arthropoda", "Ecdysozoa", "Protostomia", "Bilateria",
    "Metazoa", "Eumetazoa", "Opisthokonta", "Eukaryota", "Mandibulata",
    "Pancrustacea", "Crustacea", "Hexapoda", "Chelicerata", "Myriapoda",
    "Multicrustacea", "Malacostraca", "Insecta", "Arachnida", "Branchiopoda",
    "Maxillopoda", "Collembola", "Entognatha", "Pterygota", "Neoptera",
    "Endopterygota", "Holometabola", "Dicondylia", "Eumalacostraca",
    "cellular organisms",
}

_KNOWN_ORDERS = {
    "Decapoda", "Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera",
    "Hemiptera", "Orthoptera", "Ixodida", "Araneae", "Isopoda", "Amphipoda",
    "Cladocera", "Blattodea", "Siphonaptera", "Scorpiones", "Mesostigmata",
    "Trombidiformes", "Sarcoptiformes", "Odonata", "Thysanoptera",
    "Phthiraptera", "Trichoptera", "Neuroptera", "Euphausiacea", "Stomatopoda",
}


def order_of(tokens: list[str]) -> str:
    """Best-effort taxonomic order.

    UniProt's plain lineage field carries no rank labels, so this is a
    heuristic: known orders first, then conventional suffixes, with
    above-order taxa explicitly rejected.
    """
    for t in tokens:
        if t in _KNOWN_ORDERS:
            return t
    for t in tokens:
        if t in _NOT_ORDERS:
            continue
        if t.endswith(("ptera", "poda", "formes", "odida", "ida", "acea")):
            return t
    return ""


def do_build(cfg, seed: int, max_query: int, max_ref: int) -> int:
    raw = cfg.p("raw", "benchmark")
    out = cfg.p("proteins", "benchmark")
    out.mkdir(parents=True, exist_ok=True)

    src = raw / "sprot_exp_go.tsv.gz"
    if not src.exists():
        print(f"missing {src} — run with --fetch first")
        return 1

    print(f"loading {src.name}")
    recs = load_records(src)
    print(f"  {len(recs):,} reviewed proteins with experimental GO, "
          f"{MIN_LEN}-{MAX_LEN} aa")

    for r in recs:
        toks = lineage_tokens(r["lineage"])
        r["tokens"] = set(toks)
        r["order"] = order_of(toks)
        r["is_query_phylum"] = QUERY_PHYLUM in r["tokens"]

    queries = [r for r in recs if r["is_query_phylum"]]
    refs_all = [r for r in recs if not r["is_query_phylum"]]
    print(f"  {QUERY_PHYLUM}: {len(queries):,}   other: {len(refs_all):,}")

    if len(queries) < 200:
        print(f"\nWARNING: only {len(queries)} {QUERY_PHYLUM} queries. "
              "Benchmark will be noisy; consider widening to Metazoa.")

    rng = random.Random(seed)
    rng.shuffle(queries)
    if max_query and len(queries) > max_query:
        queries = queries[:max_query]
        print(f"  subsampled queries to {len(queries):,}")

    # The reference pool is EVERY protein, arthropods included. It must be,
    # or the tiers are meaningless: T3 removes Arthropoda, which only tests
    # anything if Arthropoda were there to begin with. Each query's own
    # accession is excluded from its own neighbour search at scoring time
    # (31_calibrate.py), which is the only self-leak that matters.
    ref_pool = list(recs)
    if max_ref and len(ref_pool) > max_ref:
        rng.shuffle(ref_pool)
        # Guarantee every query stays in the pool so T0/T1/T2 have the
        # close relatives they are supposed to have.
        keep = {r["acc"] for r in ref_pool[:max_ref]}
        extra = [r for r in queries if r["acc"] not in keep]
        ref_pool = ref_pool[:max_ref] + extra
    n_arth = sum(1 for r in ref_pool if r["is_query_phylum"])
    print(f"  reference pool: {len(ref_pool):,} "
          f"({n_arth:,} {QUERY_PHYLUM}, {len(ref_pool)-n_arth:,} other)")

    # ---- tier masks ---------------------------------------------------
    # Stored as exclusion rules rather than materialised sets: four copies
    # of a 100k-protein reference would be wasteful, and the embedding is
    # computed once and reused across all tiers.
    tiers = {
        "T0_full":       {"exclude": "none",
                          "desc": "all references except the query itself"},
        "T1_no_species": {"exclude": "same_taxid",
                          "desc": "drop the query's own species"},
        "T2_no_order":   {"exclude": "same_order",
                          "desc": "drop the query's taxonomic order"},
        "T3_no_phylum":  {"exclude": "same_phylum",
                          "desc": f"drop all {QUERY_PHYLUM} — the orphan condition"},
    }

    def write_fasta(records, path: Path) -> None:
        with path.open("w") as fh:
            for r in records:
                fh.write(f">{r['acc']}\n")
                s = r["seq"]
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + "\n")

    write_fasta(queries, out / "query.faa")
    write_fasta(ref_pool, out / "reference.faa")
    print(f"  wrote query.faa ({len(queries):,}) and "
          f"reference.faa ({len(ref_pool):,})")

    meta = {
        "seed": seed,
        "query_phylum": QUERY_PHYLUM,
        "min_len": MIN_LEN, "max_len": MAX_LEN,
        "n_query": len(queries), "n_reference": len(ref_pool),
        "tiers": tiers,
        "query": {r["acc"]: {"taxid": r["taxid"], "order": r["order"],
                             "phylum": QUERY_PHYLUM if r["is_query_phylum"] else "",
                             "go": sorted(r["go"])} for r in queries},
        "reference": {r["acc"]: {"taxid": r["taxid"], "order": r["order"],
                                 "phylum": QUERY_PHYLUM if r["is_query_phylum"] else "",
                                 "go": sorted(r["go"])} for r in ref_pool},
    }
    meta_path = out / "benchmark.json.gz"
    with gzip.open(meta_path, "wt") as fh:
        json.dump(meta, fh)
    print(f"  wrote {meta_path.name} ({meta_path.stat().st_size/1e6:.1f} MB)")

    # ---- what the tiers actually cost ---------------------------------
    print("\n  how much reference each tier removes (median over queries):")
    sample = queries[: min(300, len(queries))]
    for name, rule in tiers.items():
        removed = []
        for q in sample:
            # Every tier also drops the query itself — counted here so the
            # percentages match what 31_calibrate.py actually does.
            if rule["exclude"] == "none":
                n = 1
            elif rule["exclude"] == "same_taxid":
                n = sum(1 for r in ref_pool if r["taxid"] == q["taxid"])
            elif rule["exclude"] == "same_order":
                n = sum(1 for r in ref_pool
                        if q["order"] and r["order"] == q["order"]) or 1
            else:
                n = sum(1 for r in ref_pool if r["is_query_phylum"])
            removed.append(n)
        removed.sort()
        med = removed[len(removed) // 2]
        print(f"    {name:<14} removes {med:>7,} of {len(ref_pool):,} "
              f"({100*med/len(ref_pool):5.1f}%)  {rule['desc']}")

    print("\n  top query orders:")
    for o, n in Counter(q["order"] or "(unresolved)" for q in queries).most_common(8):
        print(f"    {o:<18} {n:,}")

    print("\nNext:")
    print("  bash slurm/embed_fanout.sh esm2 "
          f"{out / 'reference.faa'}")
    print("  bash slurm/embed_fanout.sh esm2 "
          f"{out / 'query.faa'}")
    print("  python scripts/31_calibrate.py --model esm2")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--max-query", type=int, default=5000)
    ap.add_argument("--max-ref", type=int, default=60000)
    args = ap.parse_args()

    cfg = cfgmod.load()
    seed = args.seed if args.seed is not None else cfg["seed"]

    if args.fetch:
        rc = do_fetch(cfg)
        if rc or not args.build:
            return rc
    if args.build:
        return do_build(cfg, seed, args.max_query, args.max_ref)

    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
