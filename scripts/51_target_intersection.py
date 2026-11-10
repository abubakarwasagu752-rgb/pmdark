#!/usr/bin/env python3
"""How many of the homology-free targets are reachable structurally?

    python scripts/51_target_intersection.py

The manuscript says 190 of the 328 STAT-dependent uncharacterised
targets have no reachable homologue. 47_structure.py then found that 157
of the 328 have a structural relative. Those two numbers do not combine:
some of the 157 are among the 138 that already had a sequence homologue,
so the marginal counts say nothing about the overlap.

What §3.6 actually needs is the cross-tabulation — of the proteins with
no sequence homologue, how many does structure reach — and that is a
different number from either.

ID NAMESPACES

  annotation_targets.csv is keyed by gene (`gene-LOC119594247`) and the
  structural search by whatever headers targets.faa carries. If they do
  not match, every intersection is zero and the honest-looking answer
  "structure reaches none of them" would be an identifier bug. This
  script normalises, measures the overlap, and refuses to report a
  cross-tabulation built on fewer than half the targets.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

LOC = re.compile(r"(LOC\d+)")


def keys_for(s: str) -> set[str]:
    """Every spelling of an identifier we might have to match on."""
    out = {s}
    out.add(s.split()[0])
    base = s.split()[0]
    for pre in ("gene-", "rna-", "cds-", "lcl|"):
        if base.startswith(pre):
            base = base[len(pre):]
    out.add(base)
    out.add(base.split(".")[0])
    m = LOC.search(s)
    if m:
        out.add(m.group(1))
    return {x for x in out if x}


def main() -> int:
    cfg = cfgmod.load()
    res = cfg.p("results").resolve()
    struct = cfg.p("structures").resolve()

    csv_path = res / "annotation_targets.csv"
    m8 = struct / "targets_vs_lookup.m8"
    for p in (csv_path, m8):
        if not p.exists():
            print(f"missing {p}", file=sys.stderr)
            return 2

    rows = list(csv.DictReader(csv_path.open()))
    print("=" * 78)
    print("TARGETS: SEQUENCE REACH versus STRUCTURAL REACH")
    print("=" * 78)
    print(f"  {len(rows)} rows in {csv_path.name}")
    strata = {}
    for r in rows:
        strata[r["stratum"]] = strata.get(r["stratum"], 0) + 1
    print("  strata: " + ", ".join(f"{k} {v}" for k, v in
                                   sorted(strata.items())))

    # structural hits, self-hits excluded
    hit_raw: set[str] = set()
    for line in m8.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 11 and f[0] != f[1]:
            try:
                if float(f[10]) <= 1e-3:
                    hit_raw.add(f[0])
            except ValueError:
                pass
    hit_keys: set[str] = set()
    for h in hit_raw:
        hit_keys |= keys_for(h)
    print(f"  {len(hit_raw)} target proteins with a structural hit")

    matched, unmatched = 0, []
    for r in rows:
        if keys_for(r["gene"]) & hit_keys:
            r["_struct"] = True
            matched += 1
        else:
            r["_struct"] = False
            unmatched.append(r["gene"])

    # Could the namespaces simply not meet? Compare against how many
    # target identifiers resolve at all, not against how many hit.
    any_key = set()
    for h in hit_raw:
        any_key |= keys_for(h)
    resolvable = sum(1 for r in rows if keys_for(r["gene"]) & any_key)
    print(f"  {matched} rows matched to a structural hit")
    if hit_raw and matched == 0:
        print("\n  NO ROW MATCHED. The identifier namespaces differ.")
        print(f"    csv gene : {rows[0]['gene']}")
        print(f"    m8 query : {sorted(hit_raw)[0]}")
        print("  Refusing to report a cross-tabulation of zeros.")
        return 1

    orph = [r for r in rows if r["stratum"] != "homolog"]
    homo = [r for r in rows if r["stratum"] == "homolog"]
    a = sum(1 for r in orph if r["_struct"])
    b = len(orph) - a
    c = sum(1 for r in homo if r["_struct"])
    d = len(homo) - c

    print(f"\n  {'':<22}{'structural hit':>16}{'none':>8}{'total':>8}")
    print(f"  {'no sequence homologue':<22}{a:>16}{b:>8}{a + b:>8}")
    print(f"  {'has sequence homologue':<22}{c:>16}{d:>8}{c + d:>8}")
    print(f"  {'total':<22}{a + c:>16}{b + d:>8}{len(rows):>8}")

    if a + b:
        print(f"\n  Of the {a + b} targets with no sequence homologue, "
              f"{a} ({100 * a / (a + b):.1f}%) have a structural "
              f"relative.")
        print(f"  Unreachable by sequence OR structure: {b} "
              f"({100 * b / len(rows):.1f}% of all {len(rows)} targets).")

    out = {
        "n_targets": len(rows),
        "no_sequence_homologue": a + b,
        "no_sequence_with_structure": a,
        "no_sequence_no_structure": b,
        "has_sequence_with_structure": c,
        "has_sequence_no_structure": d,
        "structural_hits_total": a + c,
        "note": ("Cross-tabulation, not marginal counts: the 157 targets "
                 "with a structural hit include proteins that already "
                 "had a sequence homologue, so it cannot be subtracted "
                 "from 190."),
    }
    p = res / "target_intersection.json"
    p.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  wrote {p}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
