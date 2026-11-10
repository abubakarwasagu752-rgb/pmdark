#!/usr/bin/env python3
"""Statistics for the structural result, and the confound that could sink it.

    python scripts/48_structure_stats.py

Read-only. Reads what 47_structure.py produced and writes
06_results/structural_stats.json plus the identifier list of the proteins
that survive structural hardening.

THE RESULT THIS HAS TO DEFEND

  47 found the dark set clustering at 2.9% against 30.0% for a
  length-matched background — structurally HETEROGENEOUS, the opposite of
  the hypothesis that went in. That is a strong claim and it rests on 3Di
  strings that ProstT5 predicted rather than measured.

THE CONFOUND

  ProstT5 is one of the three models used for annotation transfer. If it
  represents a protein badly, its embedding is poor AND its 3Di is poor,
  together. "Proteins ProstT5 handles badly" plausibly overlaps with
  "proteins nothing can annotate", which is the dark set exactly. A
  degenerate 3Di string — one that collapses toward a couple of coil
  states — aligns to nothing and clusters with nothing, for reasons that
  have no biology in them at all.

  So the headline comparison would look identical under two stories:

    biology   these proteins genuinely share no folds
    artefact  ProstT5 produced mush for these proteins specifically

  One piece of evidence already argues against pure mush: 18 of the 70
  found a structural relative in the lookup set, which a garbage string
  cannot do. But that is 18, and says nothing about the other 52.

THE TEST

  Look at the strings. A degenerate prediction has low per-residue
  entropy and is dominated by one or two characters. If dark and
  background 3Di have comparable entropy, the clustering difference is
  not explained by prediction quality and the result stands. If dark 3Di
  is visibly degenerate, the coherence claim is withdrawn and only the
  hardening is reported.

  The verdict is computed, printed, and written to the JSON. It is not
  left to whoever reads the numbers later.
"""

from __future__ import annotations

import json
import math
import random
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

# Entropy below this, in bits per residue, is a prediction that has
# collapsed rather than a protein that is genuinely simple. Chosen as
# roughly half the entropy of a uniform draw over the 20-letter 3Di
# alphabet (log2(20) = 4.32), i.e. the string carries less than half the
# information its alphabet allows.
DEGENERATE_BITS = 2.16


def log(m: str) -> None:
    print(m, flush=True)


def find_foldseek(cfg) -> str | None:
    p = shutil.which("foldseek")
    if p:
        return p
    work = cfg.p("work")
    for c in (work / "foldseek" / "bin" / "foldseek",
              Path.home() / "pmdark_work" / "foldseek" / "bin" / "foldseek"):
        if c.is_file():
            return str(c)
    return None


def read_mmseqs_db(base: Path) -> dict[str, str]:
    """Read an mmseqs/foldseek flat database without calling foldseek.

    The GPU createdb writes SHARDED databases — `db_ss.0 .. db_ss.7`
    with `db_ss.index.0 ..` — so a check for a single `db_ss` file finds
    nothing even though the data is there. Going through
    `foldseek convert2fasta` would also have to know about the shards.
    Reading the format directly sidesteps both: an index line is
    `key<TAB>offset<TAB>length`, and the record at that offset is
    terminated by NUL.

    Returns {header name: sequence}; empty dict if nothing is readable.
    """
    def shards(stem: Path) -> list[tuple[Path, Path]]:
        out = []
        if stem.is_file() and Path(f"{stem}.index").is_file():
            out.append((stem, Path(f"{stem}.index")))
        i = 0
        while True:
            d, ix = Path(f"{stem}.{i}"), Path(f"{stem}.index.{i}")
            if not (d.is_file() and ix.is_file()):
                break
            if d.stat().st_size > 0:
                out.append((d, ix))
            i += 1
        return out

    def load(stem: Path) -> dict[int, str]:
        rec: dict[int, str] = {}
        for data, index in shards(stem):
            raw = data.read_bytes()
            for line in index.read_text().splitlines():
                f = line.split("\t")
                if len(f) < 3:
                    continue
                try:
                    k, off, ln = int(f[0]), int(f[1]), int(f[2])
                except ValueError:
                    continue
                chunk = raw[off:off + ln]
                rec[k] = chunk.split(b"\x00")[0].decode(
                    "utf-8", "replace").strip()
        return rec

    seqs = load(base)
    heads = load(Path(f"{base}_h")) if Path(f"{base}_h").is_file() or \
        Path(f"{base}_h.0").is_file() else {}
    out: dict[str, str] = {}
    for k, v in seqs.items():
        name = heads.get(k, str(k)).split()[0] if heads else str(k)
        if v:
            out[name] = v
    return out


def read_fasta(path: Path) -> dict[str, str]:
    out, cur, buf = {}, None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    out[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf)
    return out


def entropy_bits(s: str) -> float:
    if not s:
        return 0.0
    n = len(s)
    return -sum((c / n) * math.log2(c / n) for c in Counter(s).values())


def dominant_frac(s: str) -> float:
    return (Counter(s).most_common(1)[0][1] / len(s)) if s else 0.0


def mannwhitney(a: list[float], b: list[float]):
    """U test with a normal approximation and tie correction."""
    if not a or not b:
        return None, None
    merged = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks, i, n = [0.0] * len(merged), 0, len(merged)
    while i < n:
        j = i
        while j + 1 < n and merged[j + 1][0] == merged[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    r1 = sum(r for r, (_, g) in zip(ranks, merged) if g == 0)
    n1, n2 = len(a), len(b)
    u1 = r1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    tie = sum(t ** 3 - t for t in Counter(v for v, _ in merged).values())
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - tie / (n * (n - 1))))
    if sd == 0:
        return u1, None
    z = (u1 - mu) / sd
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u1, p


def cliffs_delta(a: list[float], b: list[float]) -> float | None:
    if not a or not b:
        return None
    gt = lt = 0
    bs = sorted(b)
    import bisect
    for x in a:
        # bisect_left counts b strictly BELOW x, i.e. a > b.
        # These two were transposed, which inverted the sign: the
        # artifact case reported delta +1.00 for a dark set whose
        # entropy was far LOWER than the background's.
        gt += bisect.bisect_left(bs, x)
        lt += len(bs) - bisect.bisect_right(bs, x)
    return (gt - lt) / (len(a) * len(b))


def fisher(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact on a 2x2, by enumeration."""
    def lchoose(n, k):
        return (math.lgamma(n + 1) - math.lgamma(k + 1)
                - math.lgamma(n - k + 1))

    r1, r2 = a + b, c + d
    c1, nn = a + c, a + b + c + d

    def prob(x):
        return math.exp(lchoose(r1, x) + lchoose(r2, c1 - x)
                        - lchoose(nn, c1))

    p0 = prob(a) * (1 + 1e-9)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                        if prob(x) <= p0))


def boot_diff(k1: int, n1: int, k2: int, n2: int, seed: int, R: int = 10000):
    """Percentile interval on the difference of two proportions."""
    rng = random.Random(seed)
    v1 = [1] * k1 + [0] * (n1 - k1)
    v2 = [1] * k2 + [0] * (n2 - k2)
    d = []
    for _ in range(R):
        s1 = sum(rng.choice(v1) for _ in range(n1)) / n1
        s2 = sum(rng.choice(v2) for _ in range(n2)) / n2
        d.append(s1 - s2)
    d.sort()
    return d[int(0.025 * R)], d[int(0.975 * R)]


def wilson(k: int, n: int):
    if n == 0:
        return None, None
    z, p = 1.959964, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def clustered(tsv: Path, want: set[str]) -> tuple[int, int, int]:
    """(proteins in multi-member clusters, n clusters, largest)."""
    from collections import defaultdict
    mem: dict[str, set[str]] = defaultdict(set)
    if tsv.exists():
        for line in tsv.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 2:
                mem[f[0]].add(f[1])
    sizes = [len(v & want) for v in mem.values()]
    return sum(s for s in sizes if s > 1), len(mem), max(sizes, default=0)


def main() -> int:
    cfg = cfgmod.load()
    res = cfg.p("results").resolve()
    out = cfg.p("structures").resolve()
    seed = cfg.raw.get("seed", 1)
    fs = find_foldseek(cfg)
    rep: dict = {"seed": seed}

    print("=" * 78)
    print("STRUCTURAL STATISTICS")
    print("=" * 78)

    # ---- 3Di strings ----------------------------------------------------
    tdi: dict[str, dict[str, str]] = {}
    for name in ("hardened", "background", "targets"):
        fa = out / f"{name}_3di.fasta"
        if fa.exists():
            tdi[name] = read_fasta(fa)
            continue
        # read the sharded database directly; the base db carries the
        # headers, the _ss db carries the 3Di
        got = read_mmseqs_db(out / f"{name}db_ss")
        if got:
            heads = read_mmseqs_db(out / f"{name}db")
            if len(heads) == len(got):
                got = {h: s for h, s in zip(heads.keys(), got.values())}
            tdi[name] = got
            fa.write_text("".join(f">{k}\n{v}\n" for k, v in got.items()))
            log(f"  extracted {len(got)} 3Di strings -> {fa.name}")
            continue
        if fs:
            r = subprocess.run([fs, "convert2fasta",
                                str(out / f"{name}db_ss"), str(fa)],
                               capture_output=True, text=True)
            if fa.exists():
                tdi[name] = read_fasta(fa)
            elif r.stderr:
                log(f"  convert2fasta {name}: {r.stderr.strip()[:200]}")

    if "hardened" not in tdi or "background" not in tdi:
        print("\n  looked for, under " + str(out) + ":")
        for name in ("hardened", "background"):
            found = sorted(q.name for q in out.glob(f"{name}db*"))
            print(f"    {name}db* -> " + (", ".join(found[:8]) or "nothing"))
        print(f"    foldseek: {fs or 'not found'}")
        print("\n  Could not extract 3Di strings. Needs foldseek and the")
        print("  databases 47 built. Without them the confound cannot be")
        print("  tested, and the coherence result must not be reported.")
        return 1

    print(f"\n  3Di strings: " + ", ".join(
        f"{k} {len(v)}" for k, v in tdi.items()))

    stats = {}
    for name, seqs in tdi.items():
        ent = [entropy_bits(s) for s in seqs.values() if s]
        dom = [dominant_frac(s) for s in seqs.values() if s]
        ln = [len(s) for s in seqs.values() if s]
        if not ent:
            continue
        ent_s, dom_s = sorted(ent), sorted(dom)
        stats[name] = {
            "n": len(ent),
            "entropy_median": ent_s[len(ent_s) // 2],
            "entropy_mean": sum(ent) / len(ent),
            "entropy_q1": ent_s[len(ent_s) // 4],
            "entropy_q3": ent_s[3 * len(ent_s) // 4],
            "dominant_median": dom_s[len(dom_s) // 2],
            "length_median": sorted(ln)[len(ln) // 2],
            "degenerate_n": sum(1 for e in ent if e < DEGENERATE_BITS),
            "_ent": ent,
        }

    print(f"\n  {'set':<12}{'n':>5}{'3Di entropy (bits)':>22}"
          f"{'dominant char':>15}{'len':>7}{'degenerate':>12}")
    for name, s in stats.items():
        print(f"  {name:<12}{s['n']:>5}"
              f"{s['entropy_median']:>13.2f} [{s['entropy_q1']:.2f}"
              f"-{s['entropy_q3']:.2f}]"
              f"{100 * s['dominant_median']:>14.1f}%"
              f"{s['length_median']:>7}"
              f"{s['degenerate_n']:>7}/{s['n']}")

    a_e = stats["hardened"]["_ent"]
    b_e = stats["background"]["_ent"]
    _, p_ent = mannwhitney(a_e, b_e)
    d_ent = cliffs_delta(a_e, b_e)
    gap = stats["hardened"]["entropy_median"] - \
        stats["background"]["entropy_median"]
    print(f"\n  dark minus background entropy: {gap:+.3f} bits"
          f"   Cliff's delta {d_ent:+.3f}"
          f"   p = {'<0.001' if p_ent is not None and p_ent < 1e-3 else f'{p_ent:.3f}' if p_ent is not None else 'n/a'}")

    rep["three_di"] = {k: {kk: vv for kk, vv in v.items() if kk != "_ent"}
                       for k, v in stats.items()}
    rep["three_di_test"] = {"entropy_gap_bits": gap, "cliffs_delta": d_ent,
                            "p": p_ent, "degenerate_threshold_bits":
                            DEGENERATE_BITS}

    # ---- coherence ------------------------------------------------------
    ids = {n: set(read_fasta(out / f"{n}.fasta")) for n in
           ("hardened", "background", "targets")
           if (out / f"{n}.fasta").exists()}
    coh = {}
    for name in ("hardened", "background", "targets"):
        tsv = out / f"clu_{name}.tsv"
        if name not in ids or not tsv.exists():
            continue
        k, nclu, big = clustered(tsv, ids[name])
        n = len(ids[name])
        lo, hi = wilson(k, n)
        coh[name] = {"n": n, "clustered": k, "fraction": k / n,
                     "ci": [lo, hi], "clusters": nclu, "largest": big}

    verdict = "not evaluated"
    if "hardened" in coh and "background" in coh:
        A, B = coh["hardened"], coh["background"]
        a, c = A["clustered"], B["clustered"]
        b, d = A["n"] - a, B["n"] - c
        p = fisher(a, b, c, d)
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        se = math.sqrt(1 / (a + .5) + 1 / (b + .5) + 1 / (c + .5) + 1 / (d + .5))
        lo_d, hi_d = boot_diff(a, A["n"], c, B["n"], seed)
        print(f"\n  COHERENCE")
        print(f"    dark       {a}/{A['n']} = {100 * A['fraction']:.1f}% "
              f"[{100 * A['ci'][0]:.1f}, {100 * A['ci'][1]:.1f}]  "
              f"largest cluster {A['largest']}")
        print(f"    background {c}/{B['n']} = {100 * B['fraction']:.1f}% "
              f"[{100 * B['ci'][0]:.1f}, {100 * B['ci'][1]:.1f}]  "
              f"largest cluster {B['largest']}")
        print(f"    difference {100 * (A['fraction'] - B['fraction']):+.1f} pp"
              f"  bootstrap 95% CI [{100 * lo_d:+.1f}, {100 * hi_d:+.1f}]")
        print(f"    odds ratio {orr:.3f} "
              f"[{math.exp(math.log(orr) - 1.96 * se):.3f}, "
              f"{math.exp(math.log(orr) + 1.96 * se):.3f}]"
              f"   Fisher p = {p:.3g}")
        rep["coherence"] = {
            "dark": A, "background": B, "fisher_p": p,
            "odds_ratio": orr,
            "or_ci": [math.exp(math.log(orr) - 1.96 * se),
                      math.exp(math.log(orr) + 1.96 * se)],
            "difference_pp": 100 * (A["fraction"] - B["fraction"]),
            "difference_ci_pp": [100 * lo_d, 100 * hi_d],
        }

        # ---- the verdict, computed ------------------------------------
        deg = stats["hardened"]["degenerate_n"] / stats["hardened"]["n"]
        confounded = (gap < -0.25) or (deg > 0.25) or (
            d_ent is not None and d_ent < -0.33)
        sig = hi_d < 0 or lo_d > 0
        if not sig:
            verdict = ("no difference: the interval on the difference "
                       "includes zero")
        elif confounded:
            verdict = ("CONFOUNDED: dark 3Di is measurably poorer than "
                       "background 3Di, so the clustering difference may "
                       "be prediction quality rather than biology. Report "
                       "the hardening; withhold the coherence claim.")
        else:
            verdict = ("supported: the dark set clusters less than a "
                       "length-matched background, and its 3Di is not "
                       "measurably poorer, so prediction quality does not "
                       "explain the difference")

    # ---- hardening ------------------------------------------------------
    m8 = out / "hardened_vs_lookup.m8"
    if m8.exists() and "hardened" in ids:
        hit = set()
        for line in m8.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 11 and f[0] != f[1]:
                try:
                    if float(f[10]) <= 1e-3:
                        hit.add(f[0])
                except ValueError:
                    pass
        hit &= ids["hardened"]
        survivors = sorted(ids["hardened"] - hit)
        n = len(ids["hardened"])
        lo, hi = wilson(len(hit), n)
        print(f"\n  HARDENING")
        print(f"    {len(hit)}/{n} = {100 * len(hit) / n:.1f}% "
              f"[{100 * lo:.1f}, {100 * hi:.1f}] have a structural relative")
        print(f"    orphan set: 94 -> {n} -> {len(survivors)}")
        rep["hardening"] = {"n": n, "with_structural_hit": len(hit),
                            "fraction": len(hit) / n, "ci": [lo, hi],
                            "surviving": len(survivors)}
        p_ids = res / "structure_hardened_ids.txt"
        p_ids.write_text("\n".join(survivors) + "\n")
        print(f"    wrote {p_ids.name} ({len(survivors)} identifiers) — "
              f"Fmax must be recomputed on these")

    rep["verdict"] = verdict
    print(f"\n  VERDICT  {verdict}")
    print("\n  Caveat recorded in the JSON: 3Di is predicted by ProstT5, "
          "which\n  is also one of the three annotation models, so this "
          "check is not\n  independent of that arm. It is independent of "
          "the ESM-2 result.")
    rep["caveat"] = ("3Di predicted by ProstT5, not measured. ProstT5 is "
                     "also one of the three annotation-transfer models, so "
                     "this analysis is not independent of that arm; it is "
                     "independent of the ESM-2 headline.")
    (res / "structural_stats.json").write_text(
        json.dumps(rep, indent=2) + "\n")
    print(f"\n  wrote {res / 'structural_stats.json'}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
