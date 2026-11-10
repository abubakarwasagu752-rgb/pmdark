#!/usr/bin/env python3
"""Is the dark set's low 3Di entropy a model failure or a real property?

    python scripts/49_disorder_matched.py                 # plan
    python scripts/49_disorder_matched.py --go --stage pool    # GPU
    python scripts/49_disorder_matched.py --go --stage match   # GPU
    python scripts/49_disorder_matched.py --go --stage disorder

WHAT 48 LEFT UNRESOLVED

  The dark set clusters at 2.9% against 30.0% for a length-matched
  background, but its ProstT5 3Di is measurably poorer: 2.40 bits
  against 3.05, delta -0.395, 31 of 70 below the degeneracy threshold.
  So the clustering gap has three possible explanations and 48 could not
  separate them:

    artefact   ProstT5 fails on proteins it represents poorly. The
               strings are mush and mush clusters with nothing.
    disorder   These proteins are genuinely disordered. There is no fold
               to predict and none to share, so low 3Di entropy is a
               CORRECT prediction and the low clustering is its
               consequence. This is a finding, not a failure.
    diversity  Real, well-formed, mutually unrelated folds.

  The first means withdraw the claim. The second and third mean report
  it, but they are different papers.

TWO TESTS, AND WHY THE SECOND NEEDS NO DEPENDENCY

  1. ENTROPY-MATCHED BACKGROUND. Match the background on 3Di entropy as
     well as length, then recluster. If the dark set still clusters less
     when compared against proteins whose predictions are equally poor,
     prediction quality does not explain the gap. This is the same move
     that rescued the STAT analysis, where matching on expression decile
     separated a real effect from a power artefact.

  2. AMINO-ACID SEQUENCE ENTROPY. ProstT5 cannot influence the
     complexity of the input sequence. If dark proteins have genuinely
     lower-complexity sequences, then low 3Di entropy is the right
     answer and the disorder explanation is supported. This needs no
     model, no download and no licence, which is exactly why it is the
     primary discriminator rather than a secondary one.

     Two further independent signals ride along: the Uversky
     charge-hydropathy criterion, which is arithmetic over Kyte-Doolittle
     values and a net charge, and metapredict if it happens to be
     importable. Neither is a T5.
"""

from __future__ import annotations

import argparse
import bisect
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

# Kyte & Doolittle (1982) hydropathy.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
CALIPER = 0.2      # SD units, as in 40_matched_enrichment
MAX_SMD = 0.1      # matched sets must agree this closely


def log(m: str) -> None:
    print(m, flush=True)


def find_foldseek(cfg) -> str | None:
    p = shutil.which("foldseek")
    if p:
        return p
    for c in (cfg.p("work") / "foldseek" / "bin" / "foldseek",
              Path.home() / "pmdark_work" / "foldseek" / "bin" / "foldseek"):
        if c.is_file():
            return str(c)
    return None


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


def read_mmseqs_db(base: Path) -> dict[int, str]:
    """Sharded mmseqs flat database -> {key: record}."""
    def shards(stem: Path):
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

    rec: dict[int, str] = {}
    for data, index in shards(base):
        raw = data.read_bytes()
        for line in index.read_text().splitlines():
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                k, off, ln = int(f[0]), int(f[1]), int(f[2])
            except ValueError:
                continue
            rec[k] = raw[off:off + ln].split(b"\x00")[0].decode(
                "utf-8", "replace").strip()
    return rec


def db_to_dict(base: Path) -> dict[str, str]:
    """{protein name: 3Di string}.

    Names come from `<db>_h`, NOT from `<db>`. In a foldseek database the
    base file holds amino-acid sequences and the `_h` file holds the
    headers; reading the base gave every protein a key that was its own
    sequence, and the first thing that tried to look one up by name died
    with a KeyError whose text was a 128-residue protein.
    """
    seqs = read_mmseqs_db(Path(f"{base}_ss"))
    heads = read_mmseqs_db(Path(f"{base}_h"))
    out: dict[str, str] = {}
    for k, v in seqs.items():
        if not v:
            continue
        nm = heads.get(k)
        if not nm:
            continue
        out[nm.split()[0]] = v
    return out


def entropy_bits(s: str) -> float:
    if not s:
        return 0.0
    n = len(s)
    return -sum((c / n) * math.log2(c / n) for c in Counter(s).values())


def uversky(seq: str) -> tuple[float, float, bool]:
    """Mean net charge, normalised mean hydropathy, disordered?

    Uversky's boundary R = 2.785 H - 1.151; a protein above the line has
    too much charge for its hydrophobicity to fold. Pure arithmetic over
    the sequence, so nothing a language model did can reach it.
    """
    s = [c for c in seq.upper() if c in KD]
    if not s:
        return 0.0, 0.0, False
    h = sum(KD[c] for c in s) / len(s)
    hn = (h + 4.5) / 9.0
    pos = sum(1 for c in s if c in "KR")
    neg = sum(1 for c in s if c in "DE")
    r = abs(pos - neg) / len(s)
    return r, hn, r > (2.785 * hn - 1.151)


def mannwhitney(a, b):
    if not a or not b:
        return None
    merged = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    n = len(merged)
    ranks, i = [0.0] * n, 0
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
        return None
    z = (u1 - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def cliffs_delta(a, b):
    if not a or not b:
        return None
    bs = sorted(b)
    gt = sum(bisect.bisect_left(bs, x) for x in a)
    lt = sum(len(bs) - bisect.bisect_right(bs, x) for x in a)
    return (gt - lt) / (len(a) * len(b))


def smd(a, b):
    if len(a) < 2 or len(b) < 2:
        return None
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    sp = math.sqrt((va + vb) / 2)
    return (ma - mb) / sp if sp else None


def fisher(a, b, c, d):
    def lc(n, k):
        return (math.lgamma(n + 1) - math.lgamma(k + 1)
                - math.lgamma(n - k + 1))
    r1, r2, c1, nn = a + b, c + d, a + c, a + b + c + d

    def pr(x):
        return math.exp(lc(r1, x) + lc(r2, c1 - x) - lc(nn, c1))
    p0 = pr(a) * (1 + 1e-9)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0))


def wilson(k, n):
    if not n:
        return None, None
    z, p = 1.959964, k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def match_two(dark: dict[str, tuple[float, float]],
              pool: dict[str, tuple[float, float]],
              seed: int) -> tuple[list[str], dict]:
    """Nearest-neighbour match on (length, 3Di entropy), standardised,
    without replacement, inside a caliper. Greedy in a shuffled order so
    the result does not depend on identifier sort order."""
    vals = list(dark.values()) + list(pool.values())
    out = {}
    for i in (0, 1):
        xs = [v[i] for v in vals]
        m = sum(xs) / len(xs)
        sd = math.sqrt(sum((x - m) ** 2 for x in xs) / max(1, len(xs) - 1))
        out[i] = (m, sd or 1.0)

    def z(v):
        return tuple((v[i] - out[i][0]) / out[i][1] for i in (0, 1))

    dz = {k: z(v) for k, v in dark.items()}
    pz = {k: z(v) for k, v in pool.items()}
    rng = random.Random(seed)
    order = list(dz)
    rng.shuffle(order)
    free = set(pz)
    picked, unmatched = [], 0
    for q in order:
        a = dz[q]
        best, bd = None, None
        for c in free:
            b = pz[c]
            d = math.dist(a, b)
            if bd is None or d < bd:
                best, bd = c, d
        if best is None or bd > CALIPER * math.sqrt(2):
            unmatched += 1
            continue
        picked.append(best)
        free.discard(best)
    return picked, {"unmatched": unmatched, "caliper": CALIPER}


def run(cmd, label=""):
    log(f"  $ {label or Path(str(cmd[0])).name} {' '.join(str(c) for c in cmd[1:4])} ...")
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write((r.stderr or r.stdout)[-3000:] + "\n")
    return r.returncode


def clustered(tsv: Path, want: set[str]):
    from collections import defaultdict
    mem = defaultdict(set)
    if tsv.exists():
        for line in tsv.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 2:
                mem[f[0]].add(f[1])
    sizes = [len(v & want) for v in mem.values()]
    return sum(s for s in sizes if s > 1), max(sizes, default=0)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--go", action="store_true")
    ap.add_argument("--stage", default="all",
                    choices=["pool", "match", "disorder", "all"])
    ap.add_argument("--gpu", type=int, default=1)
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results").resolve()
    out = cfg.p("structures").resolve()
    prot = cfg.p("proteins").resolve()
    seed = cfg.raw.get("seed", 1)
    fs = find_foldseek(cfg)
    weights = out / "prostt5_weights"

    qfa = prot / "benchmark" / "query.faa"
    ofa = prot / "benchmark" / "orphan_check" / "orphans.faa"
    hfa = out / "hardened.fasta"
    # Merged, not replaced. Each stage fills only its own block, and the
    # file is written whole at the end — so a `--stage disorder` run on
    # its own used to silently discard a `matched` block that a previous
    # `--stage match` had computed, and the loss showed up only as a
    # missing group in a figure much later. Anything this run recomputes
    # still overwrites its own block.
    rep: dict = {"seed": seed}
    _prev = res / "disorder_matched.json"
    if _prev.exists():
        try:
            old = json.loads(_prev.read_text())
            if isinstance(old, dict):
                rep.update({k: v for k, v in old.items() if k != "seed"})
                keep = [k for k in ("sequence", "matched") if k in old]
                if keep:
                    log(f"merging into existing {_prev.name} "
                        f"(keeping {', '.join(keep)} unless recomputed)")
        except (ValueError, OSError) as e:
            log(f"could not read existing {_prev.name} ({e}); "
                f"starting a fresh report")

    print("=" * 78)
    print("DISORDER AND ENTROPY-MATCHED COHERENCE")
    print("=" * 78)
    for lbl, p in (("foldseek", fs), ("query", qfa), ("orphans", ofa),
                   ("hardened", hfa), ("weights", weights)):
        ok = bool(p) and (p is fs or Path(p).exists())
        print(f"  {lbl:<10}{'ok   ' if ok else 'MISSING '}{p}")
    if not (fs and qfa.exists() and ofa.exists() and hfa.exists()):
        print("\n  cannot run")
        return 1

    allq = read_fasta(qfa)
    orphan = set(read_fasta(ofa))
    dark_ids = set(read_fasta(hfa))
    pool_ids = [i for i in allq if i not in orphan]
    print(f"\n  dark {len(dark_ids)}   pool {len(pool_ids)}   "
          f"query {len(allq)}")

    if not args.go:
        print("\nPLAN ONLY — nothing written. Re-run with --go.")
        return 0

    # ---- stage pool: 3Di for the whole background pool ------------------
    pool_fa = out / "pool.fasta"
    if args.stage in ("pool", "all"):
        pool_fa.write_text("".join(f">{i}\n{allq[i]}\n" for i in pool_ids))
        log(f"wrote pool.fasta ({len(pool_ids)} proteins)")
        if not (out / "pooldb_ss.0").exists() and not (out / "pooldb_ss").exists():
            cmd = [fs, "createdb", pool_fa, out / "pooldb",
                   "--prostt5-model", weights]
            if args.gpu:
                cmd += ["--gpu", "1"]
            if run(cmd, "foldseek createdb") != 0:
                return 1
        log("pool 3Di database ready")

    # ---- features -------------------------------------------------------
    dark_3di = db_to_dict(out / "hardeneddb")
    pool_3di = db_to_dict(out / "pooldb")
    if not pool_3di:
        print("\n  no pool 3Di — run --stage pool first")
        return 1
    # A name that is 40+ residues of pure amino-acid alphabet is a
    # sequence that leaked into the key, not an identifier. Catch it
    # here rather than at the first lookup that fails.
    def looks_like_sequence(x: str) -> bool:
        return len(x) > 40 and set(x.upper()) <= set("ACDEFGHIKLMNPQRSTVWYXBZU")
    for lbl, dd in (("dark", dark_3di), ("pool", pool_3di)):
        bad = [k for k in list(dd)[:50] if looks_like_sequence(k)]
        if bad:
            print(f"\n  {lbl} database keys are sequences, not names — "
                  f"headers were read from the wrong file. Refusing.")
            return 1
    missing = [i for i in list(dark_3di) + list(pool_3di) if i not in allq]
    if missing:
        print(f"\n  {len(missing)} identifiers are absent from "
              f"{qfa.name}; first: {missing[0][:60]}")
        return 1
    print(f"\n  3Di strings: dark {len(dark_3di)}, pool {len(pool_3di)}")

    feat_dark = {i: (float(len(allq[i])), entropy_bits(dark_3di[i]))
                 for i in dark_3di if i in allq}
    feat_pool = {i: (float(len(allq[i])), entropy_bits(pool_3di[i]))
                 for i in pool_3di if i in allq}

    # ---- stage disorder: independent of any language model --------------
    if args.stage in ("disorder", "all"):
        print("\n  SEQUENCE PROPERTIES (no model involved)")
        rows = {}
        for lbl, ids in (("dark", list(dark_ids)), ("pool", pool_ids)):
            seqs = [allq[i] for i in ids if i in allq]
            aa_ent = [entropy_bits(s) for s in seqs]
            uv = [uversky(s) for s in seqs]
            rows[lbl] = {
                "n": len(seqs),
                "aa_entropy": aa_ent,
                "disorder_frac": sum(1 for _, _, d in uv if d) / len(uv),
                "net_charge": [r for r, _, _ in uv],
                "hydropathy": [h for _, h, _ in uv],
            }
        for lbl in ("dark", "pool"):
            e = sorted(rows[lbl]["aa_entropy"])
            h = sorted(rows[lbl]["hydropathy"])
            print(f"    {lbl:<6}n={rows[lbl]['n']:<5}"
                  f"AA entropy {e[len(e) // 2]:.3f}  "
                  f"hydropathy {h[len(h) // 2]:.3f}  "
                  f"Uversky-disordered "
                  f"{100 * rows[lbl]['disorder_frac']:.1f}%")
        p_aa = mannwhitney(rows["dark"]["aa_entropy"],
                           rows["pool"]["aa_entropy"])
        d_aa = cliffs_delta(rows["dark"]["aa_entropy"],
                            rows["pool"]["aa_entropy"])
        # Denominators come from rows[], which counted only the proteins
        # whose sequence was actually found. The earlier version divided a
        # numerator over `i in allq` by the full len(dark_ids); identical
        # while every identifier resolves, and silently wrong the moment
        # one does not.
        n_dark, n_pool = rows["dark"]["n"], rows["pool"]["n"]
        a_d = round(rows["dark"]["disorder_frac"] * n_dark)
        p_d = round(rows["pool"]["disorder_frac"] * n_pool)
        p_uv = fisher(a_d, n_dark - a_d, p_d, n_pool - p_d)
        # Odds ratio with a Woolf interval, Haldane-Anscombe corrected.
        # Reported here rather than worked out downstream: a figure or a
        # manuscript cannot cite a number that exists in neither.
        _a, _b, _c, _d2 = a_d + .5, n_dark - a_d + .5, \
            p_d + .5, n_pool - p_d + .5
        uv_or = (_a * _d2) / (_b * _c)
        uv_se = math.sqrt(1 / _a + 1 / _b + 1 / _c + 1 / _d2)
        print(f"    AA entropy  delta {d_aa:+.3f}  "
              f"p {'<0.001' if p_aa is not None and p_aa < 1e-3 else f'{p_aa:.3f}'}")
        print(f"    Uversky     {a_d}/{n_dark} = "
              f"{100 * a_d / n_dark:.1f}%  vs  {p_d}/{n_pool} = "
              f"{100 * p_d / n_pool:.1f}%")
        print(f"                OR {uv_or:.2f} "
              f"[{math.exp(math.log(uv_or) - 1.96 * uv_se):.2f}, "
              f"{math.exp(math.log(uv_or) + 1.96 * uv_se):.2f}]"
              f"   Fisher p = {p_uv:.3g}")
        meta = None
        try:
            import metapredict as mp
            dd = [float(sum(x > 0.5 for x in mp.predict_disorder(allq[i]))
                        / max(1, len(allq[i]))) for i in dark_ids
                  if i in allq]
            pp = [float(sum(x > 0.5 for x in mp.predict_disorder(allq[i]))
                        / max(1, len(allq[i])))
                  for i in pool_ids[:300]]
            meta = {"dark_median": sorted(dd)[len(dd) // 2],
                    "pool_median": sorted(pp)[len(pp) // 2],
                    "delta": cliffs_delta(dd, pp),
                    "p": mannwhitney(dd, pp)}
            print(f"    metapredict dark {meta['dark_median']:.3f} vs "
                  f"pool {meta['pool_median']:.3f}  "
                  f"delta {meta['delta']:+.3f}")
        except ImportError:
            print("    metapredict not installed (optional): "
                  "pip install metapredict")
        rep["sequence"] = {
            "aa_entropy_delta": d_aa, "aa_entropy_p": p_aa,
            "uversky_dark": a_d / n_dark,
            "uversky_pool": p_d / n_pool, "uversky_p": p_uv,
            "uversky_dark_n": a_d, "uversky_dark_total": n_dark,
            "uversky_pool_n": p_d, "uversky_pool_total": n_pool,
            "uversky_or": uv_or,
            "uversky_or_ci": [math.exp(math.log(uv_or) - 1.96 * uv_se),
                              math.exp(math.log(uv_or) + 1.96 * uv_se)],
            "uversky_boundary": "R > 2.785*H - 1.151 (Kyte-Doolittle, "
                                "normalised 0-1)",
            "metapredict": meta,
            "dark_aa_entropy_median": sorted(
                rows["dark"]["aa_entropy"])[rows["dark"]["n"] // 2],
            "pool_aa_entropy_median": sorted(
                rows["pool"]["aa_entropy"])[rows["pool"]["n"] // 2],
            # The points themselves, so the charge-hydropathy panel is
            # drawn from the same arithmetic that produced the fraction
            # above and cannot disagree with it.
            "charge_hydropathy": {
                lbl: [[round(h, 4), round(r, 4)] for h, r in
                      zip(rows[lbl]["hydropathy"], rows[lbl]["net_charge"])]
                for lbl in ("dark", "pool")},
        }

    # ---- stage match: entropy-matched background ------------------------
    if args.stage in ("match", "all"):
        picked, info = match_two(feat_dark, feat_pool, seed)
        print(f"\n  ENTROPY-MATCHED BACKGROUND")
        print(f"    matched {len(picked)}/{len(feat_dark)}  "
              f"({info['unmatched']} outside the caliper)")
        if len(picked) < 0.8 * len(feat_dark):
            print("    REFUSING: fewer than 80% matched. The pool does "
                  "not contain enough proteins with comparably poor 3Di, "
                  "which is itself evidence that the dark set is unusual.")
            rep["matched"] = {"matched": len(picked), "refused": True}
        else:
            # Kept, not just printed. These two numbers are the entire
            # evidence that the matching worked; a reviewer asked to
            # accept a matched comparison will want them in the paper,
            # and a value that lives only in a terminal scrollback cannot
            # be put there.
            balance = {}
            for i, nm, key in ((0, "length", "length"),
                               (1, "3Di entropy", "three_di_entropy")):
                s = smd([feat_dark[k][i] for k in feat_dark],
                        [feat_pool[k][i] for k in picked])
                balance[key] = s
                flag = "" if abs(s) <= MAX_SMD else "   <-- above 0.1"
                print(f"    SMD {nm:<12}{s:+.4f}{flag}")
            mfa = out / "matched.fasta"
            mfa.write_text("".join(f">{i}\n{allq[i]}\n" for i in picked))
            if run([fs, "createdb", mfa, out / "matcheddb",
                    "--prostt5-model", weights] +
                   (["--gpu", "1"] if args.gpu else []),
                   "foldseek createdb") != 0:
                return 1
            # foldseek refuses to overwrite a cluster database and exits
            # with "clu_matched.dbtype exists already!", which on a rerun
            # killed the stage AFTER the matching had been recomputed and
            # before anything was written. These are this script's own
            # outputs of this exact step, regenerated in seconds from
            # matcheddb, so a rerun clears them rather than stopping.
            stale = sorted(out.glob("clu_matched*"))
            if stale:
                log(f"clearing {len(stale)} file(s) from a previous "
                    f"clustering of the matched set")
                for f_ in stale:
                    f_.unlink()
            t = out / "tmp_clu_matched"
            if t.exists():
                shutil.rmtree(t, ignore_errors=True)
            t.mkdir(exist_ok=True)
            if run([fs, "cluster", out / "matcheddb", out / "clu_matched",
                    t, "-c", "0.8", "-e", "0.001",
                    "--alignment-type", "2"], "foldseek cluster") != 0:
                return 1
            if run([fs, "createtsv", out / "matcheddb", out / "matcheddb",
                    out / "clu_matched", out / "clu_matched.tsv"],
                   "foldseek createtsv") != 0:
                return 1
            km, bigm = clustered(out / "clu_matched.tsv", set(picked))
            kd, bigd = clustered(out / "clu_hardened.tsv", dark_ids)
            nd, nm_ = len(dark_ids), len(picked)
            p = fisher(kd, nd - kd, km, nm_ - km)
            lo1, hi1 = wilson(kd, nd)
            lo2, hi2 = wilson(km, nm_)
            print(f"    dark            {kd}/{nd} = {100 * kd / nd:.1f}% "
                  f"[{100 * lo1:.1f}, {100 * hi1:.1f}]")
            print(f"    matched control {km}/{nm_} = "
                  f"{100 * km / nm_:.1f}% [{100 * lo2:.1f}, {100 * hi2:.1f}]")
            print(f"    Fisher p = {p:.3g}")
            _a, _b = kd + .5, nd - kd + .5
            _c, _d2 = km + .5, nm_ - km + .5
            m_or = (_a * _d2) / (_b * _c)
            m_se = math.sqrt(1 / _a + 1 / _b + 1 / _c + 1 / _d2)
            print(f"    OR {m_or:.3f} "
                  f"[{math.exp(math.log(m_or) - 1.96 * m_se):.3f}, "
                  f"{math.exp(math.log(m_or) + 1.96 * m_se):.3f}]")
            rep["matched"] = {
                "n_matched": nm_, "dark_clustered": kd,
                "matched_clustered": km,
                "dark_n": nd,
                "dark_fraction": kd / nd, "matched_fraction": km / nm_,
                "dark_ci": [lo1, hi1], "matched_ci": [lo2, hi2],
                "fisher_p": p, "refused": False,
                "odds_ratio": m_or,
                "or_ci": [math.exp(math.log(m_or) - 1.96 * m_se),
                          math.exp(math.log(m_or) + 1.96 * m_se)],
                "difference_pp": 100 * (kd / nd - km / nm_),
                "largest_dark": bigd, "largest_matched": bigm,
                "smd": balance, "caliper_sd": CALIPER, "max_smd": MAX_SMD,
                "unmatched": info["unmatched"],
                "matched_ids": sorted(picked),
            }

    # ---- verdict --------------------------------------------------------
    v = []
    m = rep.get("matched", {})
    s = rep.get("sequence", {})
    if m and not m.get("refused"):
        if m["fisher_p"] < 0.05 and m["dark_fraction"] < m["matched_fraction"]:
            v.append("the clustering deficit SURVIVES matching on 3Di "
                     "entropy, so prediction quality does not explain it")
        elif m["fisher_p"] >= 0.05:
            v.append("the clustering deficit DISAPPEARS once 3Di entropy "
                     "is matched, so it was prediction quality")
    if s:
        if s.get("aa_entropy_delta") is not None and \
                s["aa_entropy_delta"] < -0.2:
            v.append("dark amino-acid sequences are genuinely lower in "
                     "complexity, which no language model can have "
                     "caused, so the low 3Di entropy is a correct "
                     "prediction rather than a failure")
        if s.get("uversky_dark", 0) > s.get("uversky_pool", 0) + 0.1:
            v.append("the dark set is enriched for predicted disorder by "
                     "the Uversky criterion, which is arithmetic over the "
                     "sequence and independent of ProstT5")
    rep["verdict"] = v or ["inconclusive"]
    print("\n  VERDICT")
    for line in rep["verdict"]:
        print(f"    - {line}")
    (res / "disorder_matched.json").write_text(json.dumps(
        rep, indent=2, default=lambda o: None) + "\n")
    print(f"\n  wrote {res / 'disorder_matched.json'}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
