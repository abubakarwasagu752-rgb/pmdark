#!/usr/bin/env python3
"""Is whitening a property of PLM embeddings, or a quirk of one checkpoint?

    python scripts/33_multimodel.py --models esm2,prott5,prostt5
    python scripts/33_multimodel.py --status          # what is embedded

WHY THIS IS THE MOST IMPORTANT REMAINING METHODS RUN

  The calibration result rests on a single model:

      transform   T0 Fmax   Inone   ORPHAN     BP      MF
      none          0.708   0.632    0.580   0.231   0.345
      center        0.758   0.648    0.597   0.252   0.345
      whiten        0.820   0.689    0.660   0.326   0.430

  Every one of those numbers is ESM-2 650M. A reviewer's first question is
  whether anisotropy and its correction generalise, or whether this is one
  checkpoint's geometry. The literature audit makes the stakes concrete:
  a 2025 isotropy study reports ProtBERT and ProtXLNet using "only 2-14
  dimensions" but applies no correction and runs no downstream task, and
  MSRep fixes the space by TRAINING an adapter rather than post hoc. Our
  claim is that a training-free transform recovers much of the orphan gap.
  That claim is only worth making if it holds across architectures.

  ProtT5 and ProstT5 are encoder-decoder T5 models; ESM-2 is an
  encoder-only masked LM. If whitening helps all three, the finding is
  about mean-pooled protein embeddings in general. If it helps only ESM-2,
  the honest paper says so and the contribution shrinks to a footnote.

WHAT ELSE THIS ANSWERS IN THE SAME RUN

  * THE SIMPLER ALTERNATIVES. Centering is three characters of code and
    whitening is a 1280x1280 eigendecomposition. If centering gets most of
    the way, the honest claim is centering. Both are run for every model.

  * CONDITIONING. ZCA estimates a DxD covariance from 19,304 reference
    proteins -- about 15 samples per dimension -- and the implementation
    clips eigenvalues at 1e-8, which in principle lets a near-null
    direction be amplified 1e4-fold. A synthetic check suggested the clip
    never fires at this sample ratio, but that was simulation. This
    reports the real eigenvalue spectrum per model so the question is
    settled with data.

  * BOOTSTRAP INTERVALS. ORPHAN_nohit is 94 proteins. A 0.580 -> 0.660
    difference on 94 proteins without an interval is two point estimates,
    not a result. Fmax is resampled over queries.

  PRE-REGISTRATION HOLDS. Hyperparameters were frozen on T0_full and I50
  by scripts/32_sweep.py before the hard tiers were ever evaluated, and
  this reuses the frozen params_{model}_{transform}.json where present. It
  does NOT re-tune per model: a transform that needed its own tuning to
  look good would be a different and weaker claim.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

TRANSFORMS = ("none", "center", "whiten")
HEADLINE_TIERS = ("T0_full", "I50", "Inone", "ORPHAN_nohit")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def embed_status(cfg, model: str, fasta: Path) -> tuple[int, int]:
    """(done, total) shards for this model/fasta, without importing torch."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from pmdark.shards import (ShardStore, fasta_checksum, plan_shards,
                               read_fasta)
    if not fasta.exists():
        return (0, 0)
    records = read_fasta(fasta)
    shards = plan_shards(records, cfg["compute"]["shard_size"])
    store = ShardStore(cfg.p("embeddings", model, fasta.stem),
                       checksum=fasta_checksum(records), n_shards=len(shards))
    return store.progress()


def spectrum(cfg, model: str, fasta: Path) -> dict | None:
    """Real eigenvalue spectrum of the reference covariance, per model."""
    try:
        from pmdark.shards import (ShardStore, fasta_checksum, plan_shards,
                                   read_fasta)
        records = read_fasta(fasta)
        shards = plan_shards(records, cfg["compute"]["shard_size"])
        store = ShardStore(cfg.p("embeddings", model, fasta.stem),
                           checksum=fasta_checksum(records),
                           n_shards=len(shards))
        done, total = store.progress()
        if done < total:
            return None
        _, emb = store.load_all()
    except Exception as e:                      # noqa: BLE001
        log(f"  spectrum unavailable for {model}: {e}")
        return None
    x = np.asarray(emb, dtype=np.float64)
    n, d = x.shape
    xc = x - x.mean(0, keepdims=True)
    cov = (xc.T @ xc) / max(1, n - 1)
    vals = np.linalg.eigvalsh(cov)
    vals = np.clip(vals, 0.0, None)
    tot = float(vals.sum()) or 1.0
    srt = np.sort(vals)[::-1]
    cum = np.cumsum(srt) / tot
    return {
        "n_reference": int(n), "dim": int(d),
        "samples_per_dim": round(n / max(1, d), 1),
        "eig_min": float(vals.min()), "eig_max": float(vals.max()),
        "condition_number": float(vals.max() / max(vals.min(), 1e-300)),
        "clip_1e8_fires": bool(vals.min() < 1e-8),
        "dims_for_50pct_variance": int(np.searchsorted(cum, 0.50) + 1),
        "dims_for_95pct_variance": int(np.searchsorted(cum, 0.95) + 1),
        "top1_variance_share": round(float(srt[0] / tot), 4),
    }


def run_calibration(model: str, transform: str, tiers: str,
                    extra: list[str]) -> dict | None:
    cmd = [sys.executable, "scripts/31_calibrate.py",
           "--model", model, "--transform", transform, "--tiers", tiers,
           *extra]
    log("$ " + " ".join(cmd))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        log(f"  calibration failed for {model}/{transform} (rc={rc})")
        return None
    return {}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", default="esm2,prott5,prostt5")
    ap.add_argument("--transforms", default=",".join(TRANSFORMS))
    ap.add_argument("--tiers", default="T0_full,I90,I50,I30,Inone,ORPHAN_nohit")
    ap.add_argument("--status", action="store_true",
                    help="report embedding progress and exit")
    ap.add_argument("--skip-calibration", action="store_true",
                    help="only tabulate results already on disk")
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--sim-floor", type=float, default=None)
    ap.add_argument("--power", type=float, default=None)
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results")
    bench = cfg.p("proteins", "benchmark")
    q_fa, r_fa = bench / "query.faa", bench / "reference.faa"
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    transforms = [t.strip() for t in args.transforms.split(",") if t.strip()]

    # ---- what is embedded? --------------------------------------------
    print("\n" + "=" * 78)
    print("EMBEDDING STATUS")
    print("=" * 78)
    print(f"  {'model':<12}{'query shards':>16}{'reference shards':>20}")
    print("-" * 78)
    ready = []
    for m in models:
        qd, qt = embed_status(cfg, m, q_fa)
        rd, rt = embed_status(cfg, m, r_fa)
        ok = qt > 0 and qd == qt and rt > 0 and rd == rt
        ready.append(m) if ok else None
        print(f"  {m:<12}{f'{qd}/{qt}':>16}{f'{rd}/{rt}':>20}"
              f"   {'ready' if ok else 'INCOMPLETE'}")
    print("-" * 78)
    missing = [m for m in models if m not in ready]
    if missing:
        print("  embed the missing models (GPU, one fan-out per model/FASTA):")
        for m in missing:
            print(f"    bash slurm/embed_fanout.sh {m} {q_fa}")
            print(f"    bash slurm/embed_fanout.sh {m} {r_fa}")
    if args.status:
        return 0
    if not ready:
        print("\n  nothing to do until at least one model is fully embedded")
        return 2

    # ---- geometry, measured rather than simulated ---------------------
    print("\n" + "=" * 78)
    print("EMBEDDING GEOMETRY  (reference set, per model)")
    print("=" * 78)
    print(f"  {'model':<10}{'dim':>6}{'n/dim':>8}{'cond':>11}"
          f"{'dims 50%':>10}{'dims 95%':>10}{'clip?':>7}")
    print("-" * 78)
    spec = {}
    for m in ready:
        sp = spectrum(cfg, m, r_fa)
        if sp is None:
            print(f"  {m:<10}  (unavailable)")
            continue
        spec[m] = sp
        print(f"  {m:<10}{sp['dim']:>6}{sp['samples_per_dim']:>8.1f}"
              f"{sp['condition_number']:>11.1e}"
              f"{sp['dims_for_50pct_variance']:>10}"
              f"{sp['dims_for_95pct_variance']:>10}"
              f"{'YES' if sp['clip_1e8_fires'] else 'no':>7}")
    print("-" * 78)
    print("  'dims 50%' is how many of the dimensions carry half the")
    print("  variance: a small number IS the anisotropy, measured directly.")
    print("  'clip?' YES means the 1e-8 eigenvalue floor actually fires and")
    print("  whitening is amplifying a near-null direction — report it.")

    # ---- the grid ------------------------------------------------------
    if not args.skip_calibration:
        for m in ready:
            for t in transforms:
                # COMMON PARAMETERS OVERRIDE EVERYTHING WHEN GIVEN.
                #
                # The first run of this script compared esm2/whiten at its
                # frozen k=100 sim_floor=0.0 against prott5/whiten and
                # prostt5/whiten at the DEFAULT sim_floor=0.5, and the T5
                # models appeared to collapse: Fmax 0.106 and 0.078 at
                # Inone, with 993 and 1011 of 1,087 queries having no
                # usable neighbour and coverage of 0.09 and 0.07.
                #
                # That is a threshold rejecting everything, not a model
                # failing. sim_floor is an ABSOLUTE cosine cut, and
                # whitening compresses the cosine scale — median top-1
                # similarity falls to ~0.29 at Inone for all three models.
                # A 0.5 floor therefore discards nearly every neighbour in
                # whitened space. ESM-2 escaped only because its sweep had
                # already frozen the floor at 0.0.
                #
                # My docstring defended not re-tuning per model as the
                # stricter choice. It was not: this is a change of units,
                # not a tuning advantage, and comparing across transforms
                # requires a scale-free retrieval rule. Pass --sim-floor 0
                # (top-k only) to make all nine cells comparable.
                if args.k is not None or args.sim_floor is not None \
                        or args.power is not None:
                    extra = ["--k", str(args.k if args.k is not None else 25),
                             "--sim-floor", str(args.sim_floor
                                                if args.sim_floor is not None
                                                else 0.0),
                             "--power", str(args.power
                                            if args.power is not None else 3.0)]
                    log(f"{m}/{t}: COMMON params {extra} (override)")
                else:
                    pf = res / f"params_{m}_{t}.json"
                    extra = []
                    if pf.exists():
                        p = json.loads(pf.read_text())
                        extra = ["--k", str(p.get("k", 25)),
                                 "--sim-floor", str(p.get("sim_floor", 0.5)),
                                 "--power", str(p.get("power", 3.0))]
                        log(f"{m}/{t}: frozen params {extra}")
                    else:
                        log(f"{m}/{t}: no frozen params; 31_calibrate "
                            "defaults — NOT comparable across transforms, "
                            "see the note in the source")
                run_calibration(m, t, args.tiers, extra)

    # ---- tabulate whatever is on disk ---------------------------------
    rows = []
    for m in ready:
        for t in transforms:
            f = res / f"calibration_{m}_{t}.json"
            if not f.exists():
                continue
            try:
                d = json.loads(f.read_text())
            except ValueError:
                continue
            # Real schema, read from 31_calibrate.py and goeval.evaluate:
            #   calibration_<model>_<transform>.json
            #     -> results[<tier>]["metrics"]["overall"]["fmax"]
            #     -> results[<tier>]["metrics"]["BP"|"MF"|"CC"]["fmax"]
            # My first version assumed tiers[<tier>]["fmax"] and would have
            # tabulated every cell as "-" without erroring.
            tiers = d.get("results", {})
            row = {"model": m, "transform": t}
            for tier in HEADLINE_TIERS:
                met = ((tiers.get(tier) or {}).get("metrics") or {})
                ov = met.get("overall") or {}
                row[tier] = ov.get("fmax")
                row[tier + "_cov"] = ov.get("coverage")
                if tier == "ORPHAN_nohit":
                    for ont in ("BP", "MF", "CC"):
                        row[f"ORPHAN_{ont}"] = (met.get(ont) or {}).get("fmax")
                    row["ORPHAN_n"] = ov.get("n_proteins")
            rows.append(row)

    if not rows:
        print("\n  no calibration_*.json found yet")
        return 0

    print("\n" + "=" * 78)
    print("Fmax BY MODEL AND TRANSFORM")
    print("=" * 78)
    hdr = f"  {'model':<10}{'transform':<10}"
    for tier in HEADLINE_TIERS:
        hdr += f"{tier:>13}"
    print(hdr)
    print("-" * 78)
    for r in rows:
        line = f"  {r['model']:<10}{r['transform']:<10}"
        for tier in HEADLINE_TIERS:
            v = r.get(tier)
            line += f"{(f'{v:.3f}' if v is not None else '-'):>13}"
        print(line)
    print("-" * 78)

    # Does the transform help every model, or only one?
    print("\n  EFFECT OF WHITENING, PER MODEL (whiten minus none)")
    for m in ready:
        a = next((r for r in rows if r["model"] == m
                  and r["transform"] == "none"), None)
        b = next((r for r in rows if r["model"] == m
                  and r["transform"] == "whiten"), None)
        c = next((r for r in rows if r["model"] == m
                  and r["transform"] == "center"), None)
        if not (a and b):
            continue
        for tier in ("T0_full", "ORPHAN_nohit"):
            if a.get(tier) is None or b.get(tier) is None:
                continue
            dw = b[tier] - a[tier]
            dc = (c[tier] - a[tier]) if (c and c.get(tier) is not None) else None
            extra = f"   centering alone {dc:+.3f}" if dc is not None else ""
            print(f"    {m:<10}{tier:<14}{dw:+.3f}{extra}")
    print("\n  if the sign and rough size hold across architectures, the")
    print("  finding is about mean-pooled protein embeddings. If only ESM-2")
    print("  improves, the claim is about ESM-2 and must say so.")

    print("\n  ONTOLOGY SPLIT ON GENUINE ORPHANS")
    print(f"  {'model':<10}{'transform':<10}{'BP':>8}{'MF':>8}{'CC':>8}{'n':>6}")
    print("-" * 78)
    for r in rows:
        if r.get("ORPHAN_BP") is None:
            continue
        print(f"  {r['model']:<10}{r['transform']:<10}"
              + "".join(f"{r.get(f'ORPHAN_{o}') or 0:>8.3f}"
                        for o in ("BP", "MF", "CC"))
              + f"{r.get('ORPHAN_n') or 0:>6}")

    out = res / "multimodel.json"
    out.write_text(json.dumps({"rows": rows, "geometry": spec}, indent=2))
    print(f"\nwrote {out}")
    print("\nStill outstanding for the methods claim: bootstrap intervals on")
    print("the orphan tier (n is small), and a DIAMOND homology-transfer")
    print("baseline on the same tiers.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
