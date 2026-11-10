#!/usr/bin/env python3
"""Check each model tokenizes and pools correctly BEFORE spending GPU hours.

    conda activate ~/pmdark_work/envs/serprot
    python scripts/33b_verify_tokenization.py

WHY THIS EXISTS

  A tokenization or pooling mistake does not crash. It produces embeddings
  that are quietly wrong, and the multi-model comparison would then report
  "whitening does not help ProstT5" as a RESULT when it is a BUG. That
  conclusion would go in the paper.

  One such bug was already found by reading the code: _prepare prepends
  <AA2fold> for ProstT5, which that model requires in order to select the
  amino-acid input direction, and the mean-pooling mask excluded only the
  trailing </s>. The control token was therefore averaged in with the
  residues. Its weight is 1/(L+1), so it biases short and long proteins by
  DIFFERENT amounts — a length-dependent nuisance direction, which is
  exactly what this project's anisotropy analysis is trying to measure.
  Centering cannot remove it.

WHAT IT CHECKS, PER MODEL

  1. The tokenizer loads at all (T5Tokenizer without the `legacy` kwarg,
     which transformers 5.x removed).
  2. For ProstT5, <AA2fold> maps to exactly ONE token id. If it is split
     into word pieces the prefix is wrong and the model is being driven in
     an undefined direction.
  3. The token sequence, so the special-token layout is visible rather
     than assumed: where BOS/EOS/prefix actually land.
  4. The pooling mask that Embedder would build, printed as a diagram, and
     a check that every unmasked position is a residue and every residue
     is unmasked.

  Tokenizers only — no model weights are downloaded, so this runs on the
  login node in seconds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod  # noqa: E402

TEST = "MKVLAAGIVGLW"          # 12 residues, all standard
AA = set("ACDEFGHIKLMNPQRSTVWY")


def find_plm_python() -> list[str]:
    """Interpreters that can import transformers. Same approach as
    find_diamond in 30b: search known env locations instead of trusting a
    hardcoded path, because ~/pmdark_work/envs/serprot does not exist."""
    import subprocess
    seen: list[str] = []
    roots = [Path.home() / "pmdark_work" / "envs",
             Path.home() / ".conda" / "envs",
             Path.home() / "miniconda3" / "envs",
             Path.home() / "anaconda3" / "envs",
             Path("/opt/conda/envs")]
    cands = [Path(sys.executable)]
    for r in roots:
        try:
            cands += sorted(r.glob("*/bin/python"))
        except OSError:
            continue
    for c in cands:
        if not c.exists() or str(c) in seen:
            continue
        try:
            r = subprocess.run([str(c), "-c", "import transformers"],
                               capture_output=True, timeout=120)
            if r.returncode == 0:
                seen.append(str(c))
        except (OSError, subprocess.TimeoutExpired):
            continue
    return seen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", default="esm2,prott5,prostt5")
    args = ap.parse_args()

    cfg = cfgmod.load()
    from pmdark.embed import Embedder, EmbedConfig

    # transformers lives in the PLM env, which is not pmbio and is not
    # necessarily at ~/pmdark_work/envs/serprot. Locate an interpreter that
    # actually has it rather than guessing a path.
    have_tf = True
    try:
        import transformers  # noqa: F401
    except ImportError:
        have_tf = False
        print("\n" + "!" * 74)
        print("transformers is not importable in this interpreter:")
        print(f"  {sys.executable}")
        cands = find_plm_python()
        if cands:
            print("\nthese interpreters DO have it:")
            for c in cands:
                print(f"  {c}")
            print("\nrerun with one of them, e.g.")
            print(f"  {cands[0]} scripts/33b_verify_tokenization.py")
        else:
            print("\nnone found. List your environments with:")
            print("  conda env list")
            print("  ls -d ~/*/envs/* ~/.conda/envs/* 2>/dev/null")
            print("then run this script with that env's bin/python.")
        print("!" * 74)
        print("\nThe tokenizer checks are SKIPPED. The family and sequence")
        print("preparation checks below need no torch and still run.")

    ok_all = True
    for key in [m.strip() for m in args.models.split(",") if m.strip()]:
        mid = cfg["models"].get(key)
        print("\n" + "=" * 74)
        print(f"{key}   {mid}")
        print("=" * 74)
        if not mid:
            print("  not in config/config.yaml models:")
            ok_all = False
            continue

        emb = Embedder(EmbedConfig(
            model_key=key, model_id=mid, device="cpu",
            cache_dir=cfg.p("models", "hf"), dtype="float32",
            max_len=cfg["compute"]["max_len"], max_batch=1,
        ))
        print(f"  family detected: {emb._family}")

        if not have_tf:
            prepared = emb._prepare([TEST])[0]
            print(f"  prepared: {prepared!r}")
            exp_prefix = "<AA2fold> " if emb._family == "prostt5" else ""
            exp = exp_prefix + (" ".join(TEST)
                                if emb._family in ("t5", "prostt5") else TEST)
            good = prepared == exp
            print(f"  expected: {exp!r}  {'ok' if good else 'MISMATCH'}")
            ok_all = ok_all and good
            continue

        try:
            from transformers import AutoTokenizer, T5Tokenizer
            if emb._family in ("t5", "prostt5"):
                tok = T5Tokenizer.from_pretrained(
                    mid, do_lower_case=False,
                    cache_dir=str(cfg.p("models", "hf")))
            else:
                tok = AutoTokenizer.from_pretrained(
                    mid, cache_dir=str(cfg.p("models", "hf")))
        except Exception as e:                            # noqa: BLE001
            print(f"  TOKENIZER FAILED: {type(e).__name__}: {e}")
            ok_all = False
            continue
        print("  tokenizer loaded")

        # -- the ProstT5 prefix must be ONE token ------------------------
        if emb._family == "prostt5":
            ids = tok("<AA2fold>", add_special_tokens=False)["input_ids"]
            pieces = tok.convert_ids_to_tokens(ids)
            good = len(ids) == 1
            print(f"  <AA2fold> -> {len(ids)} token(s): {pieces}  "
                  f"{'ok' if good else 'WRONG — prefix is being split'}")
            if not good:
                print("    the model is then driven in an undefined "
                      "direction; do not embed with this tokenizer")
                ok_all = False

        # -- the real token layout --------------------------------------
        prepared = emb._prepare([TEST])[0]
        print(f"  prepared: {prepared!r}")
        enc = tok([prepared], add_special_tokens=True,
                  padding="longest", truncation=True,
                  max_length=cfg["compute"]["max_len"] + 2)
        ids = enc["input_ids"][0]
        toks = tok.convert_ids_to_tokens(ids)
        n = len(ids)

        drop_first = emb._family in ("esm", "saprot", "prostt5")
        keep = [True] * n
        if drop_first:
            keep[0] = False
        keep[n - 1] = False

        print(f"  {'idx':>4} {'token':<12} {'pooled?':<8}")
        for i, (t, k) in enumerate(zip(toks, keep)):
            flag = "yes" if k else "EXCLUDED"
            print(f"  {i:>4} {t:<12} {flag:<8}")

        # -- every pooled position must be a residue, and vice versa ----
        pooled_toks = [t for t, k in zip(toks, keep) if k]
        residue_like = [t.strip("▁ ").upper() for t in pooled_toks]
        bad = [t for t in residue_like if t not in AA]
        n_res = sum(1 for t in residue_like if t in AA)
        print(f"\n  pooled positions: {len(pooled_toks)}   "
              f"residues among them: {n_res}   expected: {len(TEST)}")
        if bad:
            print(f"  NON-RESIDUE TOKENS BEING POOLED: {bad}")
            ok_all = False
        elif n_res != len(TEST):
            print(f"  COUNT MISMATCH — pooling {n_res} of {len(TEST)} residues")
            ok_all = False
        else:
            print("  every pooled position is a residue, and every residue "
                  "is pooled")

    print("\n" + "=" * 74)
    if not ok_all:
        print("PROBLEMS FOUND — fix before spending GPU time")
        rc = 1
    elif not have_tf:
        # Do not say "safe to launch" when the decisive check never ran.
        print("PARTIAL — sequence preparation is correct for every model,")
        print("but the TOKENIZER checks were skipped (no transformers here).")
        print("The <AA2fold>-is-one-token check is the one that matters most")
        print("for ProstT5, and it has NOT been performed.")
        rc = 2
    else:
        print("ALL MODELS OK — safe to launch the GPU fan-outs")
        rc = 0
    print("=" * 74)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
