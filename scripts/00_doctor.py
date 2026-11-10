#!/usr/bin/env python3
"""Pre-flight check. Run this FIRST on the cluster, before anything else.

    python scripts/00_doctor.py

Verifies: config resolves, share folder is writable and large enough,
GPU is visible, required packages import, and prints a storage budget for
the databases this project needs.
"""

from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

OK, WARN, BAD = "  OK  ", " WARN ", " FAIL "

# Approximate on-disk sizes. Verify against the download page before relying
# on these for a quota request — database releases grow.
STORAGE_BUDGET = [
    ("HuggingFace model cache (ESM-2 650M, ProtT5, ProstT5, SaProt)", 25, True),
    ("ESM-2 3B (optional, L40 only)",                                 12, False),
    ("Embeddings, 28k proteins x 4 models, fp32",                      3, True),
    ("Foldseek PDB database",                                         25, True),
    ("Foldseek AlphaFold/Swiss-Prot database",                        10, True),
    ("Foldseek BFVD (viral structures)",                              20, False),
    ("Foldseek AFDB50 (full, optional - very large)",                700, False),
    ("ESMFold predicted structures (~12k dark proteins)",             40, False),
    ("Boltz-2 / co-folding outputs",                                  60, False),
    ("Public MS/MS raw files for proteogenomics",                    150, False),
    ("SRA raw reads (12 libraries, if re-assembling)",                60, False),
    ("Working / scratch headroom",                                    50, True),
]

PACKAGES = [
    ("yaml", "pyyaml"),
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("torch", "torch"),
    ("transformers", "transformers"),
    ("sentencepiece", "sentencepiece"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("Bio", "biopython"),
    ("h5py", "h5py"),
    ("tqdm", "tqdm"),
]


def line(status: str, msg: str) -> None:
    print(f"[{status}] {msg}")


def check_config():
    print("\n--- configuration ---")
    try:
        from pmdark import config as cfgmod
        cfg = cfgmod.load()
    except Exception as exc:
        line(BAD, f"config failed to load: {exc}")
        return None
    line(OK, f"loaded {cfg.path}")
    print(cfgmod.describe(cfg))
    return cfg


def check_storage(cfg) -> None:
    print("\n--- storage ---")
    share = cfg.p("share")
    try:
        share.mkdir(parents=True, exist_ok=True)
        probe = share / ".pmdark_write_test"
        probe.write_text("ok")
        probe.unlink()
        line(OK, f"share is writable: {share}")
    except Exception as exc:
        line(BAD, f"cannot write to share {share}: {exc}")
        return

    usage = shutil.disk_usage(share)
    free_gb = usage.free / 1024**3
    total_gb = usage.total / 1024**3
    line(OK, f"filesystem: {free_gb:,.0f} GB free of {total_gb:,.0f} GB")

    required = sum(gb for _, gb, req in STORAGE_BUDGET if req)
    optional = sum(gb for _, gb, req in STORAGE_BUDGET if not req)

    print(f"\n  {'component':<58} {'GB':>6}  {'tier':<9}")
    print(f"  {'-' * 58} {'-' * 6}  {'-' * 9}")
    for name, gb, req in STORAGE_BUDGET:
        print(f"  {name:<58} {gb:>6}  {'required' if req else 'optional':<9}")
    print(f"  {'-' * 58} {'-' * 6}")
    print(f"  {'REQUIRED total':<58} {required:>6}")
    print(f"  {'OPTIONAL total':<58} {optional:>6}")
    print(f"  {'FULL total':<58} {required + optional:>6}\n")

    if free_gb < required:
        line(BAD, f"share has {free_gb:,.0f} GB; core pipeline needs ~{required} GB")
    elif free_gb < required + 200:
        line(WARN, f"{free_gb:,.0f} GB free — enough for the core path, but skip "
                   "AFDB50 and keep MS/MS downloads trimmed")
    else:
        line(OK, f"{free_gb:,.0f} GB free — comfortable for the full path")

    home_usage = shutil.disk_usage(Path.home())
    line(OK if home_usage.free / 1024**3 > 20 else WARN,
         f"/home free: {home_usage.free / 1024**3:,.0f} GB "
         "(code and small tables only — never put caches here)")


def check_packages() -> None:
    print("\n--- python packages ---")
    missing = []
    for mod, pip_name in PACKAGES:
        try:
            m = importlib.import_module(mod)
            ver = getattr(m, "__version__", "?")
            line(OK, f"{pip_name:<16} {ver}")
        except ImportError:
            line(BAD, f"{pip_name:<16} MISSING")
            missing.append(pip_name)
    if missing:
        print(f"\n  pip install {' '.join(missing)}")


def check_gpu() -> None:
    print("\n--- gpu ---")
    try:
        import torch
    except ImportError:
        line(BAD, "torch not installed")
        return
    if not torch.cuda.is_available():
        line(WARN, "no CUDA device visible (expected on a login node — "
                   "re-run inside an srun/sbatch allocation)")
        return
    for i in range(torch.cuda.device_count()):
        p = torch.cuda.get_device_properties(i)
        gb = p.total_memory / 1024**3
        note = ""
        if gb < 20:
            note = "  -> use --batch 8, skip esm2_3b"
        elif gb < 40:
            note = "  -> use --batch 24, skip esm2_3b"
        else:
            note = "  -> full batch, esm2_3b OK"
        line(OK, f"cuda:{i} {p.name} {gb:.0f} GB{note}")


def check_external() -> None:
    print("\n--- external binaries (optional at this stage) ---")
    for tool, why in [
        ("foldseek", "structural homology search"),
        ("mmseqs", "sequence clustering / baseline search"),
        ("hmmscan", "Pfam domain baseline"),
        ("git", "clone the upstream annotation repo"),
        ("sbatch", "job submission"),
    ]:
        path = shutil.which(tool)
        line(OK if path else WARN, f"{tool:<10} {path or 'not found'}  ({why})")


def check_slurm() -> None:
    print("\n--- slurm ---")
    try:
        out = subprocess.run(
            ["sinfo", "-h", "-o", "%P %a %l %D %G"],
            capture_output=True, text=True, timeout=15,
        )
        if out.returncode == 0 and out.stdout.strip():
            print(out.stdout.rstrip())
        else:
            line(WARN, "sinfo returned nothing")
    except Exception as exc:
        line(WARN, f"sinfo unavailable: {exc}")


def main() -> int:
    print("=" * 74)
    print("PmDark pre-flight check")
    print("=" * 74)
    cfg = check_config()
    if cfg is not None:
        check_storage(cfg)
    check_packages()
    check_gpu()
    check_external()
    check_slurm()
    print("\nDone. Resolve every FAIL before running stage 01.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
