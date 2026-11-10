"""Configuration loading and path resolution for PmDark.

Single source of truth: config/config.yaml. Every script imports `load()`
rather than hard-coding a path, so moving the share folder is a one-line edit.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# Resolved against {share} — the archive tier, login node only.
_ARCHIVE = ("raw", "interim", "dbs", "results", "logs")

# Resolved against {work} — the only compute-visible persistent tier.
_WORK = ("models", "genome", "align", "proteins", "embeddings", "structures")

_DERIVED = _ARCHIVE + _WORK


@dataclass
class Config:
    raw: dict[str, Any]
    path: Path

    # -- convenience accessors ------------------------------------------
    @property
    def paths(self) -> dict[str, Path]:
        return self.raw["_paths"]

    def __getitem__(self, key: str) -> Any:
        return self.raw[key]

    def p(self, key: str, *parts: str) -> Path:
        """Resolve a configured directory, optionally joining sub-parts."""
        base = self.paths[key]
        return base.joinpath(*parts) if parts else base

    def mkdirs(self) -> None:
        for key in _DERIVED:
            self.paths[key].mkdir(parents=True, exist_ok=True)


def _default_config_path() -> Path:
    env = os.environ.get("PMDARK_CONFIG")
    if env:
        return Path(env).expanduser()
    # src/pmdark/config.py -> repo root -> config/config.yaml
    return Path(__file__).resolve().parents[2] / "config" / "config.yaml"


def load(path: str | Path | None = None) -> Config:
    cfg_path = Path(path).expanduser() if path else _default_config_path()
    if not cfg_path.is_file():
        raise FileNotFoundError(
            f"Config not found at {cfg_path}. "
            "Set PMDARK_CONFIG or pass an explicit path."
        )

    with cfg_path.open() as fh:
        raw = yaml.safe_load(fh)

    share = Path(str(raw["paths"]["share"])).expanduser()
    work = Path(str(raw["paths"]["work"])).expanduser()
    if "CHANGE_ME" in str(share) or "CHANGE_ME" in str(work):
        raise ValueError(
            "config/config.yaml still contains a placeholder path.\n"
            "Edit paths.share (archive) and paths.work (compute-visible)."
        )

    resolved: dict[str, Path] = {
        "share": share,
        "work": work,
        "home": Path(str(raw["paths"]["home"])).expanduser(),
    }
    for key in _DERIVED:
        tmpl = str(raw["paths"][key])
        resolved[key] = Path(
            tmpl.format(share=str(share), work=str(work))
        ).expanduser()

    raw["_paths"] = resolved
    return Config(raw=raw, path=cfg_path)


def _job_tag() -> str:
    job = os.environ.get("SLURM_JOB_ID", f"local{os.getpid()}")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    return f"pmdark_{job}" + (f"_{task}" if task else "")


def scratch(cfg: "Config", *parts: str, need_gb: float = 0.0,
            verbose: bool = True) -> Path:
    """Pick a scratch directory with enough free space, and create it.

    Slurm on this cluster does not manage local disk (sinfo reports
    TMP_DISK=0 on every node), and free space varies wildly:
        cpu01  /tmp   9.7 GB free of 79 GB   <-- nearly full
        gpu05  /tmp   64 GB free of 98 GB
    So a job cannot assume /tmp is usable. This checks, and falls back to
    the WORK tier (NFS, slower, quota'd) when local disk is too small.

    need_gb: required free space. Pass a real estimate — for alignment that
    is roughly 4x the compressed FASTQ size.
    """
    import shutil as _sh

    tag = _job_tag()
    candidates: list[tuple[str, Path]] = [
        ("node-local", Path(str(cfg.raw["paths"].get("scratch_root", "/tmp"))) / tag),
        ("work-tier", cfg.paths["work"] / "scratch" / tag),
    ]

    chosen: Path | None = None
    for label, cand in candidates:
        try:
            cand.mkdir(parents=True, exist_ok=True)
            free = _sh.disk_usage(cand).free / 1e9
        except OSError as exc:
            if verbose:
                print(f"  scratch: {label} {cand} unusable ({exc})")
            continue
        if free >= need_gb:
            if verbose:
                print(f"  scratch: using {label} {cand} ({free:.0f} GB free, "
                      f"need {need_gb:.0f} GB)")
            chosen = cand
            break
        if verbose:
            print(f"  scratch: {label} has {free:.0f} GB, need {need_gb:.0f} GB "
                  "— trying next")
        try:
            cand.rmdir()
        except OSError:
            pass

    if chosen is None:
        # Nothing met the requirement; take the roomiest and let the caller
        # decide whether to proceed.
        best, best_free = None, -1.0
        for _, cand in candidates:
            try:
                cand.mkdir(parents=True, exist_ok=True)
                f = _sh.disk_usage(cand).free / 1e9
            except OSError:
                continue
            if f > best_free:
                best, best_free = cand, f
        if best is None:
            raise RuntimeError("no usable scratch directory")
        if verbose:
            print(f"  scratch: WARNING no location has {need_gb:.0f} GB; "
                  f"using {best} with {best_free:.0f} GB")
        chosen = best

    if parts:
        chosen = chosen.joinpath(*parts)
        chosen.mkdir(parents=True, exist_ok=True)
    return chosen


def on_compute_node() -> bool:
    return "SLURM_JOB_ID" in os.environ


def describe(cfg: Config) -> str:
    """Human-readable summary, grouped by storage tier."""
    lines = [f"config : {cfg.path}", ""]

    groups = [
        ("ARCHIVE  (login node only - NOT visible to batch jobs)",
         ["share"] + list(_ARCHIVE)),
        ("WORK     (compute-visible, quota'd - keep under ~25 GB)",
         ["work"] + list(_WORK)),
        ("OTHER", ["home"]),
    ]
    for title, keys in groups:
        lines.append(title)
        for key in keys:
            val = cfg.paths.get(key)
            if val is None:
                continue
            mark = "ok" if val.exists() else "--"
            lines.append(f"  [{mark}] {key:<11} {val}")
        lines.append("")

    root = cfg.raw["paths"].get("scratch_root", "/tmp")
    lines.append(f"SCRATCH  (node-local, wiped at job end)\n  [  ] scratch     {root}/pmdark_$SLURM_JOB_ID")
    return "\n".join(lines)
