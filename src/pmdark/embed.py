"""Checkpointed protein-language-model embedding.

Supports ESM-2, ProtT5, ProstT5 (encoder) and SaProt. Each model exposes the
same interface: feed a list of sequences, get back one mean-pooled vector per
sequence (masked over real residues only — never over padding).

Designed to be launched as a Slurm job array with a hard walltime. The job
polls its own remaining time and exits cleanly before being killed, so the
in-flight shard is never half-written.
"""

from __future__ import annotations

import gc
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
# torch is imported LAZILY, inside the functions that need it.
#
# _detect_family and _prepare are pure Python — the model-family rules and
# the T5 space-separation / <AA2fold> prefix logic — and they are exactly
# what scripts/33b_verify_tokenization.py needs to check before any GPU
# time is spent. A module-level `import torch` made this file unimportable
# in the pmbio env, so the pre-flight check could not run at all without
# the GPU environment that it exists to validate.
if False:  # pragma: no cover - import-time typing only
    import torch

from .shards import Record, ShardStore, batched

# Residues outside the standard 20 are mapped to X for T5-family tokenizers
_NONSTANDARD = re.compile(r"[UZOB*]")


@dataclass
class EmbedConfig:
    model_key: str
    model_id: str
    cache_dir: Path
    max_len: int = 1022
    max_batch: int = 32
    max_tokens: int = 16384
    dtype: str = "float16"
    device: str = "cuda"


class Embedder:
    """Lazy-loading wrapper. The model is only materialised on first use."""

    def __init__(self, cfg: EmbedConfig):
        self.cfg = cfg
        self.model = None
        self.tokenizer = None
        self._family = self._detect_family(cfg.model_key, cfg.model_id)

    @staticmethod
    def _detect_family(key: str, model_id: str) -> str:
        probe = f"{key} {model_id}".lower()
        if "prostt5" in probe:
            return "prostt5"
        if "t5" in probe:
            return "t5"
        if "saprot" in probe:
            return "saprot"
        if "esm" in probe:
            return "esm"
        raise ValueError(f"Cannot infer model family from '{key}' / '{model_id}'")

    @property
    def torch_dtype(self):
        import torch
        return {"float16": torch.float16, "bfloat16": torch.bfloat16,
                "float32": torch.float32}[self.cfg.dtype]

    # -- loading --------------------------------------------------------
    def load(self) -> None:
        if self.model is not None:
            return

        import torch  # noqa: F401  (used via torch_dtype below)
        from transformers import AutoTokenizer, AutoModel, T5EncoderModel, T5Tokenizer

        cache = str(self.cfg.cache_dir)
        mid = self.cfg.model_id

        if self._family in ("t5", "prostt5"):
            # transformers 5.x removed the `legacy` kwarg from T5Tokenizer.
            # Verified on this cluster: transformers 5.22.1 reports
            # "T5Tokenizer legacy kwarg: False". Passing it raises TypeError.
            self.tokenizer = T5Tokenizer.from_pretrained(
                mid, do_lower_case=False, cache_dir=cache
            )
            self.model = T5EncoderModel.from_pretrained(
                mid, cache_dir=cache, torch_dtype=self.torch_dtype
            )
        else:
            self.tokenizer = AutoTokenizer.from_pretrained(mid, cache_dir=cache)
            self.model = AutoModel.from_pretrained(
                mid, cache_dir=cache, torch_dtype=self.torch_dtype
            )

        self.model = self.model.to(self.cfg.device).eval()
        for p in self.model.parameters():
            p.requires_grad_(False)

    # -- preprocessing ---------------------------------------------------
    def _prepare(self, seqs: list[str]) -> list[str]:
        out = []
        for s in seqs:
            s = s[: self.cfg.max_len].upper()
            if self._family in ("t5", "prostt5"):
                s = _NONSTANDARD.sub("X", s)
                s = " ".join(s)                      # T5 expects space-separated
                if self._family == "prostt5":
                    s = "<AA2fold> " + s             # AA->structure direction
            out.append(s)
        return out

    # -- inference -------------------------------------------------------
    def embed(self, seqs: list[str]) -> np.ndarray:
        import torch
        with torch.inference_mode():
            return self._embed_impl(seqs)

    def _embed_impl(self, seqs: list[str]) -> np.ndarray:
        import torch  # noqa: F401
        self.load()
        prepared = self._prepare(seqs)

        enc = self.tokenizer(
            prepared,
            add_special_tokens=True,
            padding="longest",
            truncation=True,
            max_length=self.cfg.max_len + 2,
            return_tensors="pt",
        ).to(self.cfg.device)

        out = self.model(**enc)
        hidden = out.last_hidden_state                       # (B, L, D)
        mask = enc["attention_mask"].unsqueeze(-1).to(hidden.dtype)

        # Exclude non-residue tokens from the mean. Getting this wrong
        # injects a LENGTH-DEPENDENT bias rather than a constant one: a
        # control token contributes with weight 1/(L+1), so short and long
        # proteins are biased by different amounts. Centering cannot remove
        # that, and it is precisely the kind of nuisance direction the
        # anisotropy analysis is about — so it must not be introduced here.
        #
        #   ESM / SaProt : [CLS] residues... [EOS]   -> drop first and last
        #   ProtT5       : residues... </s>          -> drop last
        #   ProstT5      : <AA2fold> residues... </s>-> drop FIRST and last
        #
        # The ProstT5 case was wrong: _prepare prepends <AA2fold> (required
        # by that model to select the amino-acid input direction) and the
        # pooling then averaged it in with the residues.
        lengths = enc["attention_mask"].sum(dim=1)
        drop_first = self._family in ("esm", "saprot", "prostt5")
        if drop_first:
            mask[:, 0, :] = 0
        for i, n in enumerate(lengths.tolist()):
            mask[i, n - 1, :] = 0

        pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
        return pooled.float().cpu().numpy().astype(np.float32)


# ----------------------------------------------------------------------
# Walltime guard
# ----------------------------------------------------------------------
class Deadline:
    """Stop cleanly before Slurm sends SIGKILL.

    Slurm exposes no reliable in-process remaining-time variable, so the
    submitting script passes --walltime-min and we reserve a safety margin.
    """

    def __init__(self, walltime_min: float, margin_min: float = 12.0):
        self.start = time.time()
        self.budget = max(0.0, (walltime_min - margin_min) * 60.0)

    @property
    def elapsed(self) -> float:
        return time.time() - self.start

    @property
    def remaining(self) -> float:
        return self.budget - self.elapsed

    def expired(self, next_step_estimate: float = 0.0) -> bool:
        return self.remaining <= next_step_estimate


def run_shards(
    records_by_shard: list[list[Record]],
    store: ShardStore,
    embedder: Embedder,
    shard_indices: list[int],
    deadline: Deadline,
    log=print,
) -> dict:
    import torch
    """Process assigned shards until done or out of time. Idempotent."""
    stats = {"completed": [], "skipped": [], "stopped_early": False}
    durations: list[float] = []

    for idx in shard_indices:
        if store.is_done(idx):
            stats["skipped"].append(idx)
            continue

        # Estimate next shard cost from observed median; be pessimistic early
        est = float(np.median(durations)) * 1.3 if durations else 0.0
        if deadline.expired(est):
            log(f"[deadline] {deadline.remaining:.0f}s left, est {est:.0f}s — stopping")
            stats["stopped_early"] = True
            break

        t0 = time.time()
        recs = records_by_shard[idx]
        ids: list[str] = []
        vecs: list[np.ndarray] = []

        for batch in batched(recs, embedder.cfg.max_tokens, embedder.cfg.max_batch):
            try:
                vecs.append(embedder.embed([r.seq for r in batch]))
                ids.extend(r.pid for r in batch)
            except torch.cuda.OutOfMemoryError:
                # Fall back to one-at-a-time for this batch rather than losing
                # the whole shard to a single pathological sequence.
                torch.cuda.empty_cache()
                log(f"  OOM on batch of {len(batch)} — retrying singly")
                for rec in batch:
                    vecs.append(embedder.embed([rec.seq]))
                    ids.append(rec.pid)

        store.save(idx, ids, np.vstack(vecs))
        dt = time.time() - t0
        durations.append(dt)
        stats["completed"].append(idx)

        done, total = store.progress()
        log(f"[shard {idx:05d}] {len(ids)} seqs in {dt:.1f}s — {done}/{total} done")

        del vecs
        gc.collect()
        torch.cuda.empty_cache()

    return stats


def slice_for_task(n_shards: int, task_id: int, n_tasks: int) -> list[int]:
    """Interleaved assignment so every array task gets a mix of long and
    short shards (shards are length-sorted, so contiguous blocks would be
    wildly unbalanced)."""
    return list(range(task_id, n_shards, n_tasks))


def gpu_report() -> str:
    import torch
    if not torch.cuda.is_available():
        return "no CUDA device visible"
    i = torch.cuda.current_device()
    props = torch.cuda.get_device_properties(i)
    return f"{props.name} ({props.total_memory / 1024**3:.0f} GB), cuda:{i}"


def autoscale_batch(default: int) -> int:
    import torch
    """Shrink batch size on small cards so one config works across partitions."""
    if not torch.cuda.is_available():
        return 1
    gb = torch.cuda.get_device_properties(0).total_memory / 1024**3
    if gb >= 44:      # L40 48 GB
        return default
    if gb >= 22:      # RTX 4090 / A4000-ish
        return max(1, default // 3)
    return max(1, default // 6)
