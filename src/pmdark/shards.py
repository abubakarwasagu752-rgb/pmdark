"""Deterministic sharding + resume logic.

Every GPU partition on this cluster caps at 12 h walltime, so no long job may
assume it will finish. The contract here:

  * Proteins are sorted by DESCENDING length, then split into fixed-size
    shards. Length-sorting keeps each batch's padding tight, which is the
    single biggest throughput win for transformer inference.
  * Each shard writes exactly one .npy + one .done marker, atomically.
  * A restarted job skips any shard whose .done marker exists.

Result: `sbatch --array` the same script as many times as needed; it converges.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np


@dataclass(frozen=True)
class Record:
    pid: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


def read_fasta(path: str | Path) -> list[Record]:
    """Minimal FASTA reader. Takes the first whitespace token as the ID."""
    records: list[Record] = []
    pid: str | None = None
    chunks: list[str] = []

    with Path(path).open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    records.append(Record(pid, "".join(chunks)))
                pid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if pid is not None:
        records.append(Record(pid, "".join(chunks)))

    if not records:
        raise ValueError(f"No sequences parsed from {path}")
    return records


def write_fasta(records: list[Record], path: str | Path, width: int = 60) -> None:
    with Path(path).open("w") as fh:
        for rec in records:
            fh.write(f">{rec.pid}\n")
            for i in range(0, len(rec.seq), width):
                fh.write(rec.seq[i : i + width] + "\n")


def fasta_checksum(records: list[Record]) -> str:
    """Stable digest of (id, seq) pairs — guards against silent input drift.

    If this changes between runs, previously written shards are invalid and
    the resume logic must not reuse them.
    """
    h = hashlib.sha256()
    for rec in records:
        h.update(rec.pid.encode())
        h.update(b"\0")
        h.update(rec.seq.encode())
        h.update(b"\n")
    return h.hexdigest()[:16]


def plan_shards(records: list[Record], shard_size: int) -> list[list[Record]]:
    """Length-descending sort, then fixed-size chunks.

    Sorting is by (-length, pid) so the plan is fully deterministic and
    reproducible across machines regardless of input file order.
    """
    ordered = sorted(records, key=lambda r: (-r.length, r.pid))
    return [ordered[i : i + shard_size] for i in range(0, len(ordered), shard_size)]


class ShardStore:
    """Filesystem-backed, resumable shard output directory."""

    def __init__(self, root: str | Path, checksum: str, n_shards: int):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.checksum = checksum
        self.n_shards = n_shards
        self._verify_manifest()

    # -- manifest -------------------------------------------------------
    @property
    def manifest_path(self) -> Path:
        return self.root / "manifest.json"

    def _verify_manifest(self) -> None:
        payload = {"checksum": self.checksum, "n_shards": self.n_shards}
        if self.manifest_path.exists():
            existing = json.loads(self.manifest_path.read_text())
            if existing != payload:
                raise RuntimeError(
                    f"Shard directory {self.root} was built from different input.\n"
                    f"  on disk : {existing}\n"
                    f"  current : {payload}\n"
                    "Delete the directory or point --out somewhere new."
                )
        else:
            self.manifest_path.write_text(json.dumps(payload, indent=2))

    # -- per-shard paths ------------------------------------------------
    def emb_path(self, idx: int) -> Path:
        return self.root / f"shard_{idx:05d}.npy"

    def ids_path(self, idx: int) -> Path:
        return self.root / f"shard_{idx:05d}.ids.txt"

    def done_path(self, idx: int) -> Path:
        return self.root / f"shard_{idx:05d}.done"

    def is_done(self, idx: int) -> bool:
        return self.done_path(idx).exists()

    def pending(self, indices: list[int] | None = None) -> list[int]:
        candidates = indices if indices is not None else range(self.n_shards)
        return [i for i in candidates if not self.is_done(i)]

    # -- atomic write ---------------------------------------------------
    def save(self, idx: int, ids: list[str], matrix: np.ndarray) -> None:
        if matrix.shape[0] != len(ids):
            raise ValueError(
                f"shard {idx}: {matrix.shape[0]} vectors vs {len(ids)} ids"
            )

        tmp_emb = self.emb_path(idx).with_suffix(".npy.tmp")
        tmp_ids = self.ids_path(idx).with_suffix(".txt.tmp")

        # np.save APPENDS '.npy' unless the path already ends with it, so
        # np.save("x.npy.tmp", a) silently writes "x.npy.tmp.npy" and the
        # rename below then fails with FileNotFoundError. Passing an open
        # file object bypasses that behaviour entirely.
        with tmp_emb.open("wb") as fh:
            np.save(fh, matrix)
        tmp_ids.write_text("\n".join(ids) + "\n")

        # rename is atomic within a filesystem; .done is written last so a
        # crash mid-write never leaves a shard that looks complete
        tmp_emb.replace(self.emb_path(idx))
        tmp_ids.replace(self.ids_path(idx))
        self.done_path(idx).write_text(f"{len(ids)}\t{matrix.shape[1]}\n")

    # -- assembly -------------------------------------------------------
    def load_all(self) -> tuple[list[str], np.ndarray]:
        missing = self.pending()
        if missing:
            raise RuntimeError(
                f"{len(missing)} of {self.n_shards} shards incomplete "
                f"(first: {missing[:5]}). Re-run the job array."
            )
        ids: list[str] = []
        mats: list[np.ndarray] = []
        for i in range(self.n_shards):
            ids.extend(self.ids_path(i).read_text().split())
            mats.append(np.load(self.emb_path(i)))
        return ids, np.vstack(mats)

    def progress(self) -> tuple[int, int]:
        done = sum(1 for i in range(self.n_shards) if self.is_done(i))
        return done, self.n_shards


def batched(records: list[Record], max_tokens: int, max_batch: int) -> Iterator[list[Record]]:
    """Token-budgeted batching.

    Records arrive length-descending, so the first record in a batch sets the
    padded width. We cap on padded tokens rather than sequence count, which
    prevents OOM on the long-sequence shards while keeping short-sequence
    shards fast.
    """
    batch: list[Record] = []
    width = 0
    for rec in records:
        new_width = max(width, rec.length)
        if batch and (new_width * (len(batch) + 1) > max_tokens or len(batch) >= max_batch):
            yield batch
            batch, width = [rec], rec.length
        else:
            batch.append(rec)
            width = new_width
    if batch:
        yield batch
