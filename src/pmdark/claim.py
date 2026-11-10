"""NFS-safe shard claiming for heterogeneous multi-partition fan-out.

The problem this solves:

  We want to submit independent job arrays to gpul40, gpu4090, gpu3060,
  gpu3080, gpu4070ti and gpua4000 at the same time, all chewing through one
  shared shard directory on /shared/a2-public. Those jobs cannot see each
  other, start at unpredictable times, and are killed at unpredictable times
  by the 12 h cap. Without coordination, two GPUs redo the same shard.

Why mkdir and not a lock file:

  O_CREAT|O_EXCL is not reliably atomic on NFSv3, and silently succeeding
  twice is exactly the failure we cannot tolerate. mkdir() IS atomic on NFS
  — it is the classic NFS mutual-exclusion primitive. A directory either
  gets created by exactly one caller or raises FileExistsError.

Staleness:

  A job killed by the walltime cap leaves its claim behind. Each claim
  carries a heartbeat file; a claim whose heartbeat has not advanced for
  `stale_after` seconds is reclaimable. Heartbeats are refreshed from the
  worker while it processes a shard.
"""

from __future__ import annotations

import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path


def _identity() -> dict:
    return {
        "host": socket.gethostname(),
        "pid": os.getpid(),
        "job": os.environ.get("SLURM_JOB_ID", "local"),
        "task": os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
        "partition": os.environ.get("SLURM_JOB_PARTITION", "?"),
        "gpu": os.environ.get("CUDA_VISIBLE_DEVICES", "?"),
        "claimed_at": time.time(),
    }


@dataclass
class Claim:
    """A held claim on one shard. Refresh it; release it when done."""

    path: Path
    owner: dict

    @property
    def heartbeat(self) -> Path:
        return self.path / "heartbeat"

    def refresh(self) -> None:
        try:
            self.heartbeat.write_text(str(time.time()))
        except OSError:
            pass  # a failed heartbeat is not fatal; staleness just ticks up

    def release(self) -> None:
        """Release without marking done — another worker may retake it."""
        try:
            for child in self.path.iterdir():
                child.unlink(missing_ok=True)
            self.path.rmdir()
        except OSError:
            pass


class ClaimRegistry:
    """Atomic shard claiming inside a shard store directory."""

    def __init__(self, root: str | Path, stale_after: float = 3600.0):
        self.root = Path(root) / ".claims"
        self.root.mkdir(parents=True, exist_ok=True)
        self.stale_after = stale_after

    def _dir(self, idx: int) -> Path:
        return self.root / f"shard_{idx:05d}"

    # -- inspection -----------------------------------------------------
    def _age(self, claim_dir: Path) -> float:
        """Seconds since this claim last showed life."""
        hb = claim_dir / "heartbeat"
        try:
            if hb.exists():
                return time.time() - float(hb.read_text().strip())
            meta = claim_dir / "owner.json"
            if meta.exists():
                return time.time() - json.loads(meta.read_text())["claimed_at"]
        except (OSError, ValueError, KeyError):
            pass
        return float("inf")  # unreadable claim is treated as stale

    def is_claimed(self, idx: int) -> bool:
        d = self._dir(idx)
        return d.exists() and self._age(d) < self.stale_after

    def owner_of(self, idx: int) -> dict | None:
        meta = self._dir(idx) / "owner.json"
        try:
            return json.loads(meta.read_text())
        except (OSError, ValueError):
            return None

    # -- acquisition ----------------------------------------------------
    def try_claim(self, idx: int) -> Claim | None:
        """Atomically claim shard `idx`, or return None if someone else has it."""
        d = self._dir(idx)

        try:
            d.mkdir()                      # <-- the atomic operation
        except FileExistsError:
            age = self._age(d)
            if age < self.stale_after:
                return None                # live claim, respect it
            # Stale: the holder was almost certainly killed by the 12 h cap.
            # Steal it by replacing the owner record. Worst case two workers
            # redo one shard, which is wasteful but never corrupting, because
            # shard writes are atomic and content-identical.
            try:
                (d / "owner.json").write_text(
                    json.dumps({**_identity(), "stole_after_s": round(age)}, indent=2)
                )
            except OSError:
                return None
            claim = Claim(d, _identity())
            claim.refresh()
            return claim
        except OSError:
            return None

        owner = _identity()
        try:
            (d / "owner.json").write_text(json.dumps(owner, indent=2))
        except OSError:
            pass
        claim = Claim(d, owner)
        claim.refresh()
        return claim

    def claim_next(self, candidates: list[int], is_done) -> Claim | None:
        """Walk candidates and claim the first that is neither done nor held."""
        for idx in candidates:
            if is_done(idx):
                continue
            claim = self.try_claim(idx)
            if claim is not None:
                if is_done(idx):           # finished between our two checks
                    claim.release()
                    continue
                claim.owner["shard"] = idx
                return claim
        return None

    # -- reporting ------------------------------------------------------
    def status(self, n_shards: int, is_done) -> dict:
        done = live = stale = free = 0
        by_partition: dict[str, int] = {}
        for i in range(n_shards):
            if is_done(i):
                done += 1
                continue
            d = self._dir(i)
            if not d.exists():
                free += 1
                continue
            if self._age(d) < self.stale_after:
                live += 1
                owner = self.owner_of(i) or {}
                key = f"{owner.get('partition', '?')}/{owner.get('host', '?')}"
                by_partition[key] = by_partition.get(key, 0) + 1
            else:
                stale += 1
        return {
            "done": done, "in_progress": live,
            "stale": stale, "unclaimed": free,
            "total": n_shards, "workers": by_partition,
        }

    def clear_stale(self) -> int:
        """Remove every stale claim. Safe to run while jobs are live."""
        removed = 0
        for d in sorted(self.root.glob("shard_*")):
            if self._age(d) >= self.stale_after:
                try:
                    for child in d.iterdir():
                        child.unlink(missing_ok=True)
                    d.rmdir()
                    removed += 1
                except OSError:
                    pass
        return removed
