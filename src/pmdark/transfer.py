"""Embedding-based GO transfer — the method the calibration benchmark grades.

Given embeddings for a reference set with known GO terms and a query set,
transfer annotations by weighted k-nearest-neighbour vote in embedding space.
This is the core of FANTASIA-style PLM annotation, reimplemented here so the
scoring is transparent and auditable rather than hidden behind a pipeline.

Scoring. For query q with neighbours n_1..n_k at cosine similarities s_i:

    score(q, term) = sum_i [ w_i * 1(term in GO(n_i)) ] / sum_i w_i

with w_i = max(0, s_i - sim_floor) ** power. The floor discards neighbours
that are merely "not opposite", and the exponent sharpens the vote so a
single close hit outweighs several distant ones. Scores land in [0, 1] and
are directly thresholdable, which is what Fmax needs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def l2_normalise(x: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(n, 1e-12, None)


def center(x: np.ndarray, mean: np.ndarray | None = None
           ) -> tuple[np.ndarray, np.ndarray]:
    """Subtract the dataset mean before computing cosine similarity.

    WHY THIS IS NOT OPTIONAL
    ------------------------
    Mean-pooled transformer embeddings are ANISOTROPIC: they occupy a narrow
    cone rather than spreading over the sphere, so raw cosine similarity
    between ANY two proteins is ~0.99. Measured on this benchmark:

        median top-1 similarity, all references available : 0.993
        median top-1 similarity, every homolog excluded   : 0.987

    A metric that cannot distinguish "closest relative in Swiss-Prot" from
    "no relative at all" is not measuring homology. It also makes a
    similarity floor useless — nothing ever falls below 0.5, which is why
    the benchmark reported 0 queries with no usable neighbour even in the
    orphan condition.

    Subtracting the mean vector re-centres the cloud on the origin and
    restores the usable dynamic range. The mean MUST be computed on the
    reference set and reused for the queries, or the two sets are projected
    into different spaces.
    """
    if mean is None:
        mean = x.mean(axis=0, keepdims=True)
    return x - mean, mean


def whiten(x: np.ndarray, stats: tuple | None = None
           ) -> tuple[np.ndarray, tuple]:
    """Centre and decorrelate (ZCA-style). Stronger than centring alone.

    Removes the dominant directions that carry length/composition signal
    rather than function. Slower and can over-correct on small reference
    sets, so centring is the default and this is opt-in.
    """
    if stats is None:
        mu = x.mean(axis=0, keepdims=True)
        xc = x - mu
        cov = (xc.T @ xc) / max(1, xc.shape[0] - 1)
        vals, vecs = np.linalg.eigh(cov)
        vals = np.clip(vals, 1e-8, None)
        w = vecs @ np.diag(1.0 / np.sqrt(vals)) @ vecs.T
        stats = (mu, w)
    mu, w = stats
    return (x - mu) @ w, stats


@dataclass
class TransferConfig:
    k: int = 25
    sim_floor: float = 0.5      # ignore neighbours below this cosine similarity
    power: float = 3.0          # sharpen the weighting
    block: int = 2048           # query rows per similarity block
    min_score: float = 0.01     # drop near-zero scores to keep output sparse


def knn(
    query: np.ndarray,
    reference: np.ndarray,
    k: int,
    block: int = 2048,
) -> tuple[np.ndarray, np.ndarray]:
    """Blocked exact cosine kNN.

    Exact rather than approximate: the reference set here is tens of
    thousands of proteins, so a full matmul is affordable and removes ANN
    recall as a confound in the benchmark. Blocked to bound memory.
    """
    q = l2_normalise(query.astype(np.float32))
    r = l2_normalise(reference.astype(np.float32))
    k = min(k, r.shape[0])

    idx_out = np.empty((q.shape[0], k), dtype=np.int32)
    sim_out = np.empty((q.shape[0], k), dtype=np.float32)

    for start in range(0, q.shape[0], block):
        stop = min(start + block, q.shape[0])
        sims = q[start:stop] @ r.T                       # (b, R)
        part = np.argpartition(-sims, kth=k - 1, axis=1)[:, :k]
        rows = np.arange(stop - start)[:, None]
        part_sims = sims[rows, part]
        order = np.argsort(-part_sims, axis=1)
        idx_out[start:stop] = part[rows, order]
        sim_out[start:stop] = part_sims[rows, order]

    return idx_out, sim_out


def transfer_go(
    query_ids: list[str],
    query_emb: np.ndarray,
    ref_ids: list[str],
    ref_emb: np.ndarray,
    ref_go: dict[str, set[str]],
    cfg: TransferConfig | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, dict]]:
    """Return (scored annotations, per-query neighbour diagnostics)."""
    cfg = cfg or TransferConfig()
    idx, sim = knn(query_emb, ref_emb, cfg.k, cfg.block)

    scored: dict[str, dict[str, float]] = {}
    diag: dict[str, dict] = {}

    for i, qid in enumerate(query_ids):
        neigh_idx, neigh_sim = idx[i], sim[i]
        weights = np.clip(neigh_sim - cfg.sim_floor, 0.0, None) ** cfg.power
        total = float(weights.sum())

        diag[qid] = {
            "top_sim": float(neigh_sim[0]),
            "mean_sim": float(neigh_sim.mean()),
            "n_above_floor": int((neigh_sim > cfg.sim_floor).sum()),
            "top_neighbour": ref_ids[int(neigh_idx[0])],
        }

        if total <= 0:
            # Every neighbour is below the floor: this query has no usable
            # evidence. Emitting nothing is the honest outcome, and the
            # benchmark counts it against recall rather than precision.
            scored[qid] = {}
            continue

        votes: dict[str, float] = {}
        for j, w in zip(neigh_idx, weights):
            if w <= 0:
                continue
            for term in ref_go.get(ref_ids[int(j)], ()):
                votes[term] = votes.get(term, 0.0) + float(w)

        scored[qid] = {t: v / total for t, v in votes.items()
                       if v / total >= cfg.min_score}

    return scored, diag


def propagate_scores(
    scored: dict[str, dict[str, float]], graph
) -> dict[str, dict[str, float]]:
    """Push scores up the GO DAG, taking the max over descendants.

    Required before evaluation: a parent term is at least as likely as its
    most likely child. Skipping this makes recall look far worse than it is.
    """
    out: dict[str, dict[str, float]] = {}
    for qid, terms in scored.items():
        prop: dict[str, float] = {}
        for term, score in terms.items():
            for anc in graph.ancestors(term):
                if anc in graph.__class__.__dict__.get("ROOTS", ()):  # noqa
                    continue
                if score > prop.get(anc, 0.0):
                    prop[anc] = score
        from .goeval import ROOTS
        out[qid] = {t: s for t, s in prop.items() if t not in ROOTS}
    return out
