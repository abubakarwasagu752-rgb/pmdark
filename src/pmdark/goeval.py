"""GO ontology handling and CAFA-style evaluation.

Implements the metrics the calibration benchmark reports. Two details matter
and are routinely got wrong:

1. ANCESTOR PROPAGATION. GO is a DAG. A protein annotated with a specific
   term is implicitly annotated with every ancestor of that term. Both the
   prediction and the ground truth must be propagated before comparison, or
   precision is silently inflated.

2. Fmax IS PROTEIN-CENTRIC. Precision and recall are computed per protein at
   each threshold and THEN averaged — not pooled over all (protein, term)
   pairs. Pooling lets a few heavily-annotated proteins dominate. Precision
   is averaged only over proteins with at least one prediction at that
   threshold, which is the CAFA convention.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

NAMESPACES = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

# Terms so general they carry no information; excluded from evaluation.
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}


@dataclass
class GOGraph:
    parents: dict[str, set[str]] = field(default_factory=dict)
    namespace: dict[str, str] = field(default_factory=dict)
    alt_to_main: dict[str, str] = field(default_factory=dict)
    obsolete: set[str] = field(default_factory=set)
    _anc_cache: dict[str, set[str]] = field(default_factory=dict, repr=False)

    # -- construction ---------------------------------------------------
    @classmethod
    def from_obo(cls, path: str | Path) -> "GOGraph":
        g = cls()
        opener = gzip.open if str(path).endswith(".gz") else open
        term: dict | None = None

        def flush(t: dict | None) -> None:
            if not t or "id" not in t:
                return
            tid = t["id"]
            if t.get("is_obsolete"):
                g.obsolete.add(tid)
                return
            g.parents[tid] = t.get("parents", set())
            if "namespace" in t:
                g.namespace[tid] = NAMESPACES.get(t["namespace"], t["namespace"])
            for alt in t.get("alt_id", []):
                g.alt_to_main[alt] = tid

        with opener(path, "rt") as fh:
            in_term = False
            for line in fh:
                line = line.rstrip()
                if line == "[Term]":
                    flush(term)
                    term = {"parents": set(), "alt_id": []}
                    in_term = True
                    continue
                if line.startswith("[") and line != "[Term]":
                    flush(term)
                    term, in_term = None, False
                    continue
                if not in_term or term is None or not line:
                    continue

                if line.startswith("id: "):
                    term["id"] = line[4:]
                elif line.startswith("alt_id: "):
                    term["alt_id"].append(line[8:])
                elif line.startswith("namespace: "):
                    term["namespace"] = line[11:]
                elif line.startswith("is_a: "):
                    term["parents"].add(line[6:].split(" ! ")[0].strip())
                elif line.startswith("relationship: part_of "):
                    # part_of is transitive for annotation propagation;
                    # other relationship types are not and are ignored.
                    term["parents"].add(
                        line[len("relationship: part_of "):].split(" ! ")[0].strip()
                    )
                elif line.startswith("is_obsolete: true"):
                    term["is_obsolete"] = True
        flush(term)
        return g

    # -- queries --------------------------------------------------------
    def canonical(self, term: str) -> str | None:
        if term in self.parents:
            return term
        main = self.alt_to_main.get(term)
        return main if main in self.parents else None

    def ancestors(self, term: str) -> set[str]:
        """All ancestors of `term`, inclusive. Iterative to avoid deep
        recursion on the BP subgraph."""
        t = self.canonical(term)
        if t is None:
            return set()
        if t in self._anc_cache:
            return self._anc_cache[t]

        seen: set[str] = set()
        stack = [t]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            if cur in self._anc_cache:
                seen |= self._anc_cache[cur]
                continue
            stack.extend(self.parents.get(cur, ()))
        self._anc_cache[t] = seen
        return seen

    def propagate(self, terms, keep_roots: bool = False) -> set[str]:
        out: set[str] = set()
        for t in terms:
            out |= self.ancestors(t)
        return out if keep_roots else out - ROOTS

    def filter_ns(self, terms, ns: str) -> set[str]:
        return {t for t in terms if self.namespace.get(t) == ns}


def information_content(
    truth: dict[str, set[str]], graph: GOGraph
) -> dict[str, float]:
    """IC = -log2(freq) over the propagated ground truth. Used to weight
    terms so that predicting 'binding' counts for less than predicting a
    specific catalytic activity."""
    counts: dict[str, int] = defaultdict(int)
    for terms in truth.values():
        for t in terms:
            counts[t] += 1
    n = max(1, len(truth))
    return {t: -np.log2(c / n) for t, c in counts.items() if c > 0}


def fmax(
    truth: dict[str, set[str]],
    scored: dict[str, dict[str, float]],
    thresholds: np.ndarray | None = None,
) -> dict:
    """Protein-centric Fmax, plus the threshold and coverage that achieve it.

    truth  : protein -> propagated true terms
    scored : protein -> {term: score in [0,1]}, also propagated
    """
    if thresholds is None:
        thresholds = np.arange(0.01, 1.00, 0.01)

    proteins = [p for p in truth if truth[p]]
    if not proteins:
        return {"fmax": 0.0, "threshold": 0.0, "coverage": 0.0,
                "precision": 0.0, "recall": 0.0, "n": 0}

    best = {"fmax": 0.0, "threshold": 0.0, "coverage": 0.0,
            "precision": 0.0, "recall": 0.0, "n": len(proteins)}

    for tau in thresholds:
        precisions, recalls = [], []
        for p in proteins:
            pred = {t for t, s in scored.get(p, {}).items() if s >= tau}
            true = truth[p]
            if pred:
                precisions.append(len(pred & true) / len(pred))
            recalls.append(len(pred & true) / len(true))

        if not precisions:
            continue
        pr = float(np.mean(precisions))     # over predicted proteins only
        rc = float(np.mean(recalls))        # over all proteins
        if pr + rc == 0:
            continue
        f = 2 * pr * rc / (pr + rc)
        if f > best["fmax"]:
            best = {"fmax": f, "threshold": float(tau),
                    "coverage": len(precisions) / len(proteins),
                    "precision": pr, "recall": rc, "n": len(proteins)}
    return best


def smin(
    truth: dict[str, set[str]],
    scored: dict[str, dict[str, float]],
    ic: dict[str, float],
    thresholds: np.ndarray | None = None,
) -> dict:
    """Semantic distance S-min (lower is better). Complements Fmax by
    weighting errors on specific terms more heavily than general ones."""
    if thresholds is None:
        thresholds = np.arange(0.01, 1.00, 0.01)
    proteins = [p for p in truth if truth[p]]
    best = {"smin": float("inf"), "threshold": 0.0}

    for tau in thresholds:
        ru, mi = [], []
        for p in proteins:
            pred = {t for t, s in scored.get(p, {}).items() if s >= tau}
            true = truth[p]
            ru.append(sum(ic.get(t, 0.0) for t in true - pred))   # remaining uncert.
            mi.append(sum(ic.get(t, 0.0) for t in pred - true))   # misinformation
        if not ru:
            continue
        s = float(np.sqrt(np.mean(ru) ** 2 + np.mean(mi) ** 2))
        if s < best["smin"]:
            best = {"smin": s, "threshold": float(tau)}
    return best


def precision_at_recall(
    truth: dict[str, set[str]],
    scored: dict[str, dict[str, float]],
    target_recall: float = 0.5,
) -> float:
    """Precision at the lowest threshold reaching `target_recall`.

    More interpretable than Fmax for our purpose: it answers "if we want to
    annotate half the true function, how often are we right?"
    """
    proteins = [p for p in truth if truth[p]]
    best_p = 0.0
    for tau in np.arange(0.99, 0.0, -0.01):
        precisions, recalls = [], []
        for p in proteins:
            pred = {t for t, s in scored.get(p, {}).items() if s >= tau}
            true = truth[p]
            if pred:
                precisions.append(len(pred & true) / len(pred))
            recalls.append(len(pred & true) / len(true))
        if not precisions:
            continue
        if float(np.mean(recalls)) >= target_recall:
            return float(np.mean(precisions))
        best_p = float(np.mean(precisions))
    return best_p


def evaluate(
    truth: dict[str, set[str]],
    scored: dict[str, dict[str, float]],
    graph: GOGraph,
    by_namespace: bool = True,
) -> dict:
    """Full metric panel, overall and per GO namespace."""
    ic = information_content(truth, graph)
    out = {
        "overall": {
            **fmax(truth, scored),
            **smin(truth, scored, ic),
            "precision_at_recall_0.5": precision_at_recall(truth, scored, 0.5),
            "n_proteins": len(truth),
        }
    }
    if by_namespace:
        for ns in ("BP", "MF", "CC"):
            t_ns = {p: graph.filter_ns(v, ns) for p, v in truth.items()}
            t_ns = {p: v for p, v in t_ns.items() if v}
            if not t_ns:
                continue
            s_ns = {
                p: {t: s for t, s in scored.get(p, {}).items()
                    if graph.namespace.get(t) == ns}
                for p in t_ns
            }
            out[ns] = {**fmax(t_ns, s_ns),
                       "precision_at_recall_0.5": precision_at_recall(t_ns, s_ns, 0.5)}
    return out
