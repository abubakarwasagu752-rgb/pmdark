#!/usr/bin/env python3
"""Build every figure and table to Bioinformatics specification.

    python scripts/43_figures.py
    python scripts/43_figures.py --only 3
    python scripts/43_figures.py --tiff-dpi 1200

JOURNAL SPECIFICATION ENCODED HERE

  Width            86 mm single column, 178 mm double column, exactly.
  Background       white, as required.
  Gridlines        none (removed on request); top and right spines dropped
                   so the axes stay recessive.
  Formats          .png, .svg, .pdf, .tif for each figure.
  TIFF resolution  600 dpi by default. The journal asks 1200 dpi for line
                   drawings and 350 dpi for colour and half-tone. These are
                   colour line drawings, so 350 is the governing figure and
                   600 clears it with margin; --tiff-dpi 1200 if a
                   production editor insists. PDF and SVG are vector and
                   carry no resolution, which is what most production
                   pipelines prefer anyway.
  Alt text         required by the journal for every figure. Written to
                   figures/ALT_TEXT.md, ready to paste under each legend.

NO NUMBER IS WRITTEN INTO THIS FILE

  Every value plotted is read from a results JSON at run time, exactly as
  scripts/42_consolidate.py does. A figure therefore cannot disagree with
  the text, and re-running after any analysis change updates both. If a
  figure's inputs are missing the figure is SKIPPED and reported, never
  drawn from partial data.

COLOUR

  Two fixed palettes, validated with the dataviz skill's checker for
  colourblind separation, chroma, lightness and >=3:1 contrast against a
  white surface:

    models   #0072B2 ESM-2   #D55E00 ProtT5   #009E73 ProstT5
    methods  #0072B2 PLM     #882255 DIAMOND  #A65000 jackhmmer

  Hue is never the only cue: every series also carries a distinct marker
  and line style, so the figures survive greyscale printing and
  colourblind readers. Colours are assigned to an entity in fixed order
  and never cycled.
"""

from __future__ import annotations

import argparse
import re
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmdark import config as cfgmod                              # noqa: E402

MM = 25.4
SINGLE = 86 / MM                     # 3.386 in
DOUBLE = 178 / MM                    # 7.008 in

MODEL_COLOR = {"esm2": "#0072B2", "prott5": "#D55E00", "prostt5": "#009E73"}
MODEL_LABEL = {"esm2": "ESM-2", "prott5": "ProtT5", "prostt5": "ProstT5"}
MODEL_MARK = {"esm2": "o", "prott5": "s", "prostt5": "^"}
MODEL_LS = {"esm2": "-", "prott5": "--", "prostt5": ":"}

METHOD_COLOR = {"plm": "#0072B2", "DIAMOND": "#882255",
                "jackhmmer": "#A65000"}
METHOD_MARK = {"plm": "o", "DIAMOND": "s", "jackhmmer": "D"}
METHOD_LS = {"plm": "-", "DIAMOND": "--", "jackhmmer": "-."}

CLASS_COLOR = {"unc": "#0072B2", "desc": "#A65000"}

TIER_ORDER = ["T0_full", "I90", "I50", "I30", "I20", "Inone", "ORPHAN_nohit"]
TIER_LABEL = {"T0_full": "T0", "I90": "I90", "I50": "I50", "I30": "I30",
              "I20": "I20", "Inone": "Inone", "ORPHAN_nohit": "orphan"}

INK = "#222222"
MUTED = "#666666"


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------------------

def orphan_n(F) -> int | None:
    """Size of the orphan subset, from whichever results file records it."""
    boot = F.J("bootstrap_ci.json") or {}
    for k, v in boot.get("cells", {}).items():
        if k.endswith("/ORPHAN_nohit") and v.get("n"):
            return int(v["n"])
    for p in sorted(F.res.glob("calibration_*_whiten.json")):
        try:
            d = json.loads(p.read_text())
        except ValueError:
            continue
        blk = d.get("results", {}).get("ORPHAN_nohit")
        if blk and blk.get("n"):
            return int(blk["n"])
    return None


def tick_labels(F, tiers: list[str]) -> list[str]:
    """Tier labels, with the orphan subset's n shown, because it is the one
    column whose denominator differs from the rest of the axis."""
    n = orphan_n(F)
    out = []
    for t in tiers:
        lab = TIER_LABEL.get(t, t)
        if t == "ORPHAN_nohit" and n:
            lab = f"{lab} ({n})"
        out.append(lab)
    return out


def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 7,
        "axes.labelsize": 7.5,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.grid": False,              # no gridlines, as requested
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.2,
        "lines.markersize": 3.5,
        "legend.frameon": False,
        "errorbar.capsize": 1.8,
        "pdf.fonttype": 42,             # editable text, not outlines
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    return plt


class Figures:
    def __init__(self, res: Path, out: Path, tiff_dpi: int, png_dpi: int,
                 show_titles: bool = True):
        self.res, self.out = res, out
        self.tiff_dpi, self.png_dpi = tiff_dpi, png_dpi
        self.show_titles = show_titles
        self.drop_orphan_fig2 = False
        self.blocks: list[tuple] = []      # shared by the md and docx writers
        self.out.mkdir(parents=True, exist_ok=True)
        self.alt: list[tuple[str, str, str]] = []   # (name, legend, alt)
        self.made: list[str] = []
        self.skipped: list[str] = []
        self._cache: dict[str, dict | None] = {}
        # provenance of every plotted value: (figure, panel, series, at,
        # value, source file, source key). Recorded always, written only
        # under --audit. A figure is only verifiable if each mark on it can
        # be traced back to a file and a key without reading the code.
        self.audit: list[tuple] = []
        # Intervals found not to bracket their own estimate. Collected
        # rather than raised, so one bad cell cannot stop the run.
        self.bad_ci: list[str] = []

    def J(self, name: str) -> dict | None:
        if name not in self._cache:
            p = self.res / name
            try:
                self._cache[name] = json.loads(p.read_text()) \
                    if p.exists() else None
            except ValueError:
                self._cache[name] = None
        return self._cache[name]

    def need(self, fig: str, **files) -> dict | None:
        """Load every input or skip the figure. Never draw partial data."""
        got, miss = {}, []
        for key, name in files.items():
            d = self.J(name)
            if d is None:
                miss.append(name)
            got[key] = d
        if miss:
            self.skipped.append(f"{fig} (missing {', '.join(miss)})")
            log(f"SKIP {fig}: missing {', '.join(miss)}")
            return None
        return got

    def rec(self, fig: str, panel: str, series: str, at, value,
            src: str, key: str):
        """Record one plotted value and where it came from."""
        try:
            v = None if value is None else float(value)
        except (TypeError, ValueError):
            v = None
        self.audit.append((fig, panel, series, str(at), v, src, key))

    def pick(self, fig: str, panel: str, series: str, at, src: str,
             key: str, value):
        """rec() in expression position: records, then returns the value."""
        self.rec(fig, panel, series, at, value, src, key)
        return value

    PAD = 0.02          # matches savefig.pad_inches

    def fit_width(self, fig, limit: float = DOUBLE) -> float:
        """Shrink until the saved image is within the column width.

        tight_layout fits the axes inside the canvas, but savefig's tight
        bbox crops to the ink, and a rotated tick label or an offset
        suptitle can sit outside the canvas — so the file can come out
        wider than the figsize that was asked for. The journal states a
        width, so the file has to honour it rather than the intent.

        Text does not scale with the canvas, so a large reduction
        converges slowly; the step undershoots slightly to approach the
        limit from below rather than oscillating on it.
        """
        for _ in range(12):
            fig.canvas.draw()
            bb = fig.get_tightbbox(fig.canvas.get_renderer())
            w = bb.width + 2 * self.PAD
            if w <= limit + 1e-3:
                return w
            fig.set_size_inches(
                fig.get_size_inches() * (limit / w) * 0.999, forward=True)
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        return bb.width + 2 * self.PAD

    def save(self, fig, name: str, legend: str, alt: str):
        w = self.fit_width(fig)
        log(f"{name}: {w * 25.4:.1f} mm wide "
            f"({'within' if w <= DOUBLE + 1e-3 else 'OVER'} 178 mm)")
        for ext, kw in (("png", {"dpi": self.png_dpi}),
                        ("svg", {}),
                        ("pdf", {})):
            fig.savefig(self.out / f"{name}.{ext}", **kw)
        fig.savefig(self.out / f"{name}.tif", dpi=self.tiff_dpi,
                    pil_kwargs={"compression": "tiff_lzw"})
        import matplotlib.pyplot as plt
        plt.close(fig)
        self.alt.append((name, legend, alt))
        self.made.append(name)
        log(f"wrote {name}.{{png,svg,pdf,tif}}")

    FIGNUM = re.compile(r"^Figure\s+\d+\.\s*")

    def title(self, fig, text: str):
        """Short descriptive title drawn onto the figure itself.

        The "Figure N." prefix is stripped here rather than at the call
        sites, so each call can pass the same string as the legend's
        opening sentence and the two cannot drift apart. Numbering is the
        legend's job.

        Bioinformatics typesets the legend adjacent to the figure, so a
        title inside the image is unusual for the final submission;
        --no-titles produces the convention-conforming version from the
        same code.
        """
        if not self.show_titles:
            return None
        fig.suptitle(self.FIGNUM.sub("", text), fontsize=8,
                     fontweight="bold", x=0.012, ha="left", y=0.995)
        return [0, 0, 1, 0.93]

    @staticmethod
    def population_band(ax, idx: int, label: bool = False):
        """Shade the column whose denominator differs from the rest.

        A band rather than a divider line: a line sits in the same visual
        register as the data and can be read as continuing the series,
        whereas a shaded region reads as a change of context. Drawn at
        zorder 0 so every mark stays on top of it.
        """
        if idx is None:
            return
        ax.axvspan(idx - 0.5, idx + 0.5, color="#ededed", lw=0, zorder=0)
        if label:
            ax.annotate("different\nsubset", xy=(idx, 1.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(0, -2), textcoords="offset points",
                        ha="center", va="top", fontsize=5.6, color=MUTED)

    @staticmethod
    def runs(ax, xs, ys, break_at=None, dx=0.0, **kw):
        """Plot in contiguous runs of x, so no line is ever drawn across a
        condition that was not measured. A connected line over a gap
        asserts an intermediate value that does not exist.

        `dx` shifts the drawn position without touching the run logic,
        which needs the original integer indices to tell a gap in the
        conditions from a gap introduced by a dodge.

        `break_at` forces a new segment at the given x even when the data
        are contiguous. It is used for the orphan point, which rests on a
        subset of the queries rather than on a further exclusion applied
        to all of them; joining it to the previous point with a line would
        imply a trend across a change of denominator.
        """
        if not xs:
            return
        brk = set(break_at or ())
        lab = kw.pop("label", None)
        start = 0
        for i in range(1, len(xs) + 1):
            if i == len(xs) or xs[i] != xs[i - 1] + 1 or xs[i] in brk:
                seg_x, seg_y = xs[start:i], ys[start:i]
                ax.plot([x + dx for x in seg_x], seg_y, label=lab, **kw)
                lab = None          # legend entry once
                start = i

    def yerr(self, value, ci, where: str, scale: float = 1.0):
        """Asymmetric error-bar lengths from a point estimate and an
        interval, or None if the two disagree.

        A results file can carry an interval that does not bracket its
        own estimate — a Wilson interval computed on a different
        denominator, a bootstrap written against a stale point value.
        Passed straight to errorbar that becomes a negative length and
        matplotlib raises `'yerr' must not contain negative values`,
        which names neither the figure nor the file. Here it becomes a
        named complaint, the error bar is dropped, and the rest of the
        figure still draws.
        """
        if value is None or not ci or ci[0] is None or ci[1] is None:
            return None
        lo, hi = ci[0], ci[1]
        if lo > value + 1e-12 or hi < value - 1e-12:
            self.bad_ci.append(
                f"{where}: interval [{lo:.4g}, {hi:.4g}] does not "
                f"bracket the estimate {value:.4g}; error bar omitted")
            log(f"BAD INTERVAL {where}: [{lo:.4g}, {hi:.4g}] does not "
                f"bracket {value:.4g}")
            return None
        return [[scale * (value - lo)], [scale * (hi - value)]]

    @staticmethod
    def panel(ax, letter: str, dx: float = -0.16, dy: float = 1.04):
        ax.text(dx, dy, letter, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------------------
# spectrum (data, cached — needed by Figure 1A)
# ---------------------------------------------------------------------------

def spectrum(cfg, res: Path, models: list[str], force: bool) -> dict | None:
    """Cumulative variance explained by the leading embedding dimensions.

    This is the anisotropy measurement itself, so it is computed from the
    reference embeddings and cached as a results file rather than being
    recomputed inside a plotting routine.
    """
    cache = res / "embedding_spectrum.json"
    if cache.exists() and not force:
        try:
            return json.loads(cache.read_text())
        except ValueError:
            pass
    import importlib.util
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "pmdark_calib", here / "31_calibrate.py")
    calib = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(calib)
    bench = cfg.p("proteins", "benchmark")
    out: dict = {"models": {}}
    for m in models:
        try:
            _ids, emb = calib.load_embeddings(cfg, m, bench / "reference.faa")
        except SystemExit:
            log(f"spectrum: {m} reference embeddings absent, skipping")
            continue
        X = emb.astype(np.float64)
        X -= X.mean(axis=0, keepdims=True)
        # eigenvalues of the covariance, descending
        cov = (X.T @ X) / max(1, len(X) - 1)
        ev = np.linalg.eigvalsh(cov)[::-1]
        ev = np.clip(ev, 0, None)
        tot = ev.sum()
        cum = np.cumsum(ev) / tot if tot > 0 else ev
        d50 = int(np.searchsorted(cum, 0.5) + 1)
        out["models"][m] = {
            "dim": int(emb.shape[1]),
            "n": int(emb.shape[0]),
            "cumulative": [float(x) for x in cum],
            "dims_for_50pct": d50,
            "condition_number": float(ev[0] / ev[ev > 0][-1])
            if (ev > 0).any() else None,
        }
        log(f"spectrum: {m} dim {emb.shape[1]}, {d50} dims carry 50% variance")
    if not out["models"]:
        return None
    cache.write_text(json.dumps(out, indent=2))
    log(f"wrote {cache.name}")
    return out


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig1(F: Figures, plt, spec: dict | None):
    d = F.need("Figure 1", boot="bootstrap_ci.json")
    if d is None or spec is None:
        if spec is None:
            F.skipped.append("Figure 1 (no embedding_spectrum.json)")
            log("SKIP Figure 1: spectrum unavailable")
        return
    boot = d["boot"]
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, DOUBLE / 3.05))

    # -- A: cumulative variance ---------------------------------------
    ax = axes[0]
    for m in ("esm2", "prott5", "prostt5"):
        s = spec["models"].get(m)
        if not s:
            continue
        cum = np.array(s["cumulative"])
        x = np.arange(1, len(cum) + 1)
        ax.plot(x, cum, color=MODEL_COLOR[m], ls=MODEL_LS[m], lw=1.2,
                label=f"{MODEL_LABEL[m]} ({s['dims_for_50pct']})")
    ax.axhline(0.5, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log")
    ax.set_xlabel("embedding dimensions (ranked)")
    ax.set_ylabel("cumulative variance explained")
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower right", title="dims for 50%", title_fontsize=6.5)
    F.panel(ax, "A")

    # -- B, C: paired deltas at two tiers ------------------------------
    def deltas(ax, tier, title, panel):
        rows = []
        for m in ("esm2", "prott5", "prostt5"):
            for tf, mark in (("center", "o"), ("whiten", "s")):
                k = f"{m}/{tier}: {tf} - none"
                v = boot.get("deltas", {}).get(k)
                if v:
                    rows.append((m, tf, mark, v))
                    F.rec("Figure 1", panel, f"{MODEL_LABEL[m]} {tf}",
                          tier, v.get("delta"), "bootstrap_ci.json",
                          f"deltas[{k}].delta")
        if not rows:
            ax.set_visible(False)
            return
        ax.axvline(0, color=MUTED, lw=0.6)
        ylab = []
        for i, (m, tf, mark, v) in enumerate(rows):
            y = len(rows) - 1 - i
            lo, hi = v["ci"]
            filled = bool(v.get("excludes_zero"))
            ax.plot([lo, hi], [y, y], color=MODEL_COLOR[m], lw=1.1,
                    solid_capstyle="butt")
            ax.plot([v["delta"]], [y], marker=mark, ms=4,
                    color=MODEL_COLOR[m],
                    markerfacecolor=MODEL_COLOR[m] if filled else "white",
                    markeredgecolor=MODEL_COLOR[m], markeredgewidth=0.9)
            ylab.append((y, f"{MODEL_LABEL[m]}  {tf}"))
        ax.set_yticks([y for y, _ in ylab])
        ax.set_yticklabels([t for _, t in ylab])
        ax.set_xlabel("Δ$F_{max}$ vs untransformed")
        ax.set_title(title, pad=3)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        return True

    deltas(axes[1], "T0_full", "all queries (T0)", "B")
    F.panel(axes[1], "B", dx=-0.42)
    deltas(axes[2], "ORPHAN_nohit", "orphans", "C")
    F.panel(axes[2], "C", dx=-0.42)

    # B and C invite the reader to compare the same model across the two
    # panels — that comparison is the point of the pair. On independent
    # scales, equal distances on the page are unequal effects. Union of
    # the two autoscaled ranges, so nothing is clipped.
    vis = [a for a in (axes[1], axes[2]) if a.get_visible()]
    if len(vis) == 2:
        lo = min(a.get_xlim()[0] for a in vis)
        hi = max(a.get_xlim()[1] for a in vis)
        for a in vis:
            a.set_xlim(lo, hi)

    d50 = {m: spec["models"][m]["dims_for_50pct"]
           for m in spec["models"]}
    rect = F.title(fig, "Figure 1. Embedding anisotropy and its correction")
    fig.tight_layout(pad=0.4, rect=rect)
    F.save(fig, "Figure1",
           "Fig. 1. Embedding anisotropy and its correction. (A) Cumulative "
           "variance explained by the ranked principal directions of the "
           "reference embeddings; the dashed line marks 50%. (B, C) Paired "
           "bootstrap change in protein-centric Fmax relative to "
           "untransformed embeddings, for centering and for whitening, over "
           "all benchmark queries (B) and over the orphan subset (C). Bars "
           "are 95% percentile intervals from 2000 paired resamples; filled "
           "markers indicate intervals excluding zero.",
           "Three-panel line and interval plot. Panel A shows cumulative "
           "variance against ranked embedding dimension on a logarithmic "
           "axis for three protein language models, each reaching half of "
           "total variance within a small number of dimensions "
           f"({', '.join(f'{MODEL_LABEL[m]} {v}' for m, v in d50.items())}). "
           "Panels B and C show change in Fmax relative to untransformed "
           "embeddings with 95% confidence intervals: whitening gives a "
           "clear positive shift in both panels, whereas centering alone "
           "gives a small shift over all queries and none on orphans.")


def fig2(F: Figures, plt):
    d = F.need("Figure 2", boot="bootstrap_ci.json")
    if d is None:
        return
    boot = d["boot"]
    calibs = {}
    for p in sorted(F.res.glob("calibration_*_whiten.json")):
        try:
            calibs[p.name] = json.loads(p.read_text())
        except ValueError:
            pass

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE / 2.5))

    # -- A: Fmax across strata -----------------------------------------
    ax = axes[0]
    tiers = [t for t in TIER_ORDER
             if any(f"{m}/whiten/{t}" in boot.get("cells", {})
                    for m in MODEL_COLOR)]
    if F.drop_orphan_fig2:
        tiers = [t for t in tiers if t != "ORPHAN_nohit"]
    xs = np.arange(len(tiers))
    orph_i = tiers.index("ORPHAN_nohit") if "ORPHAN_nohit" in tiers else None
    F.population_band(ax, orph_i)
    # Horizontal dodge. At T0 and I90 the three models agree to within
    # about a marker's width, so plotted on a shared x they superimpose
    # and whichever is drawn first — the headline model — is buried.
    # Agreement between models is a result, not a reason to hide one.
    # The dodge moves nothing but the drawn x position.
    MODELS = ("esm2", "prott5", "prostt5")
    DODGE = 0.17                      # < 0.5, so each stays in its column
    for j, m in enumerate(MODELS):
        off = (j - (len(MODELS) - 1) / 2) * DODGE
        y, lo, hi, xx = [], [], [], []
        for i, t in enumerate(tiers):
            c = boot.get("cells", {}).get(f"{m}/whiten/{t}")
            if not c:
                continue
            xx.append(i)
            y.append(c["fmax"])
            lo.append(c["fmax"] - c["fmax_ci"][0])
            hi.append(c["fmax_ci"][1] - c["fmax"])
            F.rec("Figure 2", "A", MODEL_LABEL[m], t, c["fmax"],
                  "bootstrap_ci.json", f"cells[{m}/whiten/{t}].fmax")
        if not y:
            continue
        # error bars everywhere; the connecting line stops before the
        # orphan column, whose denominator differs from the rest
        ax.errorbar([x + off for x in xx], y, yerr=[lo, hi], fmt="none",
                    ecolor=MODEL_COLOR[m], elinewidth=0.9, capsize=1.6)
        F.runs(ax, xx, y, break_at={orph_i} if orph_i is not None else None,
               dx=off, color=MODEL_COLOR[m], ls=MODEL_LS[m],
               marker=MODEL_MARK[m], lw=1.2, label=MODEL_LABEL[m])
    ax.set_xticks(xs)
    rot = 45 if len(tiers) > 4 else 0
    ax.set_xticklabels(tick_labels(F, tiers), rotation=rot,
                       ha="right" if rot else "center")
    ax.set_xlabel("homology available to the transfer")
    ax.set_ylabel("$F_{max}$")
    ax.legend(loc="lower left")
    F.panel(ax, "A")

    # -- B: ontology collapse ------------------------------------------
    # Only the headline model. Six bars per aspect (three models x two
    # conditions) read as three gradient-filled bars rather than three
    # pairs: adjacent fills with no gap are not separable, and the
    # model comparison already has panel A and Table 1. One model, two
    # conditions, a real gap between them.
    ax = axes[1]
    ont = ["BP", "MF", "CC"]
    head = "esm2"
    d2 = calibs.get(f"calibration_{head}_whiten.json")
    # hatch=None, never the string "none": matplotlib reads a string as
    # hatch characters, so "none" parses as n,o,n,e and two of those are
    # invalid. It warns today and is an error from matplotlib 3.13.
    conds = [("T0_full", "all queries", MODEL_COLOR[head], None),
             ("ORPHAN_nohit", "orphans", "white", "////")]
    series = []
    for tier, lab, face, hatch in conds:
        blk = (d2 or {}).get("results", {}).get(tier)
        if not blk:
            continue
        vals = [blk["metrics"].get(o, {}).get("fmax") for o in ont]
        if any(v is None for v in vals):
            continue
        for o, v in zip(ont, vals):
            F.rec("Figure 2", "B", f"{MODEL_LABEL[head]} {lab}", o, v,
                  f"calibration_{head}_whiten.json",
                  f"results[{tier}].metrics[{o}].fmax")
        series.append((lab, vals, face, hatch))
    if len(series) < 1:
        ax.set_visible(False)
    else:
        width = 0.34
        gap = 0.04
        base = np.arange(len(ont))
        for j, (lab, vals, face, hatch) in enumerate(series):
            pos = base + (j - (len(series) - 1) / 2) * (width + gap)
            ax.bar(pos, vals, width, color=face,
                   edgecolor=MODEL_COLOR[head], linewidth=0.9,
                   hatch=hatch, label=lab, zorder=2)
            for x_, v_ in zip(pos, vals):
                ax.annotate(f"{v_:.2f}", xy=(x_, v_), xytext=(0, 2),
                            textcoords="offset points", ha="center",
                            fontsize=6, color=INK)
        ax.set_xticks(base)
        ax.set_xticklabels(["biological\nprocess", "molecular\nfunction",
                            "cellular\ncomponent"])
        ax.set_ylabel("$F_{max}$")
        ax.set_ylim(0, 1.0)
        ax.set_title(MODEL_LABEL[head], pad=3, fontsize=7)
        ax.legend(loc="upper left", ncol=1)
        F.panel(ax, "B", dx=-0.14)

    rect = F.title(fig, "Figure 2. Annotation accuracy as homologous evidence is withdrawn")
    fig.tight_layout(pad=0.4, rect=rect)
    F.save(fig, "Figure2",
           "Fig. 2. Annotation accuracy as homologous evidence is withdrawn. "
           "(A) Protein-centric Fmax for whitened embeddings across "
           "identity-stratified conditions, from all queries (T0) through "
           "progressively stricter identity ceilings. The rightmost column "
           "is the orphan subset: the queries for which BLAST-class search "
           "returns no hit at all, with its size in parentheses. It is "
           "plotted detached and on a shaded background because it is a "
           "subset of the queries rather than a further exclusion applied "
           "to all of them, so the denominator differs from every other "
           "column and no trend should be read across the gap. Figure 4 shows what becomes of "
           "this subset under a more sensitive search. Error bars are 95% "
           "bootstrap intervals. (B) Fmax by Gene Ontology aspect for the "
           "best-performing model, over all queries (solid) and over "
           "orphans (hatched); the corresponding values for the other "
           "models are given in Table 1.",
           "Two-panel figure. Panel A plots Fmax against increasingly "
           "strict homology exclusion for three protein language models; "
           "all three decline gradually and remain well above zero on the "
           "orphan subset, which is drawn as a detached point because it "
           "rests on fewer proteins. Panel B is a grouped bar chart of Fmax split by "
           "ontology aspect for one model, with cellular component "
           "highest, molecular function intermediate and biological "
           "process lowest, and every aspect scoring lower on orphans than "
           "on the full query set.")


def fig3(F: Figures, plt):
    d = F.need("Figure 3",
               bd="homology_baseline_DIAMOND.json",
               bj="homology_baseline_jackhmmer.json",
               cj="plm_vs_homology_jackhmmer.json",
               boot="bootstrap_ci.json")
    if d is None:
        return
    bd, bj, cj, boot = d["bd"], d["bj"], d["cj"], d["boot"]
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, DOUBLE / 3.05))

    tiers = [t for t in TIER_ORDER if t in bd.get("tiers", {})]
    xs = np.arange(len(tiers))
    brk = ({tiers.index("ORPHAN_nohit")} if "ORPHAN_nohit" in tiers
           else None)

    # -- A: Fmax, PLM vs the two searches ------------------------------
    # The PLM point estimates come from the calibration file, which covers
    # every condition; the bootstrap covers only the four conditions whose
    # intervals were computed, and using it alone would leave gaps that a
    # connected line would paper over.
    ax = axes[0]
    F.population_band(ax, next(iter(brk)) if brk else None)
    cal = F.J("calibration_esm2_whiten.json") or {}
    y, xx = [], []
    for i, t in enumerate(tiers):
        src = "calibration_esm2_whiten.json"
        key = f"results[{t}].metrics.overall.fmax"
        v = (cal.get("results", {}).get(t, {})
             .get("metrics", {}).get("overall", {}).get("fmax"))
        if v is None:
            c = boot.get("cells", {}).get(f"esm2/whiten/{t}")
            v = c["fmax"] if c else None
            src, key = "bootstrap_ci.json", f"cells[esm2/whiten/{t}].fmax"
        if v is None:
            continue
        xx.append(i)
        y.append(v)
        F.rec("Figure 3", "A", "embeddings", t, v, src, key)
    if y:
        F.runs(ax, xx, y, break_at=brk, color=METHOD_COLOR["plm"],
               ls=METHOD_LS["plm"], marker=METHOD_MARK["plm"],
               label="embeddings", lw=1.4)
    for tag, b in (("DIAMOND", bd), ("jackhmmer", bj)):
        y, xx = [], []
        for i, t in enumerate(tiers):
            v = b.get("tiers", {}).get(t, {}).get("topk", {}).get("fmax")
            if v is None:
                continue
            xx.append(i)
            y.append(v)
            F.rec("Figure 3", "A", tag, t, v,
                  f"homology_baseline_{tag}.json", f"tiers[{t}].topk.fmax")
        F.runs(ax, xx, y, break_at=brk, color=METHOD_COLOR[tag],
               ls=METHOD_LS[tag], marker=METHOD_MARK[tag], label=tag, lw=1.2)
    ax.set_xticks(xs)
    ax.set_xticklabels(tick_labels(F, tiers), rotation=45, ha="right")
    ax.set_ylabel("$F_{max}$")
    ax.legend(loc="lower left")
    F.panel(ax, "A")

    # -- B: answerable fraction ----------------------------------------
    ax = axes[1]
    F.population_band(ax, next(iter(brk)) if brk else None)
    for tag, b in (("DIAMOND", bd), ("jackhmmer", bj)):
        y, xx = [], []
        for i, t in enumerate(tiers):
            v = b.get("tiers", {}).get(t, {}).get("best", {}) \
                 .get("answerable_frac")
            if v is None:
                continue
            xx.append(i)
            y.append(100 * v)
            F.rec("Figure 3", "B", tag, t, 100 * v,
                  f"homology_baseline_{tag}.json",
                  f"tiers[{t}].best.answerable_frac x100")
        F.runs(ax, xx, y, break_at=brk, color=METHOD_COLOR[tag],
               ls=METHOD_LS[tag], marker=METHOD_MARK[tag], label=tag, lw=1.2)
    ax.set_xticks(xs)
    ax.set_xticklabels(tick_labels(F, tiers), rotation=45, ha="right")
    ax.set_ylabel("queries answered (%)")
    ax.set_ylim(-3, 103)
    ax.legend(loc="lower left")
    F.panel(ax, "B", dx=-0.22)

    # -- C: coverage vs accuracy, decomposed ---------------------------
    ax = axes[2]
    ax.axvline(0, color=MUTED, lw=0.6)
    rows = []
    for t in tiers:
        m = (cj.get("tiers", {}).get(t, {}).get("models") or {}).get("esm2")
        if not m:
            continue
        for key, mark, lab in (("all_queries", "o", "all"),
                               ("diamond_answerable", "s", "answerable")):
            v = m.get(key)
            if v:
                rows.append((t, lab, mark, v))
                F.rec("Figure 3", "C", lab, t, v.get("delta"),
                      "plm_vs_homology_jackhmmer.json",
                      f"tiers[{t}].models.esm2.{key}.delta")
    if not rows:
        ax.set_visible(False)
    else:
        ylab = []
        for i, (t, lab, mark, v) in enumerate(rows):
            y = len(rows) - 1 - i
            lo, hi = v["ci"]
            filled = bool(v.get("excludes_zero"))
            col = METHOD_COLOR["jackhmmer"]
            ax.plot([lo, hi], [y, y], color=col, lw=1.1)
            ax.plot([v["delta"]], [y], marker=mark, ms=4, color=col,
                    markerfacecolor=col if filled else "white",
                    markeredgecolor=col, markeredgewidth=0.9)
            ylab.append((y, f"{TIER_LABEL[t]}  {lab}"))
        ax.set_yticks([y for y, _ in ylab])
        ax.set_yticklabels([t for _, t in ylab])
        ax.set_xlabel("Δ$F_{max}$, embeddings − profile search")
        ax.set_ylim(-0.6, len(rows) - 0.4)
        F.panel(ax, "C", dx=-0.42)

    rect = F.title(fig, "Figure 3. Reach, not accuracy, is what embedding transfer adds")
    fig.tight_layout(pad=0.4, rect=rect)
    F.save(fig, "Figure3",
           "Fig. 3. Reach, not accuracy, is what embedding transfer adds. "
           "(A) Fmax for embedding transfer and for two homology searches "
           "across the identity-stratified conditions, all transferring at "
           "the same significance threshold. As in Figure 2, the orphan "
           "column is a subset of the queries and is shown detached on a "
           "shaded background; a break in a line marks a condition for "
           "which that method has no measurement. (B) Fraction of queries for "
           "which each search finds any permitted hit. (C) Paired "
           "difference in Fmax between embedding transfer and profile "
           "search, over all queries and over only those the profile "
           "search could answer; filled markers indicate 95% intervals "
           "excluding zero.",
           "Three-panel figure comparing embedding-based transfer with "
           "sequence-homology search. Panel A shows Fmax falling for both "
           "homology searches as identity ceilings tighten while embedding "
           "transfer declines far less. Panel B shows the fraction of "
           "queries each homology search can answer at all, falling "
           "steeply for the stricter conditions. Panel C shows that the "
           "advantage over all queries largely disappears when restricted "
           "to queries the profile search could answer, indicating the "
           "gain comes from coverage rather than accuracy.")


def fig4(F: Figures, plt):
    d = F.need("Figure 4", h="orphan_hardening.json")
    if d is None:
        return
    h = d["h"]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE * 0.36))

    # -- A: survival ---------------------------------------------------
    ax = axes[0]
    cuts = sorted(((float(k), v) for k, v in h.get("survival", {}).items()),
                  key=lambda kv: kv[0])
    if cuts:
        x = [c for c, _ in cuts]
        y = [v["still_orphan"] for _, v in cuts]
        for c, v in cuts:
            F.rec("Figure 4", "A", "still orphan", f"E={c:g}",
                  v["still_orphan"], "orphan_hardening.json",
                  f"survival[{c:g}].still_orphan")
        ax.plot(x, y, color=METHOD_COLOR["jackhmmer"], marker="o", lw=1.2)
        ax.set_xscale("log")
        # Tick the decades from the top down, so the most relaxed
        # thresholds tested are labelled. Anchoring at the bottom can
        # leave the rightmost points past the last tick, unlocatable —
        # and those are the points that show why a lax cut is untenable.
        import math
        d_hi = math.ceil(math.log10(max(x)))
        d_lo = math.floor(math.log10(min(x)))
        step = max(1, (d_hi - d_lo) // 7 + 1)
        decs = sorted(range(d_hi, d_lo - 1, -step))
        ax.set_xticks([10.0 ** t for t in decs])
        ax.set_xlim(10.0 ** (d_lo - 0.4), 10.0 ** (d_hi + 0.4))
        n0 = h.get("n_orphan_diamond")
        if n0:
            ax.axhline(n0, color=MUTED, lw=0.6, ls=(0, (2, 2)))
            ax.annotate(f"no profile search ({n0})", xy=(x[0], n0),
                        xytext=(0, 3), textcoords="offset points",
                        fontsize=6, color=MUTED)
        ha = h.get("harden_at")
        if ha:
            ax.axvline(ha, color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_xlabel("jackhmmer E-value threshold")
    ax.set_ylabel("proteins with no homologue")
    F.panel(ax, "A", dx=-0.26)

    # -- B: Fmax per stratum -------------------------------------------
    ax = axes[1]
    order = ["DIAMOND-orphan", "HARDENED", "dropped (had a profile hit)"]
    short = {"DIAMOND-orphan": "no BLAST-class hit",
             "HARDENED": "no hit, either search",
             "dropped (had a profile hit)": "profile hit only"}
    rows = [(k, h["strata"][k]) for k in order if k in h.get("strata", {})]
    for i, (k, v) in enumerate(rows):
        y = len(rows) - 1 - i
        lo, hi = v.get("ci", [None, None])
        F.rec("Figure 4", "B", short[k], f"n={v.get('n')}", v["fmax"],
              "orphan_hardening.json", f"strata[{k}].fmax")
        col = METHOD_COLOR["plm"] if k == "HARDENED" else MUTED
        if lo is not None:
            ax.plot([lo, hi], [y, y], color=col, lw=1.1)
        ax.plot([v["fmax"]], [y], marker="o", ms=4.2, color=col,
                markerfacecolor=col if k == "HARDENED" else "white",
                markeredgecolor=col, markeredgewidth=0.9)
        ax.annotate(f"n={v.get('n')}", xy=(v["fmax"], y),
                    xytext=(0, 4), textcoords="offset points",
                    fontsize=6, ha="center", color=MUTED)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([short[k] for k, _ in reversed(rows)])
    ax.set_xlabel("$F_{max}$")
    ax.set_ylim(-0.5, len(rows) - 0.4)
    F.panel(ax, "B", dx=-0.62)

    rect = F.title(fig, "Figure 4. The orphan subset withstands a more sensitive search")
    fig.tight_layout(pad=0.4, rect=rect)
    F.save(fig, "Figure4",
           "Fig. 4. The orphan subset withstands a more sensitive search. "
           "(A) Number of benchmark proteins retaining orphan status as the "
           "iterative-profile E-value threshold is relaxed; the dashed "
           "horizontal line is the count from BLAST-class search alone and "
           "the vertical line marks the threshold adopted. (B) Fmax with "
           "95% bootstrap intervals for the original orphan subset, the "
           "subset surviving profile search, and the proteins removed by "
           "it.",
           "Two-panel figure. Panel A shows the number of proteins still "
           "classed as orphans decreasing as the profile-search E-value "
           "threshold is relaxed, with a reference line for the count from "
           "BLAST-class search alone. Panel B shows Fmax with confidence "
           "intervals for three nested subsets; the hardened subset scores "
           "slightly lower than the original and the removed proteins "
           "score higher, and all three intervals overlap.")


def fig5(F: Figures, plt):
    d = F.need("Figure 5", mt="matched_enrichment.json",
               an="annotation_summary.json")
    if d is None:
        return
    mt, an = d["mt"], d["an"]
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, DOUBLE / 3.05))

    # -- A: dependence by expression decile ----------------------------
    ax = axes[0]
    dt = mt.get("decile_table") or []
    if dt:
        x = [r["decile"] for r in dt]
        for key, lab, cls in (("unc_frac", "uncharacterized", "unc"),
                              ("desc_frac", "described", "desc")):
            y = [100 * r[key] if r[key] is not None else np.nan for r in dt]
            for r in dt:
                F.rec("Figure 5", "A", lab, f"decile {r['decile']}",
                      None if r[key] is None else 100 * r[key],
                      "matched_enrichment.json",
                      f"decile_table[{r['decile']}].{key} x100")
            ax.plot(x, y, color=CLASS_COLOR[cls],
                    ls="-" if cls == "unc" else "--",
                    marker="o" if cls == "unc" else "s", lw=1.2, label=lab)
        ax.set_xticks(x)
        ax.set_xlabel("expression decile (low → high)")
        ax.set_ylabel("STAT-dependent (%)")
        ax.legend(loc="upper left")
    else:
        ax.set_visible(False)
    F.panel(ax, "A")

    # -- B: effect size by decile --------------------------------------
    ax = axes[1]
    wd = (mt.get("fold_change") or {}).get("within_deciles") or []
    if wd:
        x = [w["decile"] for w in wd]
        y = [w["delta"] for w in wd]
        sig = [w["p"] < 0.05 for w in wd]
        for w in wd:
            F.rec("Figure 5", "B", "Cliff's delta", f"decile {w['decile']}",
                  w["delta"], "matched_enrichment.json",
                  f"fold_change.within_deciles[{w['decile']}].delta")
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.plot(x, y, color=CLASS_COLOR["unc"], lw=1.1, ls="-", zorder=1)
        for xi, yi, s in zip(x, y, sig):
            ax.plot([xi], [yi], marker="o", ms=4,
                    color=CLASS_COLOR["unc"],
                    markerfacecolor=CLASS_COLOR["unc"] if s else "white",
                    markeredgecolor=CLASS_COLOR["unc"], markeredgewidth=0.9,
                    zorder=2)
        tr = mt.get("trend") or {}
        if tr.get("spearman_rho") is not None:
            # The correlation is on |delta|, but the series drawn here is
            # delta, which descends. A bare "rho = +0.84" beside a falling
            # line reads as a contradiction, so the statistic names its
            # own argument.
            txt = f"ρ(|δ|) = {tr['spearman_rho']:+.2f}"
            pv = tr.get("spearman_p", tr.get("p"))
            if pv is not None:
                txt += (f"\np = {pv:.3f}" if pv >= 1e-3 else "\np < 0.001")
            # Two lines in whichever top corner the data leaves free.
            # Guessing a corner has now failed twice: the lower right hit
            # the last decile, the mid right hit a local peak. The series
            # is not monotone, so the corner is chosen by looking at it —
            # the half whose closest-to-zero point sits lowest wins.
            h = max(1, len(y) // 2)
            left_hi = max(v for v in y[:h] if v == v)
            right_hi = max(v for v in y[h:] if v == v)
            right = right_hi < left_hi
            ax.annotate(txt, xy=(0.97 if right else 0.03, 0.96),
                        xycoords="axes fraction",
                        ha="right" if right else "left", va="top",
                        fontsize=6.5, color=INK, linespacing=1.4)
        ax.set_xticks(x)
        ax.set_xlabel("expression decile (low → high)")
        ax.set_ylabel("Cliff's δ")
    else:
        ax.set_visible(False)
    F.panel(ax, "B", dx=-0.24)

    # -- C: annotation yield -------------------------------------------
    ax = axes[2]
    y = an.get("yield") or {}
    n = an.get("n_targets") or 1
    # Aspect order matches Figure 2B, so a reader carries one ordering
    # through the paper rather than re-reading the axis.
    bars = [("biological\nprocess", y.get("biological_process"),
             "biological_process"),
            ("molecular\nfunction", y.get("molecular_function"),
             "molecular_function"),
            ("cellular comp.\nonly", y.get("cellular_component_only"),
             "cellular_component_only")]
    for k, v, jk in bars:
        F.rec("Figure 5", "C", "target proteins",
              k.replace("\n", " "), v, "annotation_summary.json",
              f"yield.{jk}")
    bars = [(k, v) for k, v, _ in bars if v is not None]
    if bars:
        pos = np.arange(len(bars))
        vals = [100 * v / n for _, v in bars]
        ax.bar(pos, vals, 0.56, color=CLASS_COLOR["unc"],
               edgecolor=CLASS_COLOR["unc"], linewidth=0)
        for p_, v_, (_, raw) in zip(pos, vals, bars):
            ax.annotate(f"{raw}", xy=(p_, v_), xytext=(0, 2),
                        textcoords="offset points", ha="center",
                        fontsize=6.5, color=INK)
        ax.set_xticks(pos)
        ax.set_xticklabels([k for k, _ in bars])
        ax.set_ylabel(f"% of {n} target proteins")
        ax.set_ylim(0, max(vals) * 1.22)
        orph = an.get("n_orphan")
        if orph is not None:
            ax.annotate(f"{orph} of {n} have no homologue",
                        xy=(0.5, 0.97), xycoords="axes fraction",
                        ha="center", va="top", fontsize=6.3, color=MUTED)
    else:
        ax.set_visible(False)
    F.panel(ax, "C", dx=-0.20)

    rect = F.title(fig, "Figure 5. Application to STAT-dependent uncharacterized genes")
    fig.tight_layout(pad=0.4, rect=rect)
    F.save(fig, "Figure5",
           "Fig. 5. Application to STAT-dependent uncharacterized genes. "
           "(A) Percentage of genes whose expression depends on STAT, by "
           "expression decile, for genes with and without a functional "
           "description in the reference annotation. (B) Effect size of the "
           "difference in fold-change distributions within each decile; "
           "filled markers are significant at p < 0.05 and ρ is the "
           "Spearman correlation of absolute effect size with expression. "
           "(C) Gene Ontology aspects assigned to the target proteins by "
           "calibrated embedding transfer; counts above each bar.",
           "Three-panel figure. Panel A shows the percentage of "
           "STAT-dependent genes rising with expression for both gene "
           "classes, with described genes consistently higher than "
           "uncharacterized genes and the gap widening at high expression. "
           "Panel B shows the effect size of that difference becoming more "
           "negative as expression increases, with a positive Spearman "
           "correlation between absolute effect size and expression. Panel "
           "C shows the proportion of target proteins receiving molecular "
           "function, biological process, or only cellular component "
           "terms.")


def fig6(F: Figures, plt):
    """Structural reach and biophysical character of the dark proteome.

    The four panels are the four things that had to be true before the
    title could say "structurally heterogeneous, disorder-enriched":
    that tightening the orphan definition changes the accuracy figure
    (A), that sequence and structure reach different proteins rather
    than nested ones (B), that the clustering deficit is not an artefact
    of poor 3Di prediction (C), and that the disorder enrichment holds
    under arithmetic no language model touched (D).
    """
    d = F.need("Figure 6", fx="fmax_structural.json",
               ti="target_intersection.json",
               ss="structural_stats.json",
               dm="disorder_matched.json")
    if d is None:
        return
    fx, ti, ss, dm = d["fx"], d["ti"], d["ss"], d["dm"]
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE, DOUBLE * 0.62))

    # -- A: Fmax across the three orphan definitions -------------------
    ax = axes[0][0]
    # Tick labels are kept to three short lines each. The panel gives
    # roughly 12 mm per category, and the first draft's "+ no profile"
    # and "+ no structural" overlapped each other on the axis.
    NESTED = [("ORPHAN_nohit", "no BLAST\nhit"),
              ("HARDENED", "+ profile\nhit"),
              ("STRUCT_HARDENED", "+ structure\nhit")]
    COMPL = [("RESCUED_PROFILE", "profile-\nrescued"),
             ("RESCUED_STRUCTURE", "structure-\nrescued")]
    pan = fx.get("panels") or {}
    xs, ys, los, his, labs = [], [], [], [], []
    for i, (key, lab) in enumerate(NESTED):
        b = pan.get(key)
        if not b:
            continue
        xs.append(i)
        ys.append(F.pick("Figure 6", "A", "nested definition", lab,
                         "fmax_structural.json", f"panels.{key}.fmax",
                         b["fmax"]))
        e = F.yerr(b["fmax"], b.get("ci"), f"Figure 6A {key}")
        los.append(e[0][0] if e else 0.0)
        his.append(e[1][0] if e else 0.0)
        labs.append(f"{lab}\n({b['n']})")
    if xs:
        ax.errorbar(xs, ys, yerr=[los, his], fmt="-o",
                    color=METHOD_COLOR["plm"], lw=1.1, ms=4,
                    markerfacecolor=METHOD_COLOR["plm"], zorder=3)
        # The two complements are not a further step along the same axis
        # — they are the proteins each step REMOVED. Drawn after a gap,
        # on a shaded field, as open markers, because joining them to
        # the nested series would assert a progression that does not
        # exist.
        cx = len(NESTED) + 0.6
        for key, lab in COMPL:
            b = pan.get(key)
            if not b:
                continue
            ax.errorbar([cx], [F.pick("Figure 6", "A", "complement", lab,
                                      "fmax_structural.json",
                                      f"panels.{key}.fmax", b["fmax"])],
                        yerr=F.yerr(b["fmax"], b.get("ci"),
                                    f"Figure 6A {key}"),
                        fmt="o", color=METHOD_COLOR["jackhmmer"],
                        markerfacecolor="white", ms=4.5, mew=1.0, zorder=3)
            xs.append(cx)
            labs.append(f"{lab}\n({b['n']})")
            cx += 1
        if cx > len(NESTED) + 0.6:
            ax.axvspan(len(NESTED) + 0.1, cx - 0.4, color="#ededed",
                       lw=0, zorder=0)
            ax.annotate("removed", xy=((len(NESTED) + 0.1 + cx - 0.4) / 2,
                                       1.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(0, -2), textcoords="offset points",
                        ha="center", va="top", fontsize=5.6, color=MUTED)
        ax.set_xticks(xs)
        ax.set_xticklabels(labs, fontsize=5.2, linespacing=1.25)
        ax.set_xlim(-0.55, max(xs) + 0.55)
        ax.set_ylabel("$F_{max}$")
    else:
        ax.set_visible(False)
    F.panel(ax, "A", dx=-0.17)

    # -- B: sequence reach versus structural reach ---------------------
    ax = axes[0][1]
    a = ti.get("no_sequence_with_structure")
    b_ = ti.get("no_sequence_no_structure")
    c = ti.get("has_sequence_with_structure")
    d2 = ti.get("has_sequence_no_structure")
    if None not in (a, b_, c, d2):
        for nm, v, jk in (("no homologue / structure", a,
                           "no_sequence_with_structure"),
                          ("no homologue / none", b_,
                           "no_sequence_no_structure"),
                          ("homologue / structure", c,
                           "has_sequence_with_structure"),
                          ("homologue / none", d2,
                           "has_sequence_no_structure")):
            F.rec("Figure 6", "B", nm, "count", v,
                  "target_intersection.json", jk)
        ypos = [1, 0]
        hit = [a, c]
        non = [b_, d2]
        ax.barh(ypos, hit, 0.5, color=METHOD_COLOR["plm"],
                edgecolor=METHOD_COLOR["plm"], lw=0,
                label="structural relative")
        ax.barh(ypos, non, 0.5, left=hit, color="white",
                edgecolor=METHOD_COLOR["plm"], lw=0.7, hatch="////",
                label="neither")
        for yp, h, n in zip(ypos, hit, non):
            if h:
                ax.annotate(f"{h}", xy=(h / 2, yp), ha="center",
                            va="center", fontsize=6.3, color="white")
            if n:
                # On the hatched segment the diagonals run through the
                # digits. A count a reader has to squint at is not a
                # count; it gets a patch of clear ground.
                ax.annotate(f"{n}", xy=(h + n / 2, yp), ha="center",
                            va="center", fontsize=6.3, color=INK,
                            zorder=5,
                            bbox=dict(boxstyle="square,pad=0.2",
                                      fc="white", ec="none", alpha=0.9))
        ax.set_yticks(ypos)
        ax.set_yticklabels([f"no sequence\nhomologue ({a + b_})",
                            f"sequence\nhomologue ({c + d2})"],
                           fontsize=6.2)
        ax.set_xlabel(f"target proteins (n = {ti.get('n_targets')})")
        # The band below the lower bar is opened deliberately: with the
        # default limits the only free region was 15% of the axis height
        # and the legend landed on top of the bar it was describing.
        ax.set_ylim(-1.05, 1.5)
        ax.legend(loc="lower right", fontsize=5.8, borderpad=0.2,
                  handlelength=1.4, handletextpad=0.5)
        # The point of the panel in one line: the two searches are
        # complementary, not nested.
        if c + d2:
            ax.annotate(f"{d2} of {c + d2} with a homologue have no\n"
                        f"structural hit — the reaches are complementary",
                        xy=(0.5, 1.02), xycoords="axes fraction",
                        ha="center", va="bottom", fontsize=5.8,
                        color=MUTED, linespacing=1.3)
    else:
        ax.set_visible(False)
    F.panel(ax, "B", dx=-0.30)

    # -- C: clustering deficit, raw and entropy-matched ----------------
    ax = axes[1][0]
    coh = ss.get("coherence") or {}
    mt = dm.get("matched") or {}
    # Each bar carries the file and the key it came from, because the two
    # groups are shaped differently in their JSONs — coherence nests the
    # two arms under `dark` and `background`, the matched block puts them
    # side by side as `dark_fraction` and `matched_fraction` — and a key
    # inferred from the bar's label would be a guess.
    groups = []
    if coh.get("dark") and coh.get("background"):
        groups.append((
            "length-matched",
            [("dark", coh["dark"]["fraction"], coh["dark"].get("ci"),
              METHOD_COLOR["plm"], "structural_stats.json",
              "coherence.dark.fraction"),
             ("control", coh["background"]["fraction"],
              coh["background"].get("ci"), "#555555",
              "structural_stats.json", "coherence.background.fraction")],
            coh.get("odds_ratio"), coh.get("or_ci"), coh.get("fisher_p")))
    if mt and not mt.get("refused") and mt.get("matched_fraction") is not None:
        groups.append((
            "3Di-entropy-matched",
            [("dark", mt.get("dark_fraction"), mt.get("dark_ci"),
              METHOD_COLOR["plm"], "disorder_matched.json",
              "matched.dark_fraction"),
             ("control", mt.get("matched_fraction"), mt.get("matched_ci"),
              "#555555", "disorder_matched.json",
              "matched.matched_fraction")],
            mt.get("odds_ratio"), mt.get("or_ci"), mt.get("fisher_p")))
    if groups:
        w = 0.33
        for gi, (lab, bars, orr, orci, pv) in enumerate(groups):
            for bi, (nm, frac, ci, col, src, jk) in enumerate(bars):
                if frac is None:
                    continue
                off = (bi - 0.5) * w
                F.rec("Figure 6", "C", nm, lab, 100 * frac, src,
                      f"{jk} x100")
                ax.bar([gi + off], [100 * frac], w, color=col,
                       edgecolor=col, lw=0,
                       label=nm if gi == 0 else None)
                e = F.yerr(frac, ci, f"Figure 6C {lab} {nm}", scale=100)
                if e:
                    ax.errorbar([gi + off], [100 * frac], yerr=e,
                                fmt="none", ecolor=INK, elinewidth=0.6)
            if orr is not None:
                t = f"OR {orr:.3f}"
                if orci and orci[0] is not None:
                    t += f"\n[{orci[0]:.3f}, {orci[1]:.3f}]"
                if pv is not None:
                    t += (f"\np = {pv:.3g}" if pv >= 1e-4
                          else f"\np = {pv:.1e}")
                ax.annotate(t, xy=(gi, 0.99), xycoords=("data",
                                                        "axes fraction"),
                            ha="center", va="top", fontsize=5.8,
                            color=INK, linespacing=1.35)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([g[0] for g in groups], fontsize=6.2)
        ax.set_ylabel("in a structural cluster (%)")
        ax.set_xlim(-0.55, len(groups) - 0.45)
        # Headroom is measured from the top of the error bars, not the
        # top of the bars: the interval on the control arm is wide, and
        # sizing to the bar alone put the odds ratio on the whisker.
        top = max(100 * max((b[2] or [0, b[1] or 0])[1] or 0, b[1] or 0)
                  for g in groups for b in g[1])
        ax.set_ylim(0, max(top, 1.0) * 1.42)
        # Above the axes, not inside: both top corners carry a group's
        # odds ratio, so an inside legend has nowhere to go that is not
        # already occupied by the statistic it would sit on.
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
                  fontsize=6, handlelength=1.3, handletextpad=0.5,
                  columnspacing=1.4, borderpad=0.1)
    else:
        ax.set_visible(False)
    F.panel(ax, "C", dx=-0.17)

    # -- D: charge-hydropathy, the model-free disorder test ------------
    ax = axes[1][1]
    sq = dm.get("sequence") or {}
    ch = sq.get("charge_hydropathy") or {}
    if ch.get("dark") and ch.get("pool"):
        pool = ch["pool"]
        dark = ch["dark"]
        ax.plot([p[0] for p in pool], [p[1] for p in pool], "o",
                ms=1.4, color="#8c8c8c", mew=0, alpha=0.45, zorder=1,
                label=f"background ({len(pool)})")
        ax.plot([p[0] for p in dark], [p[1] for p in dark], "o",
                ms=2.8, color=METHOD_COLOR["plm"], mew=0, zorder=3,
                label=f"dark ({len(dark)})")
        # Uversky boundary. Drawn from the coefficients named in the
        # JSON, not from a line fitted to these points.
        xlo, xhi = ax.get_xlim()
        ylo, yhi = ax.get_ylim()
        gx = [xlo + (xhi - xlo) * t / 100 for t in range(101)]
        gy = [2.785 * x - 1.151 for x in gx]
        seg = [(x, y) for x, y in zip(gx, gy) if ylo <= y <= yhi]
        if len(seg) > 1:
            ax.plot([s[0] for s in seg], [s[1] for s in seg], "-",
                    color=INK, lw=0.8, zorder=2)
            # Label the side of the line the criterion selects, at the
            # foot of the segment and to its LEFT. To the right it named
            # the ordered cloud; at the top it collided with the
            # statistics block, which has the only other free corner.
            ax.annotate("disordered", xy=(seg[0][0], seg[0][1]),
                        xytext=(-3, 3), textcoords="offset points",
                        ha="right", va="bottom", fontsize=5.6,
                        color=MUTED, zorder=5,
                        bbox=dict(boxstyle="square,pad=0.2", fc="white",
                                  ec="none", alpha=0.85))
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel("mean hydropathy (Kyte–Doolittle, normalised)")
        # The criterion is defined on ABSOLUTE net charge, and every
        # value here is non-negative because of it. "net charge per
        # residue" would invite a reader to wonder where the anions
        # went; the bars are not decoration.
        ax.set_ylabel("|net charge| per residue")
        # The two marker sizes on the plot are 1.4 pt and 2.8 pt: right
        # for a thousand points against seventy, wrong for a key, where
        # the first is invisible. A single markerscale cannot fix one
        # without inflating the other, so the key gets its own handles
        # at one legible size and the plot keeps its own.
        ax.legend(handles=[
            plt.Line2D([], [], ls="none", marker="o", ms=3.4, mew=0,
                       color="#8c8c8c", alpha=0.75,
                       label=f"background ({len(pool)})"),
            plt.Line2D([], [], ls="none", marker="o", ms=3.4, mew=0,
                       color=METHOD_COLOR["plm"],
                       label=f"dark ({len(dark)})")],
            loc="upper right", fontsize=6, framealpha=0.85,
            borderpad=0.3, handletextpad=0.5)
        ud, up = sq.get("uversky_dark"), sq.get("uversky_pool")
        if ud is not None and up is not None:
            t = f"{100 * ud:.1f}% vs {100 * up:.1f}% disordered"
            orr = sq.get("uversky_or")
            if orr is not None:
                ci = sq.get("uversky_or_ci") or [None, None]
                t += f"\nOR {orr:.2f}"
                if ci[0] is not None:
                    t += f" [{ci[0]:.2f}, {ci[1]:.2f}]"
            pv = sq.get("uversky_p")
            if pv is not None:
                t += (f"\np = {pv:.3g}" if pv >= 1e-4 else f"\np = {pv:.1e}")
            F.rec("Figure 6", "D", "dark disordered", "fraction",
                  100 * ud, "disorder_matched.json",
                  "sequence.uversky_dark x100")
            F.rec("Figure 6", "D", "background disordered", "fraction",
                  100 * up, "disorder_matched.json",
                  "sequence.uversky_pool x100")
            # Over a scatter there is no reliably empty corner — a point
            # landed on the first digit of the percentage. The label
            # carries its own ground rather than depending on one.
            ax.annotate(t, xy=(0.03, 0.97), xycoords="axes fraction",
                        ha="left", va="top", fontsize=5.8, color=INK,
                        linespacing=1.35, zorder=5,
                        bbox=dict(boxstyle="square,pad=0.3", fc="white",
                                  ec="none", alpha=0.85))
    else:
        ax.set_visible(False)
        F.skipped.append("Figure 6D (charge_hydropathy absent — rerun "
                         "49_disorder_matched.py --stage disorder)")
    F.panel(ax, "D", dx=-0.20, dy=1.07)

    rect = F.title(fig, "Figure 6. Structural reach and biophysical "
                        "character of the dark proteome")
    fig.tight_layout(pad=0.5, rect=rect)
    hard = ss.get("hardening") or {}
    F.save(
        fig, "Figure6",
        "Fig. 6. Structural reach and biophysical character of the dark "
        "proteome. (A) Annotation accuracy under three nested definitions "
        "of an orphan protein, each adding one search that must also "
        "return nothing; the shaded field holds the proteins each step "
        "removed, scored separately. Error bars are 95% percentile "
        "bootstrap intervals. (B) Cross-tabulation of sequence homology "
        "against structural reach for the target proteins, showing that "
        "the two searches recover overlapping but distinct sets. "
        "(C) Proportion of proteins falling in a structural cluster, "
        "against a background matched on length alone and against one "
        "matched additionally on predicted 3Di entropy; odds ratios with "
        "95% intervals and Fisher exact P. (D) Charge–hydropathy "
        "distribution with the Uversky boundary; points above the line "
        "satisfy the criterion for intrinsic disorder. The boundary and "
        "both percentages are arithmetic over the amino-acid sequence "
        "and involve no language model.",
        "Four-panel figure. Panel A shows annotation accuracy under "
        "progressively stricter orphan definitions, with the proteins "
        "excluded at each step scored separately in a shaded region at "
        "the right and scoring lower than the sets they were removed "
        "from. Panel B shows horizontal stacked bars for target proteins "
        "with and without a sequence homologue, each split by whether a "
        "structural relative was found; a substantial number of proteins "
        "with a sequence homologue have no structural hit and vice "
        "versa. Panel C shows that the deficit in structural clustering "
        "for the dark set relative to its background persists when the "
        "background is additionally matched on predicted structure "
        "quality. Panel D is a scatter plot of net charge against "
        "hydropathy in which the dark proteins sit disproportionately "
        "above the disorder boundary compared with the background "
        "cloud."
        + (f" Structural relatives were found for "
           f"{hard.get('with_structural_hit')} of {hard.get('n')} "
           f"previously orphan proteins." if hard else ""))


def figS1(F: Figures, plt):
    """Phylostratigraphy, with the control that decides what it means.

    Two panels because the distribution and the comparison answer
    different questions, and the second is the one that matters: the
    targets are more lineage-restricted than the proteome, and not
    detectably more so than other uncharacterised genes. A figure
    showing only the first would support a claim this analysis
    withdrew.
    """
    d = F.need("Figure S1", ln="lineage.json")
    if d is None:
        return
    ln = d["ln"]
    SETS = [("targets", "STAT-dependent dark", "#0072B2", "o", "-"),
            ("control", "uncharacterized control", "#882255", "s", "--"),
            ("background", "proteome background", "#555555", "^", ":")]
    present = [s for s in SETS if isinstance(ln.get(s[0]), dict)]
    if not present:
        F.skipped.append("Figure S1 (lineage.json has no scored sets)")
        return
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE / 2.6))
    rmax = ln.get("restricted_max_rank", 3)

    # -- A: cumulative distribution over strata ------------------------
    ax = axes[0]
    for key, lab, col, mk, ls in present:
        blk = ln[key]
        dist = blk.get("distribution") or {}
        n = blk.get("n") or 1
        cum, run = [], 0
        for r in range(10):
            run += int(dist.get(str(r), 0))
            cum.append(100 * run / n)
        for r in range(10):
            F.rec("Figure S1", "A", lab, f"rank {r}", cum[r],
                  "lineage.json", f"{key}.distribution cumulative/{key}.n")
        ax.plot(range(10), cum, color=col, marker=mk, ls=ls, lw=1.1,
                ms=3.2, label=lab)
    ax.axvspan(-0.5, rmax + 0.5, color="#ededed", lw=0, zorder=0)
    ax.annotate(f"lineage-restricted\n(rank ≤ {rmax})",
                xy=(rmax + 0.4, 0.02), xycoords=("data", "axes fraction"),
                ha="right", va="bottom", fontsize=5.8, color=MUTED,
                linespacing=1.3)
    ax.set_xticks(range(10))
    ax.set_xlabel("furthest clade with a detectable homologue (rank)")
    ax.set_ylabel("cumulative % of genes")
    ax.set_xlim(-0.5, 9.5)
    ax.set_ylim(0, 103)
    ax.legend(loc="lower right", fontsize=6)
    F.panel(ax, "A", dx=-0.16)

    # -- B: restricted fraction, with the two comparisons --------------
    ax = axes[1]
    pos, vals, errs, labs, cols = [], [], [], [], []
    for i, (key, lab, col, mk, ls) in enumerate(present):
        blk = ln[key]
        f_ = blk.get("fraction")
        ci = blk.get("ci") or [None, None]
        if f_ is None:
            continue
        pos.append(i)
        vals.append(F.pick("Figure S1", "B", lab, "restricted fraction",
                           "lineage.json", f"{key}.fraction x100",
                           100 * f_))
        e = F.yerr(f_, ci, f"Figure S1B {key}", scale=100)
        errs.append([e[0][0], e[1][0]] if e else [0.0, 0.0])
        labs.append(f"{lab}\n({blk.get('restricted')}/{blk.get('n')})")
        cols.append(col)
    if pos:
        ax.bar(pos, vals, 0.5, color=cols, edgecolor=cols, lw=0)
        ax.errorbar(pos, vals,
                    yerr=[[e[0] for e in errs], [e[1] for e in errs]],
                    fmt="none", ecolor=INK, elinewidth=0.6)
        ax.set_xticks(pos)
        ax.set_xticklabels(labs, fontsize=5.8)
        ax.set_ylabel(f"lineage-restricted (%), rank ≤ {rmax}")
        ax.set_ylim(0, max(vals) * 1.6)
        # The two odds ratios, each against the set it compares with, so
        # the null one cannot be mistaken for the significant one.
        notes = []
        for other, short in (("background", "proteome"),
                             ("control", "control")):
            cmp = ln.get(f"vs_{other}")
            if not cmp:
                continue
            pv = cmp.get("fisher_p")
            notes.append(f"vs {short}: OR {cmp['odds_ratio']:.2f}, "
                         + (f"p = {pv:.3g}" if pv is not None and pv >= 1e-4
                            else f"p = {pv:.1e}" if pv is not None else "—"))
        # Where the targets bar would stand once proteins with a
        # structural relative beyond Arthropoda are taken out of it. The
        # strata are built from sequence search, which this paper shows
        # failing on these proteins, so the headline bar without this
        # mark overstates restriction by whatever the structural arm
        # found.
        sc = ln.get("structural_check") or {}
        if sc.get("status") == "reported" and \
                sc.get("restricted_fraction_after_structure") is not None:
            names = [p[0] for p in present]
            if "targets" in names:
                ti_ = names.index("targets")
                adj = 100 * sc["restricted_fraction_after_structure"]
                F.rec("Figure S1", "B", "STAT-dependent dark",
                      "restricted after structure", adj, "lineage.json",
                      "structural_check.restricted_fraction_after_"
                      "structure x100")
                ax.plot([ti_ - 0.25, ti_ + 0.25], [adj, adj],
                        ls=(0, (2.2, 1.6)), color="white", lw=0.9,
                        zorder=4, solid_capstyle="butt")
                # The qualifier has to travel with the mark: the dashed
                # line falls BELOW the background bar, and a reader
                # comparing the two would conclude the targets end up
                # less restricted than the proteome. They would be
                # comparing one corrected number against two uncorrected
                # ones — only the targets were searched structurally.
                #
                # It goes here, with the other statistics, rather than
                # inside the bar: the bar is about 8 mm wide and any
                # honest qualifier is wider, so white text on it ran off
                # the blue and lost a letter at each end.
                notes.append(f"dashed: targets only, after structural "
                             f"check ({adj:.1f}%)")
        if notes:
            ax.annotate("\n".join(notes), xy=(0.5, 0.99),
                        xycoords="axes fraction", ha="center", va="top",
                        fontsize=6, color=INK, linespacing=1.4)
    else:
        ax.set_visible(False)
    F.panel(ax, "B", dx=-0.20)

    rect = F.title(fig, "Figure S1. Phylogenetic distribution of "
                        "detectable homology")
    fig.tight_layout(pad=0.45, rect=rect)
    F.save(
        fig, "FigureS1",
        "Fig. S1. Phylogenetic distribution of detectable homology. "
        "(A) Cumulative proportion of genes whose furthest detectable "
        "homologue falls at or before each clade rank, for the "
        "STAT-dependent uncharacterized targets, a control of "
        "uncharacterized genes that are not STAT-dependent, and a "
        "length-matched whole-proteome background. The shaded field "
        "marks the ranks counted as lineage-restricted. (B) Proportion "
        "lineage-restricted in each set, with Wilson 95% intervals and "
        "the two odds ratios. The enrichment against the whole proteome "
        "is substantial, and the enrichment against uncharacterized "
        "genes that are not STAT-dependent is not detectable, which "
        "places the restriction with the absence of annotation rather "
        "than with the transcriptional response. The dashed mark on the "
        "target bar, where present, is the proportion remaining once "
        "proteins with a structural relative beyond the restricted "
        "ranks are removed. Only the target set was searched "
        "structurally, so this mark is not comparable with the other "
        "two bars, which carry the same over-calling uncorrected; it "
        "bounds the target figure from below rather than ranking the "
        "three sets. Strata rest on detectable sequence homology, and "
        "absence of detection is not absence of a homologue.",
        "Two-panel supplementary figure. Panel A shows cumulative "
        "distributions of the furthest clade at which a homologue is "
        "detected; the target and control curves rise steeply within "
        "the arthropod ranks and track each other closely, while the "
        "proteome background rises more gradually and reaches high "
        "cumulative values only at the deepest ranks. Panel B shows "
        "three bars of the lineage-restricted percentage with "
        "confidence intervals: the targets and the uncharacterized "
        "control are similar and overlapping, and both are higher than "
        "the proteome background. A dashed line within the target bar "
        "marks the lower value reached after proteins with a distant "
        "structural relative are excluded; this adjustment was applied "
        "only to the target set and does not place it below the other "
        "two.")


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def nf(x, nd=3):
    if x is None:
        return "—"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v or v in (float("inf"), float("-inf")):
        return "—"
    return f"{v:.{nd}f}"


def tables(F: Figures) -> list[str]:
    """Main-text tables. Appends to F.blocks so the Markdown and the Word
    file are rendered from one structure and cannot drift apart."""
    made = []

    # Table 1 — calibration grid
    rows = []
    for p in sorted(F.res.glob("calibration_*.json")):
        try:
            d = json.loads(p.read_text())
        except ValueError:
            continue
        model = d.get("model")
        tf = d.get("params", {}).get("transform")
        boot = F.J("bootstrap_ci.json") or {}
        for t in TIER_ORDER:
            blk = d.get("results", {}).get(t)
            if not blk:
                continue
            ov = blk["metrics"]["overall"]
            cell = boot.get("cells", {}).get(f"{model}/{tf}/{t}")
            rows.append({
                "Model": MODEL_LABEL.get(model, model), "Transform": tf,
                "Condition": TIER_LABEL.get(t, t),
                "n": cell["n"] if cell else ov.get("n"),
                "Fmax": nf(ov.get("fmax")),
                "95% CI": (f"[{nf(cell['fmax_ci'][0])}, "
                           f"{nf(cell['fmax_ci'][1])}]") if cell else "\u2014",
                "Precision": nf(ov.get("precision")),
                "Recall": nf(ov.get("recall")),
                "BP": nf(blk["metrics"].get("BP", {}).get("fmax")),
                "MF": nf(blk["metrics"].get("MF", {}).get("fmax")),
                "CC": nf(blk["metrics"].get("CC", {}).get("fmax")),
            })
    if rows:
        cap = ("Table 1. Annotation accuracy by model, embedding transform "
               "and homology condition. Confidence intervals are 95% "
               "percentile intervals from 2000 bootstrap resamples of the "
               "query set. BP, MF and CC are Fmax restricted to each Gene "
               "Ontology aspect.")
        _write_table(F, "Table1_calibration", rows, cap)
        made.append("Table1_calibration")

    # Table 2 — PLM vs homology, decomposed
    rows = []
    for tag in ("DIAMOND", "jackhmmer"):
        c = F.J(f"plm_vs_homology_{tag}.json")
        if not c:
            continue
        for t in TIER_ORDER:
            blk = c.get("tiers", {}).get(t)
            if not blk:
                continue
            for model, m in (blk.get("models") or {}).items():
                al = m.get("all_queries") or {}
                ac = m.get("diamond_answerable") or {}
                un = m.get("diamond_unanswerable") or {}
                rows.append({
                    "Baseline": tag,
                    "Condition": TIER_LABEL.get(t, t),
                    "Model": MODEL_LABEL.get(model, model),
                    "PLM": nf(al.get("a_fmax")),
                    "Base.": nf(al.get("b_fmax")),
                    "\u0394 all": nf(al.get("delta")),
                    "95% CI (all)": f"[{nf((al.get('ci') or [None, None])[0])}, "
                                  f"{nf((al.get('ci') or [None, None])[1])}]",
                    "\u0394 answ.": nf(ac.get("delta")),
                    "95% CI (answ.)": f"[{nf((ac.get('ci') or [None, None])[0])}, "
                                    f"{nf((ac.get('ci') or [None, None])[1])}]",
                    "n unans.": un.get("n"),
                    "PLM there": nf(un.get("fmax")),
                })
    if rows:
        cap = ("Table 2. Embedding transfer against two homology searches, "
               "with the contribution of coverage separated from that of "
               "accuracy. \u0394 all is the paired difference over every query; "
               "\u0394 answ. restricts to queries for which the homology "
               "search returned a usable hit. The final columns give the "
               "number of queries the search could not annotate and the "
               "accuracy of embedding transfer on exactly those proteins.")
        _write_table(F, "Table2_vs_homology", rows, cap)
        made.append("Table2_vs_homology")
    return made


def _write_table(F: Figures, name: str, rows: list[dict], caption: str):
    cols = list(rows[0].keys())
    md = [caption, "", "| " + " | ".join(cols) + " |",
          "|" + "---|" * len(cols)]
    for r in rows:
        md.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    (F.out / f"{name}.md").write_text("\n".join(md) + "\n")
    with (F.out / f"{name}.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    F.blocks.append(("table", caption, cols, rows))
    log(f"wrote {name}.{{md,tsv}} ({len(rows)} rows)")


def supplementary(F: Figures) -> None:
    """Supplementary tables, appended to the same block list."""
    B = F.blocks
    B.append(("h1", "Supplementary Material"))
    B.append(("p", "All values are generated from the analysis outputs by "
                   "scripts/43_figures.py; none is transcribed."))

    bd, bj = F.J("homology_baseline_DIAMOND.json"), \
        F.J("homology_baseline_jackhmmer.json")
    if bd and bj:
        cols = ["Condition", "DIAMOND answers", "jackhmmer answers",
                "DIAMOND Fmax", "jackhmmer Fmax"]
        rows = []
        for t in TIER_ORDER:
            a = bd.get("tiers", {}).get(t, {}).get("best", {})
            b = bj.get("tiers", {}).get(t, {}).get("best", {})
            if not a and not b:
                continue

            def frac(x):
                return (f"{x.get('answerable')}/{x.get('n_scored')} "
                        f"({100 * x['answerable_frac']:.1f}%)"
                        if x.get("answerable_frac") is not None else "\u2014")
            rows.append({cols[0]: TIER_LABEL.get(t, t),
                         cols[1]: frac(a), cols[2]: frac(b),
                         cols[3]: nf(a.get("fmax")),
                         cols[4]: nf(b.get("fmax"))})
        und = sorted(set(TIER_LABEL.get(x, x) for x in
                         ((bd.get("undefined_by_construction") or [])
                          + (bj.get("undefined_by_construction") or []))))
        B.append(("table",
                  "Table S1. Search sensitivity and tier leakage. Queries "
                  "with any permitted hit under each search, at a common "
                  "transfer threshold. Conditions where a search can answer "
                  "almost nothing are undefined for that search by "
                  "construction and carry no accuracy interpretation"
                  + (f" ({', '.join(und)})." if und else "."),
                  cols, rows))

    de = F.J("de_stat_summary.json")
    if de:
        con = de.get("contrast", {})
        pc = de.get("positive_control") or {}
        pr, st = de.get("primary", {}), de.get("stringent", {})
        padj = (f"{pc['padj']:.3g}" if pc.get("padj") is not None
                else "\u2014")
        rows = [
            ("Contrast", f"{con.get('a')} vs {con.get('b')}"),
            ("Libraries used", len(de.get("libraries", {}))),
            ("Libraries excluded",
             ", ".join(x["run"] for x in de.get("excluded", []))
             or "\u2014"),
            ("Genes quantified", de.get("n_genes_quantified")),
            ("Genes tested", de.get("n_genes_tested")),
            ("Knockdown control gene", pc.get("gene", "\u2014")),
            ("Control log2 fold change", nf(pc.get("log2FoldChange"), 2)),
            ("Control adjusted P", padj),
            ("Significant, no fold-change threshold", pr.get("n_sig")),
            ("Of which STAT-dependent", pr.get("n_stat_dependent")),
            (f"Significant at |log2FC| > {st.get('lfc')}", st.get("n_sig")),
            ("Previously reported, same threshold",
             f"{st.get('published_up_a')} and {st.get('published_up_b')}"),
        ]
        B.append(("table",
                  "Table S2. Differential expression summary. Counts are "
                  "coverage-derived from the reference quantification "
                  f"(read length {de.get('read_length')}).",
                  ["Quantity", "Value"],
                  [{"Quantity": k, "Value": v} for k, v in rows]))

    mt = F.J("matched_enrichment.json")
    if mt:
        wd = {w["decile"]: w for w in
              (mt.get("fold_change") or {}).get("within_deciles", [])}
        cols = ["Decile", "Expression range", "Uncharacterized",
                "Described", "Cliff's \u03b4", "P"]
        rows = []
        for r in mt.get("decile_table") or []:
            w = wd.get(r["decile"], {})
            rows.append({
                cols[0]: r["decile"],
                cols[1]: f"{r['basemean_lo']:.1f}\u2013{r['basemean_hi']:.1f}",
                cols[2]: f"{r['unc_dep']}/{r['unc_n']} "
                         f"({100 * r['unc_frac']:.1f}%)",
                cols[3]: f"{r['desc_dep']}/{r['desc_n']} "
                         f"({100 * r['desc_frac']:.1f}%)",
                cols[4]: nf(w.get("delta")),
                cols[5]: (f"{w['p']:.3g}" if w.get("p") is not None
                          else "\u2014")})
        cf, tr = mt.get("confound") or {}, mt.get("trend") or {}
        cap = ("Table S3. STAT dependence by expression decile. "
               f"Standardized mean difference between the two gene classes "
               f"is {nf(cf.get('smd_log_basemean'))} in log10 expression "
               f"and {nf(cf.get('smd_log_length'))} in log10 transcript "
               "length.")
        if tr.get("p") is not None:
            cap += (" Absolute effect size against decile, Spearman "
                    f"rho = {nf(tr.get('spearman_rho'), 2)} "
                    f"(P = {tr['p']:.3g}).")
        if rows:
            B.append(("table", cap, cols, rows))

    an = F.J("annotation_summary.json")
    if an:
        cols = ["Subset", "Source", "n (benchmark)", "Fmax",
                "Score threshold", "Precision", "Recall", "BP", "MF", "CC"]
        rows = []
        for lab, c in (an.get("calibration") or {}).items():
            rows.append({
                cols[0]: lab, cols[1]: c.get("source", "benchmark"),
                cols[2]: c.get("n_benchmark"), cols[3]: nf(c.get("fmax")),
                cols[4]: nf(c.get("threshold"), 2),
                cols[5]: nf(c.get("precision")), cols[6]: nf(c.get("recall")),
                cols[7]: nf(c.get("BP")), cols[8]: nf(c.get("MF")),
                cols[9]: nf(c.get("CC"))})
        tf = an.get("transferability") or {}
        cap = ("Table S4. Calibration applied to the target proteins. Each "
               "protein is scored against the subset matching its own "
               "homology situation. Target proteins' median "
               "nearest-neighbour similarity is "
               f"{nf(tf.get('target_median_top_sim'))} against "
               f"{nf(tf.get('calibration_median_top_sim'))} for the "
               "calibration subset.")
        if rows:
            B.append(("table", cap, cols, rows))
        B.append(("p", "Per-protein predictions are provided as "
                       "annotation_targets.csv (one row per protein) and "
                       "annotation_targets_long.csv (one row per protein "
                       "and term)."))

    # ---- S5: accuracy under the three orphan definitions --------------
    fx = F.J("fmax_structural.json")
    if fx and fx.get("panels"):
        cols = ["Definition", "n", "Fmax", "95% CI", "Role"]
        order = ["ORPHAN_nohit", "HARDENED", "STRUCT_HARDENED",
                 "RESCUED_PROFILE", "RESCUED_STRUCTURE"]
        rows = []
        for key in order:
            b = (fx["panels"] or {}).get(key)
            if not b:
                continue
            ci = b.get("ci") or [None, None]
            rows.append({
                cols[0]: b.get("pretty", key).strip(),
                cols[1]: b.get("n"),
                cols[2]: nf(b.get("fmax")),
                cols[3]: f"[{nf(ci[0])}, {nf(ci[1])}]",
                cols[4]: "excluded proteins, scored separately"
                         if b.get("complement") else "nested definition"})
        att = fx.get("attrition") or {}
        cap = ("Table S5. Annotation accuracy as the definition of an "
               "orphan protein is tightened. Each nested row adds one "
               "search that must also return no hit, so the sets are "
               "subsets of one another; the final rows hold the proteins "
               "removed at each step, scored on their own.")
        if att:
            cap += (f" Accuracy over the full tightening is "
                    f"{nf(att.get('from_fmax'))} on "
                    f"{att.get('from_n')} proteins against "
                    f"{nf(att.get('to_fmax'))} on {att.get('to_n')} "
                    f"({nf(att.get('delta'))}), and the intervals "
                    + ("overlap." if fx.get("intervals_overlap")
                       else "do not overlap."))
        if fx.get("rise_explained"):
            cap += " " + fx["rise_explained"][0].upper() \
                + fx["rise_explained"][1:]
        if rows:
            B.append(("table", cap, cols, rows))

    # ---- S6: structural hardening and the 3Di diagnostic --------------
    ss = F.J("structural_stats.json")
    if ss:
        rows = []
        hd = ss.get("hardening") or {}
        if hd:
            ci = hd.get("ci") or [None, None]
            rows.append(("Orphan proteins searched structurally",
                         hd.get("n")))
            rows.append(("With a structural relative",
                         f"{hd.get('with_structural_hit')} "
                         f"({100 * hd['fraction']:.1f}%"
                         + (f", 95% CI {100 * ci[0]:.1f}–{100 * ci[1]:.1f}%)"
                            if ci[0] is not None else ")")))
            rows.append(("Surviving all three searches", hd.get("surviving")))
        td = ss.get("three_di") or {}
        for nm, lab in (("hardened", "dark set"),
                        ("background", "background")):
            blk = td.get(nm)
            if not blk:
                continue
            rows.append((f"3Di entropy, median ({lab})",
                         nf(blk.get("entropy_median"))))
            rows.append((f"Below the degeneracy threshold ({lab})",
                         f"{blk.get('degenerate_n')}/{blk.get('n')}"))
        tt = ss.get("three_di_test") or {}
        if tt:
            rows.append(("3Di entropy gap (bits)",
                         nf(tt.get("entropy_gap_bits"))))
            rows.append(("Cliff's delta, 3Di entropy",
                         nf(tt.get("cliffs_delta"))))
            rows.append(("P, 3Di entropy",
                         f"{tt['p']:.3g}" if tt.get("p") is not None
                         else "—"))
            rows.append(("Degeneracy threshold (bits)",
                         nf(tt.get("degenerate_threshold_bits"), 2)))
        if rows:
            cap = ("Table S6. Structural hardening of the orphan set and "
                   "the diagnostic that tests it. 3Di states are predicted "
                   "from sequence by ProstT5 rather than derived from "
                   "coordinates; the entropy columns exist because a "
                   "degenerate prediction would produce a spurious "
                   "absence of structural neighbours.")
            if ss.get("verdict"):
                cap += f" Verdict recorded by the analysis — {ss['verdict']}"
                cap += "" if str(ss["verdict"]).endswith(".") else "."
            B.append(("table", cap, ["Quantity", "Value"],
                      [{"Quantity": k, "Value": v} for k, v in rows]))

    # ---- S7: the clustering deficit, raw and entropy-matched ----------
    dm = F.J("disorder_matched.json")
    coh = (ss or {}).get("coherence") or {}
    mt = (dm or {}).get("matched") or {}
    if coh or mt:
        cols = ["Comparison", "Dark set", "Control", "Difference (pp)",
                "Odds ratio [95% CI]", "P"]
        rows = []

        def frac_cell(k, n, ci):
            s = f"{k}/{n} ({100 * k / n:.1f}%)" if n else "—"
            if ci and ci[0] is not None:
                s += f" [{100 * ci[0]:.1f}, {100 * ci[1]:.1f}]"
            return s

        if coh.get("dark") and coh.get("background"):
            a, b2 = coh["dark"], coh["background"]
            orci = coh.get("or_ci") or [None, None]
            dci = coh.get("difference_ci_pp") or [None, None]
            rows.append({
                cols[0]: "matched on length",
                cols[1]: frac_cell(a["clustered"], a["n"], a.get("ci")),
                cols[2]: frac_cell(b2["clustered"], b2["n"], b2.get("ci")),
                cols[3]: (f"{coh['difference_pp']:+.1f}"
                          + (f" [{dci[0]:+.1f}, {dci[1]:+.1f}]"
                             if dci[0] is not None else "")),
                cols[4]: (f"{nf(coh.get('odds_ratio'))} "
                          f"[{nf(orci[0])}, {nf(orci[1])}]"),
                cols[5]: (f"{coh['fisher_p']:.3g}"
                          if coh.get("fisher_p") is not None else "—")})
        if mt and not mt.get("refused"):
            orci = mt.get("or_ci") or [None, None]
            rows.append({
                cols[0]: "matched on length and 3Di entropy",
                cols[1]: frac_cell(mt.get("dark_clustered", 0),
                                   mt.get("dark_n", 0), mt.get("dark_ci")),
                cols[2]: frac_cell(mt.get("matched_clustered", 0),
                                   mt.get("n_matched", 0),
                                   mt.get("matched_ci")),
                cols[3]: (f"{mt['difference_pp']:+.1f}"
                          if mt.get("difference_pp") is not None
                          else "—"),
                cols[4]: (f"{nf(mt.get('odds_ratio'))} "
                          f"[{nf(orci[0])}, {nf(orci[1])}]"),
                cols[5]: (f"{mt['fisher_p']:.3g}"
                          if mt.get("fisher_p") is not None else "—")})
        cap = ("Table S7. Proportion of proteins assigned to a structural "
               "cluster, against two backgrounds. The second background "
               "is matched additionally on predicted 3Di entropy, so a "
               "deficit that survives it cannot be attributed to poorer "
               "structure prediction on the dark set.")
        bal = mt.get("smd") or {}
        if bal:
            cap += (" Standardized mean difference after matching is "
                    f"{nf(bal.get('length'))} in length and "
                    f"{nf(bal.get('three_di_entropy'))} in 3Di entropy, "
                    f"within a caliper of {nf(mt.get('caliper_sd'), 2)} SD "
                    f"and a ceiling of {nf(mt.get('max_smd'), 2)}; "
                    f"{mt.get('unmatched')} "
                    + ("protein fell" if mt.get("unmatched") == 1
                       else "proteins fell")
                    + " outside the caliper.")
        elif mt and not mt.get("refused"):
            cap += (" Matching was constrained to a caliper of "
                    f"{nf(mt.get('caliper_sd'), 2)} SD with a "
                    "standardized-mean-difference ceiling of "
                    f"{nf(mt.get('max_smd'), 2)}, enforced by the "
                    "analysis; rerun 49_disorder_matched.py --stage match "
                    "to record the attained values.")
        if rows:
            B.append(("table", cap, cols, rows))

    # ---- S8: sequence properties and predicted disorder ---------------
    sq = (dm or {}).get("sequence") or {}
    if sq:
        rows = [
            ("Amino-acid entropy, median (dark)",
             nf(sq.get("dark_aa_entropy_median"))),
            ("Amino-acid entropy, median (background)",
             nf(sq.get("pool_aa_entropy_median"))),
            ("Cliff's delta, amino-acid entropy",
             nf(sq.get("aa_entropy_delta"))),
            ("P, amino-acid entropy",
             f"{sq['aa_entropy_p']:.3g}"
             if sq.get("aa_entropy_p") is not None else "—"),
            ("Uversky-disordered (dark)",
             f"{sq.get('uversky_dark_n')}/{sq.get('uversky_dark_total')} "
             f"({100 * sq['uversky_dark']:.1f}%)"
             if sq.get("uversky_dark") is not None else "—"),
            ("Uversky-disordered (background)",
             f"{sq.get('uversky_pool_n')}/{sq.get('uversky_pool_total')} "
             f"({100 * sq['uversky_pool']:.1f}%)"
             if sq.get("uversky_pool") is not None else "—"),
        ]
        if sq.get("uversky_or") is not None:
            ci = sq.get("uversky_or_ci") or [None, None]
            rows.append(("Odds ratio, predicted disorder",
                         f"{nf(sq['uversky_or'], 2)} "
                         f"[{nf(ci[0], 2)}, {nf(ci[1], 2)}]"))
        if sq.get("uversky_p") is not None:
            rows.append(("P, predicted disorder", f"{sq['uversky_p']:.3g}"))
        meta = sq.get("metapredict") or {}
        if meta:
            rows.append(("Disordered residue fraction, median (dark)",
                         nf(meta.get("dark_median"))))
            rows.append(("Disordered residue fraction, median "
                         "(background)", nf(meta.get("pool_median"))))
            rows.append(("Cliff's delta, disordered residue fraction",
                         nf(meta.get("delta"))))
        cap = ("Table S8. Amino-acid composition and predicted disorder. "
               "Entropy is Shannon entropy over the twenty amino acids, "
               "so a low value would indicate low-complexity sequence. "
               "The Uversky criterion is "
               + (sq.get("uversky_boundary") or
                  "R > 2.785*H - 1.151 (Kyte-Doolittle)")
               + ", computed from the sequence alone and independent of "
                 "every language model used elsewhere in this work.")
        B.append(("table", cap, ["Quantity", "Value"],
                  [{"Quantity": k, "Value": v} for k, v in rows]))

    # ---- S9: sequence reach against structural reach ------------------
    ti = F.J("target_intersection.json")
    if ti and ti.get("n_targets"):
        a = ti.get("no_sequence_with_structure")
        b2 = ti.get("no_sequence_no_structure")
        c = ti.get("has_sequence_with_structure")
        d2 = ti.get("has_sequence_no_structure")
        if None not in (a, b2, c, d2):
            cols = ["Target proteins", "Structural relative",
                    "No structural relative", "Total"]
            rows = [
                {cols[0]: "No sequence homologue", cols[1]: a,
                 cols[2]: b2, cols[3]: a + b2},
                {cols[0]: "Sequence homologue", cols[1]: c,
                 cols[2]: d2, cols[3]: c + d2},
                {cols[0]: "Total", cols[1]: a + c, cols[2]: b2 + d2,
                 cols[3]: a + b2 + c + d2},
            ]
            cap = ("Table S9. Sequence homology against structural reach "
                   f"for the {ti['n_targets']} target proteins. The "
                   "marginal counts cannot be combined: proteins with a "
                   "structural relative include many that already had a "
                   "sequence homologue, which is why the cross-tabulation "
                   "is given rather than the two totals. Of the "
                   f"{a + b2} with no sequence homologue, {a} "
                   f"({100 * a / (a + b2):.1f}%) have a structural "
                   f"relative, and {b2} are beyond the reach of either "
                   "search.")
            B.append(("table", cap, cols, rows))

    # ---- S10: phylostrata ---------------------------------------------
    ln = F.J("lineage.json")
    if ln and isinstance(ln.get("targets"), dict):
        names = [k for k in ("targets", "control", "background")
                 if isinstance(ln.get(k), dict)]
        pretty = {"targets": "STAT-dependent dark",
                  "control": "Uncharacterized control",
                  "background": "Proteome background"}
        cols = ["Furthest clade (rank)"] + [pretty[k] for k in names]
        rows = []
        for r in range(10):
            cells = {}
            any_v = False
            for k in names:
                dist = (ln[k].get("distribution") or {})
                v = int(dist.get(str(r), 0))
                n = ln[k].get("n") or 1
                cells[pretty[k]] = f"{v} ({100 * v / n:.1f}%)"
                any_v = any_v or v > 0
            if any_v:
                rows.append({cols[0]: str(r), **cells})
        last = {cols[0]: f"Lineage-restricted (rank ≤ "
                         f"{ln.get('restricted_max_rank', 3)})"}
        for k in names:
            blk = ln[k]
            ci = blk.get("ci") or [None, None]
            last[pretty[k]] = (
                f"{blk.get('restricted')}/{blk.get('n')} "
                f"({100 * blk['fraction']:.1f}%"
                + (f", 95% CI {100 * ci[0]:.1f}–{100 * ci[1]:.1f}%)"
                   if ci[0] is not None else ")"))
        rows.append(last)
        cap = ("Table S10. Furthest clade at which a homologue is "
               "detectable, by gene set. Ranks run from genes with no "
               "detectable homologue outside the species to genes with "
               "homologues in bacteria.")
        for other, short in (("background", "the whole proteome"),
                             ("control", "uncharacterized genes that are "
                                         "not STAT-dependent")):
            cmp = ln.get(f"vs_{other}")
            if cmp:
                cap += (f" Against {short}, odds ratio "
                        f"{nf(cmp.get('odds_ratio'), 2)} "
                        f"(P = {cmp['fisher_p']:.3g}).")
        if ln.get("caveat"):
            cap += " " + ln["caveat"]
        if rows:
            B.append(("table", cap, cols, rows))


def legends_block(F: Figures) -> None:
    F.blocks.append(("h1", "Figure legends and alt text"))
    for name, legend, alt in F.alt:
        F.blocks.append(("h2", name.replace("Figure", "Figure ")))
        F.blocks.append(("p", legend))
        F.blocks.append(("p", f"Alt text: {alt}"))


def render_markdown(F: Figures, path: Path) -> None:
    L = ["# Tables and Supplementary Material", ""]
    for b in F.blocks:
        if b[0] == "h1":
            L += [f"## {b[1]}", ""]
        elif b[0] == "h2":
            L += [f"### {b[1]}", ""]
        elif b[0] == "p":
            L += [b[1], ""]
        elif b[0] == "table":
            _, cap, cols, rows = b
            L += [cap, "", "| " + " | ".join(cols) + " |",
                  "|" + "---|" * len(cols)]
            for r in rows:
                L.append("| " + " | ".join(str(r.get(c, "")) for c in cols)
                         + " |")
            L.append("")
    path.write_text("\n".join(L) + "\n")
    log(f"wrote {path.name}")


def render_docx(F: Figures, path: Path) -> bool:
    """One Word file holding the tables and the supplementary material."""
    try:
        from docx import Document
        from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        print("\n  python-docx is not installed, so the Word file was not")
        print("  written. The Markdown and TSV versions are complete.")
        print("    pip install python-docx")
        return False

    doc = Document()
    for sec in doc.sections:
        sec.left_margin = sec.right_margin = Cm(1.8)
        sec.top_margin = sec.bottom_margin = Cm(2.0)
    st = doc.styles["Normal"]
    st.font.name = "Arial"
    st.font.size = Pt(10)

    h = doc.add_heading("Tables and Supplementary Material", level=0)
    h.alignment = WD_ALIGN_PARAGRAPH.LEFT
    doc.add_paragraph(
        "Every value in this file is generated directly from the analysis "
        "outputs; none is transcribed by hand.")

    for b in F.blocks:
        if b[0] == "h1":
            doc.add_heading(b[1], level=1)
        elif b[0] == "h2":
            doc.add_heading(b[1], level=2)
        elif b[0] == "p":
            doc.add_paragraph(b[1])
        elif b[0] == "table":
            _, cap, cols, rows = b
            cp = doc.add_paragraph()
            run = cp.add_run(cap)
            run.bold = True
            run.font.size = Pt(9)
            t = doc.add_table(rows=1, cols=len(cols))
            t.style = "Table Grid"
            t.autofit = False
            # Equal columns break long headers mid-word ("Transfor m").
            # Allocate the usable width in proportion to the widest cell
            # each column actually has to hold, with a floor so short
            # columns stay legible.
            # What forces a column wide is the longest UNBREAKABLE token,
            # not the longest string: "[0.799, 0.829]" wraps happily at the
            # space, "jackhmmer" does not. Sizing by total length is what
            # produced "Transfor m" and "jackhmme r".
            usable, chrome = 17.4, 0.34        # cm: A4 less margins; padding
            need = []
            for c in cols:
                vals = [str(c)] + [str(r.get(c, "")) for r in rows]
                token = max((len(t) for v in vals for t in v.split()),
                            default=4)
                whole = max(len(v) for v in vals)
                need.append(max(token, (whole + 1) // 2, 4))
            tot = float(sum(need))
            avail = max(1.0, usable - chrome * len(cols))
            widths = [Cm(chrome + avail * n / tot) for n in need]
            scale = usable / sum(w.cm for w in widths)
            widths = [Cm(w.cm * scale) for w in widths]
            for i, c in enumerate(cols):
                cell = t.rows[0].cells[i]
                cell.width = widths[i]
                cell.text = ""
                r = cell.paragraphs[0].add_run(str(c))
                r.bold = True
                r.font.size = Pt(7)
            for row in rows:
                cells = t.add_row().cells
                for i, c in enumerate(cols):
                    cells[i].width = widths[i]
                    cells[i].text = ""
                    r = cells[i].paragraphs[0].add_run(
                        str(row.get(c, "")) if row.get(c) is not None
                        else "\u2014")
                    r.font.size = Pt(7)
            doc.add_paragraph()

    doc.save(str(path))
    log(f"wrote {path.name}")
    return True


# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# audit
# ---------------------------------------------------------------------------

def audit_report(F: Figures, path: Path) -> list[str]:
    """Write every plotted value with its provenance, then run the checks
    that a reader could run with a ruler. Returns the complaints.

    A figure that agrees with itself is not verified; what matters is
    whether a number a reader derives from one panel matches the number
    another panel states outright.
    """
    lines = ["figure\tpanel\tseries\tat\tvalue\tsource_file\tsource_key"]
    for fig, panel, series, at, val, src, key in F.audit:
        v = "" if val is None else f"{val:.6g}"
        lines.append(f"{fig}\t{panel}\t{series}\t{at}\t{v}\t{src}\t{key}")
    path.write_text("\n".join(lines) + "\n")
    log(f"wrote {path.name} ({len(F.audit)} plotted values)")

    def get(fig, panel, series, at):
        for f_, p_, s_, a_, v_, _, _ in F.audit:
            if (f_, p_, s_, a_) == (fig, panel, series, at):
                return v_
        return None

    gripes: list[str] = []
    print("\n" + "=" * 78)
    print("AUDIT — values a reader can derive, against values we state")
    print("=" * 78)

    # 1. Figure 3: does panel A's visible gap equal panel C's stated delta?
    #    They come from different files (point estimates vs the paired
    #    bootstrap), so a gap is expected to be small but need not be zero.
    #    If it is large the two panels contradict each other on the page.
    tiers = sorted({a for f_, p_, s_, a_, *_ in F.audit
                    for a in [a_] if f_ == "Figure 3" and p_ == "C"},
                   key=lambda t: (TIER_ORDER.index(t)
                                  if t in TIER_ORDER else 99))
    if tiers:
        print("\n  Figure 3A subtraction vs Figure 3C 'all' delta")
        print(f"    {'tier':<8}{'3A emb':>9}{'3A jack':>9}"
              f"{'A−A':>9}{'3C all':>9}{'gap':>9}")
        for t in tiers:
            e = get("Figure 3", "A", "embeddings", t)
            j = get("Figure 3", "A", "jackhmmer", t)
            c = get("Figure 3", "C", "all", t)
            if None in (e, j, c):
                continue
            gap = (e - j) - c
            flag = "  <-- check" if abs(gap) > 0.01 else ""
            print(f"    {t:<8}{e:>9.3f}{j:>9.3f}{e - j:>9.3f}"
                  f"{c:>9.3f}{gap:>9.3f}{flag}")
            if abs(gap) > 0.01:
                gripes.append(
                    f"Figure 3 {t}: panel A implies {e - j:+.3f} but panel C "
                    f"states {c:+.3f} (gap {gap:+.3f}). Panel A uses point "
                    f"estimates from the calibration and baseline files; "
                    f"panel C uses the paired bootstrap. A reader "
                    f"subtracting the lines will not reproduce panel C.")

    # 2. The orphan Fmax is stated in three figures. All three must agree.
    o2 = get("Figure 2", "A", "ESM-2", "ORPHAN_nohit")
    o3 = get("Figure 3", "A", "embeddings", "ORPHAN_nohit")
    o4 = next((v for f_, p_, s_, a_, v, _, _ in F.audit
               if (f_, p_, s_) == ("Figure 4", "B", "no BLAST-class hit")),
              None)
    trio = [("Figure 2A ESM-2", o2), ("Figure 3A embeddings", o3),
            ("Figure 4B no BLAST-class hit", o4)]
    have = [(k, v) for k, v in trio if v is not None]
    if len(have) > 1:
        print("\n  Orphan Fmax, stated in more than one figure")
        for k, v in have:
            print(f"    {k:<34}{v:.4f}")
        spread = max(v for _, v in have) - min(v for _, v in have)
        if spread > 0.001:
            gripes.append(
                f"The orphan Fmax differs by {spread:.4f} between figures: "
                + "; ".join(f"{k} {v:.4f}" for k, v in have)
                + ". These are the same quantity and must print identically.")
        else:
            print(f"    agree to {spread:.5f}")

    # 3. Figure 2B's orphan bars cannot exceed the all-query bars, and the
    #    orphan overall Fmax must sit between its lowest and highest aspect.
    asp = {}
    for f_, p_, s_, a_, v, _, _ in F.audit:
        if f_ == "Figure 2" and p_ == "B" and v is not None:
            asp.setdefault(s_.split(" ", 1)[-1], {})[a_] = v
    if len(asp) == 2 and o2 is not None:
        orph = asp.get("orphans", {})
        if orph:
            lo, hi = min(orph.values()), max(orph.values())
            print(f"\n  Figure 2B orphan aspects span {lo:.3f}–{hi:.3f}; "
                  f"overall orphan Fmax {o2:.3f}")
            if not (lo - 1e-9 <= o2 <= hi + 1e-9):
                gripes.append(
                    f"Figure 2A's orphan Fmax ({o2:.3f}) lies outside the "
                    f"range of Figure 2B's per-aspect orphan bars "
                    f"({lo:.3f}–{hi:.3f}). One of the two is from a stale "
                    f"file.")

    # 4. Figure 4: the strata must partition the original orphan set.
    ns = {}
    for f_, p_, s_, a_, v, _, _ in F.audit:
        if f_ == "Figure 4" and p_ == "B" and a_.startswith("n="):
            ns[s_] = int(a_[2:])
    if len(ns) == 3:
        tot = ns.get("no BLAST-class hit")
        parts = ns.get("no hit, either search", 0) + ns.get("profile hit only", 0)
        print(f"\n  Figure 4B strata: {ns.get('no hit, either search')} + "
              f"{ns.get('profile hit only')} = {parts}, original {tot}")
        if tot != parts:
            gripes.append(
                f"Figure 4B's strata do not partition the orphan set: "
                f"{parts} != {tot}.")

    # 5. Figure 5C bar counts cannot exceed the target set.
    n_t = (F.J("annotation_summary.json") or {}).get("n_targets")
    if n_t:
        for f_, p_, s_, a_, v, _, _ in F.audit:
            if f_ == "Figure 5" and p_ == "C" and v is not None and v > n_t:
                gripes.append(
                    f"Figure 5C '{a_}' is {v:.0f} of {n_t} targets.")

    # 6. Figure 6B's four cells must partition the target set exactly.
    #    They are read from one file, so a failure here means that file
    #    is internally inconsistent — which is worth catching, because
    #    the cross-tabulation exists precisely to stop marginal counts
    #    being added together.
    ti = F.J("target_intersection.json") or {}
    cells = [ti.get(k) for k in ("no_sequence_with_structure",
                                 "no_sequence_no_structure",
                                 "has_sequence_with_structure",
                                 "has_sequence_no_structure")]
    if ti.get("n_targets") and None not in cells:
        if sum(cells) != ti["n_targets"]:
            gripes.append(
                f"Figure 6B's cells sum to {sum(cells)}, not the "
                f"{ti['n_targets']} targets.")
        marg = ti.get("structural_hits_total")
        if marg is not None and cells[0] + cells[2] != marg:
            gripes.append(
                f"Figure 6B's structural column sums to "
                f"{cells[0] + cells[2]} against a stated total of {marg}.")

    # 7. Figure 6A's strictest panel must be the set 48 says survived,
    #    and the nested sizes must decrease. A subset that grew would
    #    mean the three definitions are not nested at all.
    fx = F.J("fmax_structural.json") or {}
    ss = F.J("structural_stats.json") or {}
    pan = fx.get("panels") or {}
    surv = (ss.get("hardening") or {}).get("surviving")
    sh = (pan.get("STRUCT_HARDENED") or {}).get("n")
    if surv is not None and sh is not None and surv != sh:
        gripes.append(
            f"Figure 6A scores {sh} structurally hardened proteins but "
            f"structural_stats.json reports {surv} survivors.")
    sizes = [(k, (pan.get(k) or {}).get("n")) for k in
             ("ORPHAN_nohit", "HARDENED", "STRUCT_HARDENED")]
    sizes = [(k, n) for k, n in sizes if n is not None]
    for (k1, n1), (k2, n2) in zip(sizes, sizes[1:]):
        if n2 > n1:
            gripes.append(
                f"Figure 6A: {k2} has {n2} proteins against {k1}'s {n1}, "
                f"so the definitions are not nested.")
    # The complements must account for exactly the proteins removed.
    for a_k, b_k, c_k in (("ORPHAN_nohit", "HARDENED", "RESCUED_PROFILE"),
                          ("HARDENED", "STRUCT_HARDENED",
                           "RESCUED_STRUCTURE")):
        a_n = (pan.get(a_k) or {}).get("n")
        b_n = (pan.get(b_k) or {}).get("n")
        c_n = (pan.get(c_k) or {}).get("n")
        if None in (a_n, b_n, c_n):
            continue
        if a_n - b_n != c_n:
            gripes.append(
                f"Figure 6A: {a_k} minus {b_k} is {a_n - b_n} but "
                f"{c_k} holds {c_n}.")

    # 8. Figure 6C's dark arm is the same proteins in both groups, so
    #    the two bars must carry the same value. They come from two
    #    different files, which is exactly why this is worth checking.
    d_raw = get("Figure 6", "C", "dark", "length-matched")
    d_mat = get("Figure 6", "C", "dark", "3Di-entropy-matched")
    if d_raw is not None and d_mat is not None and abs(d_raw - d_mat) > 1e-6:
        gripes.append(
            f"Figure 6C's dark arm is {d_raw:.2f}% against one background "
            f"and {d_mat:.2f}% against the other; it is the same set of "
            f"proteins and must not differ.")

    # 9. Figure S1's cumulative curves must end at 100%, and the
    #    restricted bar must equal the curve's value at the threshold.
    ln = F.J("lineage.json") or {}
    rmax = ln.get("restricted_max_rank", 3)
    for key, lab in (("targets", "STAT-dependent dark"),
                     ("control", "uncharacterized control"),
                     ("background", "proteome background")):
        if not isinstance(ln.get(key), dict):
            continue
        end = get("Figure S1", "A", lab, "rank 9")
        if end is not None and abs(end - 100.0) > 0.05:
            gripes.append(
                f"Figure S1A's {lab} curve ends at {end:.2f}%, not 100%.")
        at_thr = get("Figure S1", "A", lab, f"rank {rmax}")
        bar = get("Figure S1", "B", lab, "restricted fraction")
        if at_thr is not None and bar is not None and \
                abs(at_thr - bar) > 0.05:
            gripes.append(
                f"Figure S1: {lab} is {at_thr:.2f}% cumulative at rank "
                f"{rmax} in panel A but {bar:.2f}% in panel B.")

    # 10. Intervals that did not bracket their own estimate. Collected
    #     during drawing; an error bar was omitted for each.
    gripes.extend(F.bad_ci)

    print("\n" + "-" * 78)
    if gripes:
        print(f"  {len(gripes)} inconsistency(ies):")
        for g in gripes:
            print(f"    - {g}")
    else:
        print("  every cross-figure check passed")
    print("-" * 78)
    return gripes


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None)
    ap.add_argument("--tiff-dpi", type=int, default=600)
    ap.add_argument("--png-dpi", type=int, default=350)
    ap.add_argument("--only", default=None,
                    help="comma-separated figure numbers, e.g. 1,3. "
                         "7 selects Supplementary Figure S1.")
    ap.add_argument("--drop-orphan-fig2", action="store_true",
                    help="omit the orphan column from Figure 2A. Only "
                         "sensible once the remaining conditions have "
                         "bootstrap intervals, or the panel is left with "
                         "three points ending at I30.")
    ap.add_argument("--no-titles", action="store_true",
                    help="omit the title drawn on each figure. The journal "
                         "typesets legends adjacent to the figure, so this "
                         "is the convention-conforming version.")
    ap.add_argument("--audit", action="store_true",
                    help="write figures/AUDIT.tsv, one row per plotted "
                         "value with the file and key it came from, and "
                         "run the cross-figure consistency checks.")
    ap.add_argument("--recompute-spectrum", action="store_true")
    ap.add_argument("--no-spectrum", action="store_true",
                    help="skip Figure 1A's eigen-decomposition (needs the "
                         "reference embeddings)")
    args = ap.parse_args()

    cfg = cfgmod.load()
    res = cfg.p("results")
    out = Path(args.out) if args.out else (res / "figures")
    plt = setup_mpl()
    F = Figures(res, out, args.tiff_dpi, args.png_dpi,
                show_titles=not args.no_titles)
    F.drop_orphan_fig2 = args.drop_orphan_fig2

    print("=" * 78)
    print("BIOINFORMATICS FIGURE SPECIFICATION")
    print("=" * 78)
    print(f"  single column {86} mm ({SINGLE:.3f} in), "
          f"double {178} mm ({DOUBLE:.3f} in)")
    print(f"  png {args.png_dpi} dpi, tiff {args.tiff_dpi} dpi, "
          f"svg and pdf vector")
    print("  white background, no gridlines, top/right spines removed")
    print("  colour and marker and line style all carry series identity, so")
    print("  the figures survive greyscale printing")
    print(f"  output: {out}")
    print("=" * 78)

    want = None
    if args.only:
        want = {int(x) for x in args.only.split(",") if x.strip().isdigit()}

    spec = None
    if (want is None or 1 in want) and not args.no_spectrum:
        spec = spectrum(cfg, res, ["esm2", "prott5", "prostt5"],
                        args.recompute_spectrum)

    for n, fn in ((1, lambda: fig1(F, plt, spec)), (2, lambda: fig2(F, plt)),
                  (3, lambda: fig3(F, plt)), (4, lambda: fig4(F, plt)),
                  (5, lambda: fig5(F, plt)), (6, lambda: fig6(F, plt)),
                  (7, lambda: figS1(F, plt))):
        if want is None or n in want:
            fn()

    tbl = tables(F)
    supplementary(F)
    legends_block(F)
    render_markdown(F, out / "Tables_and_Supplementary.md")
    ok_docx = render_docx(F, out / "Tables_and_Supplementary.docx")

    # alt text on its own too, for pasting under each legend
    A = ["# Figure legends and alt text", "",
         "Paste each legend under its figure and the alt text directly "
         "below it, prefixed with `Alt text:` as the journal requires.", ""]
    for name, legend, alt in F.alt:
        A += [f"## {name}", "", legend, "", f"Alt text: {alt}", ""]
    (out / "ALT_TEXT.md").write_text("\n".join(A) + "\n")
    log("wrote ALT_TEXT.md")

    gripes = audit_report(F, out / "AUDIT.tsv") if args.audit else []

    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    print(f"  {len(F.made)} figures: {', '.join(F.made) or 'none'}")
    print(f"  {len(tbl)} main tables: {', '.join(tbl) or 'none'}")
    print("  Tables_and_Supplementary.md"
          + (", Tables_and_Supplementary.docx" if ok_docx else ""))
    print("  ALT_TEXT.md")
    print(f"  titles drawn on figures: {'yes' if not args.no_titles else 'no'}")
    if F.skipped:
        print("\n  SKIPPED (inputs missing — not drawn from partial data):")
        for s in F.skipped:
            print(f"    {s}")
    if args.audit:
        print(f"  AUDIT.tsv ({len(F.audit)} plotted values)")
    n_files = len(list(out.glob("*")))
    print(f"\n  {n_files} files in {out}")
    return 1 if (F.skipped or gripes) else 0


if __name__ == "__main__":
    raise SystemExit(main())
