#!/usr/bin/env python3
"""
04_merge_and_figures.py -- combine every ddG result and produce the figures.

===========================================================================
TL;DR
===========================================================================
    in   results/02_foldx/fam<N>[_<tag>].csv     (from 02, optional)
         results/03_flexddg/fam<N>[_<tag>].csv   (from 03, optional)
         results/02_foldx/work_fam<N>*/repair_energies.csv
    out  results/04_analysis/all_ddg.csv        every result, long format
         results/04_analysis/correlations.csv   every number in the figures
         results/04_analysis/fig1_scatter.pdf
         results/04_analysis/fig2_summary.pdf
         results/04_analysis/fig3_agreement.pdf   (if both methods present)
         results/04_analysis/fig4_repair.pdf      (if 02 was run)

    run  python 04_merge_and_figures.py
         python 04_merge_and_figures.py --families 363
         python 04_merge_and_figures.py --flexddg-tag barlow
         python 04_merge_and_figures.py --x ddg_folding

02 and 03 can write tagged outputs (fam363_prelim.csv) so runs at
different settings coexist. Name the ones to analyse with --foldx-tag and
--flexddg-tag; both default to the untagged files. The tag is recorded in
the `protocol` column of all_ddg.csv, so a merged file stays interpretable
even when it mixes settings.

Runs on whatever exists. With only FoldX done it produces the FoldX
figures and skips the comparison.

all_ddg.csv is the file to hand to a collaborator: one row per variant per
method, in long format, with the measured value carried alongside so it
plots standalone. Failed variants keep their row with a blank ddG, so the
row count always matches the library. No licensed software is needed to
read or replot it.

===========================================================================
WHY THE CORRELATIONS ARE REPORTED SPLIT, NOT POOLED
===========================================================================
A single correlation over everything can look strong for reasons that have
nothing to do with ranking variants, so three splits are reported and the
pooled number is shown alongside rather than instead:

  by substitution count.  Singles and doubles occupy different regions of
  both axes, so pooling them measures the gap between the two groups as
  much as the ordering inside either. FoldX also repacks on a fixed
  backbone, which has least room on multi-point variants, so the split is
  the informative comparison against flex ddG.

  by position, within singles.  A method can rank POSITIONS correctly while
  ranking substitutions at a given position no better than chance. Pooled
  correlation rewards the former; choosing which substitution to make needs
  the latter. Subtracting each position's mean removes the between-position
  component and leaves the part that matters.

  by family, and pooled with each family z-scored.  Families differ in
  baseline interaction score and in ddG scale, so raw pooling mixes those
  offsets into the correlation.

===========================================================================
WHAT IS PLOTTED
===========================================================================
  fig1  scatter, methods down the rows, families across plus an all-family
        panel. Points coloured by substitution count.
  fig2  Spearman rho as bars: family x method x {all, singles, doubles},
        with a pooled panel that also shows the family-centred value.
  fig3  the two methods against each other, one panel per family.
  fig4  FoldX repair convergence: total energy per round, and the
        round-to-round change, which is what justifies the round count.

Spearman is the headline statistic throughout: the relationship between a
predicted energy and a normalised interaction score is monotonic at best,
not linear, so rank correlation is the honest measure. Pearson is reported
in the tables for completeness.

===========================================================================
REQUIRES
===========================================================================
pandas, numpy, matplotlib. No external tools; runs from the CSVs alone.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas not found.  pip install pandas")

METHOD_COLOUR = {"foldx": "#dd6b20", "flexddg": "#2b6cb0"}
NMUT_COLOUR = {1: "#2b6cb0", 2: "#dd6b20", 3: "#38a169"}


# ---------------------------------------------------------------- config --
def load_env(path="paths.env"):
    """Read KEY=value lines. The file is also valid bash."""
    cfg, p = {}, Path(path)
    if not p.is_file():
        sys.exit(f"ERROR: {p} not found. cp paths.template.env paths.env")
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        k, _, v = line.partition("=")
        cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg


# ------------------------------------------------------------ statistics --
def rho(x, y):
    """Spearman, on the pairs where both values are finite."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan
    return float(np.corrcoef(pd.Series(x[ok]).rank(),
                             pd.Series(y[ok]).rank())[0, 1])


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def boot_ci(x1, x2, y, n=2000, seed=0):
    """Bootstrap CI on rho(x1,y) - rho(x2,y), resampling variants.

    With a few hundred variants a rho gap of 0.1 is comfortably inside the
    noise, so a difference between methods needs an interval, not a
    point estimate.
    """
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y))
    diffs = []
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        d = rho(x1[s], y[s]) - rho(x2[s], y[s])
        if np.isfinite(d):
            diffs.append(d)
    if not diffs:
        return np.nan, np.nan, np.nan
    d = np.array(diffs)
    return (float(d.mean()), float(np.percentile(d, 2.5)),
            float(np.percentile(d, 97.5)))


# ----------------------------------------------------------------- input --
def load_results(cfg, families, args):
    """Concatenate every matching fam<N>.csv, from both methods.

    The tag lets one results directory hold several settings; whichever is
    named here is what gets analysed, and the name is carried into the
    output so the provenance of every row is visible.
    """
    root = Path(cfg.get("RESULTS_DIR", "results"))
    tags = {"foldx": args.foldx_tag, "flexddg": args.flexddg_tag}
    frames, found = [], []
    for method, sub in (("foldx", "02_foldx"), ("flexddg", "03_flexddg")):
        tag = tags[method]
        suffix = f"_{tag}" if tag else ""
        for fam in families:
            path = root / sub / f"fam{fam}{suffix}.csv"
            if not path.is_file():
                # a helpful nudge when only tagged files exist
                alts = sorted(p.name for p in (root / sub).glob(
                    f"fam{fam}_*.csv") if not p.name.endswith("_perjob.csv"))
                if alts and not tag:
                    print(f"  {method:<8} {fam}: fam{fam}.csv not found, but "
                          f"{', '.join(alts)} exist -- name one with "
                          f"--{method}-tag")
                continue
            df = pd.read_csv(path)
            df["protocol"] = tag or "default"
            for col in ("ddg_binding", "ddg_folding", "ppi_exp",
                        "ddg_binding_sd", "ddg_folding_sd", "n_mut"):
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")
            df["family"] = df["family"].astype(str)
            df["method"] = method
            frames.append(df)
            ok = (df["status"] == "ok").sum()
            found.append(f"  {method:<8} {fam}: {len(df)} row(s), "
                         f"{ok} ok, {len(df) - ok} failed   ({path})")
    for line in found:
        print(line)
    if not frames:
        sys.exit(f"\nno results found under {root}/02_foldx or "
                 f"{root}/03_flexddg -- run 02 or 03 first")
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------- tables ---
def correlation_table(data, xcol, ycol):
    """Spearman and Pearson for every family x method x subset."""
    recs = []
    subsets = [("all", None), ("singles", 1), ("doubles", 2)]
    for family in sorted(data["family"].unique()) + ["pooled"]:
        for method in sorted(data["method"].unique()):
            base = data if family == "pooled" else \
                data[data["family"] == family]
            base = base[base["method"] == method]
            for name, nm in subsets:
                g = base if nm is None else base[base["n_mut"] == nm]
                g = g.dropna(subset=[xcol, ycol])
                if len(g) < 3:
                    continue
                recs.append({"family": family, "method": method,
                             "subset": name, "n": len(g),
                             "spearman": round(rho(g[xcol], g[ycol]), 4),
                             "pearson": round(pearson(g[xcol], g[ycol]), 4)})
                # pooling families mixes their different baselines into the
                # correlation; z-scoring each family first removes that
                if family == "pooled":
                    z = g.copy()
                    for c in (xcol, ycol):
                        z[c] = z.groupby("family")[c].transform(
                            lambda v: (v - v.mean()) /
                            (v.std() if v.std() else 1))
                    recs.append({"family": "pooled (family-centred)",
                                 "method": method, "subset": name,
                                 "n": len(z),
                                 "spearman": round(rho(z[xcol], z[ycol]), 4),
                                 "pearson": round(pearson(z[xcol], z[ycol]),
                                                  4)})
    return pd.DataFrame(recs)


def position_table(data, xcol, ycol):
    """Within-position correlation for single mutants.

    Removes the between-position component, leaving the ranking of
    substitutions at a given position -- the harder and more useful task.
    """
    recs = []
    singles = data[data["n_mut"] == 1].copy()
    if singles.empty:
        return pd.DataFrame()
    singles["position"] = pd.to_numeric(
        singles["mutations"].astype(str).str.extract(r"(\d+)")[0],
        errors="coerce")
    for family in sorted(singles["family"].unique()):
        for method in sorted(singles["method"].unique()):
            g = singles[(singles["family"] == family) &
                        (singles["method"] == method)]
            g = g.dropna(subset=[xcol, ycol, "position"])
            if len(g) < 10:
                continue
            pooled = rho(g[xcol], g[ycol])
            c = g.copy()
            for col in (xcol, ycol):
                c[col] = c[col] - c.groupby("position")[col].transform("mean")
            recs.append({"family": family, "method": method, "n": len(g),
                         "spearman_pooled": round(pooled, 4),
                         "spearman_within_position": round(
                             rho(c[xcol], c[ycol]), 4)})
    return pd.DataFrame(recs)


# --------------------------------------------------------------- figures --
def scatter_panel(ax, sub, xcol, ycol, title, xlabel, show_family=False,
                  xlim=None, ylim=None):
    if sub.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="#999999")
        ax.set_xticks([])
        ax.set_yticks([])
        return
    if show_family:
        marks = dict(zip(sorted(sub["family"].unique()),
                         ["o", "s", "^", "D", "v"]))
        for fam in sorted(sub["family"].unique()):
            for nm in sorted(sub["n_mut"].dropna().unique()):
                s = sub[(sub["family"] == fam) & (sub["n_mut"] == nm)]
                if not s.empty:
                    ax.scatter(s[xcol], s[ycol], s=13, alpha=0.6,
                               marker=marks[fam],
                               color=NMUT_COLOUR.get(int(nm), "#805ad5"),
                               edgecolors="none",
                               label=f"{fam} {int(nm)}-point")
    else:
        for nm in sorted(sub["n_mut"].dropna().unique()):
            s = sub[sub["n_mut"] == nm]
            ax.scatter(s[xcol], s[ycol], s=16, alpha=0.7,
                       color=NMUT_COLOUR.get(int(nm), "#805ad5"),
                       edgecolors="none",
                       label=f"{int(nm)}-point (n={len(s)})")
    x = sub[xcol].to_numpy(float)
    y = sub[ycol].to_numpy(float)
    if len(x) >= 3:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, np.polyval(coef, xs), color="#333333", linewidth=1.1)
    ax.set_title(f"{title}\nn={len(sub)}  rho={rho(x, y):+.3f}  "
                 f"r={pearson(x, y):+.3f}", fontsize=10)
    ax.set_xlabel(xlabel, fontsize=9)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.axvline(0, color="#e8e8e8", linewidth=0.8, zorder=0)
    ax.legend(fontsize=5.5 if show_family else 6.5, frameon=False,
              ncol=2 if show_family else 1)
    ax.tick_params(labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def shared_limits(values, pad=0.04):
    """Axis limits covering every panel, so they are visually comparable.

    Panels that scale themselves invite the eye to compare shapes that are
    drawn at different magnifications; a point at the right edge of one
    panel and the right edge of another then mean different numbers. The
    limits come from the pooled data and are applied to every panel.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        lo, hi = lo - 0.5, hi + 0.5
    m = (hi - lo) * pad
    return lo - m, hi + m


def fig_scatter(data, xcol, ycol, families, methods, out, share_x=True,
                share_y=True, marginals=True):
    """Scatter grid, optionally with marginal distributions.

    The marginals answer questions the scatter alone hides: whether the
    predictions pile up against a bound, whether the measured values are
    bimodal, and whether the apparent spread comes from a handful of points
    at one edge. A correlation computed on a target with two separated
    clusters is largely reporting cluster membership, and that is visible in
    the histogram long before it is visible in the cloud.

    Each panel becomes a 2x2 block: the scatter, a histogram of x below it,
    and a histogram of y to its right, sharing the scatter's axes.
    """
    cols = list(families) + (["all"] if len(families) > 1 else [])
    ylim = shared_limits(data[ycol]) if share_y else None
    xlims = {m: shared_limits(data[data["method"] == m][xcol])
             for m in methods} if share_x else {m: None for m in methods}

    if not marginals:
        fig, axes = plt.subplots(len(methods), len(cols),
                                 figsize=(4.4 * len(cols),
                                          4.2 * len(methods)),
                                 squeeze=False)
        for i, method in enumerate(methods):
            for j, fam in enumerate(cols):
                g = data[data["method"] == method]
                if fam != "all":
                    g = g[g["family"] == fam]
                g = g.dropna(subset=[xcol, ycol])
                scatter_panel(axes[i][j], g, xcol, ycol,
                              f"{'all families' if fam == 'all' else fam}"
                              f" -- {method}",
                              f"{method} {xcol.replace('_', ' ')}",
                              show_family=(fam == "all"),
                              xlim=xlims.get(method), ylim=ylim)
                if j == 0:
                    axes[i][j].set_ylabel(ycol, fontsize=9)
        fig.tight_layout()
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    fig = plt.figure(figsize=(4.9 * len(cols), 4.8 * len(methods)))
    outer = fig.add_gridspec(len(methods), len(cols), hspace=0.35,
                             wspace=0.3)

    for i, method in enumerate(methods):
        for j, fam in enumerate(cols):
            g = data[data["method"] == method]
            if fam != "all":
                g = g[g["family"] == fam]
            g = g.dropna(subset=[xcol, ycol])

            inner = outer[i, j].subgridspec(
                2, 2, width_ratios=(5, 1.3), height_ratios=(1.3, 5),
                hspace=0.05, wspace=0.05)
            ax = fig.add_subplot(inner[1, 0])
            ax_top = fig.add_subplot(inner[0, 0], sharex=ax)
            ax_right = fig.add_subplot(inner[1, 1], sharey=ax)

            scatter_panel(ax, g, xcol, ycol,
                          f"{'all families' if fam == 'all' else fam}"
                          f" -- {method}",
                          f"{method} {xcol.replace('_', ' ')}",
                          show_family=(fam == "all"),
                          xlim=xlims.get(method), ylim=ylim)
            if j == 0:
                ax.set_ylabel(ycol, fontsize=9)

            # the title belongs above the marginal, not on top of the
            # scatter, or the two overlap
            ax_top.set_title(ax.get_title(), fontsize=10)
            ax.set_title("")

            if not g.empty:
                colour = METHOD_COLOUR.get(method, "#666666")
                x = g[xcol].to_numpy(float)
                y = g[ycol].to_numpy(float)
                ax_top.hist(x[np.isfinite(x)], bins=30, color=colour,
                            alpha=0.65, edgecolor="none")
                ax_right.hist(y[np.isfinite(y)], bins=30, color=colour,
                              alpha=0.65, edgecolor="none",
                              orientation="horizontal")
                # the median marks where the bulk sits, which is what makes
                # a skewed or bimodal distribution obvious
                ax_top.axvline(np.nanmedian(x), color="#333333",
                               linewidth=0.9)
                ax_right.axhline(np.nanmedian(y), color="#333333",
                                 linewidth=0.9)

            # sharex/sharey propagates tick visibility, so the marginals
            # need their LABELS hidden rather than their ticks removed --
            # calling set_xticks([]) on a shared axis strips them from the
            # scatter as well
            for a in (ax_top, ax_right):
                for side in ("top", "right", "left", "bottom"):
                    a.spines[side].set_visible(False)
            ax_top.tick_params(axis="x", labelbottom=False, length=0)
            ax_top.tick_params(axis="y", labelleft=False, length=0)
            ax_right.tick_params(axis="y", labelleft=False, length=0)
            ax_right.tick_params(axis="x", labelbottom=False, length=0)
            ax.tick_params(axis="both", labelsize=8, length=3)

    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_summary(corr, families, methods, out):
    panels = list(families) + ["pooled"]
    labels = ["all", "singles", "doubles"]
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(4.2 * len(panels), 4.4), squeeze=False)
    width = 0.8 / max(1, len(methods))
    for ax, fam in zip(axes[0], panels):
        for k, method in enumerate(methods):
            g = corr[(corr["family"] == fam) & (corr["method"] == method)]
            vals, pos = [], []
            for xi, lab in enumerate(labels):
                r = g[g["subset"] == lab]
                if len(r):
                    vals.append(float(r["spearman"].iloc[0]))
                    pos.append(xi + (k - (len(methods) - 1) / 2) * width)
            if vals:
                ax.bar(pos, vals, width=width, label=method,
                       color=METHOD_COLOUR.get(method, "#666666"))
            if fam == "pooled":
                gc = corr[(corr["family"] == "pooled (family-centred)") &
                          (corr["method"] == method)]
                cv, cp = [], []
                for xi, lab in enumerate(labels):
                    r = gc[gc["subset"] == lab]
                    if len(r):
                        cv.append(float(r["spearman"].iloc[0]))
                        cp.append(xi + (k - (len(methods) - 1) / 2) * width)
                if cv:
                    ax.bar(cp, cv, width=width, facecolor="none",
                           edgecolor=METHOD_COLOUR.get(method, "#666666"),
                           hatch="///", linewidth=1.0,
                           label=f"{method}, family-centred")
        ax.axhline(0, color="#555555", linewidth=0.9)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_ylim(-1, 1)
        ax.set_title("all families" if fam == "pooled" else f"family {fam}",
                     fontsize=10)
        if fam == panels[0]:
            ax.set_ylabel("Spearman rho", fontsize=9)
        ax.legend(fontsize=6.5, frameon=False)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_agreement(data, xcol, families, out):
    wide = data.pivot_table(index=["family", "sequence", "n_mut"],
                            columns="method", values=xcol).reset_index()
    if not {"foldx", "flexddg"} <= set(wide.columns):
        return None
    wide = wide.dropna(subset=["foldx", "flexddg"])
    if wide.empty:
        return None
    fams = [f for f in families if f in set(wide["family"])]
    fig, axes = plt.subplots(1, max(1, len(fams)),
                             figsize=(4.2 * max(1, len(fams)), 4.2),
                             squeeze=False)
    # one square range on both axes, so the diagonal means the same thing
    # in every panel
    both = shared_limits(np.concatenate([wide["foldx"].to_numpy(float),
                                         wide["flexddg"].to_numpy(float)]))
    for ax, fam in zip(axes[0], fams):
        g = wide[wide["family"] == fam]
        for nm in sorted(g["n_mut"].dropna().unique()):
            s = g[g["n_mut"] == nm]
            ax.scatter(s["flexddg"], s["foldx"], s=16, alpha=0.7,
                       color=NMUT_COLOUR.get(int(nm), "#805ad5"),
                       edgecolors="none", label=f"{int(nm)}-point")
        lim = list(both) if both else [
            min(g["flexddg"].min(), g["foldx"].min()),
            max(g["flexddg"].max(), g["foldx"].max())]
        ax.plot(lim, lim, "--", color="#aaaaaa", linewidth=0.9)
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"family {fam}\nrho={rho(g['flexddg'], g['foldx']):+.3f}"
                     f"  n={len(g)}", fontsize=10)
        ax.set_xlabel("flex ddG (kcal/mol)", fontsize=9)
        if fam == fams[0]:
            ax.set_ylabel("FoldX (kcal/mol)", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return wide


def fig_repair(cfg, families, out):
    """FoldX repair convergence, the justification for the round count."""
    root = Path(cfg.get("RESULTS_DIR", "results")) / "02_foldx"
    traces = {}
    for fam in families:
        hits = sorted(root.glob(f"work_fam{fam}*/repair_energies.csv"))
        if not hits:
            continue
        for f in hits[:1]:
            d = pd.read_csv(f)
            d = d.dropna(subset=["total_energy"])
            if len(d):
                traces[fam] = d
    if not traces:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    for fam, d in sorted(traces.items()):
        axes[0].plot(d["round"], d["total_energy"], "o-", markersize=4,
                     linewidth=1.3, label=f"family {fam}")
        delta = d["total_energy"].diff().abs()
        axes[1].plot(d["round"][1:], delta[1:], "o-", markersize=4,
                     linewidth=1.3, label=f"family {fam}")
    axes[0].set_ylabel("FoldX total energy (kcal/mol)", fontsize=9)
    axes[0].set_title("repair convergence", fontsize=10)
    axes[1].set_ylabel("|change from previous round| (kcal/mol)", fontsize=9)
    axes[1].set_yscale("log")
    axes[1].set_title("round-to-round change", fontsize=10)
    rounds = sorted({int(r) for d in traces.values() for r in d["round"]})
    for ax in axes:
        ax.set_xlabel("RepairPDB round", fontsize=9)
        ax.set_xticks(rounds)
        ax.set_xlim(min(rounds) - 0.2, max(rounds) + 0.2)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(axis="y", color="#eeeeee", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return traces


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paths", default="paths.env")
    ap.add_argument("--families", default=None)
    ap.add_argument("--out-dir", default=None,
                    help="default: RESULTS_DIR/04_analysis")
    ap.add_argument("--x", dest="xcol", default="ddg_binding",
                    help="which prediction to correlate. ddg_binding "
                         "(default) or ddg_folding, or the name of any "
                         "column in the metrics table, e.g. "
                         "'fold::Backbone Hbond'. A partial name is "
                         "accepted when it matches exactly one metric")
    ap.add_argument("--metrics-csv", default=None,
                    help="a metrics table from 05_list_foldx_metrics.py, "
                         "merged in so any FoldX energy term can be plotted "
                         "with the same figures (default: "
                         "RESULTS_DIR/05_metrics/metrics.csv when it exists)")
    ap.add_argument("--y", dest="ycol", default="ppi_exp")
    ap.add_argument("--no-marginals", action="store_true",
                    help="omit the distribution histograms beside each "
                         "scatter panel")
    ap.add_argument("--free-axes", action="store_true",
                    help="let every panel scale itself. Off by default: "
                         "shared limits keep the panels comparable by eye")
    ap.add_argument("--foldx-tag", default=None,
                    help="which 02 output to use, e.g. --foldx-tag vdw0 "
                         "reads fam<N>_vdw0.csv (default: fam<N>.csv)")
    ap.add_argument("--flexddg-tag", default=None,
                    help="which 03 output to use, e.g. --flexddg-tag barlow "
                         "reads fam<N>_barlow.csv (default: fam<N>.csv)")
    ap.add_argument("--no-bootstrap", action="store_true")
    args = ap.parse_args()

    cfg = load_env(args.paths)
    families = [f.strip() for f in
                (args.families or cfg.get("FAMILIES", "")).split(",")
                if f.strip()]
    out_dir = Path(args.out_dir or
                   Path(cfg.get("RESULTS_DIR", "results")) / "04_analysis")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("inputs")
    data = load_results(cfg, families, args)

    # Optionally bring in the FoldX energy decomposition, so a single term
    # can be plotted with the same figures rather than a parallel script.
    # Terms only exist for FoldX, so flex ddG rows carry NaN and drop out.
    mpath = Path(args.metrics_csv) if args.metrics_csv else \
        Path(cfg.get("RESULTS_DIR", "results")) / "05_metrics" / "metrics.csv"
    if mpath.is_file():
        mt = pd.read_csv(mpath)
        cols = [c for c in mt.columns if "::" in c]
        if cols:
            mt["family"] = mt["family"].astype(str)
            mt["sequence"] = mt["sequence"].astype(str).str.strip().str.upper()
            data["sequence"] = data["sequence"].astype(str).str.strip().str.upper()
            data = data.merge(mt[["family", "sequence"] + cols],
                              on=["family", "sequence"], how="left")
            print(f"  metrics  : {len(cols)} term(s) from {mpath.name}")

    # resolve a partial metric name, rather than failing on a long one
    if args.xcol not in data.columns:
        hits = [c for c in data.columns
                if args.xcol.lower() in c.lower()]
        if len(hits) == 1:
            print(f"  --x '{args.xcol}' -> '{hits[0]}'")
            args.xcol = hits[0]
        elif len(hits) > 1:
            sys.exit(f"ERROR: '{args.xcol}' matches {len(hits)} columns:\n  "
                     + "\n  ".join(hits))
        else:
            avail = [c for c in data.columns if "::" in c]
            sys.exit(f"ERROR: '{args.xcol}' is not a column.\n"
                     f"  ddg_binding, ddg_folding"
                     + (f", or one of {len(avail)} metrics:\n  "
                        + "\n  ".join(avail[:20]) if avail else "")
                     + ("\n  ..." if len(avail) > 20 else ""))
        for c in (args.xcol,):
            data[c] = to_num(data[c])

    # keep only rows that carry both a prediction and a measurement
    usable = data.dropna(subset=[args.xcol, args.ycol])
    dropped = len(data) - len(usable)
    print(f"\n{len(data)} row(s) loaded, {len(usable)} with both "
          f"{args.xcol} and {args.ycol}"
          + (f", {dropped} dropped" if dropped else ""))

    # sort before numbering, so `n` is reproducible rather than an artefact
    # of the order the input files happened to be read in
    sort_keys = [c for c in ("family", "method", "sequence")
                 if c in data.columns]
    data = data.sort_values(sort_keys).reset_index(drop=True)
    data.insert(0, "n", range(1, len(data) + 1))
    data.to_csv(out_dir / "all_ddg.csv", index=False)
    ok = (data["status"] == "ok").sum() if "status" in data.columns else \
        len(data)
    print(f"wrote {out_dir / 'all_ddg.csv'}   {len(data)} row(s), "
          f"{ok} ok, {len(data) - ok} failed")
    if "protocol" in data.columns:
        seen = data.groupby("method")["protocol"].unique()
        for method, tags in seen.items():
            print(f"  {method}: protocol {', '.join(sorted(tags))}")

    present = [f for f in families if f in set(usable["family"])]
    methods = sorted(usable["method"].unique())

    corr = correlation_table(usable, args.xcol, args.ycol)
    corr.to_csv(out_dir / "correlations.csv", index=False)
    print("\n" + corr.to_string(index=False))

    pos = position_table(usable, args.xcol, args.ycol)
    if not pos.empty:
        pos.to_csv(out_dir / "within_position.csv", index=False)
        print("\nsingle mutants, position-centred:")
        print(pos.to_string(index=False))
        print("  A pooled value that collapses when centred means the "
              "method ranks\n  positions rather than substitutions.")

    # method comparison, on the variants both methods computed
    if len(methods) == 2:
        wide = usable.pivot_table(
            index=["family", "sequence", "n_mut"], columns="method",
            values=[args.xcol, args.ycol]).reset_index()
        wide.columns = ["_".join(c).strip("_") if isinstance(c, tuple) else c
                        for c in wide.columns]
        xa, xb = f"{args.xcol}_foldx", f"{args.xcol}_flexddg"
        ya = f"{args.ycol}_foldx"
        if {xa, xb, ya} <= set(wide.columns):
            both = wide.dropna(subset=[xa, xb, ya])
            if len(both) >= 10 and not args.no_bootstrap:
                m, lo, hi = boot_ci(both[xb].to_numpy(float),
                                    both[xa].to_numpy(float),
                                    both[ya].to_numpy(float))
                print(f"\nhead to head on the {len(both)} variant(s) both "
                      f"methods computed:")
                print(f"  flex ddG rho {rho(both[xb], both[ya]):+.3f}   "
                      f"FoldX rho {rho(both[xa], both[ya]):+.3f}")
                print(f"  difference {m:+.3f}  95% CI "
                      f"[{lo:+.3f}, {hi:+.3f}]")
                print("  " + ("the interval spans zero: on this many "
                              "variants the two are\n  not distinguishable"
                              if lo < 0 < hi else
                              "the interval excludes zero: the difference "
                              "survives resampling"))

    print("\nfigures")
    fig_scatter(usable, args.xcol, args.ycol, present, methods,
                out_dir / "fig1_scatter.pdf",
                share_x=not args.free_axes, share_y=not args.free_axes,
                marginals=not args.no_marginals)
    print(f"  {out_dir / 'fig1_scatter.pdf'}")
    fig_summary(corr, present, methods, out_dir / "fig2_summary.pdf")
    print(f"  {out_dir / 'fig2_summary.pdf'}")
    if len(methods) == 2:
        if fig_agreement(usable, args.xcol, present,
                         out_dir / "fig3_agreement.pdf") is not None:
            print(f"  {out_dir / 'fig3_agreement.pdf'}")
    if fig_repair(cfg, families, out_dir / "fig4_repair.pdf") is not None:
        print(f"  {out_dir / 'fig4_repair.pdf'}")

    print(f"\nall_ddg.csv is the file to share: one row per variant per "
          f"method,\nno licensed software needed to read or replot it.")


if __name__ == "__main__":
    main()
