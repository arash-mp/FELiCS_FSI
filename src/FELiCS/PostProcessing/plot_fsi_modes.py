#!/usr/bin/env python3
"""
plot_fsi_modes.py -- post-processing / visualisation for FELiCS rigid-body FSI.

Reads ``fsi_structural_modes.csv`` (written by the FSI branch of run_modal.py) and
produces a multi-panel figure with everything you need to read off a run:

  A. Spectrum            eigenvalues in the complex plane, coloured by which DOF
                         dominates, sized by structural content, with the stability
                         line and the unstable half-plane shaded.
  B. Fluid vs structure  struct_ratio_energy per mode -- separates wake modes from
                         structural branches. THIS IS NOT the DOF split.
  C. DOF participation   stacked bars: within the structural part, which DOFs move.
  D. Inter-DOF phase     phase of each DOF relative to DOF 0, with the +/-90 deg
                         flutter-favourable band marked.
  E. Trust panel         solver residual + the exact consistency identity
                         1j*phi_k = lambda*eta_k.  Filter on this BEFORE believing
                         anything else.
  F. Verdict table       per-mode classification in words.

Panels A, B and E cover EVERY mode in the file. Panels C, D and F are per-mode and
would be unreadable for hundreds of modes, so they default to the 3 most unstable
modes; change with --n-modes.

Usage
-----
    python plot_fsi_modes.py fsi_structural_modes.csv
    python plot_fsi_modes.py fsi_structural_modes.csv --n-modes 6
    python plot_fsi_modes.py out/fsi_structural_modes.csv --n-modes 5 \
        --error-tol 1e-8 --hide-drift -o modes.png

Options
-------
    --n-modes N     how many modes in the per-mode panels (default 3)
    --error-tol X   grey out / exclude modes with residual > X (default: keep all,
                    but always flag them)
    --hide-drift    hide the |lambda| ~ 0 rigid-body drift modes (K=0 free body)
    --sort-by       'growth' (default) | 'struct' | 'freq'
    -o FILE         save to FILE instead of showing interactively
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# --- locate fsi_common whether we are run as a script, via a symlink on PATH,
# --- or as "python -m FELiCS.PostProcessing.plot_fsi_modes"
try:
    from .fsi_common import (parse_csv, is_rotation, consistency_residual,
                             classify, report_health)
except ImportError:                     # run directly, or via a symlink on PATH
    import os, sys
    _here = os.path.dirname(os.path.realpath(__file__))   # realpath resolves symlinks
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from fsi_common import (parse_csv, is_rotation, consistency_residual,
                            classify, report_health)

# ----------------------------------------------------------------------
# CSV parsing: discover the DOFs from the column names
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Physics helpers
# ----------------------------------------------------------------------




# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def make_figure(d, n_show=3, error_tol=None, hide_drift=False,
                sort_by="growth", title=None):
    n_dof  = d["n_dof"]
    labels = d["labels"]
    lam    = d["omega"]
    grow   = lam.imag
    freq   = lam.real
    resid, is_noise = consistency_residual(d)

    # ---- mode selection for the per-mode panels --------------------------
    mask = np.ones(d["n_modes"], dtype=bool)
    if hide_drift:
        mask &= np.abs(lam) > 1e-8
    if error_tol is not None:
        mask &= np.nan_to_num(d["error"], nan=0.0) <= error_tol

    idx_all = np.where(mask)[0]
    if len(idx_all) == 0:
        sys.exit("No modes left after filtering. Loosen --error-tol / --hide-drift.")

    if sort_by == "growth":
        order = idx_all[np.argsort(-grow[idx_all])]
        sel_desc = "most unstable"
    elif sort_by == "struct":
        order = idx_all[np.argsort(-np.nan_to_num(d["struct_energy"][idx_all]))]
        sel_desc = "most structural"
    else:
        order = idx_all[np.argsort(freq[idx_all])]
        sel_desc = "lowest frequency"
    sel = order[:n_show]

    # ---- figure layout ---------------------------------------------------
    fig = plt.figure(figsize=(16.5, 11))
    gs  = GridSpec(3, 3, figure=fig, hspace=0.42, wspace=0.28,
                   height_ratios=[1.25, 1.0, 0.85])

    ax_spec = fig.add_subplot(gs[0, 0:2])
    ax_sre  = fig.add_subplot(gs[0, 2])
    ax_part = fig.add_subplot(gs[1, 0])
    ax_phase= fig.add_subplot(gs[1, 1])
    ax_trust= fig.add_subplot(gs[1, 2])
    ax_tab  = fig.add_subplot(gs[2, :])

    dof_colors = plt.cm.tab10(np.arange(n_dof) % 10)

    # ================= A. spectrum ========================================
    ax = ax_spec
    xlo, xhi = freq.min() - 0.15, freq.max() + 0.15
    ax.axhspan(0, max(grow.max()*1.35, 0.05), color="red", alpha=0.05, zorder=0)
    ax.axhline(0, color="k", lw=1.1, ls="--", zorder=1)
    ax.text(xhi, 0.002, " unstable", va="bottom", ha="right",
            color="red", fontsize=9, alpha=0.8)

    # size ~ structural content, colour = dominant DOF
    sre = np.nan_to_num(d["struct_energy"], nan=0.0)
    size = 40 + 420 * sre
    for k in range(n_dof):
        m = (d["dominant"] == k) & mask
        ax.scatter(freq[m], grow[m], s=size[m], color=dof_colors[k],
                   edgecolor="k", linewidth=0.6, alpha=0.85,
                   label=f"{labels[k]}-dominant", zorder=3)
    # drift modes
    dm = (np.abs(lam) < 1e-8)
    if dm.any() and not hide_drift:
        ax.scatter(freq[dm], grow[dm], s=110, facecolor="none",
                   edgecolor="grey", linewidth=1.6, marker="s",
                   label="drift mode (K=0)", zorder=4)
    # highlight selected
    ax.scatter(freq[sel], grow[sel], s=size[sel] + 190, facecolor="none",
               edgecolor="red", linewidth=1.7, zorder=5, label=f"top {len(sel)}")
    for rank, i in enumerate(sel):
        ax.annotate(f"#{rank+1}", (freq[i], grow[i]),
                    textcoords="offset points", xytext=(9, 8),
                    fontsize=10, fontweight="bold", color="red")

    ax.set_xlabel(r"$\omega_r=\mathrm{Re}(\lambda)$   (frequency)")
    ax.set_ylabel(r"$\omega_i=\mathrm{Im}(\lambda)$   (growth rate)")
    ax.set_title("A. Spectrum — colour = dominant DOF, size = structural content",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.grid(alpha=0.25)

    # ================= B. fluid vs structure ==============================
    ax = ax_sre
    sre_pos = sre[mask]
    sre_pos = sre_pos[sre_pos > 0]
    # If every mode is fluid-dominated (as happens when the eigensolver never
    # went near the structural branches), a linear axis squashes everything onto
    # zero and the panel is useless.  Switch to log so the modes stay separable.
    use_log = (sre.max() < 0.05) and (len(sre_pos) > 0)

    if use_log:
        lo = max(sre_pos.min() * 0.5, 1e-12)
        ax.set_yscale("log")
        ax.set_ylim(lo, 1.5)
        ax.axhspan(lo, 0.05, color="tab:blue",   alpha=0.10)
    else:
        ax.set_ylim(-0.02, 1.03)
        ax.axhspan(0, 0.05, color="tab:blue",   alpha=0.10)
    ax.axhspan(0.05, 0.30, color="tab:purple", alpha=0.10)
    ax.axhspan(0.30, 1.50 if use_log else 1.02, color="tab:orange", alpha=0.10)

    ax.scatter(freq[mask], np.maximum(sre[mask], 1e-30) if use_log else sre[mask],
               s=34, c="k", alpha=0.65, zorder=3)
    ax.scatter(freq[sel], np.maximum(sre[sel], 1e-30) if use_log else sre[sel],
               s=95, facecolor="none", edgecolor="red", linewidth=1.5, zorder=4)
    # label the selected modes, matching the #n numbering of panel A / table F
    for rank, i in enumerate(sel):
        y = max(sre[i], 1e-30) if use_log else sre[i]
        ax.annotate(f"#{rank+1}", (freq[i], y),
                    textcoords="offset points", xytext=(8, 6),
                    fontsize=9.5, fontweight="bold", color="red", zorder=6)
    for txt, y in [("fluid / wake", 0.02), ("coupled", 0.16), ("structural", 0.6)]:
        ax.text(0.98, y, txt, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8, alpha=0.75)
    if use_log:
        ax.text(0.02, 0.97, "log scale: ALL modes are fluid-dominated",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color="red", fontweight="bold")
    ax.set_xlabel(r"$\omega_r$")
    ax.set_ylabel("struct_ratio_energy")
    ax.set_title("B. Fluid vs structure\n(NOT the DOF split)",
                 fontsize=11, fontweight="bold", loc="left")
    ax.grid(alpha=0.25, which="both" if use_log else "major")

    # ================= C. DOF participation ===============================
    ax = ax_part
    x = np.arange(len(sel))
    bottom = np.zeros(len(sel))
    for k in range(n_dof):
        vals = np.array([d["frac"][k][i] for i in sel])
        ax.bar(x, vals, bottom=bottom, color=dof_colors[k],
               edgecolor="k", linewidth=0.5, label=labels[k])
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0.07:
                ax.text(xi, b + v/2, f"{v:.2f}", ha="center", va="center",
                        fontsize=8.5, fontweight="bold",
                        color="w" if v > 0.25 else "k")
        bottom += vals
    equi = 1.0 / n_dof
    ax.axhline(equi, color="k", ls=":", lw=1.3)
    ax.text(len(sel) - 0.45, equi, f" equipartition ({equi:.2f})",
            fontsize=8, va="bottom", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([f"#{r+1}" for r in range(len(sel))])
    ax.set_ylim(0, 1.16)
    ax.set_ylabel("energy participation $p_k$")
    ax.set_title(f"C. Which DOFs move ({sel_desc})",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8, ncol=min(n_dof, 3), loc="upper center", framealpha=0.9)

    # ================= D. inter-DOF phase =================================
    ax = ax_phase
    if n_dof > 1:
        for band in (90, -90):
            ax.axhspan(band - 35, band + 35, color="red", alpha=0.10)
        ax.axhline(90, color="red", ls="--", lw=1.1)
        ax.axhline(-90, color="red", ls="--", lw=1.1)
        ax.axhline(0, color="k", lw=0.8, alpha=0.5)
        w = 0.8 / max(n_dof - 1, 1)
        for k in range(1, n_dof):
            ph = np.array([d["phase"][k][i] for i in sel])
            ax.bar(x + (k-1)*w - 0.4 + w/2, ph, width=w*0.9,
                   color=dof_colors[k], edgecolor="k", linewidth=0.5,
                   label=f"{labels[k]} vs {labels[0]}")
            for xi, v in enumerate(ph):
                ax.text(xi + (k-1)*w - 0.4 + w/2,
                        v + (7 if v >= 0 else -12),
                        f"{v:.0f}°", ha="center", fontsize=8)
        ax.text(0.99, 0.97, "±90° = flutter-favourable",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="red")
        ax.set_ylim(-195, 195)
        ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_ylabel("phase rel. to DOF 0  [deg]")
        ax.legend(fontsize=8, loc="lower right")
    else:
        ax.text(0.5, 0.5, "single DOF:\nno inter-DOF phase",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_yticks([])
    ax.set_xticks(x)
    ax.set_xticklabels([f"#{r+1}" for r in range(len(sel))])
    ax.set_title("D. Inter-DOF phase", fontsize=11, fontweight="bold", loc="left")
    ax.grid(alpha=0.2, axis="y")

    # ================= E. trust panel =====================================
    ax = ax_trust
    err = np.nan_to_num(d["error"], nan=0.0)
    ax.scatter(np.maximum(err[mask], 1e-16), np.maximum(resid[mask], 1e-16),
               s=34, c="k", alpha=0.6, zorder=3, label="all modes")
    ax.scatter(np.maximum(err[sel], 1e-16), np.maximum(resid[sel], 1e-16),
               s=95, facecolor="none", edgecolor="red", linewidth=1.5,
               zorder=4, label=f"top {len(sel)}")
    ax.axvline(1e-6, color="red", ls="--", lw=1.1)
    ax.axhline(1e-6, color="red", ls="--", lw=1.1)
    ax.fill_betweenx([1e-16, 1e-6], 1e-16, 1e-6, color="green", alpha=0.07)
    ax.text(0.03, 0.05, "trustworthy", transform=ax.transAxes,
            fontsize=8.5, color="green")
    ax.text(0.97, 0.95, "DO NOT TRUST", transform=ax.transAxes,
            fontsize=8.5, color="red", ha="right", va="top")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("SLEPc residual (error)")
    ax.set_ylabel(r"$|i\hat\varphi_k-\lambda\hat\eta_k| \,/\, |\lambda\hat\eta_k|$")
    ax.set_title("E. Trust: convergence + bordering identity",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, which="both")

    # ================= F. verdict table ===================================
    ax = ax_tab
    ax.axis("off")
    hdr = ["#", "omega_r", "omega_i", "|lambda|"]
    for k in range(n_dof):
        hdr.append(f"p[{labels[k]}]")
    for k in range(1, n_dof):
        hdr.append(f"phase[{labels[k]}]")
    hdr += ["struct_E", "residual", "verdict"]

    rows, colors = [], []
    for rank, i in enumerate(sel):
        verdict, col = classify(d, i, noise=bool(is_noise[i]))
        r = [f"#{rank+1}", f"{freq[i]:+.4f}", f"{grow[i]:+.4f}", f"{abs(lam[i]):.4f}"]
        for k in range(n_dof):
            r.append(f"{d['frac'][k][i]:.3f}")
        for k in range(1, n_dof):
            r.append(f"{d['phase'][k][i]:+.1f}°")
        r += [f"{sre[i]:.3f}", f"{max(resid[i],err[i]):.1e}", verdict]
        rows.append(r)
        colors.append(col)

    tab = ax.table(cellText=rows, colLabels=hdr, cellLoc="center", loc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1, 1.55)
    for j in range(len(hdr)):
        tab[0, j].set_facecolor("#e8e8e8")
        tab[0, j].set_text_props(fontweight="bold")
    for r_i, col in enumerate(colors):
        cell = tab[r_i + 1, len(hdr) - 1]
        cell.set_text_props(color=col, fontweight="bold")
    ax.set_title(f"F. Verdict — {len(sel)} {sel_desc} mode(s) "
                 f"(of {int(mask.sum())} shown / {d['n_modes']} in file)",
                 fontsize=11, fontweight="bold", loc="left", pad=16)

    ttl = title or "FELiCS rigid-body FSI — modal post-processing"
    fig.suptitle(f"{ttl}    [{n_dof} DOF: {', '.join(labels)}]",
                 fontsize=13.5, fontweight="bold", y=0.985)
    return fig


# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Visualise the FELiCS FSI structural-modes CSV.")
    p.add_argument("csv", help="path to fsi_structural_modes.csv")
    p.add_argument("--n-modes", type=int, default=3,
                   help="how many modes in the per-mode panels (default: 3)")
    p.add_argument("--error-tol", type=float, default=None,
                   help="drop modes whose SLEPc residual exceeds this")
    p.add_argument("--hide-drift", action="store_true",
                   help="hide the |lambda|~0 rigid-body drift modes (K=0)")
    p.add_argument("--sort-by", choices=["growth", "struct", "freq"],
                   default="growth",
                   help="ranking for the per-mode panels (default: growth)")
    p.add_argument("-o", "--output", default=None,
                   help="save the figure instead of showing it")
    p.add_argument("--title", default=None, help="figure title")
    args = p.parse_args()

    if args.output:
        matplotlib.use("Agg")

    d = parse_csv(args.csv)
    print(f"Read {d['n_modes']} modes, {d['n_dof']} DOF: {', '.join(d['labels'])}")

    resid, is_noise = report_health(d)

    fig = make_figure(d, n_show=args.n_modes, error_tol=args.error_tol,
                      hide_drift=args.hide_drift, sort_by=args.sort_by,
                      title=args.title)
    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Figure written to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()