#!/usr/bin/env python3
"""
plot_fsi_growth.py -- reconstruct the TIME EVOLUTION of the rigid-body DOFs from a
single FSI eigenmode, so it can be overlaid on a DNS signal.

Physics
-------
FELiCS solves  A q = lambda B q  with the convention  q(t) ~ exp(-i*lambda*t).
Writing  lambda = omega_r + i*sigma :

    exp(-i*lambda*t) = exp(sigma*t) * exp(-i*omega_r*t)

so every rigid-body DOF of a single mode evolves as

    eta_k(t) = |eta_k| * exp(sigma*t) * cos(omega_r*t - arg(eta_k))
    phi_k(t) = d eta_k / dt          (identically, since 1j*phi_k = lambda*eta_k)

KEY PREDICTION FOR DNS:
    Every DOF of a given mode shares the SAME growth rate sigma and the SAME
    frequency omega_r.  They differ ONLY in amplitude and phase, both fixed by the
    eigenvector.  In the linear regime the DNS must reproduce:

      * two parallel straight lines for log|envelope| vs t, slope = sigma
      * a CONSTANT amplitude ratio  |eta_1| / |eta_0|
      * a CONSTANT phase lag        arg(eta_1) - arg(eta_0)

    Any drift in the ratio or the lag means you are NOT yet in the asymptotic
    linear regime (or nonlinearity has kicked in).

IMPORTANT -- amplitude is arbitrary
-----------------------------------
An eigenvector has no intrinsic scale: the eigensolver normalises it however it
likes (here |eta| ~ 1e-6).  To compare with DNS you must PICK a scale.  Use
--amp0 to set the initial amplitude of the reference DOF; everything else follows
from the eigenvector.  Match --amp0 to your DNS initial perturbation.

Usage
-----
    # leading (most unstable) mode, default scaling, out to t = 100
    python plot_fsi_growth.py fsi_structural_modes.csv --t-end 100

    # a specific mode by row index, with a chosen initial plunge amplitude
    python plot_fsi_growth.py fsi_structural_modes.csv --mode 0 --amp0 1e-3 --t-end 200

    # pick the mode nearest a given eigenvalue
    python plot_fsi_growth.py fsi_structural_modes.csv --near 0.81+0.025j

    # export the predicted signal so you can overlay DNS yourself
    python plot_fsi_growth.py fsi_structural_modes.csv --save-csv lsa_prediction.csv

    # overlay a DNS time series (CSV with columns: t, eta0, eta1, ...)
    python plot_fsi_growth.py fsi_structural_modes.csv --dns dns_motion.csv
"""

import argparse
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# --- locate fsi_common whether we are run as a script, via a symlink on PATH,
# --- or as "python -m FELiCS.PostProcessing.plot_fsi_..."
try:
    from .fsi_common import (parse_csv, is_rotation, consistency_residual,
                             classify, report_health)
except ImportError:                     # run directly / via symlink
    import os, sys
    _here = os.path.dirname(os.path.realpath(__file__))   # realpath: resolve symlink
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from fsi_common import (parse_csv, is_rotation, consistency_residual,
                            classify, report_health)


# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
def reconstruct(d, i, t, scale):
    """eta_k(t) and phi_k(t) for mode i, scaled by `scale`."""
    lam = d["omega"][i]
    out_eta, out_phi = [], []
    for k in range(d["n_dof"]):
        eh = d["eta"][k][i] * scale
        ph = d["phi"][k][i] * scale
        osc = np.exp(-1j * lam * t)          # = exp(sigma t) * exp(-i omega_r t)
        out_eta.append(np.real(eh * osc))
        out_phi.append(np.real(ph * osc))
    return out_eta, out_phi


def main():
    p = argparse.ArgumentParser(
        description="Reconstruct DOF growth in time from an FSI eigenmode "
                    "(for comparison with DNS).")
    p.add_argument("csv", help="fsi_structural_modes.csv")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--mode", type=int, default=None,
                   help="row index of the mode to use (default: most unstable)")
    g.add_argument("--near", type=str, default=None,
                   help="use the mode nearest this eigenvalue, e.g. '0.81+0.025j'")
    p.add_argument("--t-end", type=float, default=100.0,
                   help="final time (default: 100)")
    p.add_argument("--t-start", type=float, default=0.0, help="initial time")
    p.add_argument("--npts", type=int, default=4000, help="samples (default: 4000)")
    p.add_argument("--amp0", type=float, default=1e-3,
                   help="amplitude |eta| of the reference DOF at t=t_start "
                        "(default: 1e-3). MATCH THIS TO YOUR DNS.")
    p.add_argument("--ref-dof", type=int, default=0,
                   help="which DOF --amp0 refers to (default: 0)")
    p.add_argument("--degrees", action="store_true",
                   help="plot rotation DOFs in degrees instead of radians")
    p.add_argument("--linear-limit", type=float, default=None,
                   help="mark where |eta| of any DOF exceeds this (linear-regime "
                        "validity bound; for rotation, given in the plotted unit)")
    p.add_argument("--dns", default=None,
                   help="CSV with DNS signal: columns t, eta0, eta1, ... to overlay")
    p.add_argument("--save-csv", default=None,
                   help="write the predicted time series to this CSV")
    p.add_argument("-o", "--output", default=None, help="save figure to file")
    args = p.parse_args()

    if args.output:
        matplotlib.use("Agg")

    d = parse_csv(args.csv)
    n_dof, labels = d["n_dof"], d["labels"]

    # ---- pick the mode ---------------------------------------------------
    if args.mode is not None:
        i = args.mode
        if not (0 <= i < d["n_modes"]):
            sys.exit(f"--mode {i} out of range (file has {d['n_modes']} modes).")
        how = f"row {i}"
    elif args.near is not None:
        target = complex(args.near.replace("i", "j"))
        i = int(np.argmin(np.abs(d["omega"] - target)))
        how = f"nearest to {target}"
    else:
        i = int(np.argmax(d["omega"].imag))
        how = "most unstable"

    lam   = d["omega"][i]
    sigma = lam.imag
    omr   = lam.real

    # ---- scale: eigenvectors have no intrinsic amplitude -----------------
    ref = args.ref_dof
    a_ref = abs(d["eta"][ref][i])
    if a_ref == 0:
        sys.exit(f"Reference DOF {ref} has zero amplitude in this mode; "
                 f"pick another with --ref-dof.")
    # scale so that |eta_ref| = amp0 AT t = t_start
    scale = args.amp0 / (a_ref * np.exp(sigma * args.t_start))

    t = np.linspace(args.t_start, args.t_end, args.npts)
    eta, phi = reconstruct(d, i, t, scale)

    # unit conversion for rotation DOFs
    unit, conv = [], []
    for k in range(n_dof):
        if is_rotation(labels[k]) and args.degrees:
            unit.append("deg"); conv.append(180.0 / np.pi)
        elif is_rotation(labels[k]):
            unit.append("rad"); conv.append(1.0)
        else:
            unit.append("L");   conv.append(1.0)
    eta_p = [eta[k] * conv[k] for k in range(n_dof)]

    # ---- report ----------------------------------------------------------
    print("=" * 72)
    print(f"Mode {i} ({how}):  lambda = {omr:+.6f} {sigma:+.6f}j")
    print("-" * 72)
    print(f"  growth rate   sigma   = {sigma:+.6f}")
    print(f"  frequency     omega_r = {omr:+.6f}")
    if abs(omr) > 1e-12:
        T = 2 * np.pi / abs(omr)
        print(f"  period        T       = {T:.4f}   (Strouhal = {omr/(2*np.pi):.4f})")
    else:
        T = np.nan
        print("  period        T       = inf   (NON-OSCILLATORY: divergence-type mode)")
    if abs(sigma) > 1e-12:
        word = "doubling" if sigma > 0 else "halving "
        print(f"  e-folding time  1/|sigma|    = {1/abs(sigma):.3f}"
              f"  ({'GROWTH' if sigma>0 else 'decay'})")
        print(f"  {word} time   ln2/|sigma|  = {np.log(2)/abs(sigma):.3f}")
        if abs(omr) > 1e-12:
            print(f"  growth per cycle exp(sigma*T) = {np.exp(sigma*T):.4f}")
    else:
        print("  neutral mode (sigma ~ 0)")
    print("  SLEPc residual = %.2e" % d["error"][i])
    print("-" * 72)
    print(f"  scaling: |eta[{labels[ref]}]| = {args.amp0:g} at t = {args.t_start:g}"
          f"   (eigenvector scale factor {scale:.4g})")
    print("-" * 72)
    a0 = abs(d["eta"][ref][i])
    for k in range(n_dof):
        ak = abs(d["eta"][k][i])
        ph = np.degrees(np.angle(d["eta"][k][i]) - np.angle(d["eta"][ref][i]))
        ph = (ph + 180) % 360 - 180
        lag = ph / 360.0 * T if np.isfinite(T) else np.nan
        amp0 = ak * scale * conv[k]
        line = (f"  {labels[k]:<16} |eta|/|eta_ref| = {ak/a0:9.4e}   "
                f"amp(t0) = {amp0:10.4g} {unit[k]:<3}  "
                f"phase = {ph:+7.1f} deg")
        if k != ref and np.isfinite(lag):
            line += f"  (lags by {lag:+.3f} t)"
        print(line)
    print("=" * 72)
    print("DNS CHECK: all DOFs must show the SAME sigma and SAME omega_r;")
    print("           the amplitude ratio and phase above must stay CONSTANT in time.")
    print("=" * 72)

    # ---- linear-validity horizon ----------------------------------------
    t_lim = None
    if args.linear_limit is not None:
        env = np.zeros_like(t)
        for k in range(n_dof):
            env = np.maximum(env, abs(d["eta"][k][i]) * scale * conv[k]
                             * np.exp(sigma * t))
        over = np.where(env > args.linear_limit)[0]
        if len(over):
            t_lim = t[over[0]]
            print(f"\nLINEAR-REGIME LIMIT: |eta| exceeds {args.linear_limit:g} "
                  f"at t = {t_lim:.2f}. Beyond that the LSA prediction (and the "
                  f"transpiration BC) is no longer valid.\n")

    # ---- DNS overlay -----------------------------------------------------
    dns = None
    if args.dns:
        import pandas as pd
        dd = pd.read_csv(args.dns)
        dd.columns = [c.strip() for c in dd.columns]
        tcol = dd.columns[0]
        dns = {"t": dd[tcol].values, "cols": []}
        for k in range(n_dof):
            for cand in (f"eta{k}", labels[k], dd.columns[k+1] if k+1 < len(dd.columns) else None):
                if cand and cand in dd.columns:
                    dns["cols"].append(dd[cand].values)
                    break
            else:
                dns["cols"].append(None)

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 8.6))
    fig.subplots_adjust(hspace=0.33, wspace=0.27)
    colors = plt.cm.tab10(np.arange(n_dof) % 10)

    # A: time series (twin axis: displacement vs angle have different units)
    ax = axes[0, 0]
    ax.plot(t, eta_p[0], color=colors[0], lw=1.3, label=f"{labels[0]} [{unit[0]}]")
    ax.plot(t,  abs(d["eta"][0][i])*scale*conv[0]*np.exp(sigma*t),
            color=colors[0], ls="--", lw=0.9, alpha=0.6)
    ax.plot(t, -abs(d["eta"][0][i])*scale*conv[0]*np.exp(sigma*t),
            color=colors[0], ls="--", lw=0.9, alpha=0.6)
    ax.set_ylabel(f"{labels[0]}  [{unit[0]}]", color=colors[0])
    ax.tick_params(axis="y", labelcolor=colors[0])
    if n_dof > 1:
        ax2 = ax.twinx()
        for k in range(1, n_dof):
            ax2.plot(t, eta_p[k], color=colors[k], lw=1.3,
                     label=f"{labels[k]} [{unit[k]}]")
        ax2.set_ylabel(" / ".join(f"{labels[k]} [{unit[k]}]" for k in range(1, n_dof)),
                       color=colors[1])
        ax2.tick_params(axis="y", labelcolor=colors[1])
    if dns:
        for k in range(n_dof):
            if dns["cols"][k] is not None:
                tgt = ax if k == 0 else ax2
                tgt.plot(dns["t"], dns["cols"][k] * conv[k], color=colors[k],
                         ls=":", lw=2.0, alpha=0.8)
        ax.plot([], [], "k:", lw=2, label="DNS")
    if t_lim:
        ax.axvline(t_lim, color="red", ls="-.", lw=1.2)
        ax.text(t_lim, ax.get_ylim()[1], " linear limit", color="red",
                fontsize=8, va="top")
    ax.axhline(0, color="k", lw=0.6, alpha=0.4)
    ax.set_xlabel("t")
    ax.set_title("A. Reconstructed motion (dashed = envelope)",
                 fontsize=11, fontweight="bold", loc="left")
    h1, l1 = ax.get_legend_handles_labels()
    if n_dof > 1:
        h2, l2 = ax2.get_legend_handles_labels()
        h1 += h2; l1 += l2
    ax.legend(h1, l1, fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)

    # B: log envelope -> THE growth-rate comparison
    ax = axes[0, 1]
    for k in range(n_dof):
        env = abs(d["eta"][k][i]) * scale * conv[k] * np.exp(sigma * t)
        ax.plot(t, env, color=colors[k], lw=1.6, label=f"{labels[k]} [{unit[k]}]")
    if dns:
        for k in range(n_dof):
            if dns["cols"][k] is not None:
                ax.plot(dns["t"], np.abs(dns["cols"][k] * conv[k]),
                        color=colors[k], ls=":", lw=1.4, alpha=0.7)
        ax.plot([], [], "k:", lw=2, label="DNS |signal|")
    ax.set_yscale("log")
    ax.set_xlabel("t")
    ax.set_ylabel("|envelope|")
    ax.set_title(f"B. Growth check — slope = sigma = {sigma:+.5f}\n"
                 f"(parallel straight lines in DNS = linear regime)",
                 fontsize=11, fontweight="bold", loc="left")
    if t_lim:
        ax.axvline(t_lim, color="red", ls="-.", lw=1.2)
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.25, which="both")

    # C: phase portrait / Lissajous  (rotation vs translation)
    ax = axes[1, 0]
    if n_dof > 1:
        ax.plot(eta_p[0], eta_p[1], color="tab:purple", lw=1.0, alpha=0.85)
        ax.plot(eta_p[0][0], eta_p[1][0], "o", color="green", ms=7, label="t start")
        ax.plot(eta_p[0][-1], eta_p[1][-1], "s", color="red", ms=7, label="t end")
        ax.set_xlabel(f"{labels[0]}  [{unit[0]}]")
        ax.set_ylabel(f"{labels[1]}  [{unit[1]}]")
        phdeg = np.degrees(np.angle(d["eta"][1][i]) - np.angle(d["eta"][0][i]))
        phdeg = (phdeg + 180) % 360 - 180
        ax.set_title(f"C. Orbit — phase {phdeg:+.1f}° "
                     f"({'open ellipse' if abs(abs(phdeg)-90)<45 else 'near-linear'})",
                     fontsize=11, fontweight="bold", loc="left")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.axvline(0, color="k", lw=0.5, alpha=0.4)
    else:
        ax.text(0.5, 0.5, "single DOF:\nno orbit", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")

    # D: summary text
    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        r"$q(t)=\mathrm{Re}\left[\hat q\, e^{-i\lambda t}\right]$,"
        r"   $\lambda=\omega_r+i\sigma$",
        "",
        rf"$\lambda = {omr:+.5f} {sigma:+.5f}\,i$",
        rf"$\sigma$ (growth)  $= {sigma:+.5f}$",
        rf"$\omega_r$ (freq)   $= {omr:+.5f}$",
    ]
    if np.isfinite(T):
        lines += [rf"$T = 2\pi/\omega_r = {T:.3f}$",
                  rf"growth/cycle $= e^{{\sigma T}} = {np.exp(sigma*T):.4f}$"]
    else:
        lines += ["NON-OSCILLATORY (divergence-type)"]
    if abs(sigma) > 1e-12:
        lines += [rf"e-folding $= {1/abs(sigma):.2f}$",
                  rf"doubling  $= {np.log(2)/abs(sigma):.2f}$"]
    lines += ["", r"$\bf{DNS\ must\ reproduce:}$"]
    for k in range(n_dof):
        ak = abs(d["eta"][k][i]) / a0
        ph = np.degrees(np.angle(d["eta"][k][i]) - np.angle(d["eta"][ref][i]))
        ph = (ph + 180) % 360 - 180
        lines.append(rf"  {labels[k]}:  ratio $={ak:.4g}$,  phase $={ph:+.1f}^\circ$")
    lines += ["", "Same $\\sigma$ and $\\omega_r$ for EVERY DOF.",
              "Ratio & phase CONSTANT in time."]
    ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes,
            va="top", ha="left", fontsize=10.5, family="sans-serif")
    ax.set_title("D. What to compare against DNS",
                 fontsize=11, fontweight="bold", loc="left")

    fig.suptitle(f"FSI mode {i} — predicted linear growth "
                 f"({n_dof} DOF: {', '.join(labels)})",
                 fontsize=13, fontweight="bold", y=0.98)

    # ---- optional CSV export --------------------------------------------
    if args.save_csv:
        import pandas as pd
        out = {"t": t}
        for k in range(n_dof):
            out[f"eta_{labels[k]}"] = eta[k]      # raw units (rad for rotation)
            out[f"phi_{labels[k]}"] = phi[k]
        pd.DataFrame(out).to_csv(args.save_csv, index=False)
        print(f"Predicted time series written to {args.save_csv}")

    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Figure written to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
