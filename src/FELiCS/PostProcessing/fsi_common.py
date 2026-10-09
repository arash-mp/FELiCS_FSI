#  ___________________________________   _______________________________________________________
# /-----------------------------------\ /-------------------------------------------------------\
# |   (         (                (     |  This source code is part of FELiCS                     |
# |   )\ )     ) )        (      )\ )  |  (F)inite (E)lement (Li)nearized (C)ombustion (S)olver  |
# |  (()/(  (  (()/( (    )\   (()/(   |                                                         |
# |  /(_)) )\  /(_)))\  (((_)  /(_))   |  Licensed under the GNU GPLv3                           |
# |  (_)_)((_) (_)) ((_) )\___ (_))    |                                                         |
# |  | __|| __|| |   (_)((/ __|/ __|   |  (C) 2018-2025: The FELiCS Developers (www.felics.eu)   |
# |  | _| | _| | |__ | | | (__ \__ \   |  Visit          www.felics.eu                           |
# |  |_|  |___||____||_|  \___||___/   |  Contact        info@felics.eu                          |
# \___________________________________/ \_______________________________________________________/
#
"""
fsi_common.py -- shared helpers for the rigid-body FSI post-processing scripts.

This is the SINGLE place where the layout of ``fsi_structural_modes.csv`` is
understood.  Both ``plot_fsi_modes.py`` and ``plot_fsi_growth.py`` import from
here, so if the CSV format ever changes, only this file needs editing.

Deliberately depends on NOTHING from FELiCS itself (no dolfinx / PETSc / SLEPc):
these are pure post-processing utilities that only need numpy / pandas, so they
run on a laptop against a CSV copied off the HPC.
"""

import re
import sys

import numpy as np


# ----------------------------------------------------------------------
# CSV layout
# ----------------------------------------------------------------------
def parse_csv(path):
    """Read fsi_structural_modes.csv and discover the DOFs from the header.

    Column convention written by run_modal.py, per DOF k with label `tag`:
        eta{k}_{tag}_r/_i, phi{k}_{tag}_r/_i,
        eta{k}_{tag}_abs, dof{k}_{tag}_frac, dof{k}_{tag}_phase_deg
    plus the globals: omega_r, omega_i, dominant_dof, struct_ratio,
    struct_ratio_energy, error.

    Labels never contain a comma (coordinates use ';'), so the header is safe to
    split.  Older files without the abs/frac/phase columns are handled by deriving
    an approximate split from the amplitudes.

    Returns a dict with keys:
        omega, eta, phi, frac, phase, struct_ratio, struct_energy, error,
        dominant, labels, n_dof, n_modes
    """
    import pandas as pd

    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]

    # discover DOF index + label from the eta{k}_{label}_r columns
    dofs = {}
    for col in df.columns:
        m = re.match(r"^eta(\d+)_(.+)_r$", col)
        if m:
            dofs[int(m.group(1))] = m.group(2)
    if not dofs:                                   # very old, untagged file
        for col in df.columns:
            m = re.match(r"^eta(\d+)_r$", col)
            if m:
                dofs[int(m.group(1))] = f"DOF{m.group(1)}"
    if not dofs:
        sys.exit(f"ERROR: no eta columns found in {path}. Is this an FSI CSV?")

    n_dof  = max(dofs) + 1
    labels = [dofs.get(k, f"DOF{k}") for k in range(n_dof)]

    d = {}
    d["omega"] = df["omega_r"].values + 1j * df["omega_i"].values
    d["eta"], d["phi"], d["frac"], d["phase"] = [], [], [], []

    for k in range(n_dof):
        lab = labels[k]
        e_r, e_i = f"eta{k}_{lab}_r", f"eta{k}_{lab}_i"
        p_r, p_i = f"phi{k}_{lab}_r", f"phi{k}_{lab}_i"
        if e_r not in df.columns:                  # untagged fallback
            e_r, e_i, p_r, p_i = (f"eta{k}_r", f"eta{k}_i",
                                  f"phi{k}_r", f"phi{k}_i")
        d["eta"].append(df[e_r].values + 1j * df[e_i].values)
        d["phi"].append(df[p_r].values + 1j * df[p_i].values)

        fcol, pcol = f"dof{k}_{lab}_frac", f"dof{k}_{lab}_phase_deg"
        d["frac"].append(df[fcol].values if fcol in df.columns else None)
        d["phase"].append(df[pcol].values if pcol in df.columns else None)

    # graceful degradation for CSVs written before the frac/phase columns existed
    if any(f is None for f in d["frac"]):
        print("NOTE: this CSV predates the participation columns; deriving an "
              "approximate split from |eta| (no M/K weighting).")
        mags = np.array([np.abs(e) for e in d["eta"]])
        tot  = mags.sum(axis=0)
        tot[tot == 0] = 1.0
        d["frac"] = [mags[k] / tot for k in range(n_dof)]
    if any(p is None for p in d["phase"]):
        ref = np.angle(d["eta"][0])
        d["phase"] = [
            np.degrees((np.angle(d["eta"][k]) - ref + np.pi) % (2 * np.pi) - np.pi)
            for k in range(n_dof)
        ]

    nan = np.full(len(df), np.nan)
    d["struct_ratio"]  = df["struct_ratio"].values \
        if "struct_ratio" in df.columns else nan
    d["struct_energy"] = df["struct_ratio_energy"].values \
        if "struct_ratio_energy" in df.columns else nan
    d["error"]         = df["error"].values \
        if "error" in df.columns else np.zeros(len(df))
    d["dominant"] = (df["dominant_dof"].values.astype(int)
                     if "dominant_dof" in df.columns
                     else np.argmax(np.array(d["frac"]), axis=0))

    d["labels"]  = labels
    d["n_dof"]   = n_dof
    d["n_modes"] = len(df)
    return d


def is_rotation(label):
    """True if a DOF label denotes a rotation (used for rad/deg handling)."""
    return str(label).lower().startswith("rot")


# ----------------------------------------------------------------------
# Physics helpers
# ----------------------------------------------------------------------
def consistency_residual(d, noise_rel=1e-3):
    """Check the bordering identity  1j*phi_k = lambda*eta_k.

    It must hold for every converged mode.  Reported as a RELATIVE residual, but
    ONLY for modes whose structural amplitude is above the noise floor: modes in
    which the body barely moves have |eta| ~ 1e-10 while the eigenvector is
    normalised to 1, so dividing by that amplitude turns pure round-off into a
    huge apparent "relative error".  For those, the residual is scaled by the
    largest structural amplitude in the file instead, which is the honest measure.

    Returns (residual, is_noise).  Where is_noise is True, the participation
    fractions and phases are MEANINGLESS and must not be interpreted.
    """
    lam = d["omega"]
    amp = np.zeros(d["n_modes"])
    for k in range(d["n_dof"]):
        amp = np.maximum(amp, np.abs(d["eta"][k]))
    amp_max  = amp.max() if amp.max() > 0 else 1.0
    is_noise = amp < noise_rel * amp_max

    worst = np.zeros(d["n_modes"])
    for k in range(d["n_dof"]):
        num = np.abs(1j * d["phi"][k] - lam * d["eta"][k])
        den = np.where(is_noise, np.abs(lam) * amp_max,
                       np.abs(lam * d["eta"][k]))
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(den > 0, num / den, 0.0)
        worst = np.maximum(worst, np.nan_to_num(r))
    return worst, is_noise


def classify(d, i, noise=False, drift_tol=1e-8, fluid_cut=0.05,
             phase_band=35.0, hybrid_rel=0.6):
    """Return (verdict, colour) for mode i.

    Applies the reading rules of the integration notes:
      * lambda ~ 0                      -> rigid-body drift mode (K = 0)
      * structural amplitude at noise    -> fluid mode, body essentially still
      * struct_ratio_energy < fluid_cut  -> fluid / wake mode
      * every DOF carries a decent share -> hybrid; + growth + phase ~ +/-90 deg
                                            -> FLUTTER candidate
    NOTE these are heuristics for triage, not a verdict: confirm with a sweep.
    """
    lam  = d["omega"][i]
    grow = lam.imag
    sre  = d["struct_energy"][i]
    fr   = [d["frac"][k][i] for k in range(d["n_dof"])]
    dom  = int(np.argmax(fr))
    lo   = min(fr) if d["n_dof"] > 1 else 1.0

    if abs(lam) < drift_tol:
        return "drift mode (K=0)", "grey"

    if noise:
        return ("fluid mode (body ~ still) | UNSTABLE" if grow > 0
                else "fluid mode (body ~ still)"), "tab:blue"

    equi   = 1.0 / d["n_dof"]
    hybrid = (d["n_dof"] > 1) and (lo >= hybrid_rel * equi)
    near90 = any(abs(abs(d["phase"][k][i]) - 90.0) <= phase_band
                 for k in range(1, d["n_dof"]))

    if not np.isnan(sre) and sre < fluid_cut:
        return ("fluid / wake mode | UNSTABLE" if grow > 0
                else "fluid / wake mode | stable"), "tab:blue"

    if hybrid and near90 and grow > 0:
        return "FLUTTER candidate", "red"

    if hybrid:
        base, col = "hybrid (coupled)", "tab:purple"
    else:
        base = f"{d['labels'][dom]}-dominated"
        col  = "tab:green" if dom == 0 else "tab:orange"

    base += " | UNSTABLE" if grow > 0 else " | stable"
    return base, col


def report_health(d):
    """Print the checks that must pass before ANY mode is interpreted.
    Returns (residual, is_noise) for reuse by the caller."""
    resid, is_noise = consistency_residual(d)

    bad_r = int(np.sum(resid > 1e-6))
    bad_e = int(np.sum(np.nan_to_num(d["error"], nan=0.0) > 1e-6))
    if bad_r:
        print(f"WARNING: {bad_r} mode(s) violate 1j*phi = lambda*eta by >1e-6 "
              f"-- the bordering may be mis-wired, or those modes are unconverged.")
    if bad_e:
        print(f"WARNING: {bad_e} mode(s) have SLEPc residual > 1e-6 "
              f"-- do not interpret them.")

    n_unstable = int(np.sum(d["omega"].imag > 0))
    print(f"{n_unstable} unstable mode(s) (omega_i > 0).")

    # the diagnostic that matters most: is there a structural branch at all?
    sre  = np.nan_to_num(d["struct_energy"], nan=0.0)
    absl = np.abs(d["omega"])
    if sre.max() < 1e-3:
        print("\n" + "!" * 74)
        print("NO STRUCTURAL BRANCH IN THIS SPECTRUM.")
        print(f"  Every mode has struct_ratio_energy < {sre.max():.1e} -- i.e. the "
              f"body is")
        print("  essentially STILL in all of them. These are all fluid/wake modes.")
        print(f"  |lambda| spans {absl.min():.3f} .. {absl.max():.3f}; the solver "
              f"never looked near the origin.")
        print("  With K = 0 (free body) the structural branches sit at SMALL "
              "|lambda|.")
        print("  -> Re-run with a small-magnitude EigenValueGuess, e.g. [0.02], "
              "[0.05], [0.1],")
        print("     in ADDITION to your shedding-frequency guess.")
        print("!" * 74 + "\n")

    n_noise = int(is_noise.sum())
    if n_noise:
        print(f"NOTE: {n_noise} mode(s) have structural amplitudes at the noise "
              f"floor; their participation/phase columns are meaningless.")

    return resid, is_noise
