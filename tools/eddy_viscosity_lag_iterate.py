#!/usr/bin/env python3
"""Iterate a modal FSI run with Case.EddyViscosityLag until the eigenvalue in the lag factor
equals the eigenvalue found.

usage
    python tools/eddy_viscosity_lag_iterate.py case.json --start 0.70 -0.03 [--tolerance 2e-4]
        [--max-runs 4] [--felics path/to/src/main.py] [--python python]

Each run n writes <case>_lag<n>.json with Eigenvalue set to the current value and ExportFolder
<export>_lag<n>, runs FELiCS on it, and reads the structural mode with the largest structural
energy fraction from <export>_lag<n>/fsi_structural_modes.csv.  Stops when the eigenvalue changes
by less than --tolerance.  The eigenvalue convention is FELiCS's (omega + i sigma, sigma > 0 growing).
"""
# Standard library
import argparse
import csv
import json
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__.split("\n\nusage")[0])
parser.add_argument("case")
parser.add_argument("--start", nargs=2, type=float, required=True, metavar=("REAL", "IMAG"))
parser.add_argument("--tolerance", type=float, default=2e-4)
parser.add_argument("--max-runs", type=int, default=4)
parser.add_argument("--felics", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "main.py"))
parser.add_argument("--python", default=sys.executable)
arguments = parser.parse_args()

with open(arguments.case) as handle:
    case = json.load(handle)
lag = case["Case"].get("EddyViscosityLag", {})
if lag.get("Model", "None") == "None":
    sys.exit("Case.EddyViscosityLag.Model is 'None': nothing to iterate.")
stem = os.path.splitext(arguments.case)[0]
export = case["Export"]["ExportFolder"]
eigenvalue = complex(*arguments.start)
history = []
for run in range(1, arguments.max_runs + 1):
    case["Case"]["EddyViscosityLag"]["Eigenvalue"] = [eigenvalue.real, eigenvalue.imag]
    case["Export"]["ExportFolder"] = f"{export}_lag{run}"
    case_file = f"{stem}_lag{run}.json"
    with open(case_file, "w") as handle:
        json.dump(case, handle, indent=4)
    log_file = f"{stem}_lag{run}.log"
    with open(log_file, "w") as handle:
        result = subprocess.run([arguments.python, arguments.felics, "-f", case_file, "-p"], stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"run {run} failed (see {log_file})")
    with open(os.path.join(case["Export"]["ExportFolder"], "fsi_structural_modes.csv")) as handle:
        rows = list(csv.DictReader(handle))
    best = max(rows, key=lambda row: float(row["struct_ratio_energy"]))
    found = complex(float(best["omega_r"]), float(best["omega_i"]))
    history.append((eigenvalue, found))
    print(f"run {run}: factor evaluated at {eigenvalue:.5f}, eigenvalue found {found:.5f}", flush=True)
    if abs(found - eigenvalue) < arguments.tolerance:
        print("converged")
        break
    eigenvalue = found
else:
    print(f"not converged after {arguments.max_runs} runs (last change {abs(history[-1][1] - history[-1][0]):.2e})")
