# Eddy-viscosity lag (delayed turbulence response)

With `TurbulenceModel: "File"` the eddy viscosity of the mean flow is frozen and the coherent Reynolds stress follows the perturbation strain instantly, $\tilde{\boldsymbol\tau} = -2\overline{\mu_t}\,\mathbf{S}'$. The option `Case.EddyViscosityLag` lets the stress respond with a delay, with one time scale $T$. The default, `{"Model": "None"}`, changes nothing, and `TimeScale: 0` recovers the frozen model exactly.

Assumptions:
- the frozen eddy viscosity $\overline{\mu_t}$ of the file is kept; only its action on the perturbation is delayed;
- the delay is linear and has a single time scale (uniform, or restricted to a region);
- the mean flow, the base-flow residual and the moving base state of the rigid-body FSI are not changed;
- 2-D Cartesian cases only.

Nomenclature:
- $\lambda = \omega + \mathrm{i}\sigma$: eigenvalue, perturbations $\sim e^{-\mathrm{i}\lambda t}$
- $T$: turbulence response time (`TimeScale`)
- $\overline{\mu_t}$: eddy viscosity of the mean-flow file; $\mathbf{S}'$: perturbation strain rate

## Models

`"Relaxation"`: the stress relaxes to its frozen value with time constant $T$,

$$
T\,\frac{\partial\tilde{\boldsymbol\tau}}{\partial t} + \tilde{\boldsymbol\tau} = -2\overline{\mu_t}\,\mathbf{S}'
\quad\Rightarrow\quad
\hat{\boldsymbol\tau} = -2\,\frac{\overline{\mu_t}}{1 - \mathrm{i}\lambda T}\,\mathbf{S}(\hat{\mathbf{u}}) .
$$

The stress lags the strain by $\arctan(\omega T)$ (at most $90^\circ$) with its amplitude reduced by the cosine of that angle.

`"Delay"`: the stress is the frozen stress evaluated a time $T$ earlier,

$$
\tilde{\boldsymbol\tau}(t) = -2\overline{\mu_t}\,\mathbf{S}'(t - T)
\quad\Rightarrow\quad
\hat{\boldsymbol\tau} = -2\,\overline{\mu_t}\,e^{\mathrm{i}\lambda T}\,\mathbf{S}(\hat{\mathbf{u}}) ,
$$

which lags by $\omega T$ without attenuation.

In both cases the eddy viscosity seen by the perturbation is $\overline{\mu_t}$ times a complex factor (`EddyViscosityLag.factor`). It replaces $\overline{\mu_t}$ in the perturbation stress of the momentum equation (`FieldProperties.tau`, through `mu_tot_perturbation`) and in the perturbation wall traction of the rigid-body force (`RigidBodyMotionFSI._build_force_vector`). Nothing else uses it.

## Regions

`Region: "All"` applies $T$ everywhere. `"Wake"` applies it for $x > $ `TrailingEdgeX`, `"Band"` for `XRange[0]` $< x <$ `XRange[1]`; both with cubic smooth ramps of half width `RampWidth`. Outside the region the eddy viscosity stays frozen.

## Solving the eigenproblem

The factor contains $\lambda$, so the eigenproblem is nonlinear in the eigenvalue. FELiCS evaluates the factor with the value given in `Eigenvalue` (`[real, imaginary]`, FELiCS convention). The run must be repeated with the eigenvalue it finds until the two agree; two to four runs are usually enough. The script `tools/eddy_viscosity_lag_iterate.py` does this loop:

```
python tools/eddy_viscosity_lag_iterate.py case.json --start 0.70 -0.03 --tolerance 2e-4 --max-runs 4
```

It writes `case_lag1.json`, `case_lag2.json`, … with the updated `Eigenvalue` and export folder, runs them one after the other and prints the sequence of eigenvalues (the structural mode with the largest structural energy fraction, from `fsi_structural_modes.csv`).

## Choosing the time scale

$T$ is a physical input, not a tuning parameter: it should be measured. For a separation bubble, project the DNS wall shear onto the structural modes and read the phase lag $\varphi$ of the reattachment point behind the motion; then $T = \tan\varphi/\omega$ (relaxation) or $T = \varphi/\omega$ (delay). The laminar separation point should respond with $T \approx 0$, which is a check of the procedure.

## Settings

```json
"EddyViscosityLag": {
    "Model": "Relaxation",
    "TimeScale": 2.2,
    "Eigenvalue": [0.70, -0.03],
    "Region": "Band",
    "XRange": [0.72, 1.01],
    "RampWidth": 0.02
}
```

| key | meaning |
|---|---|
| `Model` | `"None"` (default), `"Relaxation"`, `"Delay"` |
| `TimeScale` | $T$ in the units of the mean flow (0 = frozen) |
| `Eigenvalue` | `[real, imaginary]` used in the factor |
| `Region` | `"All"`, `"Wake"`, `"Band"` |
| `TrailingEdgeX` | (`Wake`) start of the lagged region |
| `XRange` | (`Band`) `[x_start, x_end]` of the lagged region |
| `RampWidth` | (`Wake`, `Band`) half width of the smooth ramps |

No defaults exist for the model keys: a `Model` other than `"None"` needs exactly the keys of its region.
