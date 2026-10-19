# PlasmaWaveLab — Wolfram Language (Mathematica) edition

A Wolfram Language implementation of the same models and algorithms as the
Python package in `src/plasmawavelab/`, for users who work in Mathematica.

| File | Contents |
|---|---|
| `PlasmaWaveLab.wl` | Package: Sagdeev model, critical Mach number, soliton amplitude and profile, small-amplitude KdV limit, IF-RK4 Fourier KdV solver, invariants, exact two-soliton phase shifts |
| `Verification.wls` | Reproduces the manuscript's verification results and checks them against the Python reference values (PASS/FAIL) |

Requires Mathematica / Wolfram Engine 12.0 or newer (no extra packages).

## Quick start

```mathematica
Get["path/to/PlasmaWaveLab/mathematica/PlasmaWaveLab.wl"]

CriticalMach[]                       (* 1.5852010652... *)
PeakPotential[1.2]                   (* 0.5243854849... *)
prof = SolitaryProfile[1.2];
ListLinePlot[Transpose[{prof["Xi"], prof["Phi"]}]]

(* KdV soliton propagation *)
x = Most @ Subdivide[-60., 60., 512];
{t, u} = KdVPropagate[KdVSoliton[x, 0.], "Length" -> 120., "FinalTime" -> 2.];
KdVInvariants[Last[u], 120./512]

(* small-amplitude limit: Sagdeev soliton vs ion-acoustic KdV soliton *)
IonAcousticKdVCoefficients[]         (* <|"Alpha" -> 1, "Beta" -> 1/2|> *)
SmallAmplitudeSoliton[0., 1.01]      (* 3 (M-1) = 0.03 *)
```

## Verification

From a terminal:

```bash
wolframscript -file mathematica/Verification.wls
```

or open `Verification.wls` in Mathematica and choose *Evaluation → Evaluate
Notebook*. Each check prints `PASS`/`FAIL`. Figures, a CSV table and
`summary.json` are written to `mathematica/output/`.

## Correspondence with the Python API

| Python (`plasmawavelab`) | Wolfram Language |
|---|---|
| `sagdeev(phi, mach)` | `SagdeevPotential[phi, M]` |
| `sagdeev_prime(phi, mach)` | `SagdeevPrime[phi, M]` |
| `critical_mach()` | `CriticalMach[]` |
| `peak_potential(mach)` | `PeakPotential[M]` |
| `solitary_profile(mach, xi_max, points)` | `SolitaryProfile[M, "XiMax" -> …, "Points" -> …]` |
| `small_amplitude_soliton(xi, mach)` | `SmallAmplitudeSoliton[xi, M]` |
| `ion_acoustic_kdv_coefficients()` | `IonAcousticKdVCoefficients[]` |
| `kdv_soliton(x, t, amplitude, alpha, beta, center)` | `KdVSoliton[x, t, "Amplitude" -> …, "Alpha" -> …, "Beta" -> …, "Center" -> …]` |
| `kdv_propagate(u0, length, final_time, dt, alpha, beta, save_every)` | `KdVPropagate[u0, "Length" -> …, "FinalTime" -> …, "TimeStep" -> …, "Alpha" -> …, "Beta" -> …, "SaveEvery" -> …]` |
| `kdv_invariants(u, dx, alpha, beta)` | `KdVInvariants[u, dx, "Alpha" -> …, "Beta" -> …]` |
| — | `KdVTwoSolitonPhaseShifts[A1, A2]` |

The FFT convention (`FourierParameters -> {1, -1}`), wavenumber ordering,
two-thirds mask and IF-RK4 stages are identical to the Python code, so results
agree to within round-off and solver-tolerance differences.
