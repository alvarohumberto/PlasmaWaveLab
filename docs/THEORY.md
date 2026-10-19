# Physical models and numerical methods

## Sagdeev model

Cold, positive inertial ions with Boltzmann electrons in a 1D planar,
stationary, collisionless, unmagnetized, normalized electrostatic model.
With travelling coordinate `xi`, Mach number `M > 1` and normalized
potential `phi`:

`(phi')^2 / 2 + S(phi;M) = 0`, `phi'' = -dS/dphi`.

`S(phi;M) = 1 - exp(phi) + M^2*[1 - sqrt(1 - 2*phi/M^2)]`.

The physical real-valued domain is `phi < M^2/2`. A positive-polarity
solitary wave requires a nonzero root `A` in this domain, together with
`S(phi) < 0` for `0 < phi < A`. The implementation numerically brackets
the root and constructs a symmetric profile by integrating from the
crest. Numerical tails are truncated when the potential drops to a small
fraction of `A`, because the far-field equilibrium is a saddle in the
spatial dynamical system.

Do not assume a positive-polarity solitary wave exists for every `M>1`.

## KdV model

Periodic `u_t + alpha*u*u_x + beta*u_xxx = 0` with positive coefficients.
Reference infinite-domain soliton:

`u(x,t) = A sech^2(sqrt(alpha*A/(12*beta))*(x - x0 - alpha*A*t/3))`.

Propagation uses periodic FFT spatial differentiation, an integrating-factor
four-stage Runge–Kutta method, and 2/3 filtering for the nonlinear term.
The numerical method is not claimed to conserve all invariants exactly.

## Link between the two models

For M = 1 + δ with 0 < δ ≪ 1, expanding S gives S ≈ −δφ² + φ³/3, so the
Sagdeev soliton tends to `phi = 3*delta*sech^2(sqrt(delta/2)*xi)`. This is the
soliton of the KdV equation with `alpha = 1`, `beta = 1/2` (Washimi & Taniuti
1966) moving at speed δ in the frame of the ion-acoustic speed. The package
exposes this mapping as `ion_acoustic_kdv_coefficients()` and
`small_amplitude_soliton(xi, mach)`; the difference from the full Sagdeev
soliton is O(δ) and is quantified by `examples/verification.py`.

## Critical Mach number

The soliton branch ends when the turning point reaches the ion-density limit
φ = M²/2, i.e. when 1 − exp(M²/2) + M² = 0, giving M_c ≈ 1.5852
(`critical_mach()`). For M ≥ M_c no solitary wave exists.

## Units and applicability

All inputs and outputs are dimensionless. Real laboratory quantities and
non-Maxwellian/kinetic electrons, finite ion temperature, magnetic fields,
collisions and dust species are not modeled.
