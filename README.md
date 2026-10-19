# PlasmaWaveLab

[![Python tests](https://github.com/alvarohumberto/PlasmaWaveLab/actions/workflows/tests.yml/badge.svg)](https://github.com/alvarohumberto/PlasmaWaveLab/actions/workflows/tests.yml)

**PlasmaWaveLab v0.3.0** is an open-source, tested Python package for
ion-acoustic solitary waves at two levels of description: the fully nonlinear
Sagdeev pseudopotential of the cold-ion/Boltzmann-electron plasma, and its
weakly nonlinear Korteweg–de Vries (KdV) limit, solved with a Fourier
integrating-factor RK4 scheme. The two are linked through their explicit
small-amplitude mapping, so each module checks the other.

**Author:** Alvaro H. Salas, FIZMAKO Research Group, Universidad Nacional de
Colombia, Sede Manizales. ORCID: https://orcid.org/0000-0001-9343-6062

## Features

- Sagdeev pseudopotential with domain checks, critical Mach number
  `M_c = 1.5852…`, soliton amplitude and symmetric profile.
- Small-amplitude bridge: `ion_acoustic_kdv_coefficients()` (α = 1, β = 1/2)
  and `small_amplitude_soliton()`.
- Periodic KdV solver with Fourier differentiation, 2/3 de-aliasing and
  fourth-order integrating-factor Runge–Kutta; exact solitons; spectral
  mass / L² / Hamiltonian diagnostics.
- `examples/verification.py` reproduces every figure and number in the
  SoftwareX manuscript (see [docs/VALIDATION.md](docs/VALIDATION.md)).
- 21 pytest tests, GitHub Actions CI on Python 3.10–3.13, CLI with CSV/PNG output.

## Quick start (Windows PowerShell)

Install Python 3.10+ and unzip the archive. In the extracted `PlasmaWaveLab`
folder run:

```powershell
py -m pip install -e ".[test]"
py -m pytest -q
py -m plasmawavelab sagdeev --mach 1.2 --output results/sagdeev
py -m plasmawavelab kdv --n 512 --time 2 --output results/kdv
py examples/verification.py --output results/verification
```

Linux/macOS: replace `py` with `python3` and run `python3 -m pip install -e '.[test]'`. Outputs are created under `results/`.

For the entry-point CLI you can also run:

```bash
plasmawavelab sagdeev --mach 1.2
plasmawavelab kdv --time 2
```

## API

```python
import numpy as np
from plasmawavelab import (critical_mach, peak_potential, solitary_profile,
                           kdv_soliton, kdv_propagate, kdv_invariants)

print(critical_mach())          # 1.5852010652...
print(peak_potential(1.2))      # 0.5243854849...
x = np.linspace(-60, 60, 512, endpoint=False)
u0 = kdv_soliton(x, t=0)
times, states = kdv_propagate(u0, length=120, final_time=2)
print(kdv_invariants(states[-1], dx=120/512))
```

## Repository layout

| Location | Contents |
|---|---|
| `src/plasmawavelab/` | Open-source Python package and CLI |
| `tests/` | Physics and numerical regression tests |
| `.github/workflows/` | Automatic Python tests on GitHub |
| `docs/THEORY.md` | Governing assumptions and numerics |
| `docs/VALIDATION.md` | Validation evidence and outstanding work |
| `manuscript/` | SoftwareX manuscript (LaTeX source + PDF) |
| `docs/GITHUB_UPLOAD.md` | GitHub publishing instructions |
| `figures/` | CLI outputs and `verification/` (manuscript figures) |
| `examples/` | `run_both.py` and `verification.py` |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | MIT license |

## Scientific scope and limitations

One-dimensional, unmagnetised, collisionless plasma with cold ions and
Boltzmann electrons; positive-polarity solitons only. Not a PIC, kinetic,
magnetised, multi-species or dusty-plasma code. Sagdeev far-field tails are
truncated at 10⁻⁵ φ_m. The KdV solver is periodic and conserves mass to
round-off; L² norm and Hamiltonian drift is small but non-zero.

See [theory](docs/THEORY.md) and [verification](docs/VALIDATION.md).

## Citation and license

Citation metadata: [`CITATION.cff`](CITATION.cff). MIT-licensed, see [`LICENSE`](LICENSE).
