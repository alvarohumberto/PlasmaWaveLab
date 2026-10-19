"""Verification tests backing the numbers reported in the manuscript."""
import numpy as np
import pytest
from plasmawavelab import (critical_mach, peak_potential, sagdeev,
                           solitary_profile, small_amplitude_soliton,
                           ion_acoustic_kdv_coefficients, kdv_soliton,
                           kdv_propagate, kdv_invariants)
from plasmawavelab.__main__ import main


def test_critical_mach_value():
    mc = critical_mach()
    assert abs(mc - 1.5852010652) < 1e-9
    assert abs(1 - np.exp(mc**2/2) + mc**2) < 1e-12


def test_no_soliton_above_critical_mach():
    with pytest.raises(ValueError):
        peak_potential(critical_mach() + 1e-3)


def test_amplitude_monotone_below_critical_mach():
    machs = np.linspace(1.02, 1.58, 15)
    amps = np.array([peak_potential(m) for m in machs])
    assert np.all(np.diff(amps) > 0)
    assert np.all(amps < machs**2/2)


@pytest.mark.parametrize('mach', [1.01, 1.02])
def test_small_amplitude_kdv_limit(mach):
    delta = mach - 1
    amp = peak_potential(mach)
    assert abs(amp/(3*delta) - 1) < 1.0*delta  # O(delta) correction
    x, p, _, _ = solitary_profile(mach, xi_max=150, points=3001)
    rel = np.max(np.abs(p - small_amplitude_soliton(x, mach)))/(3*delta)
    assert rel < 1.5*delta


def test_kdv_coefficients_reproduce_sagdeev_limit():
    alpha, beta = ion_acoustic_kdv_coefficients()
    delta = 0.01
    xi = np.linspace(-50, 50, 101)
    ref = 3*delta/np.cosh(np.sqrt(delta/2)*xi)**2
    assert np.allclose(kdv_soliton(xi, 0, 3*delta, alpha, beta), ref)


def test_spectral_convergence_in_space():
    L = 120.
    errs = []
    for n in (96, 192, 384):
        x = np.linspace(-L/2, L/2, n, endpoint=False)
        _, U = kdv_propagate(kdv_soliton(x, 0), length=L, final_time=1, dt=1e-3)
        errs.append(np.max(np.abs(U[-1] - kdv_soliton(x, 1))))
    assert errs[0] > 1e3*errs[1] and errs[2] < 1e-12


def test_fourth_order_in_time():
    L, n = 120., 512
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    u0, ref = kdv_soliton(x, 0, amplitude=2.), kdv_soliton(x, 2, amplitude=2.)
    e = [np.max(np.abs(kdv_propagate(u0, length=L, final_time=2, dt=dt)[1][-1] - ref))
         for dt in (0.02, 0.01)]
    assert 3.6 < np.log2(e[0]/e[1]) < 4.4


def test_long_time_invariant_drift():
    L, n = 120., 256
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    u0 = kdv_soliton(x, 0, center=-30)
    _, U = kdv_propagate(u0, length=L, final_time=20, dt=0.01, save_every=500)
    i0, i1 = kdv_invariants(U[0], L/n), kdv_invariants(U[-1], L/n)
    for key in i0:
        assert abs(i1[key] - i0[key]) <= 1e-9*abs(i0[key])


def test_two_soliton_phase_shifts():
    L, n, T = 150., 512, 100.
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    u0 = kdv_soliton(x, 0, amplitude=2.0, center=-50) + \
        kdv_soliton(x, 0, amplitude=0.5, center=-25)
    _, U = kdv_propagate(u0, length=L, final_time=T, dt=0.02, save_every=10**6)
    u = U[-1]
    k1, k2 = np.sqrt(2.0/12), np.sqrt(0.5/12)
    lg = np.log((k1 + k2)/(k1 - k2))
    xb = x[np.argmax(u)]
    xs = x[np.argmax(np.where(x < xb - 10, u, 0))]
    dx = L/n
    assert abs((xb - (-50 + 2*T/3)) - lg/k1) < dx
    assert abs((xs - (-25 + 0.5*T/3)) + lg/k2) < dx
    assert abs(u.max() - 2.0) < 1e-3


def test_cli_runs(tmp_path):
    main(['sagdeev', '--mach', '1.3', '--output', str(tmp_path/'s')])
    main(['kdv', '--n', '128', '--time', '0.1', '--dt', '0.01',
          '--output', str(tmp_path/'k')])
    assert (tmp_path/'s'/'sagdeev_profile.csv').exists()
    assert (tmp_path/'k'/'kdv_final.csv').exists()
