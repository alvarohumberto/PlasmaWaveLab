"""Reproduce every figure and number reported in the SoftwareX manuscript.

Usage (from the repository root, after ``pip install -e .``):

    python examples/verification.py --output figures/verification

Writes PNG/PDF figures, CSV tables and ``summary.json`` into ``--output``.
Typical runtime on a laptop: well under one minute.
"""
import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

import plasmawavelab as pwl
from plasmawavelab import (sagdeev, critical_mach, peak_potential,
                           solitary_profile, small_amplitude_soliton,
                           kdv_soliton, kdv_propagate, kdv_invariants)

# Fixed categorical order (colour-blind-checked), always paired with
# line style / marker so identity is never carried by colour alone.
C = ['#2a78d6', '#eb6834', '#1baf7a', '#6250d6']
INK, MUTED = '#0b0b0b', '#8a8984'
plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'legend.fontsize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': MUTED, 'xtick.color': '#52514e', 'ytick.color': '#52514e',
    'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.6,
    'lines.linewidth': 1.6, 'savefig.bbox': 'tight', 'figure.dpi': 110,
})


def save(fig, out, name):
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=300)
    plt.close(fig)


def fig_sagdeev(out, summary):
    mc = critical_mach()
    machs = [1.1, 1.2, 1.4, 1.58]
    styles = ['-', '--', '-.', ':']
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for m, c, ls in zip(machs, C, styles):
        amp = peak_potential(m)
        phi = np.linspace(0, min(amp*1.06, m*m/2*(1-1e-12)), 800)
        a.plot(phi, sagdeev(phi, m), color=c, ls=ls, label=f'M = {m:g}')
        x, p, _, _ = solitary_profile(m, xi_max=25, points=1001)
        b.plot(x, p, color=c, ls=ls, label=f'M = {m:g}')
    a.axhline(0, color=MUTED, lw=0.8)
    a.set(xlabel=r'$\phi$', ylabel=r'$S(\phi;M)$', ylim=(None, 0.06),
          title='(a) Sagdeev pseudopotential')
    b.set(xlabel=r'$\xi$', ylabel=r'$\phi(\xi)$', xlim=(-15, 15),
          title='(b) solitary-wave profiles')
    b.legend(frameon=False, loc='upper right')
    fig.tight_layout(); save(fig, out, 'fig_sagdeev')
    summary['critical_mach'] = mc
    summary['phi_m'] = {f'{m:g}': peak_potential(m) for m in machs}


def fig_existence(out, summary):
    mc = critical_mach()
    M = np.linspace(1.005, mc - 1e-4, 120)
    amp = np.array([peak_potential(m) for m in M])
    deltas = np.array([0.005, 0.01, 0.02, 0.04, 0.08, 0.16])
    amp_err, prof_err = [], []
    for d in deltas:
        m = 1 + d
        a0 = peak_potential(m)
        x, p, _, _ = solitary_profile(m, xi_max=max(60, 12/np.sqrt(d/2)), points=4001)
        amp_err.append(abs(a0/(3*d) - 1))
        prof_err.append(np.max(np.abs(p - small_amplitude_soliton(x, m)))/(3*d))
    amp_err, prof_err = np.array(amp_err), np.array(prof_err)
    slope = np.polyfit(np.log(deltas), np.log(prof_err), 1)[0]
    np.savetxt(out / 'mach_sweep.csv', np.c_[M, amp], delimiter=',',
               header='mach,phi_m', comments='')
    np.savetxt(out / 'small_amplitude.csv', np.c_[deltas, amp_err, prof_err],
               delimiter=',', header='delta,rel_amp_error,rel_profile_error', comments='')

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    a.plot(M, amp, color=C[0], label=r'$\phi_m$ (PlasmaWaveLab)')
    a.plot(M, 3*(M-1), color=C[1], ls='--', label=r'KdV limit $3(M-1)$')
    a.plot(M, M**2/2, color=MUTED, ls=':', label=r'ion limit $M^2/2$')
    a.axvline(mc, color=MUTED, lw=0.8)
    a.annotate(rf'$M_c$ = {mc:.4f}', (mc, 0.15), xytext=(-6, 0),
               textcoords='offset points', ha='right', color='#52514e', fontsize=8)
    a.set(xlabel='Mach number M', ylabel='amplitude', title='(a) existence domain')
    a.legend(frameon=False, loc='upper left')
    b.loglog(deltas, amp_err, 'o-', color=C[0], ms=4, label='amplitude')
    b.loglog(deltas, prof_err, 's--', color=C[1], ms=4, label='profile (max norm)')
    b.loglog(deltas, deltas, color=MUTED, lw=0.8, label=r'$\propto\delta$')
    b.set(xlabel=r'$\delta = M-1$', ylabel='relative error vs KdV soliton',
          title='(b) small-amplitude limit')
    b.legend(frameon=False)
    fig.tight_layout(); save(fig, out, 'fig_existence')
    summary['small_amplitude'] = {'delta': deltas.tolist(),
                                  'rel_amp_error': amp_err.tolist(),
                                  'rel_profile_error': prof_err.tolist(),
                                  'observed_order': float(slope)}


def fig_convergence(out, summary):
    L = 120.
    ns = np.array([48, 64, 96, 128, 160, 192, 256, 320, 384])
    es = []
    for n in ns:
        x = np.linspace(-L/2, L/2, n, endpoint=False)
        _, U = kdv_propagate(kdv_soliton(x, 0), length=L, final_time=2, dt=1e-3)
        es.append(np.max(np.abs(U[-1] - kdv_soliton(x, 2))))
    es = np.array(es)
    n = 512
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    u0, ref = kdv_soliton(x, 0, amplitude=2.), kdv_soliton(x, 4, amplitude=2.)
    dts = np.array([0.08, 0.04, 0.02, 0.01, 0.005])
    et = np.array([np.max(np.abs(kdv_propagate(u0, length=L, final_time=4, dt=dt)[1][-1] - ref))
                   for dt in dts])
    orders = np.log2(et[:-1]/et[1:])
    np.savetxt(out / 'convergence_space.csv', np.c_[ns, es], delimiter=',',
               header='n,max_error', comments='')
    np.savetxt(out / 'convergence_time.csv', np.c_[dts, et], delimiter=',',
               header='dt,max_error', comments='')

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    a.semilogy(ns, es, 'o-', color=C[0], ms=4)
    a.set(xlabel='grid points N', ylabel=r'$\|u_N-u_{\rm sol}\|_\infty$',
          title='(a) spatial: spectral convergence')
    b.loglog(dts, et, 'o-', color=C[0], ms=4, label='IF-RK4')
    b.loglog(dts, et[-1]*(dts/dts[-1])**4, color=MUTED, ls='--', lw=0.9,
             label=r'slope 4')
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b.set(xlabel=r'time step $\Delta t$', ylabel=r'$\|u-u_{\rm sol}\|_\infty$',
          title='(b) temporal: fourth order')
    b.legend(frameon=False)
    fig.tight_layout(); save(fig, out, 'fig_convergence')
    summary['spatial'] = dict(zip(map(int, ns), map(float, es)))
    summary['temporal'] = {'dt': dts.tolist(), 'error': et.tolist(),
                           'observed_orders': orders.tolist()}


def peak_position(x, u, mask):
    """Sub-grid crest location by Fourier (band-limited) refinement."""
    i = int(np.argmax(np.where(mask, u, -np.inf)))
    n, dx = len(x), x[1] - x[0]
    uh = np.fft.fft(u); k = 2*np.pi*np.fft.fftfreq(n, d=dx)
    xs = np.linspace(x[i] - dx, x[i] + dx, 2001)
    vals = np.real(np.exp(1j*np.outer(xs - x[0], k)) @ uh)/n
    j = int(np.argmax(vals))
    return xs[j], vals[j]


def phase_shifts(A1, A2, alpha=1., beta=1.):
    """Exact KdV two-soliton phase shifts (A1 > A2) for u_t+alpha u u_x+beta u_xxx=0."""
    k1, k2 = np.sqrt(alpha*A1/(12*beta)), np.sqrt(alpha*A2/(12*beta))
    lg = np.log((k1 + k2)/(k1 - k2))
    return lg/k1, -lg/k2


def fig_collision(out, summary):
    L, n = 150., 512
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    u0 = kdv_soliton(x, 0, amplitude=2.0, center=-50) + \
        kdv_soliton(x, 0, amplitude=0.5, center=-25)
    t0 = time.perf_counter()
    t, U = kdv_propagate(u0, length=L, final_time=100, dt=0.005, save_every=100)
    wall = time.perf_counter() - t0
    inv = [kdv_invariants(u, L/n) for u in U]
    T = t[-1]
    xb, ab = peak_position(x, U[-1], np.ones_like(x, bool))
    xs_, as_ = peak_position(x, U[-1], x < xb - 10)
    d1_num, d2_num = xb - (-50 + 2.0*T/3), xs_ - (-25 + 0.5*T/3)
    d1_th, d2_th = phase_shifts(2.0, 0.5)
    drift = {k: np.abs(np.array([i[k] for i in inv]) - inv[0][k])/abs(inv[0][k])
             for k in inv[0]}

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9),
                               gridspec_kw={'width_ratios': [1.25, 1]})
    im = a.pcolormesh(x, t, U, cmap='Blues', shading='auto', vmin=0, vmax=2.0,
                      rasterized=True)
    a.grid(False)
    a.set(xlabel='x', ylabel='t', title='(a) two-soliton collision, u(x,t)')
    fig.colorbar(im, ax=a, pad=0.02, label='u')
    for (k, d), c, ls, lab in zip(drift.items(), C, ['-', '--', ':'],
                                  ['mass', r'$L^2$ norm', 'Hamiltonian']):
        b.semilogy(t[1:], np.maximum.accumulate(np.maximum(d[1:], 1e-16)), color=c, ls=ls, label=lab)
    b.set(xlabel='t', ylabel=r'max relative drift up to $t$', ylim=(3e-17, 1e-6),
          title='(b) invariant drift')
    b.legend(frameon=False, loc='center right')
    fig.tight_layout(); save(fig, out, 'fig_collision')
    summary['collision'] = {'n': n, 'dt': 0.005, 'final_time': 100,
                            'steps': int(round(100/0.005)), 'wall_seconds': wall,
                            'max_rel_drift': {k: float(v.max()) for k, v in drift.items()},
                            'final_amplitudes': [float(ab), float(as_)],
                            'phase_shift_numerical': [float(d1_num), float(d2_num)],
                            'phase_shift_exact': [float(d1_th), float(d2_th)]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--output', default='figures/verification')
    out = Path(ap.parse_args().output); out.mkdir(parents=True, exist_ok=True)
    import scipy
    cpu = platform.processor()
    try:
        with open('/proc/cpuinfo') as fh:
            cpu = next(l.split(':', 1)[1].strip() for l in fh if l.startswith('model name'))
    except (OSError, StopIteration):
        pass
    summary = {'plasmawavelab': pwl.__version__, 'python': platform.python_version(),
               'numpy': np.__version__, 'scipy': scipy.__version__,
               'matplotlib': matplotlib.__version__, 'machine': platform.machine(),
               'system': platform.system(), 'cpu': cpu}
    t0 = time.perf_counter()
    fig_sagdeev(out, summary)
    fig_existence(out, summary)
    fig_convergence(out, summary)
    fig_collision(out, summary)
    summary['total_wall_seconds'] = time.perf_counter() - t0
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
