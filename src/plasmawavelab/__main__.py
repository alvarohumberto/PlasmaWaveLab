"""Command-line interface: ``python -m plasmawavelab {sagdeev,kdv} [options]``."""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from . import __version__
from .core import (sagdeev, solitary_profile, critical_mach, kdv_soliton,
                   kdv_propagate, kdv_invariants)


def build_parser():
    ap = argparse.ArgumentParser(
        prog='plasmawavelab',
        description='PlasmaWaveLab reproducible examples')
    ap.add_argument('model', choices=['sagdeev', 'kdv'])
    ap.add_argument('--version', action='version',
                    version=f'%(prog)s {__version__}')
    ap.add_argument('--output', default='results',
                    help='output directory for PNG and CSV files')
    g = ap.add_argument_group('Sagdeev options')
    g.add_argument('--mach', type=float, default=1.2,
                   help=f'Mach number, 1 < M < M_c ~ {critical_mach():.4f}')
    g.add_argument('--xi-max', type=float, default=45.0,
                   help='half-width of the travelling-coordinate window')
    k = ap.add_argument_group('KdV options')
    k.add_argument('--n', type=int, default=512, help='even number of grid points')
    k.add_argument('--length', type=float, default=120.0, help='periodic domain length')
    k.add_argument('--time', type=float, default=2.0, help='final time')
    k.add_argument('--dt', type=float, default=0.002, help='time step')
    k.add_argument('--amplitude', type=float, default=1.0, help='soliton amplitude A')
    k.add_argument('--alpha', type=float, default=1.0, help='nonlinear coefficient')
    k.add_argument('--beta', type=float, default=1.0, help='dispersion coefficient')
    return ap


def run_sagdeev(args, output):
    x, p, amp, end = solitary_profile(args.mach, xi_max=args.xi_max)
    phi = np.linspace(0, amp, 500)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].plot(phi, sagdeev(phi, args.mach))
    axs[0].axhline(0, color='0.6', lw=0.8)
    axs[0].set(xlabel=r'$\phi$', ylabel=r'$S(\phi;M)$', title=f'M = {args.mach:g}')
    axs[1].plot(x, p)
    axs[1].set(xlabel=r'$\xi$', ylabel=r'$\phi$', title=rf'$\phi_m$ = {amp:.6g}')
    fig.tight_layout(); fig.savefig(output / 'sagdeev.png', dpi=170); plt.close(fig)
    np.savetxt(output / 'sagdeev_profile.csv', np.c_[x, p], delimiter=',',
               header='xi,phi', comments='')
    print(f'Sagdeev amplitude = {amp:.9g}; profile terminated at xi={end:.3g}')


def run_kdv(args, output):
    L, n = args.length, args.n
    x = np.linspace(-L/2, L/2, n, endpoint=False)
    kw = dict(amplitude=args.amplitude, alpha=args.alpha, beta=args.beta)
    u0 = kdv_soliton(x, 0, **kw)
    t, U = kdv_propagate(u0, length=L, final_time=args.time, dt=args.dt,
                         alpha=args.alpha, beta=args.beta)
    exact = kdv_soliton(x, args.time, **kw)
    error = np.max(np.abs(U[-1] - exact))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, U[0], label='initial')
    ax.plot(x, U[-1], label='numerical')
    ax.plot(x, exact, '--', label='exact')
    ax.set(xlabel='x', ylabel='u', title=f'KdV soliton, max error = {error:.3g}')
    ax.legend(); fig.tight_layout(); fig.savefig(output / 'kdv.png', dpi=170); plt.close(fig)
    np.savetxt(output / 'kdv_final.csv', np.c_[x, U[-1], exact], delimiter=',',
               header='x,numerical,exact', comments='')
    inv0 = kdv_invariants(U[0], L/n, args.alpha, args.beta)
    inv1 = kdv_invariants(U[-1], L/n, args.alpha, args.beta)
    print(f'KdV soliton max error = {error:.6g}; invariants initial={inv0}; final={inv1}')


def main(argv=None):
    args = build_parser().parse_args(argv)
    output = Path(args.output); output.mkdir(exist_ok=True, parents=True)
    if args.model == 'sagdeev':
        run_sagdeev(args, output)
    else:
        run_kdv(args, output)


if __name__ == '__main__':
    main()
