"""Sagdeev electrostatic ion-acoustic structures and Fourier KdV propagation.

All quantities dimensionless. Sagdeev: cold positive ions, Boltzmann electrons,
stationary travelling wave, equilibrium densities unity. KdV is a separate
weakly nonlinear reduced model with independently configurable coefficients.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp


def sagdeev(phi, mach):
    p = np.asarray(phi, dtype=float)
    if mach <= 1:
        raise ValueError('Supersonic Mach number M > 1 required')
    if np.any(p >= mach**2 / 2):
        raise ValueError('Ion-density singularity: phi must be < M^2/2')
    return -np.expm1(p) + mach**2 * (1 - np.sqrt(1 - 2*p/mach**2))


def sagdeev_prime(phi, mach):
    p = np.asarray(phi, dtype=float)
    if np.any(p >= mach**2/2):
        raise ValueError('phi must be below M^2/2')
    return 1/np.sqrt(1-2*p/mach**2)-np.exp(p)


def critical_mach():
    """Upper Mach limit M_c of the cold-ion/Boltzmann-electron soliton branch.

    M_c solves S(M^2/2; M) = 0, i.e. 1 - exp(M^2/2) + M^2 = 0
    (M_c ~= 1.5852).  For M >= M_c the ion density becomes singular before
    the pseudopotential returns to zero, so no solitary wave exists.
    """
    return float(brentq(lambda m: 1 - np.exp(m*m/2) + m*m, 1.1, 2.0,
                        xtol=1e-15))


def ion_acoustic_kdv_coefficients():
    """Coefficients (alpha, beta) of the ion-acoustic KdV limit.

    In the present normalisation the weakly nonlinear limit M = 1 + delta,
    phi = O(delta) of the Sagdeev model is u_t + u u_x + u_xxx/2 = 0 in the
    frame moving at the ion-acoustic speed (Washimi & Taniuti 1966).
    """
    return 1.0, 0.5


def small_amplitude_soliton(xi, mach):
    """Leading-order KdV approximation of the Sagdeev solitary wave.

    phi(xi) = 3*delta*sech^2(sqrt(delta/2)*xi), delta = M - 1.  Its relative
    error with respect to the exact fluid--Poisson soliton is O(delta).
    """
    if mach <= 1:
        raise ValueError('M must exceed 1')
    delta = mach - 1
    alpha, beta = ion_acoustic_kdv_coefficients()
    # Equivalent to kdv_soliton(xi, t=0, amplitude=3*delta, alpha, beta).
    return kdv_soliton(xi, 0.0, amplitude=3*delta, alpha=alpha, beta=beta)


def peak_potential(mach, xtol=1e-13, samples=2500):
    """Locate nonzero positive Sagdeev turning point; raise if none exists.

    S is sampled at ``samples`` points below the ion-density limit M^2/2 to
    bracket the root, which is then refined by Brent's method to ``xtol``.
    """
    if mach <= 1:
        raise ValueError('M must exceed 1')
    top = mach**2 / 2 * (1-1e-11)
    grid = np.linspace(max(1e-8, top*1e-5), top, samples)
    vals = sagdeev(grid, mach)
    # Ignore roundoff at phi~0: seek negative interval followed by positive.
    neg = np.where(vals < -1e-12)[0]
    if not len(neg):
        raise ValueError('No resolved positive-polarity Sagdeev well for this M')
    last = neg[-1]
    positives = np.where(vals[last+1:] > 0)[0]
    if not len(positives):
        raise ValueError('No positive turning point before ion-density limit '
                         f'(M={mach:g} >= critical Mach {critical_mach():.6f})')
    upper = last + 1 + positives[0]
    return float(brentq(lambda p: float(sagdeev(p,mach)), grid[last], grid[upper], xtol=xtol))


def solitary_profile(mach=1.2, xi_max=45, points=1201, rtol=2e-10, atol=2e-12,
                     tail=1e-5):
    """Symmetric profile from crest by solving phi'' = -S'(phi).

    Integration (Dormand-Prince, tolerances ``rtol``/``atol``) stops when phi
    falls to ``tail*phi_m``; the profile is zero beyond that point.
    Returns (xi, phi, phi_m, xi_end).
    """
    amp = peak_potential(mach)
    x_half = np.linspace(0, xi_max, points//2+1)
    def fun(x, y):
        return [y[1], -float(sagdeev_prime(y[0], mach))]
    # Small-amplitude far-field is saddle-unstable; terminate near zero rather
    # than integrating through zero and attempting artificial negative tails.
    def event(x, y): return y[0] - amp*tail
    event.terminal = True
    event.direction = -1
    sol = solve_ivp(fun, (0, xi_max), [amp,0.0], rtol=rtol, atol=atol,
                    dense_output=True, events=event, max_step=0.05)
    if not sol.success:
        raise RuntimeError(sol.message)
    end = float(sol.t[-1])
    x = np.linspace(-xi_max,xi_max,points)
    profile = np.zeros_like(x)
    inside = np.abs(x)<=end
    profile[inside] = sol.sol(np.abs(x[inside]))[0]
    return x, profile, amp, end


def kdv_soliton(x,t,amplitude=1.,alpha=1.,beta=1.,center=0.):
    """Exact solitary wave of u_t + alpha*u*u_x + beta*u_xxx = 0."""
    if amplitude <= 0 or alpha <= 0 or beta <= 0:
        raise ValueError('amplitude, alpha, beta must be positive')
    kappa=np.sqrt(alpha*amplitude/(12*beta))
    z=kappa*(np.asarray(x)-center-alpha*amplitude*t/3)
    return amplitude / np.cosh(np.clip(z,-350,350))**2


def kdv_propagate(u0, length=120., final_time=2., dt=0.002,
                  alpha=1., beta=1., save_every=100):
    """Dealiased Fourier integrating-factor RK4 on a periodic grid.

    Returns (times, states) with shape (n_snapshots, n_grid).  Discrete
    invariants are diagnostics, not guaranteed exact under time stepping.
    """
    u=np.asarray(u0,dtype=float)
    n=len(u)
    if n<16 or n%2 or length<=0 or dt<=0 or final_time<0 or save_every<1:
        raise ValueError('Require even n>=16, positive length/dt/save_every')
    dx=length/n
    k=2*np.pi*np.fft.fftfreq(n,d=dx)
    L=1j*beta*k**3
    mask=(np.abs(np.fft.fftfreq(n)*n) <= n//3)
    U=np.fft.fft(u)
    def nonlin(W):
        a=np.fft.ifft(W).real
        return -0.5j*alpha*k*(np.fft.fft(a*a)*mask)
    steps=int(np.ceil(final_time/dt))
    h=final_time/steps if steps else 0.
    E=np.exp(L*h); Eh=np.exp(L*h/2)
    times=[0.]; snapshots=[u.copy()]
    for step in range(1,steps+1):
        a=nonlin(U)
        b=nonlin(Eh*(U+h*a/2))/Eh
        c=nonlin(Eh*(U+h*b/2))/Eh
        d=nonlin(E*(U+h*c))/E
        U=E*(U+h*(a+2*b+2*c+d)/6)
        if step%save_every==0 or step==steps:
            times.append(step*h); snapshots.append(np.fft.ifft(U).real.copy())
    return np.asarray(times),np.asarray(snapshots)


def kdv_invariants(u,dx,alpha=1.,beta=1.):
    """Approximate KdV mass, L2 and Hamiltonian on periodic grid."""
    u=np.asarray(u)
    n=len(u); k=2*np.pi*np.fft.fftfreq(n,d=dx)
    ux=np.fft.ifft(1j*k*np.fft.fft(u)).real
    return {'mass':float(dx*np.sum(u)), 'l2':float(dx*np.sum(u*u)),
            'hamiltonian':float(dx*np.sum(beta*ux*ux/2-alpha*u**3/6))}
