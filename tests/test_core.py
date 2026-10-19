import numpy as np
from plasmawavelab import (sagdeev, sagdeev_prime, peak_potential,
                           solitary_profile,kdv_soliton,kdv_propagate,kdv_invariants)

def test_sagdeev_potential_derivative():
    m=1.2; x=0.15; h=1e-6
    derivative=(sagdeev(x+h,m)-sagdeev(x-h,m))/(2*h)
    assert abs(derivative-sagdeev_prime(x,m))<1e-7

def test_turning_point():
    m=1.2; a=peak_potential(m)
    assert 0<a<m*m/2
    assert abs(sagdeev(a,m))<1e-10
    x,u,amp,end=solitary_profile(m,xi_max=25)
    assert abs(u[len(u)//2]-a)<1e-8 and np.min(u)>=-1e-8
    assert end>0

def test_kdv_exact_soliton():
    L=120.; n=512; x=np.linspace(-L/2,L/2,n,endpoint=False)
    u0=kdv_soliton(x,0)
    ts, us=kdv_propagate(u0,length=L,final_time=0.5,dt=0.002)
    err=np.sqrt(np.mean((us[-1]-kdv_soliton(x,0.5))**2))
    assert err<1e-4,err
    inv0=kdv_invariants(us[0],L/n); inv1=kdv_invariants(us[-1],L/n)
    assert abs(inv0['mass']-inv1['mass'])<1e-8
