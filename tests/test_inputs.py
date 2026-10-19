import numpy as np
import pytest
from plasmawavelab import peak_potential, sagdeev, kdv_soliton, kdv_propagate

@pytest.mark.parametrize('mach', [0.0, 0.9, 1.0])
def test_invalid_mach(mach):
    with pytest.raises(ValueError):
        peak_potential(mach)

def test_sagdeev_singularity_rejected():
    with pytest.raises(ValueError):
        sagdeev(0.72, 1.2)

def test_invalid_grid_rejected():
    with pytest.raises(ValueError):
        kdv_propagate(np.ones(15))

def test_kdv_zero_time():
    u = np.arange(16, dtype=float)
    t, v = kdv_propagate(u, final_time=0)
    assert len(t) == 1 and np.array_equal(v[0], u)

def test_soliton_speed_and_amplitude():
    amplitude = 2.0; alpha = 3.0; t = 1.5
    speed = alpha * amplitude / 3
    assert abs(float(kdv_soliton(speed * t,t,amplitude=amplitude,alpha=alpha)) - amplitude) < 1e-12
