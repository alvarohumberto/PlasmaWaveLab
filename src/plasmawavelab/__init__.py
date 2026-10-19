"""PlasmaWaveLab: reproducible electrostatic plasma-wave benchmarks."""
from .core import (sagdeev, sagdeev_prime, critical_mach, peak_potential,
                   solitary_profile, small_amplitude_soliton,
                   ion_acoustic_kdv_coefficients, kdv_soliton, kdv_propagate,
                   kdv_invariants)

__version__ = '0.3.0'

__all__ = ['sagdeev', 'sagdeev_prime', 'critical_mach', 'peak_potential',
           'solitary_profile', 'small_amplitude_soliton',
           'ion_acoustic_kdv_coefficients', 'kdv_soliton', 'kdv_propagate',
           'kdv_invariants', '__version__']
