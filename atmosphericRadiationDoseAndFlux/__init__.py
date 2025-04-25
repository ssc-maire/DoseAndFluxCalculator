"""
Atmospheric Radiation Dose and Flux Calculator
=============================================

A package for calculating radiation dose and flux rates in the atmosphere
for various particle types at different altitudes.

This package provides tools to:
- Calculate radiation doses from energy or rigidity spectra
- Support for protons and alpha particles (+ heavier ions)
- Calculate various dose metrics (effective dose, ambient dose equivalent, etc.)
- Calculate neutron flux at different energy ranges
- Perform calculations at different altitudes

The package has been optimized using Numba for high-performance computing.

Author: Space Environment and Protection Group, University of Surrey
"""

from .doseAndFluxCalculator import (
    calculate_from_energy_spec,
    calculate_from_rigidity_spec,
    calculate_from_energy_spec_array,
    calculate_from_rigidity_spec_array
)

from .particle import Particle
from .units import Distance

__version__ = "1.0.0"
