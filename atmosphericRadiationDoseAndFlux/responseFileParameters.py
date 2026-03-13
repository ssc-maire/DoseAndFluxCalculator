import numpy as np
from numba import njit, jit
from typing import Tuple, Dict, List, Union, Optional

from .particle import Particle
from .units import Distance

@njit
def calculate_altitude_layer_params(altitude_meters: float, altitude_km: float) -> Tuple[int, float]:
    """
    Calculate altitude layer index and interpolation factor using Numba optimization.
    
    Parameters
    ----------
    altitude_meters : float
        Altitude in meters
    altitude_km : float
        Altitude in kilometers
        
    Returns
    -------
    Tuple[int, float]
        Tuple containing:
        - Altitude layer index
        - Interpolation factor (f1)
        
    Raises
    ------
    Exception
        If altitude is greater than 100 km
    """
    if altitude_km > 100.0:
        # Can't use exceptions in njit, so return invalid values that will be caught
        return -1, -1.0

    if altitude_km < 0.025:
        layer = 1
        f1 = 1.0
    elif altitude_km < 1.025:
        layer = int((int(altitude_meters)-25)//50) + 1
        hr = (int(altitude_meters)-25) % 50
        f1 = 1.0 - hr/50.0
    elif altitude_km < 1.15:
        layer = 21
        f1 = 1.0 - (int(altitude_meters)-1025)/1025.0
    elif altitude_km < 5.05:
        layer = int((int(altitude_meters)-1150)//100) + 22
        hr = (int(altitude_meters)-1150) % 100
        f1 = 1.0 - hr/100.0
    elif altitude_km < 5.3:
        layer = 61
        f1 = 1.0 - (altitude_km-5.05)/0.25
    elif altitude_km < 15.1:
        layer = int((int(altitude_meters)-5300)//200) + 62
        hr = (int(altitude_meters)-5300) % 200
        f1 = 1.0 - hr/200.0
    elif altitude_km < 16.5:
        layer = 111
        f1 = 1.0 - (altitude_km-15.1)/1.4
    elif altitude_km < 38.5:
        layer = int((int(altitude_meters)-16500)//1000) + 112
        hr = (int(altitude_meters)-16500) % 1000
        f1 = 1.0 - hr/1000.0
    elif altitude_km < 40.5:
        layer = 134
        f1 = 1.0 - (altitude_km-38.5)/2.0
    elif altitude_km < 62.5:
        layer = 135
        f1 = 1.0 - (altitude_km-40.5)/22.0
    elif altitude_km < 97.5:
        layer = 136
        f1 = 1.0 - (altitude_km-40.5)/57.5
    else:
        layer = 137
        f1 = 1.0

    # Convert from 1-indexed to 0-indexed for Python arrays
    altitude_layer_index = layer - 1
    
    return altitude_layer_index, f1

@njit
def calculate_weighted_fluxes(energyBins: np.ndarray, fluxes: np.ndarray, 
                             cutoffEnergy: float) -> Tuple[np.ndarray, int]:
    """
    Calculate weighted flux values based on energy cutoff using Numba optimization.
    
    Parameters
    ----------
    energyBins : np.ndarray
        Energy bin edges
    fluxes : np.ndarray
        Flux values
    cutoffEnergy : float
        Cutoff energy in MeV
        
    Returns
    -------
    Tuple[np.ndarray, int]
        Tuple containing:
        - Weighted flux values
        - Energy index for cutoff
    """
    # Make a copy of fluxes to avoid modifying the original
    weightedFluxes = fluxes.copy()
    
    # Find energy index for cutoff
    energyIndex = 0
    while energyIndex < len(energyBins) and cutoffEnergy > energyBins[energyIndex]:
        energyIndex += 1
    
    # Apply weighting based on cutoff energy
    if energyIndex > 0:
        energyIndex = energyIndex - 1
        
        ed = energyBins[energyIndex+1] - cutoffEnergy
        et = energyBins[energyIndex+1] - energyBins[energyIndex]
        
        fe = ed/et
        weightedFluxes[energyIndex] = weightedFluxes[energyIndex] * fe
    
    return weightedFluxes, energyIndex


class ResponseFileParameters:
    """
    Class handling response file parameters for atmospheric radiation calculations.
    
    This class determines the appropriate altitude layers, interpolation factors,
    and weighted flux values for radiation response calculations.
    
    Attributes
    ----------
    altitude : Distance
        Altitude object
    particle : Particle
        Particle object
    altitudeLayerIndex : int
        Index of the altitude layer
    f1 : float
        Interpolation factor
    altIndexAbove : int
        Index of the altitude layer above
    weightedFluxes : np.ndarray
        Array of weighted flux values
        
    Methods
    -------
    __init__(altitude, energyBins, fluxes, particle)
        Initialize the response file parameters
    getAltitudeResponseLayer(altitude)
        Calculate altitude response layer and interpolation factor
    getResponseParametersFromSpectrum(energyBins, fluxes)
        Calculate response parameters from energy spectrum
    """

    def __init__(self, altitude: Distance, energyBins: np.ndarray, 
                fluxes: np.ndarray, particle: Particle):
        """
        Initialize the response file parameters.
        
        Parameters
        ----------
        altitude : Distance
            Altitude object with distance in meters and other units
        energyBins : np.ndarray
            Energy bin edges in MeV
        fluxes : np.ndarray
            Flux values
        particle : Particle
            Particle object (proton or alpha)
        """
        self.altitude = altitude
        self.particle = particle
        
        # Get altitude response layer parameters
        self.getAltitudeResponseLayer(altitude)
        
        # Get response parameters from spectrum
        self.getResponseParametersFromSpectrum(energyBins, fluxes)

    def getAltitudeResponseLayer(self, altitude: Distance):
        """
        Calculate altitude response layer and interpolation factor.
        
        Parameters
        ----------
        altitude : Distance
            Altitude object
            
        Raises
        ------
        Exception
            If altitude is greater than 100 km
        """
        # Use Numba-optimized function for layer calculation
        altitude_layer_index, f1 = calculate_altitude_layer_params(
            altitude.meters, altitude.km
        )
        
        # Check for error from Numba function (altitude > 100 km)
        if altitude_layer_index == -1:
            raise Exception("altitude.km has to be less than 100 km!")
            
        self.altitudeLayerIndex = altitude_layer_index
        self.f1 = f1

    def getResponseParametersFromSpectrum(self, energyBins: np.ndarray, fluxes: np.ndarray):
        """
        Calculate response parameters from energy spectrum.
        
        Parameters
        ----------
        energyBins : np.ndarray
            Energy bin edges in MeV
        fluxes : np.ndarray
            Flux values
        """
        # Rigidity cutoff is not relevant for these calculations
        rigidityCutoff = 0.0
        
        # Calculate altitude index above (bounded by max layer)
        altIndexAbove = self.altitudeLayerIndex + 1
        if altIndexAbove > 137:
            altIndexAbove = 137
            
        # Calculate cutoff energy
        cutoffEnergy = 1000.0 * (np.sqrt(((rigidityCutoff**2)/self.particle.atomicMass) + (0.938**2)) - 0.938)
        
        # Use Numba-optimized function for weighted flux calculation
        # Convert to correct types for Numba
        energyBins_np = np.array(energyBins, dtype=np.float64)
        fluxes_np = np.array(fluxes, dtype=np.float64)
        
        weightedFluxes, _ = calculate_weighted_fluxes(energyBins_np, fluxes_np, cutoffEnergy)
        
        # Store results
        self.altIndexAbove = altIndexAbove
        self.weightedFluxes = weightedFluxes