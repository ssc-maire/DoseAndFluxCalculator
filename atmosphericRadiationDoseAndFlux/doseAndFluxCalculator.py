import numpy as np
import pandas as pd
from numba import njit, jit, prange
import warnings

from . import particle
from . import particleResponse
from . import settings
from . import units

import ParticleRigidityCalculationTools as PRCT

import inspect

from typing import Callable, List, Union, Dict, Any, Tuple

# Display warning about alpha particle behavior
print("Warning from atmosphericRadiationDoseAndFlux module: currently using an alpha particle as input actually calculates the contribution from alpha + all simulated heavier ions, rather than just alpha particles!")

# Numba-optimized function for energy integration calculation
@njit
def calculate_energy_integrated_flux(fluxes: np.ndarray, energy_differences: np.ndarray) -> np.ndarray:
    """
    Calculate energy-integrated flux values from differential flux values.
    
    Parameters:
    -----------
    fluxes : np.ndarray
        Array of differential flux values (particles/cm²/sr/(MeV/n)/s)
    energy_differences : np.ndarray
        Array of energy bin width values (MeV/n)
        
    Returns:
    --------
    np.ndarray
        Energy-integrated flux values (particles/cm²/s)
    """
    return fluxes * energy_differences * np.pi  # units of particles / cm2 / s

@settings.allowCalculationForTotalOfParticles
def calculate_from_energy_spec(
                            inputEnergyDistributionFunctionMeV: Callable, 
                            altitudesInkm: List[float], 
                            particleName: str = "proton", 
                            inputEnergyBins: np.ndarray = 10**(0.1*(np.array(range(1,52))-1)+1),
                            verticalCutOffRigidity: float = 0.0) -> pd.DataFrame:
    """
    Calculate dose and flux rates using an energy spectrum as input.
    
    Parameters:
    -----------
    inputEnergyDistributionFunctionMeV : Callable
        Function that takes kinetic energy per nucleon (MeV/n) and returns
        differential flux (particles/cm²/sr/(MeV/n)/s)
    altitudesInkm : List[float]
        List of altitudes in kilometers
    particleName : str, default="proton"
        Particle type ("proton", "alpha", or "both")
    inputEnergyBins : np.ndarray, default=10**(0.1*(np.array(range(1,52))-1)+1)
        Energy bin edges in MeV/n. The default is the MAIRE 50-bin grid.
    verticalCutOffRigidity : float, default=0.0
        Vertical cutoff rigidity in GV
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing dose and flux rates at each altitude
        
    Notes:
    ------
    This function calculates radiation dose and flux rates based on the input energy distribution.
    The vertical cutoff rigidity feature is currently disabled pending review.
    """
    try:
        assert verticalCutOffRigidity == 0.0
    except AssertionError:
        raise Exception("Vertical cutoff rigidity is not currently supported, testing has revealed an error in previous vertical cut-off rigidity calculations, please manually set your spectrum to include vertical cut-off rigidity for now!")

    # Calculate energy bin midpoints
    energyBinMidPoints = (inputEnergyBins[1:] + inputEnergyBins[:-1])/2
    
    # Get flux values at each midpoint
    inputFluxesMeV = list(map(inputEnergyDistributionFunctionMeV, energyBinMidPoints))

    # Calculate dose and flux rates
    outputDF = calculate_from_energy_spec_array(
                            inputEnergyBins,
                            inputFluxesMeV, 
                            altitudesInkm, 
                            particleName=particleName,
                            verticalCutOffRigidity=verticalCutOffRigidity)

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_rigidity_spec(
                            inputRigidityDistributionFunctionGV: Callable, 
                            altitudesInkm: List[float], 
                            particleName: str = "proton", 
                            inputRigidityBins: np.ndarray = None,
                            verticalCutOffRigidity: float = 0.0) -> pd.DataFrame:
    """
    Calculate dose and flux rates using a rigidity spectrum as input.
    
    Parameters:
    -----------
    inputRigidityDistributionFunctionGV : Callable
        Function that takes rigidity (GV) as input and returns differential flux (particles/cm²/sr/GV/s)
    altitudesInkm : List[float]
        List of altitudes in kilometers
    particleName : str, default="proton"
        Particle type ("proton", "alpha", or "both")
    inputRigidityBins : np.ndarray, optional
        Rigidity bin edges in GV, if None will be derived from default energy bins
    verticalCutOffRigidity : float, default=0.0
        Vertical cutoff rigidity in GV
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing dose and flux rates at each altitude
        
    Notes:
    ------
    This function converts rigidity-based inputs to the MAIRE MeV/n energy grid
    internally. Default rigidity bins are generated from that MeV/n grid using
    total kinetic energy \(A \times E_n\), not by treating MeV/n as total MeV.
    """
    # Initialize particle object for calculations
    particleForCalculations = particle.Particle(particleName)

    # Generate rigidity bins from the MAIRE MeV/n energy grid if not provided
    if inputRigidityBins is None:
        inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1)
        inputRigidityBins = np.array(PRCT.convertPerNucleonEnergyToTotalRigidity(
            inputEnergyBins,
            particleMassAU=particleForCalculations.atomicMass,
            particleChargeAU=particleForCalculations.atomicCharge))

    # Calculate rigidity bin midpoints
    rigidityBinMidPoints = (inputRigidityBins[1:] + inputRigidityBins[:-1])/2
    
    # Get flux values at each midpoint
    inputFluxesGV = list(map(inputRigidityDistributionFunctionGV, rigidityBinMidPoints))

    # Calculate dose and flux rates
    outputDF = calculate_from_rigidity_spec_array(
                            inputRigidityBins,
                            inputFluxesGV, 
                            altitudesInkm, 
                            particleName=particleName,
                            verticalCutOffRigidity=verticalCutOffRigidity)

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_energy_spec_array(
                            inputEnergyBins: np.ndarray,
                            inputFluxesMeV: Union[np.ndarray, List[float]], 
                            altitudesInkm: Union[np.ndarray, List[float]], 
                            particleName: str = "proton",
                            verticalCutOffRigidity: float = 0.0) -> pd.DataFrame:
    """
    Calculate dose and flux rates using an array of energy bins and flux values.
    
    Parameters:
    -----------
    inputEnergyBins : np.ndarray
        Energy bin edges in MeV/n
    inputFluxesMeV : Union[np.ndarray, List[float]]
        Differential flux values at bin midpoints (particles/cm²/sr/(MeV/n)/s)
    altitudesInkm : Union[np.ndarray, List[float]]
        Altitudes in kilometers
    particleName : str, default="proton"
        Particle type ("proton", "alpha", or "both")
    verticalCutOffRigidity : float, default=0.0
        Vertical cutoff rigidity in GV
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing dose and flux rates at each altitude
        
    Notes:
    ------
    This is a lower-level function that handles array inputs directly.
    """
    # Format input variables and apply cutoff rigidity if specified
    particleForCalculations, inputEnergyBinsArray, altitudesInkmArray, inputFluxesArray = settings.formatInputVariables(
        inputEnergyBins, 
        inputFluxesMeV, 
        altitudesInkm, 
        particleName, 
        verticalCutOffRigidity
    )

    # Calculate energy bin widths
    inputEnergyDifferences = inputEnergyBinsArray[1:] - inputEnergyBinsArray[:-1]
    
    # Calculate integrated flux values (optimized with Numba)
    inputFluxesIntegrated = calculate_energy_integrated_flux(inputFluxesArray, inputEnergyDifferences)

    # Initialize output DataFrame
    outputDF = pd.DataFrame({
        "altitude (km)": altitudesInkmArray,
    })

    # Calculate doses for each response type
    for doseType in particleResponse.fullListOfDoseResponseTypes:
        doseResponseForParticle = particleResponse.fullDoseResponseDict[doseType](particleForCalculations, doseType)
    
        coordinateDosesList = []
        for altitudeInkm in altitudesInkmArray:
            altitude = units.Distance(altitudeInkm * 1000.0)
            totalDose = doseResponseForParticle.calculateDose(altitude, inputEnergyBinsArray, inputFluxesIntegrated)
            coordinateDosesList.append(totalDose)

        outputDF[doseType] = coordinateDosesList

    # Calculate SEU and SEL rates from thermal neutron flux
    outputDF["SEU"] = outputDF["tn2"] * 1e-13
    outputDF["SEL"] = outputDF["tn2"] * 1e-8

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_rigidity_spec_array(
                            inputRigidityBins: np.ndarray,
                            inputFluxesGV: Union[np.ndarray, List[float], Callable], 
                            altitudesInkm: Union[np.ndarray, List[float]], 
                            particleName: str = "proton",
                            verticalCutOffRigidity: float = 0.0) -> pd.DataFrame:
    """
    Calculate dose and flux rates using an array of rigidity bins and flux values.
    
    Parameters:
    -----------
    inputRigidityBins : np.ndarray
        Rigidity bin edges in GV
    inputFluxesGV : Union[np.ndarray, List[float], Callable]
        Differential flux values at bin midpoints (particles/cm²/sr/GV/s) or function to calculate them
    altitudesInkm : Union[np.ndarray, List[float]]
        Altitudes in kilometers
    particleName : str, default="proton"
        Particle type ("proton", "alpha", or "both")
    verticalCutOffRigidity : float, default=0.0
        Vertical cutoff rigidity in GV
        
    Returns:
    --------
    pd.DataFrame
        DataFrame containing dose and flux rates at each altitude
        
    Notes:
    ------
    This function converts rigidity-based inputs to the MAIRE MeV/n energy grid
    internally. Rigidity is total GV; energy is kinetic energy per nucleon.
    """
    # Initialize particle object for calculations
    particleForCalculations = particle.Particle(particleName)

    # Convert total rigidity bins to MeV/n energy bins
    inputEnergyBins = np.array(PRCT.convertTotalRigidityToPerNucleonEnergy(inputRigidityBins,
                                particleMassAU=particleForCalculations.atomicMass,
                                particleChargeAU=particleForCalculations.atomicCharge))

    # Calculate rigidity bin midpoints
    rigidityMidPoints = (inputRigidityBins[1:] + inputRigidityBins[:-1])/2

    # Process input flux values
    if inspect.isfunction(inputFluxesGV):
        inputFluxesGVarray = np.array(list(map(inputFluxesGV, rigidityMidPoints)))
    elif isinstance(inputFluxesGV, (int, float, list, np.ndarray)):
        inputFluxesGVarray = settings.convertListOrFloatToArray(inputFluxesGV)
    else:
        raise Exception("inputFluxesGV not specified as a valid type!")

    # Convert total-rigidity fluxes to per-nucleon energy fluxes
    inputFluxesMeV = PRCT.convertTotalRigiditySpecToPerNucleonEnergySpec(rigidityMidPoints,
                                                                   inputFluxesGVarray, 
                                                                   particleMassAU=particleForCalculations.atomicMass,
                                                                   particleChargeAU=particleForCalculations.atomicCharge)["Energy distribution values"]

    # Calculate dose and flux rates using energy-based method
    outputDoseFluxRates = calculate_from_energy_spec_array(
                            inputEnergyBins,
                            inputFluxesMeV, 
                            altitudesInkm, 
                            particleName=particleName,
                            verticalCutOffRigidity=verticalCutOffRigidity)

    return outputDoseFluxRates

