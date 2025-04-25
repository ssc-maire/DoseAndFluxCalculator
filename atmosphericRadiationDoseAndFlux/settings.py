import os
from functools import wraps
import numpy as np
import pandas as pd
from numba import njit, jit

from .particle import Particle

import ParticleRigidityCalculationTools as PRCT

import inspect
from typing import Union, List, Callable, Any, Dict, Tuple, Optional

@njit
def convert_array_float_numba(value: np.ndarray) -> np.ndarray:
    """
    Numba-optimized version of array conversion.
    
    Parameters
    ----------
    value : np.ndarray
        Input array to ensure as float array
        
    Returns
    -------
    np.ndarray
        Array converted to float64 type
    """
    return np.array(value, dtype=np.float64)

def convertListOrFloatToArray(value) -> np.ndarray:
    """
    Convert various input types to numpy arrays.
    
    Parameters
    ----------
    value : Union[float, int, list, np.ndarray, pd.Series]
        Input value to convert to numpy array
        
    Returns
    -------
    np.ndarray
        The input value converted to a numpy array
    """
    if isinstance(value, np.ndarray):
        outputValue = value
    elif isinstance(value, (float, int)):
        outputValue = np.array([value])
    else:
        outputValue = np.array(value)
    
    return outputValue

def runAcrossBothParticleTypes(args, kwargs, doseFluxCalcFunc):
    """
    Run a calculation for both proton and alpha particles and sum the results.
    
    Parameters
    ----------
    args : tuple
        Positional arguments to pass to the calculation function
    kwargs : dict
        Keyword arguments to pass to the calculation function
    doseFluxCalcFunc : Callable
        The dose/flux calculation function to run
        
    Returns
    -------
    pd.DataFrame
        Combined results from both particle calculations
    """
    doseFluxOutputList = []

    for particleName in ["proton", "alpha"]:
        newKwargs = kwargs.copy()  # Create a copy to avoid modifying the original
        newKwargs["particleName"] = particleName

        outputDF = doseFluxCalcFunc(*args, **newKwargs)
        doseFluxOutputList.append(outputDF)

    # Sum the results from both particle types
    totalOutputDoseFlux = sum(doseFluxOutputList)
    totalOutputDoseFlux["altitude (km)"] = outputDF["altitude (km)"]

    return totalOutputDoseFlux

def allowCalculationForTotalOfParticles(doseFluxCalcFunc):
    """
    Decorator to allow dose/flux calculations for combined particle types.
    
    This decorator enables a dose or flux calculation function to handle
    the special case where particleName='both', which calculates and combines
    results for both protons and alpha particles.
    
    Parameters
    ----------
    doseFluxCalcFunc : Callable
        The dose/flux calculation function to decorate
        
    Returns
    -------
    Callable
        Decorated function that can handle 'both' particle type
    """
    @wraps(doseFluxCalcFunc)
    def totalDoseFlux(*args, **kwargs):
        try:
            if "particleName" in kwargs and kwargs["particleName"] == "both":
                totalOutputDoseFlux = runAcrossBothParticleTypes(args, kwargs, doseFluxCalcFunc)
            else:
                totalOutputDoseFlux = doseFluxCalcFunc(*args, **kwargs)
        except IndexError:
            raise ValueError("No or incorrect input arguments supplied!")

        return totalOutputDoseFlux

    return totalDoseFlux

@njit
def apply_cutoff_rigidity_numba(energy_midpoints: np.ndarray, flux_array: np.ndarray, 
                               cutoff_energy: float) -> np.ndarray:
    """
    Apply vertical cutoff rigidity to flux values using Numba optimization.
    
    Parameters
    ----------
    energy_midpoints : np.ndarray
        Array of energy bin midpoints in MeV
    flux_array : np.ndarray
        Array of flux values to modify
    cutoff_energy : float
        Energy cutoff value in MeV
        
    Returns
    -------
    np.ndarray
        Flux array with cutoff applied
    """
    result = np.zeros_like(flux_array)
    for i in range(len(energy_midpoints)):
        if energy_midpoints[i] >= cutoff_energy:
            result[i] = flux_array[i]
    return result

def formatInputVariables(inputEnergyBins, inputFluxesMeV, altitudesInkm, particleName, verticalCutOffRigidity):
    """
    Format and validate input variables for dose/flux calculations.
    
    Parameters
    ----------
    inputEnergyBins : Union[List[float], np.ndarray]
        Energy bin edges in MeV
    inputFluxesMeV : Union[Callable, float, int, list, np.ndarray, pd.Series]
        Differential flux values at bin midpoints (particles/cm²/sr/MeV/s)
    altitudesInkm : Union[float, List[float], np.ndarray]
        Altitudes in kilometers
    particleName : str
        Particle type ('proton' or 'alpha')
    verticalCutOffRigidity : float
        Vertical cutoff rigidity in GV
        
    Returns
    -------
    Tuple[Particle, np.ndarray, np.ndarray, np.ndarray]
        Tuple containing:
        - Particle object
        - Energy bins array
        - Altitudes array
        - Flux values array with cutoff applied
        
    Raises
    ------
    Exception
        If flux values don't match the number of energy bins
    """
    # Create Particle object
    particleForCalculations = Particle(particleName)

    # Convert rigidity to energy for cutoff
    verticalCutOffRigidityInMeV = PRCT.convertParticleRigidityToEnergy(verticalCutOffRigidity,
                                  particleMassAU=particleForCalculations.atomicMass,
                                  particleChargeAU=particleForCalculations.atomicCharge)

    # Convert inputs to numpy arrays
    inputEnergyBinsArray = convertListOrFloatToArray(inputEnergyBins)
    altitudesInkmArray = convertListOrFloatToArray(altitudesInkm)

    # Calculate energy bin midpoints
    middleOfEnergyBinsArray = (inputEnergyBinsArray[1:] + inputEnergyBinsArray[:-1])/2

    # Process flux input based on type
    if inspect.isfunction(inputFluxesMeV):
        inputFluxesArrayMeV_noVcutOff = np.array(list(map(inputFluxesMeV, middleOfEnergyBinsArray)))
    elif isinstance(inputFluxesMeV, (int, float, list, np.ndarray, pd.Series)):
        inputFluxesArrayMeV_noVcutOff = convertListOrFloatToArray(inputFluxesMeV)
    else:
        raise Exception("inputFluxesMeV not specified as a valid type!")

    # Validate array lengths
    if len(inputEnergyBinsArray) != len(inputFluxesArrayMeV_noVcutOff) + 1:
        raise Exception("Number of bins does not match number of flux values!")

    # Apply cutoff rigidity using Numba-optimized function if arrays are compatible
    if verticalCutOffRigidity > 0:
        # For Numba optimization, ensure arrays are float64
        middleOfEnergyBinsArray_np = np.array(middleOfEnergyBinsArray, dtype=np.float64)
        inputFluxesArrayMeV_noVcutOff_np = np.array(inputFluxesArrayMeV_noVcutOff, dtype=np.float64)
        cutoff_energy = float(verticalCutOffRigidityInMeV)
        
        inputFluxesArray = apply_cutoff_rigidity_numba(
            middleOfEnergyBinsArray_np, 
            inputFluxesArrayMeV_noVcutOff_np,
            cutoff_energy
        )
    else:
        # If no cutoff, just use the original values
        inputFluxesArray = inputFluxesArrayMeV_noVcutOff

    return particleForCalculations, inputEnergyBinsArray, altitudesInkmArray, inputFluxesArray

# Constants and paths
directory_of_this_file = os.path.dirname(os.path.realpath(__file__))
homeDirectory = directory_of_this_file
dataFileDirectory = f"{directory_of_this_file}/data/"
