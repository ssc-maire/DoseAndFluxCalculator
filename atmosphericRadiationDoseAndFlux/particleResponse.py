import numpy as np
from numba import njit, jit
from typing import Tuple, Dict, List, Union, Optional

from .particle import Particle
from .responseFileParameters import ResponseFileParameters
from .settings import dataFileDirectory
from .units import Distance

import importlib_resources
import pkg_resources

# Numba-optimized calculation of dose response
@njit
def calculate_dose_response(weighted_fluxes: np.ndarray, 
                           response_values: np.ndarray, 
                           alt_index: int, 
                           alt_index_above: int, 
                           f1: float) -> float:
    """
    Calculate dose response using Numba optimization.
    
    Parameters
    ----------
    weighted_fluxes : np.ndarray
        Array of weighted flux values
    response_values : np.ndarray
        Array of response values
    alt_index : int
        Index for altitude layer
    alt_index_above : int
        Index for altitude layer above
    f1 : float
        Interpolation factor
        
    Returns
    -------
    float
        Calculated dose value
    """
    output_dose = 0.0
    for energy_idx in range(50):
        first_term = response_values[alt_index, energy_idx] * f1
        second_term = response_values[alt_index_above, energy_idx] * (f1 - 1)
        output_dose += weighted_fluxes[energy_idx] * (first_term - second_term)
    return output_dose

# Numba-optimized calculation of neutron flux response
@njit
def calculate_neutron_response(weighted_fluxes: np.ndarray, 
                              response_values: np.ndarray, 
                              alt_index: int, 
                              alt_index_above: int, 
                              f1: float, 
                              translation_index: int) -> float:
    """
    Calculate neutron flux response using Numba optimization.
    
    Parameters
    ----------
    weighted_fluxes : np.ndarray
        Array of weighted flux values
    response_values : np.ndarray
        Array of response values
    alt_index : int
        Index for altitude layer
    alt_index_above : int
        Index for altitude layer above
    f1 : float
        Interpolation factor
    translation_index : int
        Index translation for different neutron energy ranges
        
    Returns
    -------
    float
        Calculated neutron flux value
    """
    output_dose = 0.0
    for energy_idx in range(50):
        first_term = response_values[alt_index, energy_idx + translation_index] * f1
        second_term = response_values[alt_index_above, energy_idx + translation_index] * (f1 - 1)
        output_dose += weighted_fluxes[energy_idx] * (first_term - second_term)
    return output_dose

class ParticleResponse:
    """
    Base class for particle response calculations.
    
    This class provides the foundation for calculating particle response
    in the atmosphere, including dose rates and neutron fluxes.
    
    Attributes
    ----------
    particle : Particle
        Particle object (proton or alpha)
    doseType : str
        Type of dose/flux to calculate
    particleResponseArray : np.ndarray
        Array of response values loaded from response file
        
    Methods
    -------
    __init__(particle, doseTypeName)
        Initialize the particle response object
    getPathToResponseFile()
        Get the path to the response file (implemented by subclasses)
    calculateDose(altitude, inputEnergyBins, inputFluxesIntegrated)
        Calculate dose or flux at the specified altitude
    getDoseResponseTerms(altitudeLayerIndex, altIndexAbove, energyIndex, f1)
        Get dose response terms for calculation (implemented by subclasses)
    """

    def __init__(self, particle: Particle, doseTypeName: str):
        """
        Initialize the particle response calculator.
        
        Parameters
        ----------
        particle : Particle
            Particle object (proton or alpha)
        doseTypeName : str
            Type of dose/flux to calculate (e.g., 'edose', 'tn1')
        """
        self.particle = particle
        self.doseType = doseTypeName

        pathToRelevantResponseFile = self.getPathToResponseFile()
        self.particleResponseArray = np.genfromtxt(pathToRelevantResponseFile)

    def calculateDose(self, altitude: Distance, inputEnergyBins: np.ndarray, 
                    inputFluxesIntegrated: np.ndarray) -> float:
        """
        Calculate dose or flux at the specified altitude.
        
        Parameters
        ----------
        altitude : Distance
            Altitude object with distance in meters
        inputEnergyBins : np.ndarray
            Energy bin edges in MeV
        inputFluxesIntegrated : np.ndarray
            Integrated flux values
            
        Returns
        -------
        float
            Calculated dose or flux value
        """
        # Create response parameters object
        ResponseParameters = ResponseFileParameters(altitude, inputEnergyBins, 
                                                  inputFluxesIntegrated, self.particle)

        # Get altitude interpolation parameters
        altitudeLayerIndex = ResponseParameters.altitudeLayerIndex
        altIndexAbove = ResponseParameters.altIndexAbove
        f1 = ResponseParameters.f1
        weightedFluxes = ResponseParameters.weightedFluxes

        # Use the appropriate Numba-optimized function based on type
        if isinstance(self, DoseRateResponse):
            return calculate_dose_response(
                np.array(weightedFluxes, dtype=np.float64),
                self.particleResponseArray,
                altitudeLayerIndex,
                altIndexAbove,
                f1
            )
        elif isinstance(self, NeutronFluxResponse):
            translation_index = self.energyIndexTranslationDict[self.doseType]
            return calculate_neutron_response(
                np.array(weightedFluxes, dtype=np.float64),
                self.particleResponseArray,
                altitudeLayerIndex,
                altIndexAbove,
                f1,
                translation_index
            )
        else:
            # Fallback to non-optimized calculation if needed
            outputDose = 0.0
            for energyIndex in range(0, 50):
                (firstDoseResponseTerm, secondDoseResponseTerm) = self.getDoseResponseTerms(
                    altitudeLayerIndex, altIndexAbove, energyIndex, f1)
                outputDose += (weightedFluxes[energyIndex] * 
                              (firstDoseResponseTerm - secondDoseResponseTerm))
            return outputDose

class DoseRateResponse(ParticleResponse):
    """
    Class for calculating dose rates (edose, adose, dosee).
    
    This class handles human-relevant dose responses such as effective dose,
    ambient dose equivalent, etc.
    
    Methods
    -------
    getPathToResponseFile()
        Get the path to the dose response file
    getDoseResponseTerms(altitudeLayerIndex, altIndexAbove, energyIndex, f1)
        Get dose response terms for calculation
    """

    def getPathToResponseFile(self) -> str:
        """
        Get the path to the dose response file.
        
        Returns
        -------
        str
            Path to the response file
        """
        return pkg_resources.resource_stream(__name__, 
                                            f"data/{self.particle.particleName}/{self.doseType}.rpf")

    def getDoseResponseTerms(self, altitudeLayerIndex: int, altIndexAbove: int, 
                           energyIndex: int, f1: float) -> Tuple[float, float]:
        """
        Get dose response terms for calculation.
        
        Parameters
        ----------
        altitudeLayerIndex : int
            Index for altitude layer
        altIndexAbove : int
            Index for altitude layer above
        energyIndex : int
            Index for energy bin
        f1 : float
            Interpolation factor
            
        Returns
        -------
        Tuple[float, float]
            First and second dose response terms
        """
        firstDoseResponseTerm = self.particleResponseArray[altitudeLayerIndex, energyIndex] * f1
        secondDoseResponseTerm = self.particleResponseArray[altIndexAbove, energyIndex] * (f1 - 1)

        return (firstDoseResponseTerm, secondDoseResponseTerm)

class NeutronFluxResponse(ParticleResponse):
    """
    Class for calculating neutron fluxes (tn1, tn2, tn3).
    
    This class handles neutron flux responses for different energy ranges.
    
    Attributes
    ----------
    energyIndexTranslationDict : Dict[str, int]
        Dictionary mapping dose type to energy index offsets
        
    Methods
    -------
    getPathToResponseFile()
        Get the path to the neutron response file
    getDoseResponseTerms(altitudeLayerIndex, altIndexAbove, energyIndex, f1)
        Get neutron flux response terms for calculation
    """

    energyIndexTranslationDict = {
        "tn1": 0,    # Thermal neutrons
        "tn2": 50,   # Epithermal neutrons
        "tn3": 100,  # High-energy neutrons
    }

    def getPathToResponseFile(self) -> str:
        """
        Get the path to the neutron response file.
        
        Returns
        -------
        str
            Path to the neutron response file
        """
        return pkg_resources.resource_stream(__name__, 
                                            f"data/{self.particle.particleName}/neutron.rpf")

    def getDoseResponseTerms(self, altitudeLayerIndex: int, altIndexAbove: int, 
                           energyIndex: int, f1: float) -> Tuple[float, float]:
        """
        Get neutron flux response terms for calculation.
        
        Parameters
        ----------
        altitudeLayerIndex : int
            Index for altitude layer
        altIndexAbove : int
            Index for altitude layer above
        energyIndex : int
            Index for energy bin
        f1 : float
            Interpolation factor
            
        Returns
        -------
        Tuple[float, float]
            First and second neutron flux response terms
        """
        translationIndex = self.energyIndexTranslationDict[self.doseType]

        firstDoseResponseTerm = self.particleResponseArray[altitudeLayerIndex, 
                                                         energyIndex + translationIndex] * f1
        secondDoseResponseTerm = self.particleResponseArray[altIndexAbove, 
                                                          energyIndex + translationIndex] * (f1 - 1)

        return (firstDoseResponseTerm, secondDoseResponseTerm)

# Define dose response types
humanDoseTypes = ["edose", "adose", "dosee"]
humanDoseResponseDict = dict.fromkeys(humanDoseTypes, DoseRateResponse)

neutronDoseTypes = ["tn1", "tn2", "tn3"]
neutronDoseResponseDict = dict.fromkeys(neutronDoseTypes, NeutronFluxResponse)

# Combined dose response types and classes
fullListOfDoseResponseTypes = humanDoseTypes + neutronDoseTypes
fullDoseResponseDict = {**humanDoseResponseDict, **neutronDoseResponseDict}