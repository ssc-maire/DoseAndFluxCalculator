import os

import numpy as np
import pandas as pd

from .particle import Particle

import ParticleRigidityCalculationTools as PRCT

import inspect

def convertListOrFloatToArray(value)->np.array:
    if isinstance(value,np.ndarray):
        outputValue = value
    elif isinstance(value,float) or isinstance(value,int):
        outputValue = np.array([value])
    else:
        outputValue = np.array(value)
    
    return outputValue

def runAcrossBothParticleTypes(args, kwargs, doseFluxCalcFunc):
        doseFluxOutputList = []

        for particleName in ["proton", "alpha"]:

            newKwargs = kwargs
            newKwargs["particleName"] = particleName

            outputDF = doseFluxCalcFunc(*args, **newKwargs)

            doseFluxOutputList.append(outputDF)

        totalOutputDoseFlux = sum(doseFluxOutputList)
        totalOutputDoseFlux["altitude (km)"] = outputDF["altitude (km)"]

        return totalOutputDoseFlux

def allowCalculationForTotalOfParticles(doseFluxCalcFunc):

    def totalDoseFlux(*args, **kwargs):

        try:
            if "particleName" in kwargs.keys():
                if (kwargs["particleName"] == "both"):

                    totalOutputDoseFlux = runAcrossBothParticleTypes(args, kwargs, doseFluxCalcFunc)

                else:
                    totalOutputDoseFlux = doseFluxCalcFunc(*args, **kwargs)

            else:
                totalOutputDoseFlux = doseFluxCalcFunc(*args, **kwargs)

        except IndexError:
            raise ValueError("No or incorrect input arguments supplied!")

        return totalOutputDoseFlux

    return totalDoseFlux

def formatInputVariables(inputEnergyBins, inputFluxesMeV, altitudesInkm, particleName, verticalCutOffRigidity):
    particleForCalculations = Particle(particleName)

    verticalCutOffRigidityInMeV = PRCT.convertParticleRigidityToEnergy(verticalCutOffRigidity,
                                    particleMassAU = particleForCalculations.atomicMass,
                                    particleChargeAU = particleForCalculations.atomicCharge).iloc[0]

    inputEnergyBinsArray = convertListOrFloatToArray(inputEnergyBins)
    altitudesInkmArray = convertListOrFloatToArray(altitudesInkm)

    middleOfEnergyBinsArray = (inputEnergyBinsArray[1:] - inputEnergyBinsArray[:-1])/2

    if inspect.isfunction(inputFluxesMeV):
        inputFluxesArrayMeV_noVcutOff = np.array(list(map(inputFluxesMeV,middleOfEnergyBinsArray)))

    elif type(inputFluxesMeV) in [int, float, list, np.ndarray, pd.Series]:
        inputFluxesArrayMeV_noVcutOff = convertListOrFloatToArray(inputFluxesMeV)
    
    else:
        raise Exception("inputFluxesMeV not specified as a valid type!")

    if particleName == "alpha":
        print("Warning: currently using an alpha particle as input actually calculates the contribution from alpha + all simulated heavier ions, rather than just alpha particles!")

    if len(inputEnergyBinsArray) != len(inputFluxesArrayMeV_noVcutOff) + 1:
        raise Exception("Number of bins does not match number of flux values!")

    inputFluxesArray = inputFluxesArrayMeV_noVcutOff * (middleOfEnergyBinsArray >= verticalCutOffRigidityInMeV)

    return particleForCalculations,inputEnergyBinsArray,altitudesInkmArray,inputFluxesArray

directory_of_this_file = os.path.dirname(os.path.realpath(__file__))
homeDirectory = directory_of_this_file
dataFileDirectory = f"{directory_of_this_file}/data/"
