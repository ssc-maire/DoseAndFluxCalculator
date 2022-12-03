import numpy as np
import pandas as pd

import particle
import particleResponse
import settings
import units

import ParticleRigidityCalculationTools as PRCT

import inspect

from typing import Callable

@settings.allowCalculationForTotalOfParticles
def calculate_from_energy_spec(
                            inputEnergyDistributionFunctionMeV:Callable, 
                            altitudesInkm:list, 
                            particleName="both", 
                            inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1),
                            verticalCutOffRigidity = 0.0):

    energyBinMidPoints = (inputEnergyBins[1:] + inputEnergyBins[:-1])/2
    inputFluxesMeV = list(map(inputEnergyDistributionFunctionMeV,energyBinMidPoints))

    outputDF = calculate_from_energy_spec_array(
                            inputEnergyBins,
                            inputFluxesMeV, 
                            altitudesInkm, 
                            particleName = particleName,
                            verticalCutOffRigidity = verticalCutOffRigidity)

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_rigidity_spec(
                            inputRigidityDistributionFunctionGV:Callable, 
                            altitudesInkm:list, 
                            particleName="both", 
                            inputRigidityBins = None,
                            verticalCutOffRigidity = 0.0):

    if inputRigidityBins == None:
        particleForCalculations = particle.Particle(particleName)

        inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1)
        inputRigidityBins = np.array(PRCT.convertParticleEnergyToRigidity(inputEnergyBins,
                                    particleMassAU = particleForCalculations.atomicMass,
                                    particleChargeAU = particleForCalculations.atomicCharge))

    rigidityBinMidPoints = (inputRigidityBins[1:] + inputRigidityBins[:-1])/2
    inputFluxesGV = list(map(inputRigidityDistributionFunctionGV,
                              rigidityBinMidPoints))

    outputDF = calculate_from_rigidity_spec_array(
                            inputRigidityBins,
                            inputFluxesGV, 
                            altitudesInkm, 
                            particleName = particleName,
                            verticalCutOffRigidity = verticalCutOffRigidity)

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_energy_spec_array(
                            inputEnergyBins:list,
                            inputFluxesMeV:list, 
                            altitudesInkm:list, 
                            particleName="both", 
                            #inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1),
                            verticalCutOffRigidity = 0.0):
    
    # flux units: particles/cm2/sr/MeV/s 
    # altitudesInkm: altitude:km
    # energy units: MeV
    # verticalCutOffRigidity units: GV

    particleForCalculations, inputEnergyBinsArray, altitudesInkmArray, inputFluxesArray = settings.formatInputVariables(inputEnergyBins, 
                                                                                                               inputFluxesMeV, 
                                                                                                               altitudesInkm, 
                                                                                                               particleName, 
                                                                                                               verticalCutOffRigidity)

    ############################################################

    inputEnergyDifferences = inputEnergyBinsArray[1:] - inputEnergyBinsArray[:-1]
    inputFluxesIntegrated = inputFluxesArray * inputEnergyDifferences * np.pi # units of particles / cm2 / s

    ############################################################

    outputDF = pd.DataFrame({
        "altitude (km)":altitudesInkmArray,
    })

    for doseType in particleResponse.fullListOfDoseResponseTypes:
        doseResponseForParticle = particleResponse.fullDoseResponseDict[doseType](particleForCalculations,doseType)
    
        coordinateDosesList = []
        for altitudeInkm in altitudesInkmArray:
            altitude = units.Distance(altitudeInkm * 1000.0)
            totalDose = doseResponseForParticle.calculateDose(altitude, inputEnergyBinsArray, inputFluxesIntegrated)
            coordinateDosesList.append(totalDose)

        outputDF[doseType] = coordinateDosesList

    outputDF["SEU"] = outputDF["tn2"] * 1e-13
    outputDF["SEL"] = outputDF["tn2"] * 1e-8

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_rigidity_spec_array(
                            inputRigidityBins:list,
                            inputFluxesGV:list, 
                            altitudesInkm:list, 
                            particleName="both",
                            verticalCutOffRigidity = 0.0):

    particleForCalculations = particle.Particle(particleName)

    inputEnergyBins = np.array(PRCT.convertParticleRigidityToEnergy(inputRigidityBins,
                                particleMassAU = particleForCalculations.atomicMass,
                                particleChargeAU = particleForCalculations.atomicCharge))

    rigidityMidPoints = (inputRigidityBins[1:] + inputRigidityBins[:-1])/2

    if inspect.isfunction(inputFluxesGV):
        inputFluxesGVarray = np.array(list(map(inputFluxesGV,rigidityMidPoints)))

    elif type(inputFluxesGV) in [int, float, list, np.ndarray]:
        inputFluxesGVarray = settings.convertListOrFloatToArray(inputFluxesGV)
    
    else:
        raise Exception("inputFluxesGV not specified as a valid type!")

    inputFluxesMeV = PRCT.convertParticleRigiditySpecToEnergySpec(rigidityMidPoints,
                                                                             inputFluxesGVarray, 
                                                                             particleMassAU=particleForCalculations.atomicMass,
                                                                             particleChargeAU=particleForCalculations.atomicCharge)["Energy distribution values"]

    outputDoseFluxRates = calculate_from_energy_spec_array(
                            inputEnergyBins,
                            inputFluxesMeV, 
                            altitudesInkm, 
                            particleName = particleName,
                            verticalCutOffRigidity = verticalCutOffRigidity)

    return outputDoseFluxRates

