import numpy as np
import pandas as pd

import particle
import particleResponse
import settings
import units

import ParticleRigidityCalculationTools as PRCT

import inspect

@settings.allowCalculationForTotalOfParticles
def calculate_from_energy_spectrum(
                            inputEnergyBins,
                            inputFluxesMeV, 
                            altitudesInkm, 
                            particleName, 
                            #inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1),
                            verticalCutOffRigidity = 0.0):
    
    # flux units: particles/cm2/sr/MeV/s 
    # altitudesInkm: altitude:km
    # energy units: MeV
    # verticalCutOffRigidity units: GV

    particleForCalculations = particle.Particle(particleName)

    verticalCutOffRigidityInMeV = PRCT.convertParticleRigidityToEnergy(verticalCutOffRigidity,
                                    particleMassAU = particleForCalculations.atomicMass,
                                    particleChargeAU = particleForCalculations.atomicCharge).iloc[0]

    inputEnergyBinsArray = settings.convertListOrFloatToArray(inputEnergyBins)
    altitudesInkmArray = settings.convertListOrFloatToArray(altitudesInkm)

    middleOfEnergyBinsArray = (inputEnergyBinsArray[1:] - inputEnergyBinsArray[:-1])/2

    if inspect.isfunction(inputFluxesMeV):
        inputFluxesArrayMeV_noVcutOff = np.array(list(map(inputFluxesMeV,middleOfEnergyBinsArray)))

    elif type(inputFluxesMeV) in [int, float, list, np.ndarray, pd.Series]:
        inputFluxesArrayMeV_noVcutOff = settings.convertListOrFloatToArray(inputFluxesMeV)
    
    else:
        raise Exception("inputFluxesMeV not specified as a valid type!")

    if particleName == "alpha":
        print("Warning: currently using an alpha particle as input actually calculates the contribution from alpha + all simulated heavier ions, rather than just alpha particles!")

    if len(inputEnergyBinsArray) != len(inputFluxesArrayMeV_noVcutOff) + 1:
        raise Exception("Number of bins does not match number of flux values!")

    inputFluxesArray = inputFluxesArrayMeV_noVcutOff * (middleOfEnergyBinsArray >= verticalCutOffRigidityInMeV)

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

    return outputDF

@settings.allowCalculationForTotalOfParticles
def calculate_from_rigidity_spectrum(inputFluxesGV, 
                            altitudesInkm, 
                            particleName, 
                            inputRigidityBins = None,
                            verticalCutOffRigidity = 0.0):

    particleForCalculations = particle.Particle(particleName)

    if inputRigidityBins == None:
        inputEnergyBins = 10**(0.1*(np.array(range(1,52))-1)+1)
        inputRigidityBins = np.array(PRCT.convertParticleEnergyToRigidity(inputEnergyBins,
                                    particleMassAU = particleForCalculations.atomicMass,
                                    particleChargeAU = particleForCalculations.atomicCharge))
    else:
        inputEnergyBins = np.array(PRCT.convertParticleRigidityToEnergy(inputRigidityBins,
                                    particleMassAU = particleForCalculations.atomicMass,
                                    particleChargeAU = particleForCalculations.atomicCharge))

    rigidityMidPoints = (inputRigidityBins[1:] - inputRigidityBins[:-1])/2

    if inspect.isfunction(inputFluxesGV):
        inputFluxesGVarray = np.array(list(map(inputFluxesGV,rigidityMidPoints)))

    elif type(inputFluxesGV) in [int, float, list, np.ndarray]:
        inputFluxesGVarray = settings.convertListOrFloatToArray(inputFluxesGV)
    
    else:
        raise Exception("inputFluxesGV not specified as a valid type!")

    inputFluxesMeV = PRCT.convertParticleRigiditySpecToEnergySpec(rigidityMidPoints,
                                                                             inputFluxesGVarray, 
                                                                             particleMassAU=particleForCalculations.atomicMass,
                                                                             particleChargeAU=particleForCalculations.atomicCharge)

    outputDoseFluxRates = calculate_from_energy_spectrum(inputFluxesMeV, 
                            altitudesInkm, 
                            particleName, 
                            inputEnergyBins = inputEnergyBins,
                            verticalCutOffRigidity = verticalCutOffRigidity)

    return outputDoseFluxRates

if __name__ == "__main__":
    print(calculate_from_rigidity_spectrum(lambda x:x**-5,5.0,"proton",verticalCutOffRigidity=2.0))
    print(calculate_from_rigidity_spectrum(lambda x:x**-5,5.0,"alpha",verticalCutOffRigidity=2.0))
    print(calculate_from_rigidity_spectrum(lambda x:x**-5,5.0,"both",verticalCutOffRigidity=2.0))

