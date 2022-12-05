
import numpy as np
from .particle import Particle
from .responseFileParameters import ResponseFileParameters
from .settings import dataFileDirectory
from .units import Distance

import importlib_resources
import pkg_resources

class ParticleResponse():

    def __init__(self, particle:Particle, doseTypeName):

        self.particle = particle
        self.doseType = doseTypeName

        pathToRelevantResponseFile = self.getPathToResponseFile()

        self.particleResponseArray = np.genfromtxt(pathToRelevantResponseFile)

    def calculateDose(self, altitude:Distance, inputEnergyBins, inputFluxesIntegrated):

        ResponseParameters = ResponseFileParameters(altitude, inputEnergyBins, inputFluxesIntegrated, self.particle)

        altitudeLayerIndex = ResponseParameters.altitudeLayerIndex
        altIndexAbove = ResponseParameters.altIndexAbove
        f1 = ResponseParameters.f1
        weightedFluxes = ResponseParameters.weightedFluxes

        outputDose = 0.0
        for energyIndex in range(0,50):

            #firstDoseResponseTerm = self.particleResponseArray[altitudeLayerIndex,energyIndex]*f1
            #secondDoseResponseTerm = self.particleResponseArray[altIndexAbove,energyIndex]*(f1-1)

            (firstDoseResponseTerm, secondDoseResponseTerm) = self.getDoseResponseTerms(altitudeLayerIndex, altIndexAbove, energyIndex, f1)

            outputDose = outputDose + (weightedFluxes[energyIndex] * (firstDoseResponseTerm - secondDoseResponseTerm))

        return outputDose

class DoseRateResponse(ParticleResponse):

    def getPathToResponseFile(self):

        #return f"{dataFileDirectory}{self.particle.particleName}/{self.doseType}.rpf"

        #return importlib_resources.files(f"atmosphericRadiationDoseAndFlux.data.{self.particle.particleName}").joinpath(f"{self.doseType}.rpf")
        return pkg_resources.resource_stream(__name__,f"data/{self.particle.particleName}/{self.doseType}.rpf")

    def getDoseResponseTerms(self, altitudeLayerIndex, altIndexAbove, energyIndex, f1):

        firstDoseResponseTerm = self.particleResponseArray[altitudeLayerIndex,energyIndex]*f1
        secondDoseResponseTerm = self.particleResponseArray[altIndexAbove,energyIndex]*(f1-1)

        return (firstDoseResponseTerm, secondDoseResponseTerm)


class NeutronFluxResponse(ParticleResponse):

    energyIndexTranslationDict = {
        "tn1":0,
        "tn2":50,
        "tn3":100,
    }

    def getPathToResponseFile(self):

        #return f"{dataFileDirectory}{self.particle.particleName}/neutron.rpf"

        #return importlib_resources.files(f"atmosphericRadiationDoseAndFlux.data.{self.particle.particleName}").joinpath(f"neutron.rpf")
        return pkg_resources.resource_stream(__name__,f"data/{self.particle.particleName}/neutron.rpf")

    def getDoseResponseTerms(self, altitudeLayerIndex, altIndexAbove, energyIndex, f1):

        translationIndex = self.energyIndexTranslationDict[self.doseType]

        firstDoseResponseTerm = self.particleResponseArray[altitudeLayerIndex,energyIndex + translationIndex]*f1
        secondDoseResponseTerm = self.particleResponseArray[altIndexAbove,energyIndex + translationIndex]*(f1-1)

        return (firstDoseResponseTerm, secondDoseResponseTerm)

humanDoseTypes = ["edose","adose","dosee"]
humanDoseResponseDict = dict.fromkeys(humanDoseTypes, DoseRateResponse)

neutronDoseTypes = ["tn1","tn2","tn3"]
neutronDoseResponseDict = dict.fromkeys(neutronDoseTypes, NeutronFluxResponse)

fullListOfDoseResponseTypes = humanDoseTypes + neutronDoseTypes
fullDoseResponseDict = {**humanDoseResponseDict, **neutronDoseResponseDict}