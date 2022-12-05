import numpy as np
from .particle import Particle
from .units import Distance


class ResponseFileParameters():

    def __init__(self, altitude:Distance, energyBins:np.array, fluxes:np.array, particle:Particle):

        self.particle = particle
        self.getAltitudeResponseLayer(altitude)
        self.getResponseParametersFromSpectrum(energyBins, fluxes)

    def getAltitudeResponseLayer(self, altitude:Distance):

        self.altitude = altitude

        if altitude.km > 100.0:
            raise Exception("altitude.km has to be less than 100 km!")

        if (altitude.km < 0.025):
            layer = 1
            f1 = 1
        elif (altitude.km < 1.025):
            layer = int((int(altitude.meters)-25)/50) + 1
            hr = (int(altitude.meters)-25) % 50
            f1 = 1. - hr/50.
        elif (altitude.km < 1.15):
            layer = 21
            f1 = 1. - (int(altitude.meters)-1025)/1025.
        elif (altitude.km < 5.05):
            layer = int((int(altitude.meters)-1150)/100) + 22
            hr = (int(altitude.meters)-1150) % 100
            f1 = 1. - hr/100.
        elif (altitude.km < 5.3):
            layer = 61
            f1 = 1. - (altitude.km-5.05)/0.25
        elif (altitude.km < 15.1):
            layer = int((int(altitude.meters)-5300)/200) + 62
            hr = (int(altitude.meters)-5300) % 200
            f1 = 1. - hr/200.
        elif (altitude.km < 16.5):
            layer = 111
            f1 = 1. - (altitude.km-15.1)/1.4
        elif (altitude.km < 38.5):
            layer = int((int(altitude.meters)-16500)/1000) + 112
            hr = (int(altitude.meters)-16500) % 1000
            f1 = 1. - hr/1000.
        elif (altitude.km < 40.5):
            layer = 134
            f1 = 1. - (altitude.km-38.5)/2.
        elif (altitude.km < 62.5):
            layer = 135
            f1 = 1. - (altitude.km-40.5)/22.
        elif (altitude.km < 97.5):
            layer = 136
            f1 = 1. - (altitude.km-40.5)/57.5
        else:
            layer = 137
            f1 = 1

        self.altitudeLayerIndex = layer - 1 #needs 1 taken away from it, due to the conversion of the original algorithm for this from Fortran
        self.f1 = f1

    def getResponseParametersFromSpectrum(self, energyBins:np.array, fluxes:np.array):

        rigidityCutoff = 0.0 #rigidity cutoff is not relevant for these calculations

        altIndexAbove = self.altitudeLayerIndex + 1
        if (altIndexAbove > 137): 
            altIndexAbove = 137

        cutoffEnergy = 1000.*(np.sqrt(((rigidityCutoff**2)/self.particle.atomicMass) + (0.938**2))-0.938)

        energyIndex = 0
        while (cutoffEnergy > energyBins[energyIndex]):
            energyIndex = energyIndex + 1
        
        weightedFluxes = fluxes
        if energyIndex != 0:
            energyIndex = energyIndex-1
            
            ed = energyBins[energyIndex+1] - cutoffEnergy
            et = energyBins[energyIndex+1] - energyBins[energyIndex]
            
            fe = ed/et
            weightedFluxes[energyIndex]  = weightedFluxes[energyIndex]*fe

        self.altIndexAbove = altIndexAbove
        self.weightedFluxes = weightedFluxes