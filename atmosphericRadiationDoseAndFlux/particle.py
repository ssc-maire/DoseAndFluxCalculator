import ParticleRigidityCalculationTools

class Particle():

    def __init__(self, particleName):
        self.particleName = particleName

        if particleName == "proton":
            self.atomicCharge = 1
        elif particleName == "alpha":
            self.atomicCharge = 2
        else:
            raise Exception("Error: currently only protons and alpha particles can be used!")

        self.atomicMass = ParticleRigidityCalculationTools.getAtomicMass(self.atomicCharge)