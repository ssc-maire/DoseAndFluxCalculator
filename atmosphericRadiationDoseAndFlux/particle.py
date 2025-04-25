import ParticleRigidityCalculationTools
from typing import Dict, Union, List, Optional

class Particle:
    """
    Class representing a cosmic ray particle with defined properties.
    
    This class encapsulates properties of cosmic ray particles used in atmospheric
    radiation dose and flux calculations. Currently supports protons and alpha particles.
    
    Attributes
    ----------
    particleName : str
        Name of the particle type ('proton' or 'alpha')
    atomicCharge : int
        Atomic charge of the particle in elementary charge units
    atomicMass : float
        Atomic mass of the particle in atomic mass units (u)
    
    Methods
    -------
    __init__(particleName)
        Initialize a Particle object with the specified name
        
    Notes
    -----
    Currently only protons and alpha particles are supported. Alpha particles are
    used to represent alpha particles plus all simulated heavier ions.
    """

    def __init__(self, particleName: str):
        """
        Initialize a particle object.
        
        Parameters
        ----------
        particleName : str
            Name of the particle ('proton' or 'alpha')
            
        Raises
        ------
        Exception
            If the particle name is not 'proton' or 'alpha'
        """
        self.particleName = particleName

        if particleName == "proton":
            self.atomicCharge = 1
        elif particleName == "alpha":
            self.atomicCharge = 2
        else:
            raise Exception("Error: currently only protons and alpha particles can be used!")

        # Get atomic mass using the PRCT module
        self.atomicMass = ParticleRigidityCalculationTools.getAtomicMass(self.atomicCharge)