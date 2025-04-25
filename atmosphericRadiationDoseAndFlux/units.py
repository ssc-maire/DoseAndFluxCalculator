from numba import jit

class Distance:
    """
    Class representing distances with automatic unit conversions.
    
    This class handles distance conversions between meters, kilometers,
    feet, and kilofeet for atmospheric radiation calculations.
    
    Attributes
    ----------
    meters : float
        Distance in meters
    km : float
        Distance in kilometers
    ft : float
        Distance in feet
    kft : float
        Distance in kilofeet
        
    Class Attributes
    ----------------
    oneMeterInft : float
        Conversion factor from meters to feet (3.28084)
    """

    oneMeterInft = 3.28084

    def __init__(self, distanceInMeters: float):
        """
        Initialize a Distance object.
        
        Parameters
        ----------
        distanceInMeters : float
            Distance in meters
        """
        self.meters = distanceInMeters
        self.km = distanceInMeters/1000.0

        self.ft = self.meters * self.oneMeterInft
        self.kft = self.ft / 1000.0