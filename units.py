class Distance():

    oneMeterInft = 3.28084

    def __init__(self, distanceInMeters:float):

        self.meters = distanceInMeters
        self.km = distanceInMeters/1000.0

        self.ft = self.meters * self.oneMeterInft
        self.kft = self.ft / 1000.0