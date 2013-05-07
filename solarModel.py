import Pysolar as solar
import numpy as np
import datetime

class SolarModel:
	def __init__(self, latitude, longitude, size = 10, efficiency=0.17):
		"""latitude is defined as positive in the north, negative in the south,
		longitude is defined as negative west of greenwich, positive east.
		Use degrees for both"""
		self.latitude = latitude
		self.longitude = longitude
		self.size = size #m^2
		self.efficiency = efficiency

	def availableWattage(self, datetime):
		"""expects utc datetime"""
		altitude = solar.GetAltitude(self.latitude, self.longitude, datetime)
		radiation = solar.radiation.GetRadiationDirect(datetime, altitude)
		return radiation*self.size*self.efficiency