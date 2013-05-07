from powerLoad import yearLoad
from windModel import WindModel
from solarModel import SolarModel

import numpy as np
import datetime as dt

class MicroGrid:

	def __init__(self, latitude, longitude, year=2012, timezone=3):
		"""timezone is hours from utc"""
		self.solarModel = SolarModel(latitude, longitude)
		self.windModel = WindModel()
		self.powerLoad = yearLoad()
		self.year = year
		self.timezone = timezone

	def calculateShortages(self):
		shortages = np.array((len(self.powerLoad),len(self.powerLoad[0])))
		currentDay = 0
		for day in self.powerLoad:
			currentShortages = []
			currentDate = dt.datetime(self.year,1,1) + dt.timedelta(days=currentDay)
			for hour in day:
				solar = self.solarModel.availableWattage(currentDate-dt.timedelta(hours=self.timezone))
				wind = self.windModel.availableWattage(currentDate)
				difference = hour - (solar+wind)
				currentShortages.append(difference)
				currentDate += dt.timedelta(hours=1)
			currentDay += 1
			shortages = np.concatenate([shortages, np.array(currentShortages)])
		return shortages

test = MicroGrid(42,-71)
print test.calculateShortages()