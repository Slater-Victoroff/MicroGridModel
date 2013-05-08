from powerLoad import yearLoad
from windModel import WindModel
from solarModel import SolarModel
from batteryModel import Battery

import numpy as np
import datetime as dt
import matplotlib.pyplot as plot

from scipy.optimize import minimize

import cProfile

class MicroGrid:

	def __init__(self, latitude, longitude, year=2012, timezone=3, solarPrice=lambda x: 500*x,windPrice=lambda x: 500*x**2,
		batteryPrice = lambda capacity:capacity*0.283,shortageCostFunction = lambda params: params[1]*0.95, 
		humanCost = lambda x: 1.5*x, humanEmergencyMax = 100, humanEfficiency = 0.7, inverterEfficiency=0.85,
		chargeConverterEfficiency=0.9, solarEfficiency=0.17, windEfficiency=0.22, batterySize = np.array([6,220]), 
		solarSize = 50, windSize = 13, batterySocBounds = np.array([0.2,0.9]), meanPower = 2000):
		"""timezone is hours from utc, latitude and longitude should both be defined in degrees with standard signs.

		the cost functions should ultimately be more sensible, but for now these should do. Currently the solar
		cost function assumes that square meters of array is all that cost depends on. Wind Priced is based on the 
		rotor diameter, but should probably be based on cut in speed and whatnot later on. 

		Battery cost is probably about what it needs to be since battery cost by capacity is pretty linear.
		shortageCostFunction right now is assuming that is gets a list/vector with the first entry being the 
		number of timesteps(hours currently) that the shortage has been going on, and the size of the energy deficiency.
		
		HumanCost right now is supposed to take in the number of people-hours and return price of a person. Eventually this
		probably makes more sense as a step function with each added person taking a specific cost and outputting a
		specific power with some noise, but that can come later. humanEmergencyMax is the maximum power that a human
		can provice in an emergency in watts. humanEfficiency is the efficiency of whatever generator is being cranked
		by the person.

		The efficiencies are self explanatory, and the sizes are as defined in the specific modules with the battery sizing
		being in the form of [volts, amp hours]

		meanPower is currently the center of the powerload distribution, will include other powerLoad variables later,
		but they seem like overkill for right now."""
		self.solarModel = SolarModel(latitude, longitude, solarSize, solarEfficiency)
		self.windModel = WindModel(windSize, windEfficiency)
		self.powerLoad = yearLoad(meanPower)
		self.batteryModel = Battery(batterySize[0], batterySize[1], batterySocBounds)
		self.price = solarPrice(solarSize) + windPrice(windSize) + batteryPrice(self.batteryModel.powerCapacity)
		self.shortageCost = shortageCostFunction
		self.humanCost = humanCost
		self.humanEmergencyMax = humanEmergencyMax
		self.humanEfficiency = humanEfficiency
		self.energyWaste = 0
		self.lostHours = 0
		self.inverterEfficiency = inverterEfficiency
		self.chargeConverterEfficiency = chargeConverterEfficiency
		self.year = year
		self.timezone = timezone

	def calculateCostFunction(self):
		batteryStates = []
		currentDay = 0
		for day in self.powerLoad:
			currentDate = dt.datetime(self.year,1,1) + dt.timedelta(days=currentDay)
			for hour in day:
				self.timeStep(currentDate, hour)
				currentDate += dt.timedelta(hours=1)
				batteryStates.append(self.batteryModel.storedPower)
			currentDay += 1
		return self.price

	def timeStep(self, dateTime, demand):
		#print self.batteryModel.storedPower
		trueDemand = demand/self.inverterEfficiency #Any power is going to come through the inverter
		solar = self.solarModel.availableWattage(dateTime-dt.timedelta(hours=self.timezone))
		wind = self.windModel.availableWattage(dateTime)
		#Power lost through charge conversion
		self.energyWaste += (1-self.chargeConverterEfficiency)*(solar+wind)

		usefulSolarWind = self.chargeConverterEfficiency*(solar+wind)
		if usefulSolarWind >= trueDemand:
			#Power lost to inverter
			self.energyWaste += trueDemand-demand
			#Efficiency is 1 since it's already been through the Charge converter
			self.energyWaste += self.batteryModel.addEnergy(usefulSolarWind-trueDemand,1)
			self.lostHours = 0
		else:
			#Everything lost by shoving things through an inverter
			self.energyWaste += (1-self.inverterEfficiency)*usefulSolarWind
			remainingEnergy = trueDemand-usefulSolarWind
			#Again, the efficiency is one because the power at this point has
			#already been calculated assumed it's gone through the inverter
			batteryEnergyUsed = self.batteryModel.takeEnergy(remainingEnergy, 1)
			remainingEnergy -= batteryEnergyUsed
			self.energyWaste += (1-self.inverterEfficiency)*batteryEnergyUsed
			#If the battery wasn't enough
			if remainingEnergy > 0:
				self.lostHours += 1
				#Ignoring the inverter since humanPower would be going to either the battery or the 
				#Already inverter-adjusted energy deficit.
				maxHumanPower = self.humanEmergencyMax*self.humanEfficiency*self.chargeConverterEfficiency
				self.price += self.humanCost(1)
				self.energyWaste += self.humanEmergencyMax-maxHumanPower
				#If the person is enough
				if maxHumanPower >= remainingEnergy:
					self.energyWaste += (1-self.inverterEfficiency)*maxHumanPower
					#This shouldn't actually increase energy waste at all, but if there's something
					#Silly happening with the model this should catch it.
					self.energyWaste += self.batteryModel.addEnergy(maxHumanPower-remainingEnergy,1)
				#Worst case scenario, power went out for this hour. Not sure how to calculate
				#power waste for this.
				else:
					missedPower = remainingEnergy - maxHumanPower
					self.price += self.shortageCost(np.array([self.lostHours,missedPower]))
					#Everything you generated this hour was a waste, probably not the best way to do this
					self.energyWaste += trueDemand-missedPower
			#If we managed to meet evrything with the battery
			else:
				self.lostHours = 0

def costFunction(sizes, latitude=42.36, longitude=-71.06):
	batteryVoltage = 12
	print sizes
	test = MicroGrid(latitude, longitude, solarSize = sizes[0], windSize = sizes[1], batterySize = np.array([batteryVoltage,sizes[2]]))
	cost = test.calculateCostFunction()
	print cost
	return cost

boundaries = ((10,None),(1,None),(100,None))

res = minimize(costFunction, (100, 10, 5000), method="L-BFGS-B", bounds = boundaries, options={"maxiter":100})

#cProfile.run("costFunction([100,7,5000])")
print res