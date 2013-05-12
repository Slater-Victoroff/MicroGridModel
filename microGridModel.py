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

	def __init__(self, latitude, longitude, year=2012, timezone=3, solarPrice=lambda x: 300*x,windPrice=lambda x: 500*x,
		batteryPrice = lambda capacity:capacity*0.283,shortageCostFunction = lambda params: params[1]*11, 
		humanCost = lambda x: 1.5*x, humanEmergencyMax = 100, humanEfficiency = 0.7, inverterEfficiency=0.85,
		chargeConverterEfficiency=0.9, solarEfficiency=0.17, windEfficiency=0.22, batterySize = np.array([6,220]), 
		solarSize = 50, windSize = 13, batterySocBounds = np.array([0.2,0.9]), meanPower = 2000, granularity=1,
		cloudiness=0.33, solarCloudLoss = 0.2):
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
		but they seem like overkill for right now.

		granularity is your bin size in hours for the year of simulation"""
		self.solarModel = SolarModel(latitude, longitude, solarSize, solarEfficiency)
		self.windModel = WindModel(windSize, windEfficiency)
		self.powerLoad = yearLoad(meanPower, granularity)
		self.batteryModel = Battery(batterySize[0], batterySize[1], batterySocBounds)
		self.price = solarPrice(solarSize) + windPrice(windSize) + batteryPrice(self.batteryModel.energyCapacity)
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
		self.granularity = granularity
		self.solarCloudLoss = solarCloudLoss
		self.cloudVector = ((np.random.rand(1,365*24*self.granularity))/cloudiness).astype(int)

	def calculateCostFunction(self):
		batteryStates = []
		vectorCounter = 0
		currentDay = 0
		for day in self.powerLoad:
			currentDate = dt.datetime(self.year,1,1) + dt.timedelta(days=currentDay)
			for hour in day:
				self.timeStep(currentDate, hour, self.cloudVector[0][vectorCounter])
				currentDate += dt.timedelta(hours=self.granularity)
				batteryStates.append(self.batteryModel.storedEnergy)
				vectorCounter += 1
			currentDay += 1
		return self.price

	def timeStep(self, dateTime, demand, cloudy):
		"""Demand is coming in units of power draw in watts"""
		trueDemand = demand/self.inverterEfficiency #Any power is going to come through the inverter
		solarWatts = self.solarModel.availableWattage(dateTime-dt.timedelta(hours=self.timezone))
		windWatts = self.windModel.availableWattage(dateTime)
		if cloudy:
			solarWatts *= self.solarCloudLoss
		#Power lost through charge conversion
		self.energyWaste += (1-self.chargeConverterEfficiency)*(solarWatts+windWatts)

		usefulSolarWind = self.chargeConverterEfficiency*(solarWatts+windWatts)
		if usefulSolarWind >= trueDemand:
			#Power lost to inverter
			self.energyWaste += (trueDemand-demand)*self.granularity #Units of Wh
			#Efficiency is effectively 1 since it's already been through the Charge converter
			self.energyWaste += self.batteryModel.addEnergy((usefulSolarWind-trueDemand)*self.granularity,1)
			self.lostHours = 0
		else:
			#Everything lost by shoving things through an inverter
			self.energyWaste += (1-self.inverterEfficiency)*usefulSolarWind*self.granularity
			remainingEnergy = (trueDemand-usefulSolarWind)*self.granularity
			#Again, the efficiency is one because the power at this point has
			#already been calculated assumed it's gone through the inverter
			batteryEnergyUsed = self.batteryModel.takeEnergy(remainingEnergy, 1)
			remainingEnergy -= batteryEnergyUsed
			self.energyWaste += (1-self.inverterEfficiency)*batteryEnergyUsed
			#If the battery wasn't enough
			if remainingEnergy > 0:
				self.lostHours += self.granularity
				#Ignoring the inverter since humanPower would be going to either the battery or the 
				#Already inverter-adjusted energy deficit.
				maxHumanPower = self.humanEmergencyMax*self.humanEfficiency*self.chargeConverterEfficiency
				self.price += self.humanCost(1)*self.granularity
				self.energyWaste += (self.humanEmergencyMax-maxHumanPower)*self.granularity
				#If the person is enough
				if maxHumanPower >= remainingEnergy:
					self.energyWaste += ((1-self.inverterEfficiency)*maxHumanPower)*self.granularity
					#This shouldn't actually increase energy waste at all, but if there's something
					#Silly happening with the model this should catch it.
					self.energyWaste += self.batteryModel.addEnergy(maxHumanPower-remainingEnergy,1)*self.granularity
				#Worst case scenario, power went out for this hour. Not sure how to calculate
				#power waste for this.
				else:
					missedEnergy = remainingEnergy - (maxHumanPower*self.granularity)
					self.price += self.shortageCost(np.array([self.lostHours,missedEnergy]))
					#Everything you generated this hour was a waste, probably not the best way to do this
					self.energyWaste += trueDemand-missedEnergy
			#If we managed to meet evrything with the battery
			else:
				self.lostHours = 0

class CostFunction:
	def __init__(self, cloudiness, outputFile):
		self.cloudiness = cloudiness
		self.outputFile = outputFile

	def costFunction(self, sizes, latitude=1.1, longitude=32.4):
		with open(self.outputFile, 'a') as output:
			output.write("CLOUDINESS: " + str(self.cloudiness)+"\n")
			batteryVoltage = 12
			print sizes
			stringSizes = str(sizes)
			output.write(stringSizes)
			output.write("\n")
			test = MicroGrid(latitude, longitude, cloudiness = self.cloudiness, solarSize = sizes[0], windSize = sizes[1], batterySize = np.array([batteryVoltage,sizes[2]]))
			cost = test.calculateCostFunction()
			print cost
			stringCost = str(cost)
			output.write(stringCost)
			output.write("\n")
		return cost

def cloudyParameterSweep(minimum, maximum, steps, method, outputFile = "results.txt"):
	change = (maximum-minimum)/float(steps)
	startingGuess = (100, 10, 5000)
	boundaries = ((10,None),(1,None),(100,None))
	maxIterations = 100
	updatingFile = "resultStream.txt"
	for i in range(0,steps):
		current = CostFunction(minimum+i*change, updatingFile)
		result = minimize(current.costFunction,startingGuess,method=method,bounds=boundaries,options={"maxiter":maxIterations})
		with open(outputFile, 'a') as output:
			output.write(str(result))
			output.write("\n")

def solarSizeSweep(minimum, maximum, steps, windSize, batterySize, outputFile="solarResults.txt"):
	change = (maximum-minimum)/float(steps)
	latitude = 1.1
	longitude = 32.4
	costs = []
	sizes = [minimum+i*change for i in range(0,steps)]
	with open(outputFile, 'a') as dataDump:
		for i in range(0,steps):
			current = MicroGrid(latitude, longitude, windSize = windSize, batterySize= batterySize, solarSize = minimum+i*change)
			costs.append(current.calculateCostFunction())
			dataDump.write("Solar size: " + str(minimum+i*change)+"\n")
			dataDump.write(str(costs)+"\n")
			print costs
	plot.plot(sizes,costs)
	plot.savefig("solarSizeSweep.png")
	return sizes[costs.index(max(costs))]

def windSizeSweep(minimum, maximum, steps, solarSize, batterySize, outputFile="windResults.txt"):
	change = (maximum-minimum)/float(steps)
	latitude = 1.1
	longitude = 32.4
	costs = []
	sizes = [minimum+i*change for i in range(0,steps)]
	with open(outputFile, 'a') as dataDump:
		for i in range(0,steps):
			current = MicroGrid(latitude, longitude, solarSize = solarSize, batterySize= batterySize, windSize = minimum+i*change)
			costs.append(current.calculateCostFunction())
			dataDump.write("Wind size: " + str(minimum+i*change)+"\n")
			dataDump.write(str(costs)+"\n")
			print costs
	plot.plot(sizes,costs)
	plot.savefig("windSizeSweep.png")
	plot.show()
	return sizes[costs.index(max(costs))]

def batterySizeSweep(minimum, maximum, steps, solarSize, windSize, outputFile="batteryResults.txt"):
	change = (maximum-minimum)/float(steps)
	latitude = 1.1
	longitude = 32.4
	costs = []
	sizes = [minimum+i*change for i in range(0,steps)]
	with open(outputFile, 'a') as dataDump:
		for i in range(0,steps):
			current = MicroGrid(latitude, longitude, solarSize = solarSize, windSize= windSize, batterySize = [12,minimum+i*change])
			costs.append(current.calculateCostFunction())
			dataDump.write("Battery size: " + str(minimum+i*change)+"\n")
			dataDump.write(str(costs)+"\n")
			print costs
	plot.plot(sizes,costs)
	plot.savefig("batterySizeSweep.png")
	return sizes[costs.index(max(costs))]

#cloudyParameterSweep(0.1,0.5,50,"TNC")
windMin = windSizeSweep(1,20,100,100,[12,3600])
#solarMin = solarSizeSweep(10,500,100,windMin,[12,3600])
#batteryMin = batterySizeSweep(100,5000,100,solarMin,windMin)

print solarMin
print windMin
print batteryMin
#cProfile.run("costFunction([100,7,5000])")

#Vary Cloudiness from 0.1 to 0.5
#Vary cost of Solar Panels from 100 to 500 keep size of system
#Vary cost of Wind Turbines from 300 to 700
#Vary cost of Batteries from 0.083 to 0.417

#Make one randomized file and run off that

#solar max = 1000; wind max = 20 m; battery max = 5000 Ah, Diesel max = 100 l/hr