import numpy as np

class Battery:

	def __init__(self, voltage=6, capacity=220, stateBounds = np.array([0.2,0.9]), storedEnergy = None):
		"""Voltage is in volts, capacity is in Ah, numbers are pretty standard for
		a car battery storedPower is the number of watt hours the battery starts with"""
		self.energyCapacity = voltage*capacity
		self.storedEnergy = storedEnergy or self.energyCapacity/2.0
		self.minState = stateBounds[0]
		self.maxState = stateBounds[1]

	def addEnergy(self, wattHours, efficiency=0.9):
		"""Efficiency is the efficiency of the charge controller feeding into the battery.

		Returns amount of wasted energy"""
		effectiveEnergy = wattHours*efficiency
		#If this won't put the battery over its max charge
		if ((self.storedEnergy+effectiveEnergy)/self.energyCapacity)<self.maxState:
			self.storedEnergy += effectiveEnergy
			return 0
		else:
			addedCapacity = self.maxState-(self.storedEnergy/self.energyCapacity)
			self.storedEnergy = self.maxState*self.energyCapacity
			return (effectiveEnergy - addedCapacity*self.energyCapacity)
	
	def takeEnergy(self, wattHours, efficiency=0.85):
		"""Efficiency is the efficiency of the inverter connected to the battery.

		Returns amount of power actually gotten"""
		requestedEnergy = wattHours/efficiency
		if (self.storedEnergy-requestedEnergy)/self.energyCapacity > self.minState:
			self.storedEnergy -= requestedEnergy
			return wattHours
		else:
			energyAvailable = self.storedEnergy - (self.minState*self.energyCapacity)
			self.storedEnergy = self.minState*self.energyCapacity
			return energyAvailable

