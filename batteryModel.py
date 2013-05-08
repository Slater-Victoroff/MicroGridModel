import numpy as np

class Battery:

	def __init__(self, voltage=6, capacity=220, stateBounds = np.array([0.2,0.9]), storedPower = None):
		"""Voltage is in volts, capacity is in Ah, numbers are pretty standard for
		a car battery storedPower is the number of watt hours the battery starts with"""
		self.powerCapacity = voltage*capacity
		self.storedPower = storedPower or self.powerCapacity/2.0
		self.minState = stateBounds[0]
		self.maxState = stateBounds[1]

	def addEnergy(self, watts, efficiency=0.9):
		"""Efficiency is the efficiency of the charge controller feeding into the battery.

		Returns amount of wasted power"""
		effectiveWatts = watts*efficiency
		#If this won't put the battery over its max charge
		if ((self.storedPower+effectiveWatts)/self.powerCapacity)<self.maxState:
			self.storedPower += effectiveWatts
			return 0
		else:
			addedCapacity = self.maxState-(self.storedPower/self.powerCapacity)
			self.storedPower = self.maxState*self.powerCapacity
			return (effectiveWatts - addedCapacity*self.powerCapacity)
	
	def takeEnergy(self, watts, efficiency=0.85):
		"""Efficiency is the efficiency of the inverter connected to the battery.

		Returns amount of power actually gotten"""
		requestedPower = watts/efficiency
		if (self.storedPower-requestedPower)/self.powerCapacity > self.minState:
			self.storedPower -= requestedPower
			return watts
		else:
			powerAvailable = self.storedPower - (self.minState*self.powerCapacity)
			self.storedPower = self.minState*self.powerCapacity
			return powerAvailable

