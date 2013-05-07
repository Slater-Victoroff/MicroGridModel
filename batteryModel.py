class Battery:

	def __init__(self, voltage=6, capacity=220):
		"""Voltage is in volts, capacity is in Ah, numbers are pretty standard for
		a car battery"""
		self.powerCapacity = voltage*capacity