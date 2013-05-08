import numpy as np
import math
import random
import matplotlib.pyplot as plot

class PowerLoad:

	def __init__(self, mean=2000, peakHeight=None, peakTime=None, peakStDv=None, troughHeight=None, troughTime=None, troughStDv=None, span = 24, granularity=1, dailyNoise = [0.0,0.0,0.0,0.0,0.0,0.0,0.0], hourlyNoise = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]):
		'''Mean is the average level of power consumption
		peakHeight is the highest level of powerConsumption,
		peakTime is the time at which that peak usually occurs,
		troughHeight is the lowest level of power consumption,
		troughTime is the time at which that happens,
		span is the total time of the cycle in hours,
		granularity is the width of the time bins also in hours,
		noise is the level of randomness you  want to add to each bin also in hours
		Right now assume that this is all Gaussian
		Noise order is mean, peakStDv, peakHeight, peakTime, troughStDv, troughHeight, troughTime
		same ordering for both hourly and daily noise'''
		self.mean = mean
		self.peakHeight = peakHeight or mean*1.25
		self.peakTime = peakTime or span*(2./3)
		self.troughHeight = troughHeight or mean*0.75
		self.troughTime = troughTime or span*(1./6)
		self.span = span
		self.granularity = granularity
		self.dailyNoise = dailyNoise
		self.hourlyNoise = hourlyNoise
		self.peakStDv = peakStDv or span/10.
		self.troughStDv = troughStDv or span/10.
		length = int(float(span)/granularity)
		self.times = np.zeros([1,length])
		self.values = np.zeros([1,length])

	def constructDay(self, mean = None):
		'''Actually fills in the times and values arrays'''
		randomize = lambda x: 1+(random.random()*x)-(x/2.)
		mean = mean or self.mean
		mean += self.dailyNoise[0]
		peakStDv = self.peakStDv*randomize(self.dailyNoise[1])
		peakHeight = self.peakHeight*randomize(self.dailyNoise[2])
		peakTime = self.peakTime*randomize(self.dailyNoise[3])
		troughStDv = self.troughStDv*randomize(self.dailyNoise[4])
		troughHeight = self.troughHeight*randomize(self.dailyNoise[5])
		troughTime = self.troughTime*randomize(self.dailyNoise[6])
		length = int(float(self.span)/self.granularity)
		dayValues = np.zeros(length)

		peakConstant = 1.0/((peakStDv*randomize(self.hourlyNoise[1]))*(math.sqrt(2*math.pi)))
		troughConstant = 1.0/((troughStDv*randomize(self.hourlyNoise[4]))*(math.sqrt(2*math.pi)))
		for i in range(0,len(self.times[0])):
			self.times[0,i] = (i*self.granularity)
			x = (i+0.5)*self.granularity
			peakExponent = -((x-(peakTime*randomize(self.hourlyNoise[3])))**2)/(2.0*(peakStDv*randomize(self.hourlyNoise[1]))**2)
			troughExponent = -((x-(troughTime+randomize(self.hourlyNoise[6])))**2)/(2.0*(troughStDv*randomize(self.hourlyNoise[4]))**2)
			peakValue = abs((peakHeight*randomize(self.hourlyNoise[2]))-self.mean)*(math.exp(peakExponent))
			troughValue = abs((troughHeight*randomize(self.hourlyNoise[5]))-self.mean)*(math.exp(troughExponent))
			value = peakValue-troughValue+mean*randomize(self.hourlyNoise[0])
			dayValues[i] = (value)
		return dayValues

	def plot(self, values):
		plot.plot(values)
		plot.show()
	
def yearLoad(mean=2000):
	test = PowerLoad(mean=mean, span=365, peakTime = 250, troughTime = 100, granularity = 1)
	dayTest = PowerLoad(span=24, granularity=1)
	testDay = test.constructDay()
	allData = []
	for entry in testDay:
		allData.append(dayTest.constructDay(mean=entry))
	return np.array(allData)