import numpy as np
import math
import datetime

class WindModel:

    def __init__(self, rotorDiameter=2.1, efficiency=0.22, airDensity = 1.2, cutInSpeed=2.5, ratedSpeed=14, windFunction=None):
        """Wind function is a lambda function based on a params variable
        that includes peak, vertShift, horizShift, and x"""
        self.airDensity = airDensity
        self.efficiency = efficiency
        self.rotorDiameter = rotorDiameter
        defaultWind = lambda params: (params[0]-params[1])*2.17**((-(params[4]-params[2])**2)/params[3])+params[1]
        self.ratedSpeed = ratedSpeed
        self.cutInSpeed = cutInSpeed
        self.maxPower = 0.5*self.efficiency*self.airDensity*(ratedSpeed**3)*(np.pi/4.0)*self.rotorDiameter
        if windFunction:
            self.windFunction = windFunction
        else:
            self.windFunction = defaultWind

    def wind(self, datetime):
        """Returns wind speed on given day in m/s"""
        #80 Winter Days from Jan 1 to March 20
        #91 Spring Days from March 21 to June 20
        #92 Summer Days from June 21 to Sept 20
        #91 Fall Days from Sept 21 to Dec 2
        #11 Winter Days from Dec 21 to Dec 31
        trueDay = datetime.timetuple().tm_yday
        if (trueDay<80 or trueDay>=354):
            params = np.array([6.3,5.,12.5,30.,datetime.hour])
        elif (trueDay<171):
            params = np.array([7.2, 6., 12., 45., datetime.hour])
        elif (trueDay<263):
            params = np.array([6.1,4.9,13.,70., datetime.hour])
        elif (trueDay<354):
            params = np.array([6.4,5.1,13.5,50., datetime.hour])
        return self.windFunction(params)

    def availableWattage(self, datetime, wind=None):
        """Approximating a logistic trend in the power curve.
        The magine 0.7783 number is the distance of the standard
        logistic curve from 2% rails, but that number can be changed you like"""
        logisticCurveRail = 0.7783
        windSpeed = wind or self.wind(datetime)
        if windSpeed<self.cutInSpeed:
            return 0
        elif windSpeed>self.ratedSpeed:
            return self.maxPower
        else:
            logisticCurve = lambda x : 1/(1+math.exp(-x))
            adjustedSpeed = windSpeed-(self.cutInSpeed+logisticCurveRail)
            ratio = (2*logisticCurveRail)/(self.ratedSpeed-self.cutInSpeed)
            return logisticCurve(ratio*adjustedSpeed)

    def powerCurve(self, maximum, minimum, steps):
        changes = (maximum-minimum)/steps
    def powerTest(self, providedPower = 900., speed=12.5, rotorDiameter=2.1):
        return providedPower/(0.5*self.airDensity*(speed**3)*((np.pi*rotorDiameter**2)/4.0))