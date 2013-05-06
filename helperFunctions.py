import csv
import datetime
import numpy as np

firstLagrange = "data/longitudeFiles/firstLagrangePoint.csv"
secondLagrange = "data/longitudeFiles/secondLagrangePoint.csv"
thirdLagrange = "data/longitudeFiles/thirdLagrangePoint.csv"
fourthLagrange = "data/longitudeFiles/fourthLagrangePoint.csv"
fifthLagrange = "data/longitudeFiles/fifthLagrangePoint.csv"
sixthLagrange = "data/longitudeFiles/sixthLagrangePoint.csv"

firstLatitude = "data/latitudeFiles/firstLatitudePoint.csv"
secondLatitude = "data/latitudeFiles/secondLatitudePoint.csv"

firstRadius = "data/radiusFiles/firstRadius.csv"
secondRadius = "data/radiusFiles/secondRadius.csv"
thirdRadius = "data/radiusFiles/thirdRadius.csv"
fourthRadius = "data/radiusFiles/fourthRadius.csv"
fifthRadius = "data/radiusFiles/fifthRadius.csv"

lagrangeFiles = [firstLagrange,secondLagrange,thirdLagrange, fourthLagrange, fifthLagrange, sixthLagrange]
latitudeFiles = [firstLatitude, secondLatitude]
radiusFiles = [firstRadius, secondRadius, thirdRadius, fourthRadius, fifthRadius]
nutationFile = "data/nutation/nutation.csv"

class sunPrediction:

	def __init__(self, latitude, longitude, slope, azimuthRotation, elevation=840.0, \
		pressure=1013.25, temperature=14, powerFlux =1000, degrees = True, gregorian = True):
		self.elevation = elevation #meters
		self.pressure = pressure #mbar
		self.temperature = temperature #C
		self.degrees = degrees
		self.gregorian = gregorian
		self.powerFlux = powerFlux #W/m^2
		radians = lambda x: (np.pi/180.)*x
		if degrees == True:
			self.slope = radians(slope)
			self.latitude = radians(latitude)
			self.longitude = radians(longitude)
			self.azimuthRotation = radians(azimuthRotation)
		else:
			self.slope = slope
			self.latitude = latitude
			self.longitude = longitude
			self.azimuthRotation = azimuthRotation

	def julianDay(self, gregorianDateTime, gregorianOffset = None):
		"""Conversion of gregorian calendar to julian calendar,
		since julian gives a much more accurate stepping stone
		to calculate solar incidence. Set gregorianOffset to 
		true when passing in a gregorian day, false when Julian"""
		gregorianOffset = gregorianOffset or self.gregorian
		year = gregorianDateTime.year
		month = gregorianDateTime.month
		if month <=2:
			month = month+12
			year = year-1
		if gregorianOffset == True:
			A = int(year/100.)
			calendarOffset = 2-A+int(A/4.)
		else:
			calendarOffset = 0
		partialDay = gregorianDateTime.hour*3600.0
		partialDay += gregorianDateTime.minute*60.0
		partialDay += float(gregorianDateTime.second)
		partialDay /= 24.0*60*60
		day = gregorianDateTime.day + partialDay
		calendarOffset = calendarOffset
		firstTerm = lambda year: int(365.25*(year+4716))
		secondTerm = lambda month: int(30.6001*(month+1))
		return firstTerm(year)+secondTerm(month)+day+calendarOffset-1524.5

	def julianEphemerisConversion(self, julianDay, timeDelta):
		"""timeDelta is the difference between earth's rotation time
		and terrestrial time. It's very small, so ephemeris time
		is very nearly identical."""
		return julianDay + (timeDelta/86400.0)

	def julianCentury(self, julianDay):
		return (julianDay-2451545.0)/36525.0

	def julianMillennium(self, julianDay):
		return self.julianCentury(julianDay)/10.0

	def lagrangeContribution(self, lagrangeMatrix, julianDay):
		JME = self.julianMillennium(julianDay)
		lagrangeVector = np.zeros(len(lagrangeMatrix))
		for i in range(0,len(lagrangeMatrix)):
			lagrangeVector[i]

	def lagrangeTerm(self, julianDay, filePath):
		"""Returns summed lagrange term"""
		JME = self.julianMillennium(julianDay)
		lagrange = 0
		with open(filePath, 'rb') as lagrangeData:
			reader = csv.reader(lagrangeData, delimiter=",")
			for row in reader:
				angle = float(row[1])+float(row[2])*JME
				lagrange += float(row[0]) * np.cos(angle)
		return lagrange

	def celestialLongitude(self, gregorianDateTime, lagrangeFiles=lagrangeFiles, gregorian=None, geocentric=False):
		"""Expects lagrangeFiles to be ordered filepaths from 0
		to 5. Also, this can directly take in a gregorian datetime
		object and return the heliocentric Longitude of the earth
		in radians if gregorian is set to false, you can directly
		pass in a julian day, though this is mostly for debugging,
		will return geocentric results if it is set to true, despite the name"""
		gregorian = gregorian or self.gregorian
		if gregorian:
			JDay = self.julianDay(gregorianDateTime)
		else:
			JDay = gregorianDateTime
		JME = self.julianMillennium(JDay)
		lagrangeTerms = []
		for filePath in lagrangeFiles:
			lagrangeTerms.append(self.lagrangeTerm(JDay,filePath))
		longitude = 0
		for i in range(0,len(lagrangeFiles)):
			longitude += lagrangeTerms[i]*(JME**i)
		longitude /= (10.0**8)
		normalizedLongitude = np.arccos(np.cos(longitude))
		if (geocentric):
			return normalizedLongitude + np.pi
		else:
			return normalizedLongitude

	def celestialLatitude(self, gregorianDateTime, latitudeFiles=latitudeFiles, gregorian=None, geocentric=False):
		"""Basically the same process for a few different things,
		just breaking them up to make the function calls more
		comprehensible"""
		gregorian = gregorian or self.gregorian
		latitude = self.celestialLongitude(gregorianDateTime, latitudeFiles, gregorian)
		if (geocentric):
			return -latitude
		else:
			return latitude

	def radiusVector(self, gregorianDateTime, radiusFiles=radiusFiles, gregorian=None):
		"""See latitude comment"""
		gregorian = gregorian or self.gregorian
		return self.celestialLongitude(gregorianDateTime, radiusFiles, gregorian)

	def meanMoonElongation(self, gregorianDateTime, gregorian=None):
		"""mean elongation of the moon, again, if gregorian is false
		you can pass in a julian day"""
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		firstTerm = lambda century: 297.85036 + 445267.111480*century
		secondTerm = lambda century: 0.0019142*century**2
		thirdTerm = lambda century: (century**3)/189474.0
		degreeAnswer = firstTerm(JCE)-secondTerm(JCE)+thirdTerm(JCE)
		return self.degreesToRadians(degreeAnswer)

	def meanSunAnomaly(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		firstTerm = lambda century: 357.52772+35999.050340*century
		secondTerm = lambda century: 0.000163*century**2
		thirdTerm = lambda century: (century**3)/300000.0
		degreeAnswer = firstTerm(JCE)-secondTerm(JCE)-thirdTerm(JCE)
		return self.degreesToRadians(degreeAnswer)

	def meanMoonAnomaly(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		firstTerm = lambda century: 134.96298+477198.867398*century
		secondTerm = lambda century: 0.0086972*century**2
		thirdTerm = lambda century: (century**3)/56250.0
		degreeAnswer = firstTerm(JCE)+secondTerm(JCE)+thirdTerm(JCE)
		return self.degreesToRadians(degreeAnswer)

	def moonLatitudeArgument(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		firstTerm = lambda century: 93.27191+483202.017538*century
		secondTerm = lambda century: 0.0036825*century**2
		thirdTerm = lambda century: (century**3)/327270.0
		degreeAnswer = firstTerm(JCE)-secondTerm(JCE)+thirdTerm(JCE)
		return self.degreesToRadians(degreeAnswer)

	def ascendingNodeLongitude(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		firstTerm = lambda century: 125.04452-1934.136261*century
		secondTerm = lambda century: 0.0020708*century**2
		thirdTerm = lambda century: (century**3)/450000.0
		degreeAnswer = firstTerm(JCE)+secondTerm(JCE)+thirdTerm(JCE)
		return self.degreesToRadians(degreeAnswer)

	def nutation(self, gregorianDateTime, nutationFile=nutationFile, gregorian=None):
		"""Returns true obliquity in longitude and obliquity. Returns
		it in a numpy array"""
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		longitudeDelta = 0
		obliquityDelta = 0
		with open(nutationFile, 'rb') as csvData:
			reader = csv.reader(csvData, delimiter = ",")
			for row in reader:
				constants = [float(row[i]) for i in range(0,5)]
				sineTerm = np.sin(self.nutationAngleTerm(gregorianDateTime, constants, gregorian))
				cosineTerm = np.cos(self.nutationAngleTerm(gregorianDateTime, constants, gregorian))
				longitudeDelta += (float(row[5])+(float(row[6])*JCE))*sineTerm
				obliquityDelta += (float(row[7])+(float(row[8])*JCE))*cosineTerm
		longitudeNutation = self.degreesToRadians(longitudeDelta/36000000.0)
		obliquityNutation = self.degreesToRadians(obliquityDelta/36000000.0)
		return np.array([longitudeNutation, obliquityNutation])


	def nutationAngleTerm(self, gregorianDateTime, constants, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JCE = self.julianCentury(self.julianDay(gregorianDateTime))
		else:
			JCE = self.julianCentury(gregorianDateTime)
		angleTerm = self.meanMoonElongation(gregorianDateTime, gregorian)*constants[0]
		angleTerm += self.meanSunAnomaly(gregorianDateTime, gregorian)*constants[1]
		angleTerm += self.meanMoonAnomaly(gregorianDateTime, gregorian)*constants[2]
		angleTerm += self.moonLatitudeArgument(gregorianDateTime, gregorian)*constants[3]
		angleTerm += self.ascendingNodeLongitude(gregorianDateTime, gregorian)*constants[4]
		return angleTerm 

	def trueEclipticObliquity(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			U = self.julianMillennium(self.julianDay(gregorianDateTime))/10.
		else:
			U = self.julianMillennium(gregorianDateTime)/10.
		constants = np.array([84381.448,-4680.93,-1.55,1999.25,-51.38,-249.67,-39.05,7.12,27.87,5.79,2.45])
		meanObliquity = 0
		for i in range(0,11):
			meanObliquity += (constants[i]*U**i)
		radianMeanObliquity = self.degreesToRadians(meanObliquity/3600.)
		return radianMeanObliquity+self.nutation(gregorianDateTime, gregorian=gregorian)[1]

	def aberrationCorrection(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		radius = self.radiusVector(gregorianDateTime, gregorian=gregorian)
		degreeAberation = radius*(20.4898/3600.0)
		return self.degreesToRadians(degreeAberation)

	def apparentSunLongitude(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		longitudeNutation = self.nutation(gregorianDateTime, gregorian=gregorian)[0]
		aberration = self.aberrationCorrection(gregorianDateTime, gregorian=gregorian)
		geocentricLongitude = self.celestialLongitude(gregorianDateTime, gregorian = gregorian, geocentric=True)
		return longitudeNutation + aberration + geocentricLongitude

	def apparentGreenwichSiderealTime(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		if (gregorian):
			JD = self.julianDay(gregorianDateTime)
		else:
			JD = gregorianDateTime
		JC = self.julianCentury(JD)
		dailyTerm = lambda day: 280.46061837+360.98564736629*(day-2451545.0)
		centuryTerm = lambda century: (0.000387933*(century**2))-((century**3)/38710000.0)
		meanSiderealTime = self.degreesToRadians(dailyTerm(JD)+centuryTerm(JC))
		normalizedSiderealTime = (180/np.pi)*(np.arccos(np.cos(meanSiderealTime)))
		longitudeNutation = (180/np.pi)*self.nutation(gregorianDateTime, gregorian=gregorian)[0]
		obliquity = self.trueEclipticObliquity(gregorianDateTime, gregorian)
		return (2*np.pi)-self.degreesToRadians(normalizedSiderealTime+longitudeNutation*np.cos(obliquity))

	def geocentricSunCoordinates(self, gregorianDateTime, gregorian=None):
		gregorian = gregorian or self.gregorian
		"""returns an array containing the geocentric right ascension and declination
		of the sun"""
		apparentLongitude = self.apparentSunLongitude(gregorianDateTime, gregorian=gregorian)
		eclipticObliquity = self.trueEclipticObliquity(gregorianDateTime, gregorian=gregorian)
		earthLatitude = self.celestialLatitude(gregorianDateTime, gregorian=gregorian, geocentric=True)
		longitudeEcliptic = np.sin(apparentLongitude)*np.cos(eclipticObliquity)
		latitudeEcliptic = np.tan(earthLatitude)*np.sin(eclipticObliquity)
		rightAscension = (2*np.pi)+np.arctan2((longitudeEcliptic-latitudeEcliptic),np.cos(apparentLongitude))
		tripleProjection = np.cos(earthLatitude)*np.sin(eclipticObliquity)*np.sin(apparentLongitude)
		crossLatitudeEcliptic = np.sin(earthLatitude)*np.cos(eclipticObliquity)
		declination = np.arcsin(tripleProjection+crossLatitudeEcliptic)
		return np.array([rightAscension, declination])

	def localHourAngle(self, gregorianDateTime, longitude=None, degrees=None, gregorian=None):
		"""Assumes the input latitude is in degrees, because sadly that's what
		people use instead of radians, but converting it to radians internally.
		latitude should be positive for east of Greenwich, and negative for
		west of Greenwich return is measured westward from south"""
		gregorian = gregorian or self.gregorian
		longitude = longitude or self.longitude
		degrees = degrees or self.degrees
		greenwichSiderealTime = self.apparentGreenwichSiderealTime(gregorianDateTime, gregorian)
		geocentricRightAscension = self.geocentricSunCoordinates(gregorianDateTime, gregorian)[0]
		return np.arccos(np.cos(greenwichSiderealTime+longitude-geocentricRightAscension))

	def topocentricCoordinates(self, gregorianDateTime, latitude=None, longitude=None, elevation=None, 
		pressure=None, temperature=None, degrees=None, gregorian=None):
		"""Returns numpy array with rightAscension, declination, local hour angle,
		 zenith angle, then azimuth angle topographically.
		 Latitude and longitude can be either in degrees or radians, just change the degrees
		 flag accordingly. Pressure should be in millibars, and temperature should be in Celsius.
		 elevation should be in meters.

		 Longitude is positive east of Greenwich, and negative west of Greenwich. 
		 Latitude is positive north of the equator, and negative south of the equator"""
		gregorian = gregorian or self.gregorian
		longitude = longitude or self.longitude
		latitude = latitude or self.latitude
		elevation = elevation or self.elevation
		pressure = pressure or self.pressure
		temperature = temperature or self.temperature
		degrees = degrees or self.degrees

		equitorialHorizontalParallax = (8.794/3600)/self.radiusVector(gregorianDateTime,gregorian=gregorian)
		trueLatitude = np.arctan(0.99664719*np.tan(latitude))
		x = np.cos(trueLatitude)+(elevation/6378140.)*np.cos(latitude)
		y = 0.99664719*np.sin(trueLatitude)+(elevation/6378140.)*np.sin(latitude)
		sunCoordiantes = self.geocentricSunCoordinates(gregorianDateTime, gregorian)
		hourAngle = self.localHourAngle(gregorianDateTime, longitude,False, gregorian)
		
		#Right Ascension Calculation
		raNumerator = -x*np.sin(equitorialHorizontalParallax)*np.sin(hourAngle)
		raDenominator = np.cos(sunCoordiantes[1])-(x*np.sin(equitorialHorizontalParallax)*hourAngle)
		rightAscensionParallax = np.arctan2(raNumerator,raDenominator)
		rightAscension = sunCoordiantes[0]+rightAscensionParallax

		#Declination Calculation
		declinationNumerator = (np.sin(sunCoordiantes[1])-y*np.sin(equitorialHorizontalParallax))*np.cos(rightAscensionParallax)
		declinationDenominator = np.cos(sunCoordiantes[1])-x*np.sin(equitorialHorizontalParallax)*np.cos(hourAngle)
		declination = np.arctan2(declinationNumerator, declinationDenominator)

		#Hour angle calculation
		trueHourAngle = hourAngle-rightAscensionParallax

		#Zenith angle calculation
		latitudeDeclination = np.sin(latitude)*np.sin(declination)
		secondaryTerm = np.cos(latitude)*np.cos(declination)*np.cos(trueHourAngle)
		elevationAngle = (180./np.pi)*np.arcsin(latitudeDeclination+secondaryTerm)
		seeingApproximation = (pressure/1010.)*(283.0/(273+temperature))
		degreeElevationTerm = 1.02/(60*np.tan(self.degreesToRadians(elevationAngle+(10.3/(elevationAngle+5.11)))))
		atmosphericCorrection = self.degreesToRadians(seeingApproximation*degreeElevationTerm)
		zenithAngle = (np.pi/2)-(atmosphericCorrection+self.degreesToRadians(elevationAngle))

		#Azimuth Angle Calculation
		astronomersDenominator = np.cos(trueHourAngle)*np.sin(latitude)-np.tan(declination)*np.cos(latitude)
		azimuthAngle = np.arctan2(np.sin(trueHourAngle),astronomersDenominator) + np.pi
		#Azimuth angle here is measured eastward from North

		return np.array([rightAscension, declination, trueHourAngle, zenithAngle, azimuthAngle])

	def incidenceAngle(self, gregorianDateTime, latitude=None, longitude=None, slope=None, azimuthRotation=None,  
	 elevation=None, pressure=None, temperature=None,degrees = None, gregorian = None):
		"""aka the moneymaker. This will just directly give you the incidence angle for any given
		plane on any day, anywhere in the world. Slope is the angle in degrees or radians (degrees argument)
		measured from horizontal, azimuth rotation angle is positive if east from south, negative if west from south"""
		gregorian = gregorian or self.gregorian
		longitude = longitude or self.longitude
		latitude = latitude or self.latitude
		elevation = elevation or self.elevation
		pressure = pressure or self.pressure
		temperature = temperature or self.temperature
		degrees = degrees or self.degrees
		slope = slope or self.slope
		azimuthRotation = azimuthRotation or self.azimuthRotation
		topocentricOrientation = self.topocentricCoordinates(gregorianDateTime, latitude, longitude, elevation, 
								pressure, temperature, degrees, gregorian)
		firstTerm = np.cos(topocentricOrientation[3])*np.cos(slope)
		relativeAzimuthalAngle = (topocentricOrientation[4] - np.pi)-azimuthRotation
		secondTerm = np.sin(slope)*np.sin(topocentricOrientation[3])*np.cos(relativeAzimuthalAngle)
		return np.arccos(firstTerm+secondTerm)

	def degreesToRadians(self, degrees):
		return degrees*(np.pi/180.0)

def parseDataToCsv(filePath):
	with open(filePath, 'rb') as rawData:
		counter = 0
		fullMatrix = []
		currentRow = []
		for line in rawData:
			print line
			if counter == 4:
				fullMatrix.append(currentRow)
				counter = counter%4
				currentRow = []
			if "A-" not in line:
				if counter != 0:
					currentRow.append(float(line))
				counter += 1
	with open(filePath, 'wb') as dataDump:
		writer = csv.writer(dataDump, delimiter=",",
							quoting=csv.QUOTE_MINIMAL)
		for row in fullMatrix:
			print row
			writer.writerow(row)

#print apparentGreenwichSiderealTime(2452930.312847, False)
#print(incidenceAngle(2452930.312847,39.742476, -105.1786, 30, -10, True, gregorian=False))
longitude = -71.1786 #East of Greenwich is positive, West is negative
latitude = 42.742476 #North of equator is positive, South is negative
slope = 60.0 #Slope with respect to horizontal (angle)
azimuthRotation = -10.0 #Azimuthal angle of solar panel
pressure = 1830.14 #Pressure in mbar
elevation = 820.0 #elevation in meters
temperature = 11.0 #temperature in c
testing = sunPrediction(latitude, longitude, slope, azimuthRotation, elevation, pressure, temperature)
test = datetime.datetime.now() + datetime.timedelta(hours=5)
degrees = lambda x: (180./np.pi)*x
print degrees(testing.incidenceAngle(test))