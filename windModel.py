import numpy as np
wind_hours = []
data = []

def Gaussian_wind(peak, vert_shift, horiz_shift, width, x):
        windy_data = (peak-vert_shift)* 2.17**(-(x-horiz_shift)**2/(width)) + vert_shift
        windy_data = round(windy_data,2)
        return windy_data

#80 Winter Days from Jan 1 to March 20
for i in range(24):
        load = i
        solar = 2*i
        wind = Gaussian_wind(6.3, 5, 12.5, 30 ,i)
        data.append((load,wind,solar))

#91 Spring Days from March 21 to June 20
for i in range(91):
    for i in range(24):
        load = i
        solar = 2*i
        wind = Gaussian_wind(7.2, 6, 12, 45,i)
        data.append((load,wind,solar))

#92 Summer Days from June 21 to Sept 20
for i in range(92):
    for i in range(24):
        load = i
        solar = 2*i
        wind = Gaussian_wind(6.1, 4.9, 13, 70,i)
        data.append((load,wind,solar))

#91 Fall Days from Sept 21 to Dec 2
for i in range(91):
    for i in range(24):
        load = i
        solar = 2*i
        wind = Gaussian_wind(6.4, 5.1, 13.5, 50,i)
        data.append((load,wind,solar))

#11 Winter Days from Dec 21 to Dec 31
for i in range(11):
    for i in range(24):
        load = i
        solar = 2*i
        wind = Gaussian_wind(6.3, 5, 12.5, 30,i)
        data.append((load,wind,solar))
        
print(data)
