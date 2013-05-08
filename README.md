This library is a way of quantitatively creating and optimizing microgrids utilizing solar power, wind power, and human power. The current implementation is currently not quite as good as it could be, but with feedback I hope to make it better.

**Required libraries**
This module relies upon numpy, scipy, and pysolar. 

Currently the pysolar portion of this is by far the slowest, but I plan on speeding this up by integrating it with my previous sun prediction code and making heavier use of numpy calls.

Also, this still requires people to use the source directly and I'm going to get some setup scripts up in here by the end of the week.

There are a number of ways in which you could use this library, but there are four main objects worth interacting with:

SolarModel models an array of solar panels. It takes in the latitude and longitude at which the panels are mounted, and also the size of the array in square meters and the efficiency of the panels, currently the latitude and longitude should be passed in degrees with north and east representing positive latitudes and longitudes respectively:

import dateTime
from solarModel import SolarModel
test = SolarModel(42, -71, 100, 0.17)
#Time passed to this argument must be in utc.
print test.availableWattage(datetime.datetime.utcnow())

>>>>>>15420.7362878

***Note you will get different results depending on when you run this, but it should make sense with what is happening near you at that time. Results are in Watts

The same method can be used for the windModel by passing in rotor diameter and then efficiency. There are optional parameters for air density (default 1.2) and wind function (default is a seasonal Wiebull distribution)

The wind module has additional methods wind and powerTest. Wind will return the wind speed in meters per second rather than the available wattage (which has an equivalent for the wind method)

power test will take in the provided power of the turbine, along with the windspeed and rotor diameter in m/s and m respectively and output an efficiency of that turbine.

The battery model takes in voltage and capacity (in Ah), followed by a numpy array that specifies the minimum and maximum states the battery is allowed to reach. Additionally there is a stored power option that will let you start your batteries off with a particular amount of energy in them. Default is half the max storage.

The two methods provided are addEnergy and takeEnergy that will let you try to give and take energy respectively without violating your contraints on the maximum and minimum state of charge of the battery. They return the amount of power wasted, and the amount of power actually provided by the battery respectively.

There is also a powerload object, and a full microgrid object. The microgrid is really well documented so I won't repeat the documentation here, and the powerload is pretty complicated and the defaults should work pretty well. If you want to get the full powerload for a year you should call the yearload function from the powerload file without a mean wattage that you want the powerload to emulate.

More documentation to come after finals.