# kohls-seti
# Earth-Moon-Earth Doppler Factor Predictions
Mason Kohls, Berkeley SETI Research Center, Spring 2020

## Introduction
One method for determining Earth's technosignature is to examine the radio waves reflected by the moon back towards the Earth. Radio signals emitted from the Earth and reflected off the Moon are subject to two doppler shifts: when the Moon 'bounces' the signal and when it is recieved on Earth. In this project I generate predictions for the doppler factor for arbitrary emitter and reciever locations on Earth.
## Background
The doppler shift as described in the reciever frame is . 
Consider the a signal being emitted from the earth to the moon. Let the position and velocity of the emitter be described by A and A' respectively. Let the position and velocity of the moon be described by B and B' respectively. Then, in the moon frame (reciever frame) for the emitter we have position p = A - B,  velocity v = A' - B'. The moon's position and velocity being 0. Now cos(θ) = p.v/(|p||v|). Thus we can calculate the doppler shift. 
To get the final doppler shift we simply take the newly shifted signal and apply the reverse process, using a location on Earth as the reciever and the moon as an emitter. Notice, if we use the same location on earth twice then its as if we applied the first doppler factor twice. 

In order to get an idea of scale, lets examine the largest expected doppler shift, assuming the moon's own motion is neglible with respect to the Earth's rotation. At the equator this is approximately 460 m/s. We can see the doppler shift is on the order of ten-thousand times less than the signal frequency. For readability, most of the following graphs will use a signal frequency of 10 MHz. 
![](./images/worst_case_estimate.png)

## Example Graph
An example output graph is shown below.

![](images/ex_1.png)

The top graph displays the amount of Doppler shift at a given time. The below graph displays the altitude of the moon from the view of the reciever and emitter respectively. Notice doppler shift is only given when both altitudes are above zero. To disable this filtering use:
```python
DopplerSystemEME(..., only_visible = False)
```


## Usage
Example of calculating doppler factors for two emitters.
``` python
from (name here) import DopplerSystemEME, getEarthGrid
iers.conf.iers_auto_url = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
solar_system_ephemeris.set('jpl') 

start_string = '2020-01-01T00:00:00.00Z'
timestep = units.Quantity('15 min')
timespan = units.Quantity('4 hour')
emitters = [
  EarthLocation.from_geodetic(lat = "38°25'53.6", lon ="-79°49'08.9", height=0),
  EarthLocation.from_geodetic(lat =  "78°25'53.6", lon ="-79°49'08.9", height=0)
  ]
reciever = EarthLocation.from_geodetic(lat = "38°25'53.6", lon ="-79°49'08.9", height=0)

start_time = Time(start_string, format='isot', scale='utc')
doppler_system = DopplerSystemEME(emitters, reciever, start_time, timestep, timespan, only_visible = True)
doppler_system.update()
```
To plot all the emitters at once use 
```python
  doppler_system.save_plots('mycooldirectory/')
```

To retrieve the data or a single plot 
```python
  doppler_system.get_dopplers() #Returns list of np arrays of doppler factors
  doppler_system.get_moon_altitudes_reciever() #Returns np array of moon altitude
  doppler_system.get_moon_altitude_emitters() #Reurns list of np arrays of moon altitude
  doppler_system.get_plot(idx) #Plots, but does not save, a single figure. 
```

## Tools and Dependencies
NumPy, Astropy, jplephem, matplotlib
