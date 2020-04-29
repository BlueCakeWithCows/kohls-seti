# kohls-seti
# Earth-Moon-Earth Doppler Factor Predictions
Mason Kohls, Berkeley SETI Research Center, Spring 2020

## Introduction
One method for determining Earth's technosignature is to examine the radio waves reflected by the moon back towards the Earth. Radio signals emitted from the Earth and reflected off the Moon are subject to two doppler shifts: when the Moon 'bounces' the signal and when it is recieved on Earth. In this project I generate predictions for the doppler factor for arbitrary emitter and reciever locations on Earth.
## Background
The doppler shift as described in the reciever frame is . 
Consider the a signal being emitted from the earth to the moon. Let the position and velocity of the emitter be described by A and A' respectively. Let the position and velocity of the moon be described by B and B' respectively. Then, in the moon frame (reciever frame) for the emitter we have position p = A - B,  velocity v = A' - B'. The moon's position and velocity being 0. Now cos(θ) = p.v/(|p||v|). Thus we can calculate the doppler shift. 
To get the final doppler shift we simply take the newly shifted signal and apply the reverse process, using a location on Earth as the reciever and the moon as an emitter. Notice, if we use the same location on earth twice then its as if we applied the first doppler factor twice. 

In order to get an idea of scale, lets examine the largest expected doppler shift, assuming the moon's own motion is neglible with respect to the Earth's rotation. 
We can use this fact to calculate a worst case bounds on expected doppler shift as a function of frequency. At the equator this is approximately 460 m/s. (Graph Here)

## Tools and Dependencies
NumPy, Astropy, jplephem, matplotlib
