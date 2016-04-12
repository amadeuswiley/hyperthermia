# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 13:55:38 2015


This script models heat transfer in the system, using a second-order Finite 
Difference Method (https://en.wikipedia.org/wiki/Finite_difference) of the
heat equation (https://en.wikipedia.org/wiki/Heat_equation).

As of this writing (2016 - 04 - 12), this script contains some variables for
laser operation {laseroff,laseredge} which may be useful as optimizable 
variables in the future.


The code under the line
if __name__ == '__main__':
allows the script to be run standalone, for debugging purposes.

@author: Christopher
"""

import numpy
import matplotlib.pyplot as plt


class HeatTransfer:
    def __init__(self, r, OD, alpha, Q):
        self.r = r
        self.h = r[1]-r[0]
        self.M = len(r) - 1
        self.alpha = alpha 
        # This is the tissue or fluid thermal diffusivity
        # w/ GNRs, equal to k/(density*specific heat) units are (m^2)/s.
        self.Q = Q
        # Laser energy. Varies between 8.0*10**6 and 3.0*10**7 W/(m^3). 
        self.rho = 1000.0 
        # Tissue density. Units are kg/m^3.
        self.heatcap = 3400.0 
        # Tissue heat capacity. Units are J/(kg*K).
        self.OD = OD
        # Optical density of laser. Range: [0.0,0.5]. Unitless.


        # Laser controls:
        self.laseroff = 30.0 
        # This is the time (in seconds) when the laser is turned off.
        self.laseredge = 0.02 
        # This is the radius of the laser beam, in meters. 
        # May not be larger than the cell culture itself.
        if self.laseredge > self.r[self.M]:
            print('Error: laser edge must not extend beyond tissue sample. \
            Either increase variable "edge" or decrease variable "laseredge"')

    def TimeStep(self, oldT, time, dt):
        newT = numpy.zeros((self.M+1))
        for i in range (1,self.M):    
            if self.r[i] > self.laseredge or time > self.laseroff:
                laser = 0
            else:
                laser = self.Q*(1.0-10.0**(-self.OD))/(self.rho*self.heatcap)
            newT[i] = oldT[i] 
            newT[i] += dt*(self.alpha*(oldT[i+1]-2.0*oldT[i]+oldT[i-1])     \
                         /(self.h**2)) 
            newT[i] += dt*(self.alpha*(oldT[i+1]-oldT[i-1])                 \
                         /(self.r[i]*2.0*self.h)+laser)
        
        # Boundary conditions:
        newT[0] = oldT[0] 
        newT[0] += dt*(self.alpha*(2.0*oldT[1]-2.0*oldT[0])/(self.h**2)) 
        if time <= self.laseroff:
            newT[0] += dt*self.Q*(1.0-10.0**(-self.OD))/(self.rho*self.heatcap) 
        newT[self.M] = 38.0
        # The edge of the culture is defined to be normal body temperature.
        return(newT)
        
if __name__ == '__main__':
    M = 20 
    # Intervals in r-direction.
    N = 200 
    # Total number of time steps.
    edge = 0.08
    r = numpy.linspace(0.0,edge,M+1) 
    # Second term defines radius of sample in meters.
    start = 0.0 
    # Starting time. Value in seconds.
    stop = 150.0 
    # Final time in seconds. End of the model, not the end of the laser.
    t = numpy.linspace(start,stop,N+1)
    # Second term defines number of seconds to be modelled.
    dt = t[1]-t[0]
    T = numpy.zeros((M+1,N+1)) 
    # Note that this specifies (rows,columns), where each row is a spatial node
    # and each column is a time node.
    T[:,0] = 38.0
    # Fills the first column (time step) with 38's.
    # That's normal body temperature in degrees C.
    HTmodel = HeatTransfer(r,OD=0.5) 
    # Creates a heat transfer model object--runs init.

    # Time loop:
    for j in range (1,N+1):
        time = j*dt
        T[:,j] = HTmodel.TimeStep(T[:,j-1],time,dt)

    # Generate plots:
    plt.plot(r,T[:,0::20]) 
    # Plots certain columns (timesteps) of the matrix T (temperature).
    plt.xlabel('r steps')
    plt.ylabel('Temperature (degrees C)')
