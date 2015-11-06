# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 18:36:36 2015

@author: Christopher
"""

import numpy
import matplotlib.pyplot as plt


class Diffusion:
    def __init__(self,r,edge):
        self.r = r
        self.h = r[1]-r[0]
        self.M = len(r) - 1
        self.Cnot = 1 
        # Defines the initial dimensionless concentration of HPIs.
        self.k = 0.0001 
        # Defines binding. No clue what the value should be.
        
        # Diffusion parameters:
        self.DiffusionCoefficient = 0.000001
#TODO: figure out which units this is in, and what values are reasonable.
        self.couponedge = 0.02
        # Defines where in space the drug coupon extends to. Units are meters. 
#TODO: STILL NEEDS TO KNOW WHAT THE DIFFUSION COEFFICIENT IS, BASED ON TEMPERATURE.
        if self.couponedge > edge:
            print('error: the HPI/GNR coupon cannot be bigger than the sample.\
            Fix either "edge" or "couponedge."')


    def TimeStep(self, oldC, dt):
        newC = numpy.zeros((self.M+1))
        for i in range (1,self.M):
            newC[i] = oldC[i]
            newC[i] += (oldC[i+1]-2*oldC[i]+oldC[i-1])                        \
                        *self.DiffusionCoefficient*dt/(self.h**2)
            newC[i] += self.DiffusionCoefficient*dt/(2*self.r[i]*self.h)      \
                        *(oldC[i+1]-oldC[i-1])
            newC[i] -= self.k*dt*oldC[i]
            # Boundary conditions:
            newC[0] = newC[1]
            # This is the dC/dr=0 condition.
            newC[self.M] = 0.0
            # This is the condition that the edge of the culture
            # has a zero concentration of HPI's.
        return(newC)
            
            
            
if __name__ == '__main__':
    
    M = 50 
    # Intervals in r-direction. Becomes unstable above 175.
    N = 200
    # Total number of time steps.
    edge = 0.08
    start = 0.0
    stop = 150.0
    # Final time in seconds. End of the model, not the end of the laser.
    r = numpy.linspace(0.0,edge,M+1)
    # Second term defines radius of sample in meters.
    t = numpy.linspace(start,stop,N+1)
    dt = t[1]-t[0]
    
    DIFFmodel = Diffusion(r,edge)
    # Creates an HPI diffusion model object--runs init.
    HPIconc = numpy.zeros((M+1,N+1))
    # Creates the concentration matrix.
    for i in range (0,M):
        if r[i] < DIFFmodel.couponedge:
            HPIconc[i,0] = 1   
            # Sets initial concentrations to 1 in the coupon.

    # Time Step:        
    for j in range (1,N+1):
        HPIconc[:,j] = DIFFmodel.TimeStep(HPIconc[:,j-1],dt)

    # Generate plots:
    plt.plot(r,HPIconc[:,0::20])
    # Plots certain columns (timesteps) of the matrix HPIconc.
    plt.xlabel('r (meters)')
    plt.ylabel('HPI conc (unitless)')
    plt.show()
