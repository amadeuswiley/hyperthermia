# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 13:55:38 2015

@author: Christopher
"""

import numpy
import matplotlib.pyplot as plt
import HeatTransferModel as htm
import DiffusionModel as dfm
import CellDeathModel as cdm


# Define intervals used:
M = 20 
# Intervals in r-direction. Becomes unstable above 175.
N = 2000
# Total number of time steps.
center = 0.0 
# Center location. Keep this for expansion to 2 dimensions.
edge = 0.04
# Edge of material. Typical values: [0.03, 0.04 meters].
start = 0.0 
# Starting time. Value in seconds.
stop = 200.0 
# Final time in seconds. End of the model, not the end of the laser.

# Node locations:
r = numpy.linspace(center,edge,M+1)
t = numpy.linspace(start,stop,N+1)
dt = t[1]-t[0]

# Define temperature matrix:
T = numpy.zeros((M+1,N+1)) 
# Note that this specifies (rows,columns), where each row is a spatial node
# and each column is a time node.
T[:,0] = 38.0
# Fills the first column (time step) with 38's.
# That's normal body temperature in degrees C.
HTmodel = htm.HeatTransfer(r)
# Creates a heat transfer model object--runs init.

# Define diffusion matrix:
DIFFmodel = dfm.Diffusion(r,edge)
# Creates an HPI diffusion model object--runs init.
HPIconc = numpy.zeros((M+1,N+1))
# Creates a concentration matrix.
# Set initial concentrations to 1 in the HPI/GNR coupon:
for i in range (0,M):
    if r[i] < DIFFmodel.couponedge:
        HPIconc[i,0] = 1

# Define cell death object:
DRmodel = cdm.CellDeath(r,M,N)
# Creates a death rate model object--runs init.
kt = numpy.zeros((M+1,N+1))
# Creates an integrated injury matrix.
fd = numpy.zeros((M+1,N+1))
# Creates a fraction dead matrix.

# Time loop for calculating temperatures,
# then everything that's based on temperatures.
for j in range (1,N+1):
    time = j*dt
    # Calculate new temperatures and k values based on prior time step:
    T[:,j] = HTmodel.TimeStep(T[:,j-1],time,dt)
####Calculate new diffusivity based on current temperature:
####NEED TO ADD IN DIFFUSIVITY CALL
    # Calculate new concentration based on previous time step's
    # diffusivity and concentrations:
    HPIconc[:,j] = DIFFmodel.TimeStep(HPIconc[:,j-1],dt)

    kt[:,j] = DRmodel.integrateInjury(T[:,j],dt)  # Also easy to store k like T, if desired
    fd[:,j] = DRmodel.fractionDead()

# Generate plots:
plt.plot(r,T[:,0::20])
plt.xlabel('r (meters)')
plt.ylabel('Temperature (C)')
plt.show()

plt.plot(r,fd[:,0::20])
plt.xlabel('r (meters)')
plt.ylabel('Fraction of Cells Dead')
plt.show()

plt.plot(r,HPIconc[:,0::20])
plt.xlabel('r (meters)')
plt.ylabel('HPI Concentration (unitless)')
plt.show()
