# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 13:55:38 2015

@author: Christopher
"""

import numpy
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import HeatTransferModel as htm
import DiffusionModel as dfm
import CellDeathModel as cdm

#TODO: main() needs to be a function of a couple things for the optimizer.
def main(OD=0.5):
    
    # Define intervals used:
    M = 20 
    # Intervals in r-direction. Becomes unstable above 175.
    N = 2000
    # Total number of time steps.
    center = 0.0 
    # Center location. Keep this for expansion to 2 dimensions.
    edge = 0.04
    # Edge of material. Typical values: [0.03, 0.04 meters].
    canceredge = 0.02
    # Edge of cancerous tissue.
    start = 0.0 
    # Starting time. Value in seconds.
    stop = 200.0 
    # Final time in seconds. End of the model, not the end of the laser.
    
    #Heat Transfer Model Values:
    alpha = 3.0*10**(-7.0) 
    # This is the tissue or fluid thermal diffusivity
    # w/ GNRs, equal to k/(density*specific heat) units are (m^2)/s.
    Q = 9.0*10**6 
    # Laser energy. Varies between 8.0*10**6 and 3.0*10**7 W/(m^3). 
    
    # Diffusion parameters:
    DiffusionCoefficient = 0.000001
#TODO: figure out which units this is in, and what values are reasonable.
    couponedge = 0.02
    # Defines where in space the drug coupon extends to. Units are meters. 
#TODO: STILL NEEDS TO KNOW WHAT THE DIFFUSION COEFFICIENT IS, BASED ON TEMPERATURE.
    
    
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
    HTmodel = htm.HeatTransfer(r, OD, alpha, Q)
    # Creates a heat transfer model object--runs init.
    
    # Define diffusion matrix:
    DIFFmodel = dfm.Diffusion(r,edge,couponedge,DiffusionCoefficient)
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

#    # Generate plots:
#    plt.plot(r,T[:,0::20])
#    plt.xlabel('r (meters)')
#    plt.ylabel('Temperature (C)')
#    plt.show()
#    
#    plt.plot(r,fd[:,0::20])
#    plt.xlabel('r (meters)')
#    plt.ylabel('Fraction of Cells Dead')
#    plt.show()
#    
#    plt.plot(r,HPIconc[:,0::20])
#    plt.xlabel('r (meters)')
#    plt.ylabel('HPI Concentration (unitless)')
#    plt.show()
    # Create nodes for optimizing
    healthy0 = fd[numpy.round(M+1)*canceredge/edge+0,N]
    healthy1 = fd[numpy.round(M+1)*canceredge/edge+1,N]
    healthy2 = fd[numpy.round(M+1)*canceredge/edge+2,N]
    healthy3 = fd[numpy.round(M+1)*canceredge/edge+3,N]
    healthy4 = fd[numpy.round(M+1)*canceredge/edge+4,N]
#    print(healthy0,healthy1,healthy2,healthy3,healthy4)
    optimizervalue = healthy0 + 2*healthy1 + 4*healthy2 + 16*healthy3 + 64*healthy4
    print("optimizervalue = %s " % (optimizervalue))
    print("OD + %s" % (OD))
    return optimizervalue
    
# Optimizer
#main(OD=0.5)

optimize.minimize(main,0.5)
