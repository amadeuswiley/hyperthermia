# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 13:55:38 2015

This is the driver script for the hyperthermia model. This script defines most
of the variables utilized in the submodels, creates the Python objects to be
used in the submodels, calls them in order, and yields a matrix of fractional
cell death values, over space and time. This matrix is used as the input to
the variable "optimizervalue", which looks at the fractional cell death both
in the bulk cancerous region and the healthy regions.

Using the method of Lagrange Multipliers
(https://en.wikipedia.org/wiki/Lagrange_multiplier), which is a way of turning
a constrained optimization into an unconstrained optimization, the variable
"optimizervalue" grows with the death of healthy cells and the retention of
cancerous cells. In this way, wrapping an optimizer around the entire script
and asking it to "minimize(optimizervalue)" succeeds in achieving cancer cell
death and the retention of healthy cells.

It should be noted that the scipy.optimize() is designed to work with
non-dimensionalized values. Trust me on this; don't feed it anything with a
range outside of [0,1].

As of this writing (2016 - 04 - 12), the DiffusionModel, which models the
change in concentration of Heat Shock Protein Inhibitors (drugs which increase
the cell's affinity to dying by heat), does not impact cell death rates at all.
In talking with Dr. Heys, we hypothesized that the HSPI's would lower the
activation energy for cell death; that is, it would change the "Ea" variable
in CellDeathModel.py. However, when I started this project, the data did not
exist to characterize these changes.



@author: Christopher Wiley
"""


import numpy
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import HeatTransferModel as htm
import DiffusionModel as dfm
import CellDeathModel as cdm



# In order to wrap an optimizer around the script, the model needs to be
# formulated as a function with an input and output.
def main(Qnondim):
    
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
    
    # Heat Transfer Model Values:
    alpha = 3.0*10**(-7.0) 
    # This is the tissue or fluid thermal diffusivity
    # w/ GNRs, equal to k/(density*specific heat) units are (m^2)/s.
    Q = (8.0*10**6)*Qnondim
    # Laser energy density. Varies between 8.0*10**6 and 3.0*10**7 W/(m^3).
    # The optimizer passes in Qnondim, which should (0,1). Q redimentionalizes.
    OD = 0.5
    # Optical Density of tissues; indirect measure of GNR conc'n.
    
    # Diffusion parameters:
    DiffusionCoefficient = 0.000001
#TODO: figure out which units this is in, and what values are reasonable.
    couponedge = 0.02
    # Defines where in space the drug coupon extends to. Units are meters. 
    
    
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
        
        # Calculate new diffusivity based on current temperature:
        # NEED TO ADD IN DIFFUSIVITY CALL, once the data exists.
    
        # Calculate new concentration based on previous time step's
        # diffusivity and concentrations:
        HPIconc[:,j] = DIFFmodel.TimeStep(HPIconc[:,j-1],dt)
        
        # Calculate Cell Death Values (see CellDeathModel.py)
        kt[:,j] = DRmodel.integrateInjury(T[:,j],dt)
        #Note: Also easy to store k like T, if desired
        fd[:,j] = DRmodel.fractionDead()

    # Generate plots:
    plt.plot(r,T[:,0::20])
    plt.xlabel('r (meters)')
    plt.ylabel('Temperature (C)')
    plt.show()
    
    plt.plot(r,fd[:,0::20])
    plt.xlabel('Distance from Center of Tumor, r (meters)')
    plt.ylabel('Fraction Dead (1 is all dead, 0 is all alive)')
    plt.show()
    
    plt.plot(r,HPIconc[:,0::20])
    plt.xlabel('r (meters)')
    plt.ylabel('HPI Concentration (unitless)')
    plt.show()


    # Find number of spatial nodes with cancer
    numbercancernodes = numpy.rint((M+1)*canceredge/edge)
    numbercancernodes = int(numbercancernodes)
#    Uncomment the following line if you want to see how many cancer nodes:
#    print(numbercancernodes)
    
    # Create nodes for optimizing
    healthy0 = fd[numbercancernodes+1,N]
    healthy1 = fd[numbercancernodes+2,N]
    healthy2 = fd[numbercancernodes+3,N]
    healthy3 = fd[numbercancernodes+4,N]
    healthy4 = fd[numbercancernodes+5,N]
#    Uncomment the following for debugging purposes, if needed.
#    print(healthy0,healthy1,healthy2,healthy3,healthy4)
    healthyopt = (healthy0 + 2*healthy1 + 4*healthy2 + 16*healthy3 \
    + 64*healthy4)/87
    # This variable, 'healthyopt', contributes to the variable 'optimizervalue'
    # by penalizing the death of healthy cells. It does this more strictly
    # the further into the healthy cell region. I chose to do this in a
    # quadratic manner, but that was largely arbitrary. It's divided by 87 to
    # non-dimensionalize the variable's range.
    
    
    # Find average cancer fraction dead
    cancertotal = 0
    for i in range(0,numbercancernodes):
        cancertotal += fd[i,N]
    canceraverage = cancertotal/(1+numbercancernodes)
    # Convert to a variable    
    canceropt = 1 - canceraverage
    # Create lagrange multiplier, which scales how important killing all the
    # cancer is to the optimizer
    lagrange = 25
    optimizervalue = (healthyopt + lagrange *canceropt)/(1+lagrange)
    
    
    # Can comment out the following print() lines. They're useful for debugging
    print("healthyopt= %s " % (healthyopt))
    print("optimizervalue = %s " % (optimizervalue))
    print("canceraverage = %s " % (canceraverage))
    print("Qnondim = %s" % (Qnondim))
    return optimizervalue


#To run the optimizer, uncomment the next section:
##Optimizer:
#res=optimize.minimize(main,1.0,options={'gtol': 1e-8})
#print(res)


#To run the program without the optimizer, uncomment the next section:
#Program
main(0.75618441)
