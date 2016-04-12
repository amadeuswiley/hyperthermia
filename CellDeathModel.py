# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 14:56:54 2015

This script handles cell death. It calculates cell death rates at a given
temperature, and can also integrate those rates into an overall death fraction
for a given spatial node. Cell death follows an Arrhenius-type equation.

As of this writing (2016 - 04 - 12), this model has not been finely tuned for
a specific cell line. Nor has the impact of HSPI concentration on the death
parameters been implemented, as the data did not exist when I started this
project.


@author: Christopher
"""

import numpy


class CellDeath:
    def __init__(self,r,M,N):
        self.A = 6*10**29
        # Preexponential factor for Arrhenius equation.
        self.Ea = 200000
        # Activation energy for Arrheneus equation. Units are J/mol.
        self.R = 8.3144598
        # Universal gas constant for Arrheneus equation. Units are J/(mol*K).
        self.M = M
        self.N = N
        self.kt = numpy.zeros((M+1))
        
    def integrateInjury(self, T, dt):
    # Finds the death rate of cells for each element in the matrix T.

        # Space loop:
        for i in range (0,self.M+1):
            self.kt[i] += dt*self.A*numpy.exp(-self.Ea/(self.R*(T[i]+273.15)))
        return self.kt

    def fractionDead(self):
    # Finds fraction of dead cells for each element in matrix T.
        fd = numpy.zeros((self.M+1))
        
        # Space loop:
        for i in range (self.M+1):
            fd[i] = 1 - numpy.exp(-self.kt[i])
        return fd
