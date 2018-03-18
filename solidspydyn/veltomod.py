#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:07:19 2017

@author: casierraa
"""
from __future__ import division, print_function
import numpy as np

beta  = 1000.0
alpha = 2000.0
rho   = 1000.0

c1 = rho*beta**2
c2 = 3.0*alpha**2-4.0*beta**2
c3 = alpha**2-beta**2
E = c1*c2/c3
c1 = 2.0*rho*beta**2
nu = E/c1-1.0
print(E)
print
print(nu)


#E =2960000.0
#nu = 0.3333
#rho = 1.0
#
#shear = E/(2.0*(1+nu))
#betasq = shear/rho
#Beta = np.sqrt(betasq)
#print
#print(Beta)
#
#alphasq = (2.0*shear*(1.0-nu))/(rho*(1.0-2.0*nu))
#Alpha = np.sqrt(alphasq)
#print
#print(Alpha)