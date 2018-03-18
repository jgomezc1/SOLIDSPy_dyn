#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:23:09 2017

@author: casierraa
"""

from __future__ import division, print_function
import numpy as np
from datetime import datetime
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import matplotlib.pyplot as plt
#
folder, name = pre.initial_params()
start_time = datetime.now()
#
#%% PRE-PROCESSING
inipar , nodes, mats, elements, loads = pre.readin(folder=folder)
ninc , T , dt , ac , theta = pre.intparams(inipar)
DME , IBC , neq = ass.DME(nodes, elements)
U    = np.loadtxt(folder + 'response.txt' , ndmin=2)
#
#%% PRE-PROCESSING
salida    = np.loadtxt(folder + 'salida.txt' , ndmin = 1 , dtype=np.int)
npts = pos.sheets(salida , ninc , U , IBC , 'respuesta' , folder )
#pos.PLOTsheets( 'respuesta.txt' , folder , dt , ninc , npts ,200)

vals = np.loadtxt(folder + "respuesta.txt")
plt.close("all")
plt.figure(0)
fig = plt.figure(figsize=(10, 8))
pos.plot_sheet(vals, T, amp_signal=30, amp_shift=50)
plt.figure(1)
fig = plt.figure(figsize=(10, 8))
pos.plot_pcolor(vals, T , -1, 1 )
#
# Arrival times
#
#D = 0.5
#a = 0.3
#h = 0.2
#alpha = 1.0
#t_surface = (D+a)/alpha
#t_c_G_r = (np.sqrt( D**2 +h**2) + np.sqrt( h**2 + a**2))/alpha
#t_G_r = (np.sqrt( h**2 + a**2))/alpha
##
#print ("surface time" , t_surface)
#print ("load crack surface time" , t_c_G_r)
#print ("Crack surface time" , t_G_r)
