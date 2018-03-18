#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script para realización de un análisis con post-procesamiento.
Escribe archivos VTKs, regresa el arreglo con la histora de desplazamintos
y grafica historias (en forma de sabanas) para los puntos de la suprficie
definidos en el archivo de texto salida.txt
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from solidsPy_DYN import solidsPy_DYN
import postprocesor as pos
#
U , folder , IBC , ninc , T = solidsPy_DYN()
salida    = np.loadtxt(folder + 'salida.txt' , ndmin = 1 , dtype=np.int)
npts = pos.sheets(salida , ninc , U , IBC , 'respuesta' , folder )
#
vals = np.loadtxt(folder + "respuesta.txt")
plt.close("all")
#
amp_signal=100
amp_shift =75
plt.figure(0)
fig = plt.figure(figsize=(10, 8))
pos.plot_sheet(vals, T, amp_signal, amp_shift)
#
plt.figure(1)
fig = plt.figure(figsize=(10, 8))
pos.plot_pcolor(vals, T , -1, 1 )