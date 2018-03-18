#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def ricker(nt, Tt, tc, fc):
    """
    Computes a Ricker pulse of:
        Central frequency fc
        Time centered at tc
        Time windowed at Tt
        Samped at nt points
    """
    time = np.linspace(0, Tt, nt)
    tau = np.pi*fc*(time - tc)
    Rick = (2.0*tau**2 - 1.0) * np.exp(-tau**2)
    return Rick, time

def grafsignalG(A, dt , Ngra):
    """
     Plots the generalized time signal A[ndats] into Ngra
    """
    ndats  = len(A)
    x=np.zeros([ndats], dtype=float)
    x=np.arange(0,ndats*dt,dt)
    plt.figure(Ngra)
    plt.plot(x,A)
    plt.grid()
#
    return

def Ftrans(datos, ndats, dt, fs):
    """
    Compute the Fourier spectra of datos[] of
    length ndats and sampled at dt.
    Returns the result in Samag after smoothing by the
    smoothing factor fs.
    """
    nfr = int(ndats/2)
    df = 1.0/(ndats*dt)
    x = np.arange(df,nfr*df, df)
    A = np.fft.fft(datos)
    Aa = np.abs(A)

    # Smooth the spectrum.
    Sa = Aa[1:nfr]
    Samag = smooth(Sa , x , fs) 
    nfs = nfr-1
    return x , Samag , A , nfs

def smooth(Sa, Freq, fftfs):
    """
    Parameters
    ----------
    Sa : ndarray
        Original spectrum.
    Freq : float
        Frequency.
    fftfs : float
        Smoothing factor.
    """
    Sas  = np.zeros([len(Sa)],dtype=float)
    fia = 1
    fma = 1
    suma = Sa[0] * Sa[0]
    pot = 1./2./fftfs
    fsexpi = 2**(-pot)
    fsexpm = 2**( pot)
    Sas[0] = Sa[0]
    NNfft = len(Sa)
    for i in range(1, NNfft):
#
        fi = int((i + 1) * fsexpi)
        fm = int((i + 1) * fsexpm)
        if fi < 1:
            fi = 1
        if fm > NNfft:
            fm = NNfft

        for Num in range(fia - 1, fi - 1):
#
            suma = suma - Sa[Num] * Sa[Num]

        for Num in range(fma, fm):
#
            suma = suma + Sa[Num]*Sa[Num]

        Nf = fm - fi + 1
        fia = fi
        fma = fm
        Sas[i]=np.sqrt(suma/Nf)
    return (Sas)

def grafFourier(Sas , x , nfr , Nfig):
    """
     Plots the Fourier spectral amplitude Sas into Nfig.
     Sas : Spectrum
     x   : frecuency
     xmin,xmax,ymin,ymax
    """
#
    plt.figure(Nfig)
    plt.plot(x,Sas)
    plt.grid()
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Amplitud')
#
    return
#
