# -*- coding: utf-8 -*-
"""
gaussutil.py
------------
Weights and coordinates for Gauss-Legendre quadrature [1]_. The
values for triangles is presented in section 5.5 of Bathe book [2]_.

References
----------
.. [1] Wikipedia contributors. "Gaussian quadrature." Wikipedia,
  The Free Encyclopedia, 2 Nov.  2015. Web. 25 Dec. 2015.
  url: https://en.wikipedia.org/wiki/Gaussian_quadrature
.. [2] Bathe, Klaus-Jürgen. Finite element procedures. Prentice Hall,
   Pearson Education, 2006.
"""
from __future__ import division, print_function
import numpy as np


def gpoints2x2():
    """Gauss points for a 2 by 2 grid

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([4])
    xp = np.zeros([4, 2])
    xw[:] = 1.0
    xp[0, 0] = -0.577350269189626
    xp[1, 0] = 0.577350269189626
    xp[2, 0] = -0.577350269189626
    xp[3, 0] = 0.577350269189626

    xp[0, 1] = 0.577350269189626
    xp[1, 1] = 0.577350269189626
    xp[2, 1] = -0.577350269189626
    xp[3, 1] = -0.577350269189626

    return xw, xp

def gpoints3x3():
    """Gauss points for a 3 by 3 grid

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    XW   = np.zeros([9])
    XP   = np.zeros([9, 2])
    RLGP = np.zeros([9, 2])
#
    ZERO = 0.0
    ONE  = 1.0
    RLGP[0 , 0]=-ONE
    RLGP[1 , 0]=ZERO
    RLGP[2 , 0]= ONE
    RLGP[3 , 0]=-ONE
    RLGP[4 , 0]=ZERO
    RLGP[5 , 0]= ONE
    RLGP[6 , 0]=-ONE
    RLGP[7 , 0]=ZERO
    RLGP[8 , 0]= ONE
#
    RLGP[0 , 1]=-ONE
    RLGP[1 , 1]=-ONE
    RLGP[2 , 1]=-ONE
    RLGP[3 , 1]=ZERO
    RLGP[4 , 1]=ZERO
    RLGP[5 , 1]=ZERO
    RLGP[6 , 1]= ONE
    RLGP[7 , 1]= ONE
    RLGP[8 , 1]= ONE

    XW[0]=0.555555555555556**2
    XW[1]=0.555555555555556*0.888888888888889
    XW[2]=0.555555555555556**2
#
    XW[3]=0.555555555555556*0.888888888888889
    XW[4]=0.888888888888889**2
    XW[5]=0.555555555555556*0.888888888888889
#
    XW[6]=0.555555555555556**2
    XW[7]=0.555555555555556*0.888888888888889
    XW[8]=0.555555555555556**2    
   
    G = np.sqrt(0.60)
    for i in range(9):
        for j in range(2):
            XP[ i , j ] = G*RLGP[ i , j]    

    return XW, XP


def gpoints7():
    """Gauss points for a triangle (7 points)

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([7])
    xp = np.zeros([7, 2])
    xw[0] = 0.1259391805448
    xw[1] = 0.1259391805448
    xw[2] = 0.1259391805448
    xw[3] = 0.1323941527885
    xw[4] = 0.1323941527885
    xw[5] = 0.1323941527885
    xw[6] = 0.225

    xp[0, 0] = 0.1012865073235
    xp[1, 0] = 0.7974269853531
    xp[2, 0] = 0.1012865073235
    xp[3, 0] = 0.4701420641051
    xp[4, 0] = 0.4701420641051
    xp[5, 0] = 0.0597158717898
    xp[6, 0] = 0.3333333333333

    xp[0, 1] = 0.1012865073235
    xp[1, 1] = 0.1012865073235
    xp[2, 1] = 0.7974269853531
    xp[3, 1] = 0.0597158717898
    xp[4, 1] = 0.4701420641051
    xp[5, 1] = 0.4701420641051
    xp[6, 1] = 0.3333333333333

    return xw, xp


def gpoints3():
    """Gauss points for a triangle element (3 points)

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([3])
    xp = np.zeros([3, 2])
    xw[0] = 0.3333333333333
    xw[1] = 0.3333333333333
    xw[2] = 0.3333333333333

    xp[0, 0] = 0.1666666666667
    xp[1, 0] = 0.6666666666667
    xp[2, 0] = 0.1666666666667

    xp[0, 1] = 0.1666666666667
    xp[1, 1] = 0.1666666666667
    xp[2, 1] = 0.6666666666667

    return xw, xp

def gpoints6():
    """Gauss points for a triangle (7 points)

    Returns
    -------
    xw : ndarray
      Weights for the Gauss-Legendre quadrature.
    xp : ndarray
      Points for the Gauss-Legendre quadrature.

    """
    xw = np.zeros([6])
    xp = np.zeros([6])
    xw[0] = 0.171324492379170
    xw[1] = 0.360761573048139
    xw[2] = 0.467913934572691
    xw[5] = xw[0]
    xw[4] = xw[1]
    xw[3] = xw[2]
    
    xp[0]=-0.932469514203152
    xp[1]=-0.661209386466265
    xp[2]=-0.238619186083197
    xp[5]= 0.932469514203152
    xp[4]= 0.661209386466265
    xp[3]= 0.238619186083197   

    return xw, xp
