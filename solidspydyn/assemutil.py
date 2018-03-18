# -*- coding: utf-8 -*-
"""
assemutil.py
------------

Functions to assemble the system of equations for the Finite Element
Analysis.

"""
from __future__ import division, print_function
import numpy as np
from scipy.sparse import coo_matrix
import uelutil as ue
import femutil as fem


def eqcounter(nodes):
    """Counts active equations and creates BCs array IBC

    Parameters
    ----------
    nodes : ndarray
      Array with nodes coordinates and boundary conditions.

    Returns
    -------
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.

    """
    nnodes = nodes.shape[0]
    IBC = np.zeros([nnodes, 2], dtype=np.integer)
    for i in range(nnodes):
        for k in range(2):
            IBC[i , k] = int(nodes[i , k+3])
    neq = 0
    for i in range(nnodes):
        for j in range(2):
            if IBC[i, j] == 0:
                IBC[i, j] = neq
                neq = neq + 1

    return neq, IBC


def DME(nodes, elements):
    """Counts active equations, creates BCs array IBC[]
    and the assembly operator DME[]

    Parameters
    ----------
    nodes    : ndarray.
      Array with the nodal numbers and coordinates.
    elements : ndarray
      Array with the number for the nodes in each element.

    Returns
    -------
    DME : ndarray (int)
      Assembly operator.
    IBC : ndarray (int)
      Boundary conditions array.
    neq : int
      Number of active equations in the system.

    """
    nels = elements.shape[0]
    IELCON = np.zeros([nels, 9], dtype=np.integer)
    DME = np.zeros([nels, 18], dtype=np.integer)

    neq, IBC = eqcounter(nodes)

    for i in range(nels):
        iet = elements[i, 1]
        ndof, nnodes, ngpts = fem.eletype(iet)
        for j in range(nnodes):
            IELCON[i, j] = elements[i, j+3]
            kk = IELCON[i, j]
            for l in range(2):
                DME[i, 2*j+l] = IBC[kk, l]

    return DME , IBC , neq


def retriever(elements , mats , nodes , i, uel=None):
    """Computes the elemental stiffness matrix of element i

    Parameters
    ----------
    elements : ndarray
      Array with the number for the nodes in each element.
    mats    : ndarray.
      Array with the material profiles.
    nodes    : ndarray.
      Array with the nodal numbers and coordinates.
    i    : int.
      Identifier of the element to be assembled.

    Returns
    -------
    kloc : ndarray (float)
      Array with the local stiffness matrix.
    ndof : int.
      Number of degrees of fredom of the current element.
    """
    par = np.zeros([5], dtype=np.float)
    IELCON = np.zeros([9], dtype=np.integer)
    iet = elements[i , 1]
    ndof, nnodes, ngpts = fem.eletype(iet)
    elcoor = np.zeros([nnodes, 2])
    im = np.int(elements[i, 2])
    par[:] = mats[im, :]
    for j in range(nnodes):
        IELCON[j] = elements[i, j+3]
        elcoor[j, 0] = nodes[IELCON[j], 1]
        elcoor[j, 1] = nodes[IELCON[j], 2]
    if uel is None:
        if iet == 1:
            kloc , mloc , cloc = ue.uel4nquad(elcoor , par)
        elif iet == 2:
            kloc , mloc , cloc = ue.uel6ntrian(elcoor, par)
        elif iet == 3:
            kloc , mloc , cloc = ue.uel3ntrian(elcoor, par)
        elif iet == 5:
            kloc , mloc , cloc = ue.uelspring(elcoor,  par)
        elif iet == 6:
            kloc , mloc , cloc = ue.ueltruss2D(elcoor, par)
        elif iet == 7:
            kloc , mloc , cloc = ue.ueldashpot(elcoor, par)
        elif iet == 8:
            kloc , mloc , cloc = ue.uel9nquad(elcoor , par)
        elif iet == 9:
            kloc , mloc , cloc = ue.uel3dash(elcoor , par)
    else:
        kloc, ndof, iet = uel(elcoor, par)
    
    return kloc , mloc , cloc , ndof , iet


def assembler(elements, mats, nodes, neq, DME, sparse=True, uel=None):
    """Assembles the global stiffness matrix

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    sparse : boolean (optional)
      Boolean variable to pick sparse assembler. It is True
      by default.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
      Array with the global stiffness matrix. It might be
      dense or sparse, depending on the value of _sparse_

    """
    if sparse:
        KG , MG , CG = sparse_assem(elements, mats, nodes, neq, DME, uel=uel)
    else:
        KG , MG , CG = dense_assem(elements, mats, nodes, neq, DME, uel=uel)
    return KG , MG , CG

def effective(KG , MG , CG , ac ):
    

    KE = ac[0]*MG+ac[1]*CG+KG
    
    return KE 



def dense_assem(elements, mats, nodes, neq, DME, uel=None):
    """
    Assembles the global stiffness matrix _KG_
    using a dense storing scheme

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
      Array with the global stiffness matrix in a dense numpy
      array.
    MG : ndarray (float)
      Array with the global mass matrix in a dense numpy
      array.
    """
    KG = np.zeros((neq, neq))
    MG = np.zeros((neq, neq))
    CG = np.zeros((neq, neq))
    nels = elements.shape[0]
    for el in range(nels):
        kloc , mloc , cloc , ndof , iet  = retriever(elements , mats  , nodes , el, uel=uel)
        if iet == 6:
            dme    = np.zeros([ndof], dtype=np.integer)
            dme[0] = DME[el, 0]
            dme[1] = DME[el, 1]
            dme[2] = DME[el, 3]
            dme[3] = DME[el, 4]
        else:
            dme = DME[el, :ndof]

        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        KG[glob_row, glob_col] = KG[glob_row, glob_col] +\
                                                 kloc[row, col]
                        MG[glob_row, glob_col] = MG[glob_row, glob_col] +\
                                                 mloc[row, col]
                        CG[glob_row, glob_col] = CG[glob_row, glob_col] +\
                                                 cloc[row, col]


    return KG , MG , CG


def sparse_assem(elements, mats, nodes, neq, DME, uel=None):
    """
    Assembles the global stiffness matrix _KG_
    using a sparse storing scheme

    The scheme used to assemble is COOrdinate list (COO), and
    it converted to Compressed Sparse Row (CSR) afterward
    for the solution phase [1]_.

    Parameters
    ----------
    elements : ndarray (int)
      Array with the number for the nodes in each element.
    mats    : ndarray (float)
      Array with the material profiles.
    nodes    : ndarray (float)
      Array with the nodal numbers and coordinates.
    DME  : ndarray (int)
      Assembly operator.
    neq : int
      Number of active equations in the system.
    uel : callable function (optional)
      Python function that returns the local stiffness matrix.

    Returns
    -------
    KG : ndarray (float)
      Array with the global stiffness matrix in a sparse
      Compressed Sparse Row (CSR) format.

    References
    ----------
    .. [1] Sparse matrix. (2017, March 8). In Wikipedia,
        The Free Encyclopedia.
        https://en.wikipedia.org/wiki/Sparse_matrix

    """
    rows = []
    cols = []
    kvals = []
    mvals = []
    cvals = []
    nels = elements.shape[0]
    for el in range(nels):
        kloc , mloc , cloc , ndof , iet  = retriever(elements , mats  , nodes , el, uel=uel)
        if iet == 6:
            dme    = np.zeros([ndof], dtype=np.integer)
            dme[0] = DME[el, 0]
            dme[1] = DME[el, 1]
            dme[2] = DME[el, 3]
            dme[3] = DME[el, 4]
        else:
            dme = DME[el, :ndof]
    
        for row in range(ndof):
            glob_row = dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        rows.append(glob_row)
                        cols.append(glob_col)
                        kvals.append(kloc[row, col])
                        mvals.append(mloc[row, col])
                        cvals.append(cloc[row, col])

    stiff = coo_matrix((kvals, (rows, cols)),
                       shape=(neq, neq)).tocsr()
    mass  = coo_matrix((mvals, (rows, cols)),
                      shape=(neq, neq)).tocsr()
    damp = coo_matrix((cvals, (rows, cols)),
                      shape=(neq, neq)).tocsr()

    return stiff, mass , damp


def loadasem(loads, IBC, neq , ninc , T , Tc , fc):
    """Assembles the global Right Hand Side Vector RHSG

    Parameters
    ----------
    loads : ndarray
      Array with the loads imposed in the system.
    IBC : ndarray (int)
      Array that maps the nodes with number of equations.
    neq : int
      Number of equations in the system after removing the nodes
      with imposed displacements.

    Returns
    -------
    RHSG : ndarray
      Array with the right hand side vector.
    """
    nloads = loads.shape[0]
    RHSG = np.zeros((neq, ninc))
#
    Tt= T
    Nt=ninc
    Rick, T=ricker(Nt, Tt, Tc, fc)
#    
    for i in range(nloads):
        il = int(loads[i, 0])
        ilx = IBC[il, 0]
        ily = IBC[il, 1]
        if ilx != -1:
            for k in range(ninc):
                RHSG[ilx , k] = loads[i, 1]*Rick[k]
        if ily != -1:
            for k in range(ninc):
                RHSG[ily , k] = loads[i, 2]*Rick[k]
    return RHSG

def ricker(nt, Tt, tc, fc):
	
    Rick = np.zeros(nt)
    T    = np.zeros(nt)
    dt   = Tt/(nt-1)
	
    for i in range(nt):
        tao=np.pi*fc*(dt*i-tc)
        Rick[i]=(2.*tao**2-1.)*np.exp(-tao**2)
        T[i]= i*dt

    return Rick, T






