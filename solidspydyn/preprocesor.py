# -*- coding: utf-8 -*-
"""
This module contains functions to preprocess the input files to compute
a Finite Element Analysis.

"""
from __future__ import division, print_function
import sys
import numpy as np


def readin(folder=""):
    """Read the input files"""
    nodes    = np.loadtxt(folder + 'nodes.txt' , ndmin=2)
    mats     = np.loadtxt(folder + 'mater.txt' , ndmin=2)
    elements = np.loadtxt(folder + 'eles.txt'  , ndmin=2, dtype=np.int)
    loads    = np.loadtxt(folder + 'loads.txt' , ndmin=2)
    inipar   = np.loadtxt(folder + 'inipar.txt', ndmin=2)
    

    return inipar , nodes, mats, elements, loads

def intparams(inipar):
#
    ass = np.zeros(9, dtype = float)
    
    dt     = inipar[0 , 0]
    T      = inipar[0 , 1]
    Tc     = inipar[0 , 2]
    fc     = inipar[0 , 3]
#
    m=int(T/dt)
    theta = 1.40
#
    ass[0] = 6.0/(theta*theta*dt*dt)
    ass[1] = 3.0/(theta*dt)
    ass[2] = 2.0*ass[1]
    ass[3] = theta*dt/2.0
    ass[4] = ass[0]/theta
    ass[5] = -ass[2]/theta
    ass[6] =1.0 - (3.0/theta)
    ass[7] = dt/2.0
    ass[8] = dt*dt/6.0
#
    return  m , T , Tc , fc , dt , ass , theta

def echomod(nodes, mats, elements, loads, folder=""):
    """Create echoes of the model input files"""
    np.savetxt(folder + "KNODES.txt", nodes, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KMATES.txt", mats, fmt='%5.2f', delimiter=' ')
    np.savetxt(folder + "KELEMS.txt", elements, fmt='%d', delimiter=' ')
    np.savetxt(folder + "KLOADS.txt", loads, fmt='%5.2f', delimiter=' ')


def initial_params():
    """Read initial parameters for the simulation
    
    The parameters to be read are:

    - folder: location of the input files.
    - name: name for the output files (if echo is True).
    - echo: echo output files.
    """
    # Check Python version
    version = sys.version_info.major
    if version == 3:
        raw_input = input
    elif version == 2:
        pass
    else:
        raise ValueError("You should use Python 2.x at least!")

    # Try to run with easygui
    try:
        import easygui
        folder = easygui.diropenbox(title="Folder for the job") + "/"
        name = easygui.enterbox("Enter the job name")
#        echo = easygui.buttonbox("Do you want to echo files?",
#                                 choices=["Yes", "No"])

    except:
        folder = raw_input('Enter folder (empty for the current one): ')
        name   = raw_input('Enter the job name: ')
#        echo   = raw_input('Do you want to echo files? (y/N):')

#    if echo.upper() in ["YES", "Y"]:
#        echo = True
#    else:
#        echo = False

    return folder, name


def ele_writer(cells, cell_data ,ele_tag , phy_sur,  ele_type, mat_tag, nini):
    """
    Extracts a subset of elements from a complete mesh according to the
    physical surface  phy_sur and writes down the proper fields into an
    elements array.

    Parameters
    ----------
        cell : dictionary
            Dictionary created by meshio with cells information.
        cell_data: dictionary
            Dictionary created by meshio with cells data information.
        ele_tag : string
            Element type according to meshio convention,
            e.g., quad9 or line3.
        phy_sur : int
            Physical surface for the subset.
        ele_type: int
            Element type.
        mat_tag : int
            Material profile for the subset.
        ndof : int
            Number of degrees of freedom for the elements.
        nnode : int
            Number of nodes for the element.
        nini : int
        Element id for the first element in the set.

    Returns
    -------
        nf : int
            Element id for the last element in the set
        els_array : int
            Elemental data.

    """
    eles = cells[ele_tag]
    dict_nnode = {'triangle': 3 , 'triangle6':6 , 'quad':4 }
    nnode = dict_nnode[ele_tag]
    phy_surface = cell_data[ele_tag]['physical']
    ele_id = [cont for cont, _ in enumerate(phy_surface[:])
              if phy_surface[cont] == phy_sur]
    els_array = np.zeros([len(ele_id) , 3 + nnode], dtype=int)
    els_array[: , 0] = range(nini , len(ele_id) + nini )
    els_array[: , 1] = ele_type
    els_array[: , 2] = mat_tag
    els_array[: , 3::] = eles[ele_id, :]
    nf = nini + len(ele_id)
    return nf , els_array


def node_writer(points , point_data):
    """Write nodal data as required by SolidsPy

    Parameters
    ----------
    points : dictionary
        Nodal points
    point_data : dictionary
        Physical data associatted to the nodes.

    Returns
    -------
    nodes_array : ndarray (int)
        Array with the nodal data according to SolidsPy.

    """
    nodes_array = np.zeros([points.shape[0], 5])
    nodes_array[:, 0] = range(points.shape[0])
    nodes_array[:, 1:3] = points[:, :2]
    return nodes_array


def boundary_conditions(cells, cell_data, phy_lin, nodes_array, bc_x, bc_y):
    """Impose nodal point boundary conditions as required by SolidsPy

    Parameters
    ----------
        cell : dictionary
            Dictionary created by meshio with cells information.
        cell_data: dictionary
            Dictionary created by meshio with cells data information.
        phy_lin : int
            Physical line where BCs are to be imposed.
        nodes_array : int
            Array with the nodal data and to be modified by BCs.
        bc_x, bc_y : int
            Boundary condition flag along the x and y direction:
                * -1: restrained
                * 0: free

    Returns
    -------
        nodes_array : int
            Array with the nodal data after imposing BCs according
            to SolidsPy.

    """
    lines = cells["line"]
    # Bounds contains data corresponding to the physical line.
    phy_line = cell_data["line"]["physical"]
    id_frontera = [cont for cont in range(len(phy_line))
                   if phy_line[cont] == phy_lin]
    nodes_frontera = lines[id_frontera]
    nodes_frontera = nodes_frontera.flatten()
    nodes_frontera = list(set(nodes_frontera))
    nodes_array[nodes_frontera, 3] = bc_x
    nodes_array[nodes_frontera, 4] = bc_y
    return nodes_array


def loading(cells, cell_data, phy_lin, P_x, P_y):
    """Impose nodal boundary conditions as required by SolidsPy

    Parameters
    ----------
        cell : dictionary
            Dictionary created by meshio with cells information.
        cell_data: dictionary
            Dictionary created by meshio with cells data information.
        phy_lin : int
            Physical line where BCs are to be imposed.
        nodes_array : int
            Array with the nodal data and to be modified by BCs.
        P_x, P_y : float
            Load components in x and y directions.

    Returns
    -------
        nodes_array : int
            Array with the nodal data after imposing BCs according
            to SolidsPy.

    """
    lines = cells["line"]
    # Bounds contains data corresponding to the physical line.
    phy_line = cell_data["line"]["physical"]
    id_carga = [cont for cont in range(len(phy_line))
                if phy_line[cont] == phy_lin]
    nodes_carga = lines[id_carga]
    nodes_carga = nodes_carga.flatten()
    nodes_carga = list(set(nodes_carga))
    ncargas = len(nodes_carga)
    cargas = np.zeros((ncargas, 3))
    cargas[:, 0] = nodes_carga
    cargas[:, 1] = P_x/ncargas
    cargas[:, 2] = P_y/ncargas
    return cargas

