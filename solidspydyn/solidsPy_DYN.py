# -*- coding: utf-8 -*-
"""
Computes the displacement solution for a finite element assembly
of 2D solids under point loads using as input easy-to-create
text files containing element, nodal, materials and loads data.
The input files are created out of a Gmsh (.msh) generated file
using the Python module ``meshio``.

Created by Juan Gomez and Nicolas Guarin-Zapata as part of the courses:

- IC0283 Computational Modeling
- IC0602 Introduction to the Finite Element Method

Which are part of the Civil Engineering Department at Universidad
EAFIT.

"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import preprocesor as pre
import postprocesor as pos
import assemutil as ass
import solutil as sol


def solidsPy_DYN(write_VTK = True , folder = None):
    """
    Run a complete workflow for a Finite Element Analysis

    Parameters
    ----------
    plot_contours : Bool (optional)
        Boolean variable to plot contours of the computed variables.
        By default it is True.
    compute_strains : Bool (optional)
        Boolean variable to compute Strains and Stresses at nodes.
        By default it is False.
    folder : string (optional)
        String with the path to the input files. If not provided
        it would ask for it in a pop-up window.

    Returns
    -------
    UC : ndarray (nnodes, 2)
        Displacements at nodes.
    E_nodes : ndarray (nnodes, 3), optional
        Strains at nodes. It is returned when `compute_strains` is True.
    S_nodes : ndarray (nnodes, 3), optional
        Stresses at nodes. It is returned when `compute_strains` is True.

    """
    if folder is None:
        folder, name = pre.initial_params()
    start_time = datetime.now()

    # Pre-processing
    inipar , nodes, mats, elements, loads = pre.readin(folder=folder)
    ninc , T , Tc , fc , dt , ac , theta  = pre.intparams(inipar)
    DME , IBC , neq = ass.DME(nodes, elements)
    print("Number of nodes: {}".format(nodes.shape[0]))
    print("Number of elements: {}".format(elements.shape[0]))
    print("Number of equations: {}".format(neq))

    # System assembly
    KG , MG , CG = ass.assembler(elements, mats, nodes, neq, DME)
    RHSG    = ass.loadasem(loads , IBC , neq , ninc , T , Tc , fc)
    KE = ass.effective(KG , MG , CG , ac )
    #%% INITIAL CONDITIONS                                                  
    U , V , A = sol.initial_conds(ninc , neq, RHSG , MG , KG , CG )
    print("Finished initial conditions....: {}".format(0))
    del KG
    start_time = datetime.now()

    # System solution
    U = sol.time_implicit(neq , ninc , dt , theta , ac , U , V , A ,  RHSG , MG ,  CG  , KE)
    end_time = datetime.now()
    print('Duration for system solution: {}'.format(end_time - start_time))

    # Post-processing
    start_time = datetime.now()
    if write_VTK:
        for i in range(0, ninc, 2):
            UC = pos.complete_disp(IBC, nodes, U[:, i])
            u_vec = np.zeros((nodes.shape[0], 3))
            u_vec[:, 0:2] = UC
            field = {"displacements": u_vec}
            pos.vtk_maker_chimba3(nodes, elements ,
                                 "scatteter_{}".format(i),
                                 field=field)
    
    end_time = datetime.now()
    print('Duration for post processing: {}'.format(end_time - start_time))
    print('Analysis terminated successfully!')
    return (U , folder , IBC , ninc , T)


if __name__ == '__main__':
    displacement , folder , IBC , ninc , T = solidsPy_DYN()
#    plt.show()