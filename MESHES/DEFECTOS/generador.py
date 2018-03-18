# -*- coding: utf-8 -*-
"""
Template to generate the input files for the FEM code solids_ISO.
The script uses module meshio.py to read a GMSH mesh and produce
text files nodes.txt, eles.txt , mater.txt and loads.txt

@authors: Juan Gomez
         Nicolas Guarin-Zapata
"""
from __future__ import division, print_function
from os import sys
sys.path.append('../../MAIN/')
import meshio
import solidspy.preprocesor as msh
import numpy as np
import postprocesor as pos

#
points, cells, point_data, cell_data, field_data = \
    meshio.read("prueba.msh")
#
nodes_array    = msh.node_writer(points , point_data)
nf , els1_array = msh.ele_writer(cells , cell_data , "triangle" , 1000 , 3 , 0 , 0)
nini = nf
nf , els2_array = msh.ele_writer(cells , cell_data , "triangle" , 2000 , 3 , 1 , nini)
els_array =np.append(els1_array , els2_array , axis = 0)
#
nodes_array = msh.boundary_conditions(cells , cell_data , 100 , nodes_array , -1 , -1)

cargas      = msh.loading(cells , cell_data , 300 , 0.0 , 10.0)
#
salida = pos.respuesta(cells, cell_data, 500)
#
np.savetxt("eles.txt" , els_array , fmt="%d")
np.savetxt("loads.txt", cargas, fmt=("%d", "%.6f", "%.6f"))
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f", "%.4f", "%d", "%d"))
np.savetxt("salida.txt" , salida      , fmt="%d")