# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function
import meshio
import numpy as np
import solidspy.preprocesor as msh
#
points, cells, point_data, cell_data, field_data = \
    meshio.read("bimaterial.msh")
#
nodes_array    = msh.node_writer(points , point_data)
nf , els1_array = msh.ele_writer(cells , cell_data , \
                  "triangle" , 100 , 3 , 0 , 0)
nini = nf
nf , els2_array = msh.ele_writer(cells , cell_data , \
                  "triangle" , 200 , 3 , 1 , nini)
#
nodes_array = msh.boundary_conditions(cells , cell_data ,\
                  400 , nodes_array , -1 , -1)
els_array = np.append(els1_array, els2_array , axis = 0)
#
np.savetxt("eles.txt" , els_array   , fmt="%d")
np.savetxt("nodes.txt", nodes_array , fmt=("%d", "%.4f",\
           "%.4f", "%d", "%d"))