# -*- coding: utf-8 -*-
"""
Postprocessor subroutines
-------------------------

"""
from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import femutil as fe
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib import rcParams
import meshio
import signals as sig



rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['image.cmap'] = "YlGnBu_r"
rcParams['axes.axisbelow'] = True
rcParams['mathtext.fontset'] = "cm"

def mesh2tri(nodes, elements):
    """Generate a  matplotlib.tri.Triangulation object from the mesh
    
    Parameters
    ----------
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.
    
    Returns
    -------
    tri : Triangulation
        An unstructured triangular grid consisting of npoints points
        and ntri triangles.
    
    """
    x = nodes[:, 1]
    y = nodes[:, 2]
    triangs = []
    for el in elements:
        if el[1]==1:
            triangs.append(el[[3, 4, 5]])
            triangs.append(el[[5, 6, 3]])
        if el[1]==2:
            triangs.append(el[[3, 6, 8]])
            triangs.append(el[[6, 7, 8]])
            triangs.append(el[[6, 4, 7]])
            triangs.append(el[[7, 5, 8]])
        if el[1]==3:
            triangs.append(el[3:])
    
    tri = Triangulation(x, y, np.array(triangs))
    return tri    


def tri_plot(tri, field, title="", figtitle="", levels=12, savefigs=False,
             plt_type="contourf", filename="solution_plot.pdf"):
    
    if plt_type=="pcolor":
        disp_plot = plt.tripcolor
    elif plt_type=="contourf":
        disp_plot = plt.tricontourf

    plt.figure(figtitle)
    disp_plot(tri, field, levels, shading="gouraud")
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.axis("image")
    plt.grid()
    if savefigs:
        plt.savefig(filename)


def plot_disp(UC, nodes, elements, plt_type="contourf", levels=12,
               savefigs=False, title="Solution:"):
    """Plot the nodal displacement using a triangulation

    Parameters
    ----------
    UC : ndarray (float)
      Array with the displacements.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    """
    tri = mesh2tri(nodes, elements)
    tri_plot(tri, UC[:, 0], title=r'$u_x$',
             figtitle=title + "Horizontal displacement",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="ux_sol.pdf")
    tri_plot(tri, UC[:, 1], title=r'$u_y$',
             figtitle=title + "Vertical displacement",
             levels=levels, plt_type=plt_type, savefigs=savefigs,
             filename="uy_sol.pdf")


def grafmat(k):
    """Plot stiffness matrix sparsity

    Parameters
    ----------
    k : ndarray (int)
      Stiffness matrix of the system.

    """
    plt.figure("Stiffness matrix")
    plt.spy(k)
    plt.title("Stiffness matrix")
    plt.ylabel(r"$i$ index", size=10)
    plt.xlabel(r"$j$ index", size=10)


def complete_disp(IBC, nodes, UG):
    """
    Fill the displacement vectors with imposed and computed values.
    
    IBC : ndarray (int)
      IBC (Indicator of Boundary Conditions) indicates if the nodes
      has any type of boundary conditions applied to it.
    UG : ndarray (float)
      Array with the computed displacements.
    nodes : ndarray (float)
      Array with number and nodes coordinates:
      
    Returns
    -------
    UC : ndarray (float)
      Array with the displacements.

    """
    nn = nodes.shape[0]
    UC = np.zeros([nn, 2], dtype=np.float)
    for i in range(nn):
        for j in range(2):
            kk = IBC[i, j]
            if kk == -1:
                UC[i, j] = 0.0
            else:
                UC[i, j] = UG[kk]

    return UC


def eigvals(A, tol=1e-6):
    """Eigenvalues and eigenvectors for a 2x2 symmetric matrix/tensor
    
    Parameters
    ----------
    A : ndarray
        Symmetric matrix.
    tol : float (optional)
        Tolerance for considering a matrix diagonal.

    Returns
    -------
    eig1 : float
        First eigenvalue.
    eig2 : float
        Second eigenvalue.
    vec1 : ndarray
        First eigenvector.
    vec2 : ndarray
        Second eigenvector
    
    Examples
    --------
    
    >>> A = np.array([[5, 6],
    ...              [6, 9]])
    >>> eig1, eig2, vec1, vec2 =  eigvals(A)
    >>> np.allclose(eig1, 7 + 2*np.sqrt(10))
    True
    >>> np.allclose(eig2, 7 - 2*np.sqrt(10))
    True
    >>> np.allclose(vec1, np.array([-0.584710284663765, -0.8112421851755609]))
    True
    >>> np.allclose(vec2, np.array([-0.8112421851755609,0.584710284663765]))
    True
    
    """
    if np.abs(A).max() < tol:
        eig1 = 0.0
        eig2 = 0.0
        vec1 = np.array([np.NaN, np.NaN])
        vec2 = np.array([np.NaN, np.NaN])
    elif abs(A[0, 1])/np.abs(A).max() < tol:
        eig1 = A[0, 0]
        eig2 = A[1, 1]
        vec1 = np.array([1, 0])
        vec2 = np.array([0, 1])
    else:
        tr = A[0, 0] + A[1, 1]
        det = A[0, 0]*A[1, 1] - A[0, 1]**2
        eig1 = 0.5*(tr - np.sqrt(tr**2 - 4*det))
        eig2 = 0.5*(tr + np.sqrt(tr**2 - 4*det))
        vec1 = np.array([A[0, 0] - eig2, A[0, 1]])
        vec1 = vec1/np.sqrt(vec1[0]**2 + vec1[1]**2)
        vec2 = np.array([-vec1[1], vec1[0]])
    if abs(eig2) > abs(eig1):
        eig2, eig1 = eig1, eig2
        vec2, vec1 = vec1, vec2

    return eig1, eig2, vec1, vec2


def principal_dirs(field):
    """Compute the principal directions of a tensor field

    Parameters
    ----------
    field : ndarray
        Tensor field. The tensor is written as "vector" using
        Voigt notation.

    Returns
    -------
    eigs1 : ndarray
        Array with the first eigenvalues.
    eigs2 : ndarray
        Array with the second eigenvalues.
    vecs1 : ndarray
        Array with the first eigenvectors.
    vecs2 : ndarray
        Array with the Second eigenvector.

    """
    num = field.shape[0]
    eigs1 = np.empty((num))
    eigs2 = np.empty((num))
    vecs1 = np.empty((num, 2))
    vecs2 = np.empty((num, 2))
    A = np.zeros((2, 2))
    for cont, tensor in enumerate(field):
        A[0, 0] = tensor[0]
        A[1, 1] = tensor[1]
        A[0, 1] = tensor[2]
        eig1, eig2, vec1, vec2 = eigvals(A, tol=1e-6)
        eigs1[cont] = eig1
        eigs2[cont] = eig2
        vecs1[cont, :] = vec1
        vecs2[cont, :] = vec2

    return eigs1, eigs2, vecs1, vecs2

        
def energy(UG, KG):
    r"""
    Computes the potential energy for the current sln.

    Parameters
    ----------
    UG : ndarray (float)
      Array with the computed displacements.
    KG : ndarray (float)
      Global stiffness matrix.

    Returns
    -------
    EFE : scalar (float)
      Total energy in the system. :math:`-\frac{1}{2} U^T K U`

    """
    f = KG.dot(UG)
    EFE = -0.5*f.dot(UG)

    return EFE


def vtk_maker_chimba4(nodes, elements , fname , field= None):
    path = '../MESHES/VTKs/'
    ####
    npoints = nodes.shape[0]
    points = np.zeros((npoints, 3)) 
    points[:,0:2] = nodes[:, 1:3]
    nquads = elements.shape[0]
    quads = np.zeros((nquads, 4))
    quads[:] = elements[:, 3:].copy()
    quad_mesh = {
            'points': points,
            'cells': {'quad': quads}
            }
    point_data = field
    meshio.write(
        path + fname + ".vtk" ,
        quad_mesh["points"],
        quad_mesh["cells"],
        point_data=point_data)    
    return

def vtk_maker_chimba3(nodes, elements , fname , field= None):
    path = '../MESHES/VTKs/'
    ####
    npoints = nodes.shape[0]
    points = np.zeros((npoints, 3)) 
    points[:,0:2] = nodes[:, 1:3]
    nquads = elements.shape[0]
    quads = np.zeros((nquads, 3))
    quads[:] = elements[:, 3:].copy()
    quad_mesh = {
            'points': points,
            'cells': {'triangle': quads}
            }
    point_data = field
    meshio.write(
        path + fname + ".vtk" ,
        quad_mesh["points"],
        quad_mesh["cells"],
        point_data=point_data)    
    return


def vtk_maker_chimba9(nodes, elements , fname , field= None):
    path = '../MESHES/VTKs/'
    ####
    npoints = nodes.shape[0]
    points = np.zeros((npoints, 3)) 
    points[:,0:2] = nodes[:, 1:3]
    nquads = elements.shape[0]
    quads = np.zeros((nquads, 9))
    quads[:] = elements[:, 3:].copy()
    quad_mesh = {
            'points': points,
            'cells': {'quad9': quads}
            }
    point_data = field
    meshio.write(
        path + fname + ".vtk" ,
        quad_mesh["points"],
        quad_mesh["cells"],
        point_data=point_data)    
    return

def nodal_historyH(idnod , ninc , U , IBC , fname ):
    
    uh = np.zeros((ninc))
    idof = IBC[idnod , 0]
    uh[:] = U[idof , :]
    np.savetxt(fname , uh )
    
    return 


def nodal_historyV(idnod , ninc , U , IBC , fname ):
    
    uh = np.zeros((ninc))
    idof = IBC[idnod , 1]
    uh[:] = U[idof , :]
    np.savetxt(fname , uh )
    
    return

def sheets(idnod , ninc , U , IBC , fname , folder ):
    """
    Writes a file with the nodal history for a list of nodes
    stored in idnod
    
    idnod : ndarray (int)
      List with the nodal point names.
    ninc  : intger
      Integer indicating the number of increments.
    U     : ndarray (float)
      Array with the computed displacements.
    IBC : ndarray (integer)
      Array with the equation numbers
    fname: string.
      String with the file name.
    folder: string.
      String with the folder name.
    
      
    Returns
    -------
    nn : integer.
      Integer with the size of the idnod list.

    """
    nn = idnod.shape[0]
    uv = np.zeros([nn, ninc])
    for i in range(nn):
        idof = IBC[idnod[i] , 0]
        uv[i , :] = U[idof , :]
    np.savetxt(folder + fname + ".txt", uv)
    
    return nn

def PLOTsheets( fname , folder , dt , ninc , npts , dk):
    """
    Plots the time histories for a series of nodes with
    response stored in the file fname.
    
    fname: string.
      String with the file name.
    folder: string.
      String with the folder name.
    
    dt : Scalar (float)
      Time step.
    ninc  : intger (int)
      Integer indicating the number of increments.
    npts     : Integer (int)
      Integer with the number of poiints.
    dk : integer (integer)
      Scaling factor.
         
    """
    plt.figure(0)
    DATOS = np.loadtxt(folder+fname)
    signal=np.zeros([ninc , npts], dtype=float)
    k = 0
    for j in range(npts):
        for i in range(ninc):
            signal[i , k ] =  DATOS[j , i] + k/dk
        sig.grafsignalG(signal[: , k] , 'salida' , 'Displacement' , 'l' , 0.0 , 20.0 , dt , 0)
        k = k+1    
    return

def HV_history(lisres , nres , ninc , U , IBC , fnameh , fnamev ):
    
    uh = np.zeros((ninc))
    uv = np.zeros((ninc))
    for i in range(nres):
        idnod = lisres[i]
        idofU  = IBC[idnod , 0]
        idofV  = IBC[idnod , 1]
        uh[:] = uh[:] +  U[idofU , :]
        uv[:] = uv[:] +  U[idofV , :]
    np.savetxt(fnameh , uh )
    np.savetxt(fnamev , uv )
    
    return 
def plot_sheet(vals, tmax, amp_signal=30, amp_shift=50, ax=None):
    if ax==None:
        ax = plt.gca()
    nvals, ntimes = vals.shape
    time = np.linspace(0, tmax, ntimes)
    vert_shift = np.outer(np.linspace(0, 1, nvals), np.ones(ntimes))
    ax.plot(time, amp_signal*vals.T + amp_shift*vert_shift.T, alpha=0.8,
             color="gray", lw=1.2)


def plot_sheet_cmap(vals, tmax, amp_signal=30, amp_shift=50, ax=None,
                    cmap="Reds_r"):
    if ax==None:
        ax = plt.gca()
    cm = plt.cm.get_cmap(cmap)
    nvals, ntimes = vals.shape
    time = np.linspace(0, tmax, ntimes)
    vert_shift = amp_shift*np.linspace(0, 1, nvals)
    for cont in range(nvals):
        color = cm(cont/nvals)
        ax.plot(time, amp_signal*vals[cont,:] + vert_shift[cont],
                alpha=0.8, color=color, lw=1.2)


def plot_pcolor(vals, tmax, x_min, x_max, ax=None, cmap="Reds_r"):
    if ax==None:
        ax = plt.gca()
    nvals, ntimes = vals.shape
    Y, X = np.mgrid[x_min:x_max:nvals*1j,
                0:tmax:ntimes*1j]
    ax.pcolormesh(X, Y, vals, cmap=cmap)


def plot_surf(vals, tmax, x_min, x_max, ax=None, cmap="Reds_r",
              cstride=1, rstride=1):
    if ax==None:
        ax = plt.gca()
    nvals, ntimes = vals.shape
    Y, X = np.mgrid[x_min:x_max:nvals*1j,
                0:tmax:ntimes*1j]
    ax.plot_surface(X, Y, vals, cmap=cmap, cstride=cstride,
                    rstride=rstride)
    ax.set_zlim(-0.1, 0.1)
    ax.view_init(elev=45, azim=-60)


def respuesta(cells, cell_data, phy_lin):
    """Extracts the nodes located at the physical line
       phy_line

    Parameters
    ----------
        cell : dictionary
            Dictionary created by meshio with cells information.
        cell_data: dictionary
            Dictionary created by meshio with cells data information.
        phy_lin : int
            Physical line to print nodal histories.

    Returns
    -------
        nodes_carga : int
            Array with the nodal data corresponding to the physical
            line phy_line.

    """
    lines = cells["line"]
    phy_line = cell_data["line"]["physical"]
    id_carga = [cont for cont in range(len(phy_line))
                if phy_line[cont] == phy_lin]
    nodes_carga = lines[id_carga]
    nodes_carga = nodes_carga.flatten()
    nodes_carga = list(set(nodes_carga))
    nodes_carga.sort(reverse=False)
    
    return nodes_carga
#%%
if __name__ == "__main__":
    import doctest
    doctest.testmod()