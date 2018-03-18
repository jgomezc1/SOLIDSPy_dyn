"""
solutil.py
----------

Utilities for solution of FEM systems

"""
from __future__ import division, print_function
import numpy as np
from datetime import datetime
from numpy import ndarray
from numpy.linalg import solve
from scipy.sparse.csr import csr_matrix
from scipy.sparse.linalg import spsolve , splu


def static_sol(mat, rhs):
    """Solve a static problem [mat]{u_sol} = {rhs}
    """
    if type(mat) is csr_matrix:
        u_sol = spsolve(mat, rhs)
    elif type(mat) is ndarray:
        u_sol = solve(mat, rhs)
    else:
        raise Exception("Not supported matrix storage scheme!")

    return u_sol


def initial_conds(ninc , neq, RHSG , MG , KG , CG ):
    """
     Currently homogeneous initial conditions only
    """
    Up  =np.zeros([neq  ],dtype=np.float)
    Vp  =np.zeros([neq  ],dtype=np.float)
    Ap  =np.zeros([neq  ],dtype=np.float)
    U   =np.zeros([neq,ninc],dtype=np.float)
    V   =np.zeros([neq,2],dtype=np.float)
    A   =np.zeros([neq,2],dtype=np.float)
    F  =np.zeros([neq],dtype=np.float)
    FE  =np.zeros([neq],dtype=np.float)
    F = RHSG[:, 0]
    FS = KG.dot(Up)
    FD = CG.dot(Vp)
    FE = F - FD - FS

    Ap = static_sol(MG , FE)   
    A[:, 0] = Ap

        
        
    return U , V , A

def time_implicit(icount , m , dt , theta , ass , U , V , A ,  F , MG ,  CG  , KE):
    """Uses the Wilson theta method to perform
    implicit time integration.
    Parameters
    ----------
    icount : Integer. Number of equations
    m      : Integer. Number of time increments
    dt     : Float. Time step.
    ass    : Float array. Integration constants.
    theta  : Float. Integration parameter.
    U, V, A: Float arrays with the nodal point displacement,
             velocity and acceleration. Must be passed with the
             initial conditions.
    F      : Float array with the point load amplitudes.
    MG     : Float array. Mass matrix.
    CG     : Float array. Dammping matrix.
    KG     : Float array. Stiffness matrix.
    
    """
#
#   Integration parameters
#
    a_0 = ass[0]
    a_1 = ass[1]
    a_2 = ass[2]
    a_3 = ass[3]
    a_4 = ass[4]
    a_5 = ass[5]
    a_6 = ass[6]
    a_7 = ass[7]
    a_8 = ass[8] 
#
    LU  = splu(KE)
#
    for k in range(m-1): 
        start_time = datetime.now()
        VV  =np.zeros([icount  ],dtype=np.float)
        AA  =np.zeros([icount  ],dtype=np.float)
        RHS =np.zeros([icount  ],dtype=np.float)
        FE  =np.zeros([icount  ],dtype=np.float)
        for i in range(0,icount):
            RHS[i]=F[i,k]+theta*(F[i,1+k]- F[i,k]) 
            AA[i] = (a_0*U[i , k]+a_2*V[i , 0]+2.0*A[i , 0])
            VV[i] = (a_1*U[i , k]+2*V[i , 0]  +a_3*A[i , 0])        
        FI=MG.dot(AA)
        FD=CG.dot(VV)
        FE = RHS + FI + FD
#
        Up = LU.solve(FE)
        end_time = datetime.now()
        print("Increment number....: {}".format(k+1))
        print('Duration for this increment....: {}'.format(end_time - start_time))
        for i in range(0,icount):            
            A[i, 1]=-a_4*U[i,k]+a_5*V[i , 0]+a_6*A[i , 0]+a_4*Up[i]
            V[i, 1]= V[i,0]+a_7*A[i , 0]+a_7*A[i,1]
            U[i,k+1]=U[i,k]+dt*V[i,0]+2.0*a_8*A[i,0]+a_8*A[i,1]
        A[: , 0] = A[: , 1]
        V[: , 0] = V[: , 1]
    return U