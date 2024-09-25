#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:22:15 2021

@author: rouwane
"""

import classIgaMesh as cim
import scipy.sparse as sps 
import scipy as sp
import scipy.sparse.linalg as splalg
import numpy as np 
from sksparse.cholmod import cholesky


def CorrelateGN(f,g,m,U0,reg,niter):
    U =   U0.copy()  # Copying the initial guess
    dvc  = cim.DVCEngine() 
    dvc.GnInit(f, m) 
    print('Roi gray-level dynamic')
    print(dvc.dyn)
    for ik in range(niter):
        H,b,res = dvc.ComputeGnMembers(f, g, m, U)
        H_LU = splalg.splu(H) 
        dU=H_LU.solve(b) 
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
        if err<1e-3:
            break
    # print('----')
    # print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
    return U,res   

def CorrelateStructureEA_Basis(dvc, f, g, m, nbippe, U0, Basis, niter):
    U =   U0.copy()  # Copying the initial guess
    H    = dvc.LHS(m, nbippe)
    
    lhs = Basis.T @ H @ Basis 
    
    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta    
    
    print('Starting iterations')  
    for ik in range(niter):
        b,ssd = dvc.RHS(g, m, nbippe, U )
        rhs = Basis.T @ b
        da = rhs/lhs 
        dU = Basis * da 
        U += dU 
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break
    return U, gres  

def CorrelateStructureZNEA_Basis(dvc, f, g, m, nbippe, U0, Basis, niter):
    U =   U0.copy()  # Copying the initial guess
    H    = dvc.LHS(m, nbippe)
    
    lhs = Basis.T @ H @ Basis 
    
    elemSizeX     =  m.xxsi[0][m.pp[0]+1]-m.xxsi[0][m.pp[0]]
    elemSizeEta   =  m.xxsi[1][m.pp[1]+1]-m.xxsi[1][m.pp[1]]
    elemSizeZeta  =  m.xxsi[2][m.pp[2]+1]-m.xxsi[2][m.pp[2]]

                                   
    volE    = elemSizeX*elemSizeEta*elemSizeZeta
    
    print('Starting iterations')  
    for ik in range(niter):
        b,ssde = dvc.RHS_ZN(g, m, nbippe, U )
        rhs = Basis.T @ b
        da = rhs/lhs 
        dU = Basis * da 
        U += dU 
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        grese =  np.sqrt(ssde/volE)/dvc.feDyn
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.max(grese)*100,err))
        if err<1e-3 :
            break
    return U, grese             


def CorrelateStructuredEA(dvc,f,g,m,nbippe,U0,reg,niter):
    """ Correlates with a structured B-spline mesh using 
    the elementary assembly """ 
    U =   U0.copy()  # Copying the initial guess
    print('LHS construction')
    H    = dvc.LHS(m, nbippe )
    if reg!=0: 
        L = m.Laplacian() 
    else : 
        L =  sps.csr_matrix(H.shape)
    # A = L.T.dot(L)     
    A = L
    # print('LU decomposition')
    # H_LU  = splalg.splu(H+reg*L)  
    
    print('LLt decomposition')
    H_LLt = cholesky(H+reg*A)
 
    
  
    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta
    
    print('Starting iterations')     
    for ik in range(niter):
        b,ssd = dvc.RHS(g, m, nbippe, U)
        # dU=H_LU.solve(b-reg*L.dot(U))
        dU = H_LLt.solve_A(b-reg*A.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres 

def CorrelateFE_EA(dvc, f, g, m, N,  U0, reg, niter ):
    """ Correlates with the finite element method """ 
    U =   U0.copy()  # Copying the initial guess    
    if reg!=0: 
        m.SetGaussIntegrationRule() 
        L = m.LaplacianEA() 
    else : 
        L =  sps.csr_matrix((m.ndof,m.ndof))
    A = L 
    print('Setting integration rule')
    m.SetVoxelTetraIntegrationRule(N) 
    print('Sub-voxel evaluation')
    dvc.GetImage_Mean_Std_FE(f, m)
    print('Correlation matrix assembly')
    H = dvc.LHS_FE_EA(m)
    print('LLt decomposition')
    H_LLt = cholesky(H+reg*A)
    # H_LLt = splalg.splu(H+reg*A)
    
    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = m.ComputeVolume() 
    
    print('Starting iterations')  
    for ik in range(niter):
        b,ssd = dvc.RHS_FE_EA(g, m, U)
        dU = H_LLt.solve_A(b-reg*A.dot(U))
        # dU = H_LLt.solve(b-reg*A.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres    
     
def CorrelateFE_EA_Basis(dvc, f, g, m, N, R,  U0 , niter ):
    """ Correlates with the finite element method """ 
    U =   U0.copy()  # Copying the initial guess    

    print('Setting integration rule')
    m.SetVoxelTetraIntegrationRule(N) 
    print('Sub-voxel evaluation')
    dvc.GetImage_Mean_Std_FE(f, m)
    print('Correlation matrix assembly')
    H = dvc.LHS_FE_EA(m)
    print('LLt decomposition')
    H_LLt = cholesky(R.T.dot(H.dot(R)))
    # H_LLt = splalg.splu(H+reg*A)
    
    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = m.ComputeVolume() 
    
    print('Starting iterations')  
    for ik in range(niter):
        b,ssd = dvc.RHS_FE_EA(g, m, U)
        dU = R.dot( H_LLt.solve_A(R.T.dot(b)))
        # dU = H_LLt.solve(b-reg*A.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres      

def CorrelateFE_EA_BasisScaling(dvc, f, g, m, N, U0, Basis, niter):
    """ Correlates with the finite element method """ 
    U =   U0.copy()  # Copying the initial guess    

    print('Setting integration rule')
    m.SetVoxelTetraIntegrationRule(N) 
    print('Sub-voxel evaluation')
    dvc.GetImage_Mean_Std_FE(f, m)
    print('Correlation matrix assembly')
    H = dvc.LHS_FE_EA(m)
    
    lhs = Basis.T @ H @ Basis 
    
    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = m.ComputeVolume() 
    
    print('Starting iterations')  
    for ik in range(niter):
        b,ssd = dvc.RHS_FE_EA(g, m, U)
        rhs = Basis.T @ b
        da = rhs/lhs 
        dU = Basis * da 
        U += dU 
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break
    return U, gres  

    
 
    

 
def CorrelateStructuredEA_CloseDisp(dvc,f,g,m,nbippe,U0,ipIndices,D1,D2,reg1,reg2,niter):
    U    = U0.copy()
    print('LHS construction')  
    # H    = dvc.LHS_thrsh(m, nbippe, ipIndices, pixEvalMethod )
    H    = dvc.LHS(m, nbippe)
    L = m.Laplacian() 
    L =  sp.sparse.csc_matrix( np.real(sp.linalg.sqrtm(L.toarray())) ) # Root square of the Laplacian operator 
    print('LU decomposition')
    # Lc   = reg1*D1.T.dot(L.dot(D1)) + reg2*D2
    Lc   = reg1*L.T.dot(D1.dot(L)) + reg2*D2
    H_LU = splalg.splu(H+Lc)
    
    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta   

    print('Starting iterations')     
    for ik in range(niter):
        # b,ssd = dvc.RHS_thrsh(g, m, nbippe, ipIndices, U, pixEvalMethod )
        b,ssd = dvc.RHS(g, m, nbippe, U )
        dU=H_LU.solve(b-Lc.dot(U)+reg2*D2*U0)
        # print(np.max(dU),np.min(dU))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres      


def CorrelateStructuredEA_LM(dvc, f,g,m, nbippe, U0, ipIndices, reg, niter):
    """
    Correlates with levenberg marquardt
    Solving (H+lambdaI)=b so that H is not singular 
    """
    U = U0.copy()
    ndof = len(U)
    print('LHS construction')
    # H    = dvc.LHS_thrsh(m, nbippe, ipIndices, pixEvalMethod )
    H    = dvc.LHS(m, nbippe)
    I = sps.identity(ndof,format='csr')
    print('LU decomposition')
    H_LU = splalg.splu(H+reg*I)

    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta    

    print('Starting iterations')     
    for ik in range(niter):
        # b,ssd = dvc.RHS_thrsh(g, m, nbippe, ipIndices, U, pixEvalMethod )
        b,ssd = dvc.RHS(g, m, nbippe, U)
        dU=H_LU.solve(b)
        # print(np.max(dU),np.min(dU))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres   


def CorrelateStructuredEA_DOF_Removal(dvc, f,g,m, nbippe, U0, dofKeeped, ipIndices, reg, niter):
    """ 
    Correlates with DOF removal 
    """
    U = np.zeros(U0.shape)
    U[dofKeeped] = U0[dofKeeped] 
    dU = np.zeros(U0.shape)
    
    # U = U0.copy()
    print('LHS construction')
    H    = dvc.LHS_thrsh(m, nbippe, ipIndices )
    if reg!=0: 
        L = m.Laplacian() 
    else : 
        L =  sps.csr_matrix(H.shape)
    ind = np.ix_(dofKeeped,dofKeeped)    
    Lcr = (reg*L)[:,ind[1][0]][ind[1][0],:] # Reduced matrix 
    # Lc  = reg*L 
    Hr = H[:,ind[1][0]][ind[1][0],:]
    print('LU decomposition')
    H_LU = splalg.splu(Hr+Lcr) 
    # H_LU = splalg.splu(H+Lc)  


    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta
    
    print('Starting iterations')     
    for ik in range(niter):
        b,ssd = dvc.RHS_thrsh(g, m, nbippe, ipIndices, U )
        # print(np.max(b), np.min(b), ssd)
        br = b[dofKeeped]
        dU[dofKeeped]=H_LU.solve(br-Lcr.dot(U[dofKeeped]))
        # dU=H_LU.solve(b)
        # print(np.max(dU),np.min(dU))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres   
    

def CorrelateStructuredEA_Disk_Tikhonov(dvc, f,g,m, nbippe, U0, D1, D2, Lsr, reg1, reg2, niter):
    """ Correlates with a structured B-spline mesh using 
    the elementary assembly """ 
    U =   U0.copy()  # Copying the initial guess
    print('LHS construction')
    H    = dvc.LHS(m, nbippe)
    # L = m.Laplacian() 
    # LL = L.T.dot(L)     
    print('LU decomposition')
    # Lc = reg1*D1.T.dot(L.dot(D1))   + reg2*D2.T.dot(L.dot(D2))
    Lc = reg1*Lsr.T.dot(D1.dot(Lsr))  + reg2*Lsr.T.dot(D2.dot(Lsr))
    H_LU = splalg.splu(H + Lc  )   
    
  
    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta
    
    print('Starting iterations')     
    for ik in range(niter):
        b,ssd = dvc.RHS(g, m, nbippe, U )
        dU=H_LU.solve(b-Lc.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres     

def CorrelationStructuredEA_Tikhonov(dvc, f,g,m, nbippe, U0, A, D, reg, niter  ):
    U =   U0.copy()  # Copying the initial guess
    print('LHS construction')
    H    = dvc.LHS(m, nbippe)
    Lc = reg*A.T.dot(D.dot(A)) 
    # Lc   = reg*D.T.dot(A.dot(D))
    print('LU decomposition')
    H_LU = splalg.splu(H+Lc) 
    
    roiSizeX     =  m.xxsi[0][-1]-m.xxsi[0][0]
    roiSizeEta   =  m.xxsi[1][-1]-m.xxsi[1][0]  
    roiSizeZeta  =  m.xxsi[2][-1]-m.xxsi[2][0]   

    fdynROI    = np.max(dvc.fip)-np.min(dvc.fip)
    volROI     = roiSizeX*roiSizeEta*roiSizeZeta
    
    print('Starting iterations')     
    for ik in range(niter):
        b,ssd = dvc.RHS(g, m, nbippe, U )
        dU=H_LU.solve(b-Lc.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        gres = np.sqrt(ssd/volROI)/fdynROI # global residual 
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,gres*100,err))
        if err<1e-3 :
            break         
        
    return U, gres     
    
    
    

def CorrelateStructuredZNEA(dvc,f,g,m,nbippe,U0,reg,niter ):
    """ Correlates with a structured B-spline mesh using 
    the elementary assembly """ 
    U =   U0.copy()  # Copying the initial guess
    print('LHS construction')
    H    = dvc.LHS(m, nbippe)
    if reg!=0: 
        L = m.Laplacian() 
    else : 
        L =  sps.csr_matrix(H.shape)
    # LL = L.T.dot(L)     
    print('LU decomposition')
    H_LU = splalg.splu(H+reg*L)  
    
  
    elemSizeX     =  m.xxsi[0][m.pp[0]+1]-m.xxsi[0][m.pp[0]]
    elemSizeEta   =  m.xxsi[1][m.pp[1]+1]-m.xxsi[1][m.pp[1]]
    elemSizeZeta  =  m.xxsi[2][m.pp[2]+1]-m.xxsi[2][m.pp[2]]

                                   
    volE    = elemSizeX*elemSizeEta*elemSizeZeta
    
    print('Starting iterations')
    for ik in range(niter):
        b,ssde = dvc.RHS_ZN(g, m, nbippe, U )
        dU=H_LU.solve(b-reg*L.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        grese =  np.sqrt(ssde/volE)/dvc.feDyn
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.max(grese)*100,err))
        if err<1e-3 :
            break         
        
    return U, np.max(grese)

        
 
def Correlate(f,g,m,U0,reg,niter ):
    U =   U0.copy()  # Copying the initial guess
    dvc  = cim.DVCEngine() 
    H    = dvc.ComputeLHS(f, m )
    print('Roi gray-level dynamic')
    print(dvc.dyn)
    if reg!=0: 
        L = m.HomogeneousGaussLaplacian() 
    else : 
        L =  sps.csr_matrix(H.shape)
    # LL = L.T.dot(L)    
    # H_LU = splalg.splu(H) 
    H_LU=splalg.splu(H+reg*L)
    # H_LU = splalg.splu(H+reg*LL)
    for ik in range(niter):
        [b,res]=dvc.ComputeRHS(g,m,U)
        # dU=H_LU.solve(b)
        dU=H_LU.solve(b-reg*L.dot(U))
        # dU=H_LU.solve(b-reg*LL.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
        if err<1e-3:
            break
    # print('----')
    # print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
    return U,res,  np.std(res)/dvc.dyn


def CorrelateEquilibriumGapGN(f,g,m,U0,L,DL,regL,K,DK,regK,niter ):
    U =   U0.copy()  # Copying the initial guess
    dvc  = cim.DVCEngine() 
    dvc.GnInit(f, m) 
    print('Roi gray-level dynamic')
    print(dvc.dyn)
    Lc1  =  regK*K.T.dot(DK.dot(K))
    Lc2  =  regL*L.T.dot(DL.dot(L))
    Lc = Lc1 + Lc2 
    for ik in range(niter):
        H,b,res = dvc.ComputeGnMembers(f, g, m, U )
        H_LU=splalg.splu(H+Lc)
        # dU=H_LU.solve(b)
        dU=H_LU.solve(b-Lc.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
        if err<1e-3:
            break
    # print('----')
    # print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
    return U,res 
    

def CorrelateEquilibriumGap(f,g,m,U0,L,DL,regL,K,DK,regK,niter ):
    Uhist =[U0] 
    U =   U0.copy()  # Copying the initial guess
    dvc  = cim.DVCEngine() 
    H    = dvc.ComputeLHS(f, m )
    print('Roi gray-level dynamic')
    print(dvc.dyn)
    Lc1  =  regK*K.T.dot(DK.dot(K))
    Lc2  =  regL*L.T.dot(DL.dot(L))
    Lc = Lc1 + Lc2 
    H_LU=splalg.splu(H+Lc)
    for ik in range(niter):
        [b,res]=dvc.ComputeRHS(g,m,U)
        # dU=H_LU.solve(b)
        dU=H_LU.solve(b-Lc.dot(U))
        U+=dU
        err=np.linalg.norm(dU)/np.linalg.norm(U)
        print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
        Uhist.append(U)
        if err<1e-3:
            break
    # print('----')
    # print("Iter # %2d | disc/dyn=%2.2f %% | dU/U=%1.2e" % (ik+1,np.std(res)/dvc.dyn*100,err))
    return U,res, Uhist