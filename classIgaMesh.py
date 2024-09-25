#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 18:53:09 2020

@author: rouwane
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import nurbs as nb 
import scipy.sparse as sps 
import scipy as sp 
import scipy.sparse.linalg as splalg
from sksparse.cholmod import cholesky
from scipy.spatial import Delaunay
from evtk.hl import VtkFile,VtkPolyLine, VtkUnstructuredGrid, VtkStructuredGrid,\
VtkVertex, linesToVTK, VtkLine 
import matplotlib.collections as cols
import wide_product as wp 
import imageRoutines as ir 
import assemblyRoutines as ar 
import toolsRoutines as tr 
# import assembly_Bspline_3d as arb 
import sys 
import tools 
import gmsh 
import meshio 
import pymesh 
import sys


#%% 3D Image 
class Image3D: 
    def __init__(self):
        self.xmin = None 
    def Load(self,pix):
        if isinstance(pix,np.ndarray):
            # self.pix = pix 
            self.pix = pix.astype('double')
        else:
            raise ValueError('error while reading input data')
    def SetBoundaries(self,xmin,xmax,ymin,ymax,zmin,zmax):
        " Sets the origin of the image and the pixel size """
  
        self.xmin = xmin 
        self.ymin = ymin 
        self.xmax = xmax 
        self.ymax = ymax 
        self.zmin = zmin 
        self.zmax = zmax 
        
        # Intrinsic code projector defined as follow (see function GetGrayLevelOfPoint(x,y,z) )
        # It would have been possible to perform all computation in the (i,j,k) space 
        # x = k 
        # y = i 
        # z = -j 
        
        self.dx   = (xmax-xmin)/self.pix.shape[2] # pixel size in x direction 
        self.dy   = (ymax-ymin)/self.pix.shape[0] # pixel size in y direction  
        self.dz   = (zmax-zmin)/self.pix.shape[1] # pixel size in z direction  
        
    
        # Parameters of the L2 lumped projection (C.V Verhoosel )
    
        self.lsknotX  = 0
        self.lsknotY  = 0
        self.lsknotZ  = 0
        self.lsDegree = 0 
        self.lsNelems = 0 
        self.lsc = None 
        
        
        # Parameters of the trilinear interpolation (Knot vector of a first order B-spline volume)
        # As many elements as voxels 
        ne_xi      =  self.pix.shape[2]
        ne_eta     =  self.pix.shape[0]
        ne_zeta    =  self.pix.shape[1]
        
        # There npixel -1 elements because we start from the pixel center 

        self.knotXi1   = np.r_[self.xmin +self.dx/2, 
                         np.linspace(self.xmin +self.dx/2, self.xmax-self.dx/2, ne_xi  ),
                         self.xmax-self.dx/2] 
        self.knotEta1  = np.r_[self.ymin +self.dy/2, 
                         np.linspace(self.ymin +self.dy/2, self.ymax-self.dy/2, ne_eta ),
                         self.ymax-self.dy/2] 
        self.knotZeta1 = np.r_[self.zmin +self.dz/2, 
                         np.linspace(self.zmin +self.dz/2, self.zmax-self.dz/2, ne_zeta ),
                         self.zmax-self.dz/2] 
        
        
        # Pixel centers coordinates used for the Cardinal B-spline representation 
        
        self.pxc = self.xmin + self.dx/2 + np.arange(ne_xi)*self.dx 
        self.pyc = self.ymin + self.dy/2 + np.arange(ne_eta)*self.dy  
        self.pzc = self.zmin + self.dz/2 + np.arange(ne_zeta)*self.dz 
        
   
    def GetGrayLevelOfPoint(self, x,y,z):
        u = y
        v = self.zmin+self.zmax - z
        w = x 
        i = np.floor((u-self.ymin)/self.dy).astype('int64')
        j = np.floor((v-self.zmin)/self.dz).astype('int64')
        k = np.floor((w-self.xmin)/self.dx).astype('int64')        
        return self.pix[i,j,k] 
    
    def EvaluateL2ProjLumpledV(self,x,y,z):
        """ Method 1 : compute manually the products (very slow in Python)"""
        #phi = nb.Get3dBasisFunctionsAtPtsWd(x,y,z,self.lsknotX,self.lsknotY,self.lsknotZ,self.lsDegree[0],self.lsDegree[1],self.lsDegree[2])
        #ls = phi.dot(self.lsc)
        """ Method 2 : Compute 1d basis functions and perform Khatri–Rao product (package wide_product needed) 
            PAY ATTENTION FOR LARGE DATA SETS 
        """
        phi_xi   = nb.global_basisfunsWd(self.lsDegree[0], self.lsknotX, x)
        phi_eta  = nb.global_basisfunsWd(self.lsDegree[1], self.lsknotY, y)
        phi_zeta = nb.global_basisfunsWd(self.lsDegree[2], self.lsknotZ, z)
        phiE_phiX   = wp.wide_product(phi_eta, phi_xi)  
        phi         = wp.wide_product(phi_zeta, phiE_phiX)
        ls = phi.dot(self.lsc) 
        return ls 
    def EvaluateL2ProjLumpedStructuredV(self,x,y,z):
        phi_xi   = nb.global_basisfunsWd(self.lsDegree[0], self.lsknotX, x)
        phi_eta  = nb.global_basisfunsWd(self.lsDegree[1], self.lsknotY, y)
        phi_zeta = nb.global_basisfunsWd(self.lsDegree[2], self.lsknotZ, z)
        phiEta_phiXi   =  sps.kron(phi_eta,phi_xi, 'csc')
        phi            =  sps.kron(phi_zeta, phiEta_phiXi ,'csc')
        ls = phi.dot(self.lsc)
        return ls 
    def EvaluateL2ProjLumpled(self,x,y,z):
        """
        C++ 
        x,y,z is an unstructured cloud of points 
        loop over all points 
        """
        if self.lsc is None :
            raise ValueError('Level set coefficients must be first computed') 
        x[x<self.xmin] = self.xmin 
        x[x>self.xmax] = self.xmax 
        y[y<self.ymin] = self.ymin 
        y[y>self.ymax] = self.ymax 
        z[z<self.zmin] = self.zmin 
        z[z>self.zmax] = self.zmax 
        ls = tr.EvaluateBspline3D(x,y,z,self.lsknotX, self.lsknotY, self.lsknotZ,
                                   self.lsDegree[0], self.lsDegree[1], self.lsDegree[2],
                                   self.lsc,len(x)) 
        return ls 
    def EvaluateL2ProjLumpedAndGradient(self,x,y,z):
        """ 
        C++
        x,y,z is an unstructured cloud of points 
        loop over all points 
        """
        if self.lsc is None :
            raise ValueError('Level set coefficients must be first computed') 
        x[x<self.xmin] = self.xmin 
        x[x>self.xmax] = self.xmax 
        y[y<self.ymin] = self.ymin 
        y[y>self.ymax] = self.ymax 
        z[z<self.zmin] = self.zmin 
        z[z>self.zmax] = self.zmax          
            
        r = tr.EvaluateBsplineAndGradient3D(x,y,z,self.lsknotX, self.lsknotY, self.lsknotZ,
                                   self.lsDegree[0], self.lsDegree[1], self.lsDegree[2],
                                   self.lsc,len(x),len(x),len(x),len(x))
        ls  = r[0]
        lsx = r[1]
        lsy = r[2]
        lsz = r[3]
        return ls,lsx,lsy,lsz  
    
    def EvaluateL2ProjLumpedStructured(self,x,y,z):
        """ 
        C++
        x,y,z is a structured grid 
        loop over all points 
        """
        if self.lsc is None :
            raise ValueError('Level set coefficients must be first computed') 
        lo = len(x)*len(y)*len(z) 
        ls = tr.EvaluateBsplineStructured3D(x,y,z,self.lsknotX, self.lsknotY, self.lsknotZ,
                                   self.lsDegree[0], self.lsDegree[1], self.lsDegree[2],
                                   self.lsc,lo)
        return ls 
    def EvaluateL2ProjLumpedAndGradientStructured(self,x,y,z):
        """ 
        x,y,z is a structured grid 
        loop over all points 
        """
        if self.lsc is None :
            raise ValueError('Level set coefficients must be first computed') 
        lo = len(x)*len(y)*len(z)
        r = tr.EvaluateBsplineAndGradientStructured3D(x,y,z,self.lsknotX, self.lsknotY, self.lsknotZ,
                                   self.lsDegree[0], self.lsDegree[1], self.lsDegree[2],
                                   self.lsc,lo,lo,lo,lo)
        ls  = r[0]
        lsx = r[1]
        lsy = r[2]
        lsz = r[3]
        return ls,lsx,lsy,lsz  
    
    def EvaluateTrilinearInterpolation(self,x,y,z):
        """
        x,y,z is an unstructured cloud of points 
        Performs tri-linear interpolation by considering 
        an open knot vector and linear B-splines 
        the tri-variate functions are obtained by tensor-product 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        
        x[x<self.xmin+self.dx/2] = self.xmin +self.dx/2 
        x[x>self.xmax-self.dx/2] = self.xmax -self.dx/2
        y[y<self.ymin+self.dy/2] = self.ymin +self.dy/2 
        y[y>self.ymax-self.dy/2] = self.ymax -self.dy/2
        z[z<self.zmin+self.dz/2] = self.zmin +self.dz/2  
        z[z>self.zmax-self.dz/2] = self.zmax -self.dz/2
            
        
        v = ir.EvaluateTrilinearInterpolation(self.pix, self.knotXi1, self.knotEta1, self.knotZeta1, x,y,z, len(x) )
        
        return v  
    
    def EvaluateTrilinearInterpolationAndGradient(self,x,y,z):
        """
        x,y,z is an unstructured cloud of points 
        Performs tri-linear interpolation by considering 
        an open knot vector and linear B-splines 
        the tri-variate functions are obtained by tensor-product 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        
        x[x<self.xmin+self.dx/2] = self.xmin +self.dx/2 
        x[x>self.xmax-self.dx/2] = self.xmax -self.dx/2
        y[y<self.ymin+self.dy/2] = self.ymin +self.dy/2 
        y[y>self.ymax-self.dy/2] = self.ymax -self.dy/2
        z[z<self.zmin+self.dz/2] = self.zmin +self.dz/2  
        z[z>self.zmax-self.dz/2] = self.zmax -self.dz/2
        
 
        lo = len(x)
        r = ir.EvaluateTrilinearInterpolationAndGradient(self.pix, self.knotXi1, self.knotEta1, self.knotZeta1, x, y, z,lo,lo,lo,lo) 
        return r[0],r[1],r[2],r[3]        
        
 
    def EvaluateTrilinearInterpolationAndGradientStructured(self,x,y,z):
        """
        x,y,z is a structured cloud of points 
        Performs tri-linear interpolation by considering 
        an open knot vector and linear B-splines 
        the tri-variate functions are obtained by tensor-product 
        returns the interpolated voxels and 3d gradient vector 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        
        x[x<self.xmin+self.dx/2] = self.xmin +self.dx/2 
        x[x>self.xmax-self.dx/2] = self.xmax -self.dx/2
        y[y<self.ymin+self.dy/2] = self.ymin +self.dy/2 
        y[y>self.ymax-self.dy/2] = self.ymax -self.dy/2
        z[z<self.zmin+self.dz/2] = self.zmin +self.dz/2  
        z[z>self.zmax-self.dz/2] = self.zmax -self.dz/2
        
        lo = len(x)*len(y)*len(z)
        r = ir.EvaluateTrilinearInterpolationAndGradientStructured(self.pix, self.knotXi1, self.knotEta1, self.knotZeta1, x, y, z,lo,lo,lo,lo) 
        # Returns v, dvdx,dvdy,dvdz  
        return r[0],r[1],r[2],r[3]
    
    def EvaluateCardinalBspline(self,x,y,z,degree):
        """
        x,y,z is an unstructured cloud of points 
        Evaluate the cardinal B-spline representation of the image 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' ) 
         
        x[x<self.xmin+3*self.dx/2] = self.xmin +3*self.dx/2 
        x[x>self.xmax-3*self.dx/2] = self.xmax -3*self.dx/2
        y[y<self.ymin+3*self.dy/2] = self.ymin +3*self.dy/2 
        y[y>self.ymax-3*self.dy/2] = self.ymax -3*self.dy/2
        z[z<self.zmin+3*self.dz/2] = self.zmin +3*self.dz/2  
        z[z>self.zmax-3*self.dz/2] = self.zmax -3*self.dz/2
 
 
        lo = len(x)
        if degree ==2:
            return  ir.EvaluateCardBspline2(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo)
        elif degree==3:
            return  ir.EvaluateCardBspline3(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo) 
        else: 
            raise ValueError('Only quadratic and cubics- for tri-linear use function  EvaluateTrilinearInterpolation') 
            
    def EvaluateCardinalBsplineAndGradient(self,x,y,z,degree):
        """
        x,y,z is a structured cloud of points 
        Evaluate tthe cardinal B-spline representation of the image by exploiting the tensor-product nature 
        of B-splines 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' ) 
         
        x[x<self.xmin+3*self.dx/2] = self.xmin +3*self.dx/2 
        x[x>self.xmax-3*self.dx/2] = self.xmax -3*self.dx/2
        y[y<self.ymin+3*self.dy/2] = self.ymin +3*self.dy/2 
        y[y>self.ymax-3*self.dy/2] = self.ymax -3*self.dy/2
        z[z<self.zmin+3*self.dz/2] = self.zmin +3*self.dz/2  
        z[z>self.zmax-3*self.dz/2] = self.zmax -3*self.dz/2
 
        
        lo = len(x) 
        
        if degree ==2:
            r =  ir.EvaluateCardBsplineAndGradient2(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo,lo,lo,lo) 
            return r[0],r[1],r[2],r[3]   
        elif degree==3:
            r =  ir.EvaluateCardBsplineAndGradient3(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo,lo,lo,lo) 
            return r[0],r[1],r[2],r[3]   
        else: 
            raise ValueError('Only quadratic and cubics- for tri-linear use function  EvaluateTrilinearInterpolation')             

    def EvaluateCardinalBsplineAndGradientStructured(self,x,y,z,degree):
        """
        x,y,z is a structured cloud of points 
        Evaluate tthe cardinal B-spline representation of the image by exploiting the tensor-product nature 
        of B-splines 
        """
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' ) 
         
        x[x<self.xmin+3*self.dx/2] = self.xmin +3*self.dx/2 
        x[x>self.xmax-3*self.dx/2] = self.xmax -3*self.dx/2
        y[y<self.ymin+3*self.dy/2] = self.ymin +3*self.dy/2 
        y[y>self.ymax-3*self.dy/2] = self.ymax -3*self.dy/2
        z[z<self.zmin+3*self.dz/2] = self.zmin +3*self.dz/2  
        z[z>self.zmax-3*self.dz/2] = self.zmax -3*self.dz/2
 
        lo = len(x)*len(y)*len(z)
        
        if degree ==2:
            r =  ir.EvaluateCardBsplineAndGradient2Structured(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo,lo,lo,lo) 
            return r[0],r[1],r[2],r[3]   
        elif degree==3:
            r =  ir.EvaluateCardBsplineAndGradient3Structured(self.pix, self.xmin, self.ymin, self.zmin, self.dx, self.dy,self.dz,
                                             self.pxc,self.pyc,self.pzc,x,y,z,lo,lo,lo,lo) 
            return r[0],r[1],r[2],r[3]   
        else: 
            raise ValueError('Only quadratic and cubics- for tri-linear use function  EvaluateTrilinearInterpolation') 
 
    def LoadL2ProjLumpedParameters(self,degree,lsc):
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        self.lsDegree = degree
        oneX = np.ones(self.lsDegree[0])
        oneY = np.ones(self.lsDegree[1])
        oneZ = np.ones(self.lsDegree[2])
        
        ne_xi      =  self.pix.shape[2]
        ne_eta     =  self.pix.shape[0]
        ne_zeta    =  self.pix.shape[1]
        self.lsNelems = [ne_xi,ne_eta,ne_zeta]
        
        
        e_xi   = np.linspace(self.xmin,self.xmax, ne_xi  +1 )
        e_eta  = np.linspace(self.ymin,self.ymax, ne_eta +1 )
        e_zeta = np.linspace(self.zmin,self.zmax, ne_zeta+1 )
        
        self.lsknotX  = np.r_[ oneX*self.xmin, e_xi   , oneX*self.xmax ]
        self.lsknotY  = np.r_[ oneY*self.ymin, e_eta  , oneY*self.ymax ]
        self.lsknotZ  = np.r_[ oneZ*self.zmin, e_zeta , oneZ*self.zmax ]   
        self.lsc = lsc 
    

        
    def ComputeL2ProjLumpedCoefficients(self,degree):
        """
        Image-based goal-oriented adaptive isogeometric analysis
        with application to the micro-mechanical modeling of trabecular bone 
        C.V Verhoosel et al
        Computes the gray-level control points 
        There are as many control points as voxels  
        Computes the L2 lumped projection by sum factorization 
        thanks to the tensor product nature of B-splines only integrals of univariate 
        basis functions is computed 
        """
        
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        self.lsDegree = degree
        oneX = np.ones(self.lsDegree[0])
        oneY = np.ones(self.lsDegree[1])
        oneZ = np.ones(self.lsDegree[2])
        
        ne_xi      =  self.pix.shape[2]
        ne_eta     =  self.pix.shape[0]
        ne_zeta    =  self.pix.shape[1]
        self.lsNelems = [ne_xi,ne_eta,ne_zeta]
        
        
        e_xi   = np.linspace(self.xmin,self.xmax, ne_xi  +1 )
        e_eta  = np.linspace(self.ymin,self.ymax, ne_eta +1 )
        e_zeta = np.linspace(self.zmin,self.zmax, ne_zeta+1 )
        
        self.lsknotX  = np.r_[ oneX*self.xmin, e_xi   , oneX*self.xmax ]
        self.lsknotY  = np.r_[ oneY*self.ymin, e_eta  , oneY*self.ymax ]
        self.lsknotZ  = np.r_[ oneZ*self.zmin, e_zeta , oneZ*self.zmax ] 
        
        nbf_xi = ne_xi +self.lsDegree[0]
        nbf_eta = ne_eta + self.lsDegree[1]
        nbf_zeta = ne_zeta + self.lsDegree[2]
        
        s_lsc = nbf_xi*nbf_eta*nbf_zeta 
 
        self.lsc = ir.ComputeL2ProjLumpedCoefficientsSumFact(self.pix, self.xmin, self.ymin, self.zmin, self.zmax,
                                             self.dx, self.dy, self.dz, self.lsDegree[0], self.lsDegree[1], self.lsDegree[2],
                                             self.lsknotX, self.lsknotY, self.lsknotZ, s_lsc)

    
        
    def ComputeL2ProjLumpedCoefficientsVectorized(self,degree):
        """
        Image-based goal-oriented adaptive isogeometric analysis
        with application to the micro-mechanical modeling of trabecular bone 
        Computes the gray-level control points 
        There are as many control points as voxels  
        Vectorized Python Function 
        PAY ATTENTION WITH THE SPARSE KRONECKER FOR LARGE DATA SETS 
        """
        
        if self.xmin == None :
            raise ValueError( 'Call the function SetBoundaries first' )
        self.lsDegree = degree
        oneX = np.ones(degree[0])
        oneY = np.ones(degree[1])
        oneZ = np.ones(degree[2])
        
        ne_xi      =  self.pix.shape[2]
        ne_eta     =  self.pix.shape[0]
        ne_zeta    =  self.pix.shape[1]
        self.lsNelems = [ne_xi,ne_eta,ne_zeta]
        
        
        e_xi   = np.linspace(self.xmin,self.xmax, ne_xi  +1 )
        e_eta  = np.linspace(self.ymin,self.ymax, ne_eta +1 )
        e_zeta = np.linspace(self.zmin,self.zmax, ne_zeta+1 )
        
        self.lsknotX  = np.r_[ oneX*self.xmin, e_xi   , oneX*self.xmax ]
        self.lsknotY  = np.r_[ oneY*self.ymin, e_eta  , oneY*self.ymax ]
        self.lsknotZ  = np.r_[ oneZ*self.zmin, e_zeta , oneZ*self.zmax ]  
        
        # Perform integration by pixel assembly 
        # (p+1)*(q+1)*(r+1) Gauss point on each voxel

        nbg_xi    = degree[0]+1
        nbg_eta   = degree[1]+1  
        nbg_zeta  = degree[2]+1 
        
        Gauss_xi   =  nb.GaussLegendre(nbg_xi)
        Gauss_eta  =  nb.GaussLegendre(nbg_eta)
        Gauss_zeta =  nb.GaussLegendre(nbg_zeta)
        
        
        xi_min    =  np.kron(e_xi[:-1],np.ones(nbg_xi))    
        eta_min   =  np.kron(e_eta[:-1],np.ones(nbg_eta))
        zeta_min  =  np.kron(e_zeta[:-1],np.ones(nbg_zeta))
        
        xi_g      =  np.kron(np.ones(ne_xi)  , Gauss_xi[0])
        eta_g     =  np.kron(np.ones(ne_eta) , Gauss_eta[0])
        zeta_g    =  np.kron(np.ones(ne_zeta), Gauss_zeta[0])
        
        
        """ Going from the referance element to the parametric space  """ 
        pix        = xi_min   + 0.5*(xi_g+1)*self.dx    # Aranged gauss points in  xi direction  
        piy        = eta_min  + 0.5*(eta_g+1)*self.dy   # Aranged gauss points in  eta direction 
        piz        = zeta_min + 0.5*(zeta_g+1)*self.dz   # Aranged gauss points in  zeta direction 
         
        """ Integration points and weights """ 
        
        wg_xi       =  np.kron(np.ones(ne_xi) , Gauss_xi[1])
        wg_eta      =  np.kron(np.ones(ne_eta), Gauss_eta[1])
        wg_zeta     =  np.kron(np.ones(ne_zeta), Gauss_zeta[1])
        wg          =  np.kron(wg_zeta,np.kron(wg_eta, wg_xi))*self.dx*self.dy*self.dz/8  
        
        
        phi_xi    = nb.global_basisfunsWd(self.lsDegree[0],self.lsknotX,pix) 
        phi_eta   = nb.global_basisfunsWd(self.lsDegree[1],self.lsknotY,piy)
        phi_zeta  = nb.global_basisfunsWd(self.lsDegree[2],self.lsknotZ,piz)
        
        phiEta_phiXi   =  sps.kron(phi_eta,phi_xi, 'csc')
        phi            =  sps.kron(phi_zeta, phiEta_phiXi ,'csc')
        
        # Total integration points 
        pixT  = np.kron(np.ones(nbg_eta*ne_eta*nbg_zeta*ne_zeta), pix)
        piyT  = np.kron(np.ones(nbg_zeta*ne_zeta), np.kron(piy,np.ones(nbg_xi*ne_xi)))
        pizT  = np.kron(piz, np.ones(nbg_xi*ne_xi*nbg_eta*ne_eta))
        
        f = self.GetGrayLevelOfPoint(pixT,piyT,pizT) 

        
        self.lsc = phi.T.dot(wg*f)/phi.T.dot(wg) 
        
        
        
    def GetThresholdedC8Voxels(self,m, thrsh, pixEvalMethod):
 
        if pixEvalMethod=='trilinear':
            r = ir.GetC8MeshFromVoxelsTrilinearInterp(self.pix, 
                                                      self.knotXi1, self.knotEta1, self.knotZeta1,
                                                      thrsh, 
                                                      m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                      int(np.product(m.n_elems))  )
        return r 
        
        
 
#%% DVC Engine 
class DVCEngine:
    def __init__(self,pixEvalMethod):
        self.pixEvalMethod = pixEvalMethod
        self.f     = None 
        self.dyn   = None 
        self.mean0 = None 
        self.std0  = None 
        self.phidf = None 
        
        self.feMean = None
        self.feStd  = None 
        self.feDyn  = None  
        self.fip    = None 
        self.dfdxip = None 
        self.dfdyip = None 
        self.dfdzip = None  
        
    def GetImage_Mean_Std_FE(self, f, m ):
        nipe = m.N.shape[0] # number of integration points per element 
        ne = m.e.shape[0]
        nip = nipe * ne  # total number of integration points 
        if self.pixEvalMethod=='trilinear':
            r = ir.GetMeanImageAndStdOnFE_Mesh_TrilinearInterp(f.pix, f.knotXi1, f.knotEta1, f.knotZeta1, 
                                                               m.e, m.n, m.N, int(ne), int(ne), int(ne),
                                                            int(nip),int(nip),int(nip),int(nip)) 
        elif self.pixEvalMethod=='cBspline3':
            r = ir.GetMeanImageAndStdOnFE_Mesh_CBspline3(f.pix, f.xmin,f.ymin, f.zmin, 
                                                      f.dx, f.dy, f.dz, 
                                                      f.pxc,f.pyc,f.pzc,
                                                      m.e, m.n, m.N, int(ne), int(ne), int(ne),
                                                            int(nip),int(nip),int(nip),int(nip)) 
        elif self.pixEvalMethod=='L2LumpedProj':
            r = ir.GetMeanImageAndStdOnFE_Mesh_L2ProjLumped(f.lsc, 
                                                            f.lsknotX, f.lsknotY, f.lsknotZ, 
                                                            f.lsDegree[0], f.lsDegree[1], f.lsDegree[2],
                                                            m.e, m.n, m.N, int(ne), int(ne), int(ne),
                                                            int(nip),int(nip),int(nip),int(nip)) 
        else : 
            raise ValueError('Sub voxel method not defined')
        
        self.feMean = r[0]
        self.feStd  = r[1]
        self.feDyn  = r[2]
        self.fip    = r[3]
        self.dfdxip = r[4]
        self.dfdyip = r[5]
        self.dfdzip = r[6] 
      
        
    def GetImage_Mean_Std_Structured(self,f,m,nbippe):
        # Returns the total vector of image evaluation at integration points
        # And its gradient 
        # Returns also the mean of gray-level and std over each element of the structured mesh 
        ne   =    m.n_elems[0]*m.n_elems[1]*m.n_elems[2]  
        nipe =   nbippe[0]*nbippe[1]*nbippe[2]
        nip =  nipe*ne 
  
        if self.pixEvalMethod=='trilinear':
            # r = ir.GetMeanImageAndStdOnMesh_TrilinearInterp(f.pix, f.xmin, f.xmax, f.ymin, f.ymax, f.zmin, f.zmax, 
            #                                                 f.dx, f.dy, f.dz,  int(m.pp[0]) , int(m.pp[1]),int(m.pp[2]),
            #                                                 m.xxsi[0], m.xxsi[1], m.xxsi[2], 
            #                                                 int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
            #                                                 int(ne), int(ne), 
            #                                                 int(nip),int(nip),int(nip),int(nip)) 
            r = ir.GetMeanImageAndStdOnMesh_TrilinearInterp(f.pix, 
                                                            f.knotXi1, f.knotEta1, f.knotZeta1, 
                                                            int(m.pp[0]) , int(m.pp[1]),int(m.pp[2]),
                                                            m.xxsi[0], m.xxsi[1], m.xxsi[2], 
                                                            int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                            m.pix1d, m.piy1d, m.piz1d, 
                                                            int(ne), int(ne), int(ne),
                                                            int(nip),int(nip),int(nip),int(nip))   
        elif self.pixEvalMethod=='cBspline2': 
            r = ir.GetMeanImageAndStdOnMesh_CBspline2(f.pix, 
                                                      f.xmin,f.ymin, f.zmin, 
                                                      f.dx, f.dy, f.dz, 
                                                      f.pxc,f.pyc,f.pzc,
                                                      int(m.pp[0]) , int(m.pp[1]),int(m.pp[2]),
                                                      m.xxsi[0], m.xxsi[1], m.xxsi[2], 
                                                      int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                      m.pix1d, m.piy1d, m.piz1d, 
                                                      int(ne), int(ne), int(ne),
                                                      int(nip),int(nip),int(nip),int(nip))   
            
        elif self.pixEvalMethod=='cBspline3':
            r = ir.GetMeanImageAndStdOnMesh_CBspline3(f.pix, 
                                                      f.xmin,f.ymin, f.zmin, 
                                                      f.dx, f.dy, f.dz, 
                                                      f.pxc,f.pyc,f.pzc,
                                                      int(m.pp[0]) , int(m.pp[1]),int(m.pp[2]),
                                                      m.xxsi[0], m.xxsi[1], m.xxsi[2], 
                                                      int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                      m.pix1d, m.piy1d, m.piz1d, 
                                                      int(ne), int(ne), int(ne),
                                                      int(nip),int(nip),int(nip),int(nip)) 
        else: 
            raise ValueError('Sub-voxel evaluation method not defined') 
        
        self.feMean = r[0]
        self.feStd  = r[1]
        self.feDyn  = r[2]
        self.fip    = r[3]
        self.dfdxip = r[4]
        self.dfdyip = r[5]
        self.dfdzip = r[6] 
            
    def GetImage_Mean_Std_thrsh_Structured(self,f,m,nbippe,thrsh):
        ne   =    m.n_elems[0]*m.n_elems[1]*m.n_elems[2]  
        nipe =   nbippe[0]*nbippe[1]*nbippe[2]
        nip =  nipe*ne 
  
        
        if self.pixEvalMethod=='trilinear':
            r = ir.GetMeanImageAndStdOnMesh_thrsh_TrilinearInterp(f.pix, 
                                                            thrsh,  
                                                            f.knotXi1, f.knotEta1, f.knotZeta1, 
                                                            int(m.pp[0]) , int(m.pp[1]),int(m.pp[2]),
                                                            m.xxsi[0], m.xxsi[1], m.xxsi[2], 
                                                            int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                            m.pix1d, m.piy1d, m.piz1d, 
                                                            int(ne), int(ne), 
                                                            int(nip),int(nip),int(nip),int(nip))    
        else: 
            raise ValueError('Sub-voxel evaluation method not defined') 
            
        self.feMean = r[0]
        self.feStd  = r[1]
        self.fip    = r[2]
        self.dfdxip = r[3]
        self.dfdyip = r[4]
        self.dfdzip = r[5] 
    

        
    def GetElementalThresholdedVoxelsIndices(self,m,nbippe,thrsh):
        if self.pixEvalMethod=='trilinear': 
            r = ar.VoxelIntegrationThresholdTrilinearInterp(self.fip,thrsh, int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                            m.pp[0], m.pp[1], m.pp[2],
                                                            m.xxsi[0], m.xxsi[1], m.xxsi[2] ) 
            return np.array(r).astype('int32')
        
    def GetElementalMaskedVoxelsIndices(self,m,nbippe,maskip,pixEvalMethod):
        if self.pixEvalMethod=='trilinear': 
            r = ar.VoxelIntegrationMaskTrilinearInterp(maskip, int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                            m.pp[0], m.pp[1], m.pp[2],
                                                            m.xxsi[0], m.xxsi[1], m.xxsi[2] ) 
            return np.array(r).astype('int32')
    
    
            
    def LHS_thrsh(self, m, nbippe, ipIndices ):
        n_elems  = m.n_elems[0]*m.n_elems[1]*m.n_elems[2]
        nbf_elem = (m.pp[0]+1)*(m.pp[1]+1)*(m.pp[2]+1)
        nnz = int( 9*(nbf_elem)**2*n_elems )
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        r = ar.DVC_LHS_thrsh(self.dfdxip , self.dfdyip , self.dfdzip,
                                             ipIndices, 
                                             m.pp[0], m.pp[1], m.pp[2],
                                             m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                             int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                             m.Nxi, m.Neta, m.Nzeta, 
                                             nnz,nnz,nnz) 
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        H = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return H  
            

    
    def RHS_thrsh(self,g, m, nbippe, ipIndices, U ):
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        if self.pixEvalMethod=='trilinear':
            r = ar.DVC_RHS_thrsh_TrilinearInterp(g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, 
                                             self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                             ipIndices, 
                                             m.pp[0], m.pp[1], m.pp[2], 
                                             m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                             m.noelem, 
                                             int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                             m.pix1d, m.piy1d, m.piz1d, 
                                             m.Nxi, m.Neta, m.Nzeta, U, ndof )
            ssd = r[0]
            rhs = r[1] 
        return rhs, ssd  
    
    def RHS_ZN_thrsh(self, g, m, nbippe, ipIndices, U ):
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        if self.pixEvalMethod=='trilinear':
            r = ar.DVC_RHS_ZN_thrsh_TrilinearInterp(g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                    self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                    self.feMean, self.feStd,
                                                    ipIndices, 
                                                    m.pp[0], m.pp[1], m.pp[2], 
                                                    m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                    m.noelem, 
                                                    int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                    m.pix1d, m.piy1d, m.piz1d, 
                                                    m.Nxi, m.Neta, m.Nzeta, U, ndof ) 
            znssd = r[0]
            rhs = r[1] 
        return rhs, znssd  
   
    
    def LHS_FE_EA(self, m):
        n_elems  = m.e.shape[0] 
        nbf_elem = 4         
        nnz = int( 9*(nbf_elem)**2*n_elems ) 

        r = ar.DVC_LHS_FE(m.e, m.n, m.conn, m.N, m.dNdxi, m.dNdeta, m.dNdzeta, m.iw, self.dfdxip , self.dfdyip , self.dfdzip, nnz, nnz, nnz)
        
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        H = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (m.ndof,m.ndof )) 
        
        return H 
    
    def RHS_FE_EA(self,g,m,U):
        if self.pixEvalMethod=='trilinear':
            r = ar.DVC_RHS_FE_TrilinearInterp(g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, self.fip, self.dfdxip , self.dfdyip , self.dfdzip,
                              m.e, m.n, m.conn, m.N, m.dNdxi, m.dNdeta, m.dNdzeta, m.iw, U, m.ndof) 
        elif self.pixEvalMethod=='cBspline3':
            r = ar.DVC_RHS_FE_CBspline3(g.pix, g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                        g.dx, g.dy, g.dz, 
                                        g.pxc,g.pyc,g.pzc, 
                                        self.fip, self.dfdxip , self.dfdyip , self.dfdzip,
                                        m.e, m.n, m.conn, m.N, m.dNdxi, m.dNdeta, m.dNdzeta, m.iw, U, m.ndof) 
        elif self.pixEvalMethod=='L2LumpedProj':
            r = ar.DVC_RHS_FE_L2ProjLumped(g.lsc, g.lsknotX, g.lsknotY, g.lsknotZ, 
                                           g.lsDegree[0], g.lsDegree[1], g.lsDegree[2], 
                                           self.fip, self.dfdxip , self.dfdyip , self.dfdzip,
                                           m.e, m.n, m.conn, m.N, m.dNdxi, m.dNdeta, m.dNdzeta, m.iw, U, m.ndof) 
 
        else: 
            raise ValueError('Check input')
       
        ssd = r[0]
        rhs = r[1] 
        return rhs, ssd  
    
    
    def GLR_FE_EA(self,g,m,U):
        n_elems = m.e.shape[0] # number of element 
        nip = len(m.iw) # number of integration points (= "number of voxels" per element )
        nvoxels = int(n_elems*nip)
        if self.pixEvalMethod=='trilinear':
            res = ar.GLR_FE_TrilinearInterp(g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, self.fip, m.e, m.n, m.conn, m.N, m.iw, U, nvoxels)
        elif self.pixEvalMethod=='cBspline3':
            res = ar.GLR_FE_CBspline3(g.pix, g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                        g.dx, g.dy, g.dz, 
                                        g.pxc,g.pyc,g.pzc, 
                                        self.fip, m.e, m.n, m.conn, m.N, m.iw, U, nvoxels)
        elif self.pixEvalMethod=='L2LumpedProj':
            res = ar.GLR_FE_L2ProjLumped(g.lsc, g.lsknotX, g.lsknotY, g.lsknotZ, 
                                       g.lsDegree[0], g.lsDegree[1], g.lsDegree[2], 
                                       self.fip, m.e, m.n, m.conn, m.N, m.iw, U, nvoxels)
        else: 
            raise ValueError('Check input')
        return res 
    
    def GoIdu(self,g,m,U, xi,eta,zeta,ie):
        """
        Computes g(x+u(x)) at iso-parametric points (xi[i],eta[i],zeta[i]) located at element ie[i]
        """
        if self.pixEvalMethod=='trilinear':
            goPhi = ar.Gophi_FEMesh_TrilinearInterp(g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, xi, eta, zeta, ie, m.e, m.n, m.conn, U, len(xi))
        elif self.pixEvalMethod=='cBspline3':
            goPhi = ar.Gophi_FEMesh_CBspline3(g.pix, g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                        g.dx, g.dy, g.dz, 
                                        g.pxc,g.pyc,g.pzc, 
                                        xi, eta, zeta, ie, m.e, m.n, m.conn, U, len(xi))
        elif self.pixEvalMethod=='L2LumpedProj':
            goPhi = ar.Gophi_FEMesh_L2ProjLumped(g.lsc, g.lsknotX, g.lsknotY, g.lsknotZ, 
                                                 g.lsDegree[0], g.lsDegree[1], g.lsDegree[2],  
                                                 xi, eta, zeta, ie, m.e, m.n, m.conn, U, len(xi))
            
        else: 
            raise ValueError('Check input')
        return goPhi 
            

    def LHS(self,m,nbippe,parallel=True ):
        n_elems  = m.n_elems[0]*m.n_elems[1]*m.n_elems[2]
        nbf_elem = (m.pp[0]+1)*(m.pp[1]+1)*(m.pp[2]+1)
        nnz = int( 9*(nbf_elem)**2*n_elems )
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        if parallel==True: 
            r = ar.DVC_LHS_Structured_Parallel(self.dfdxip , self.dfdyip , self.dfdzip,
                                                      m.pp[0], m.pp[1], m.pp[2],
                                                      m.xxsi[0] , m.xxsi[1] , m.xxsi[2] ,
                                                      int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                      m.Nxi, m.Neta, m.Nzeta, 
                                                      nnz,nnz,nnz)      
        else: 
            r = ar.DVC_LHS_Structured(self.dfdxip , self.dfdyip , self.dfdzip,
                                                      m.pp[0], m.pp[1], m.pp[2],
                                                      m.xxsi[0] , m.xxsi[1] , m.xxsi[2] ,
                                                      int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                      m.Nxi, m.Neta, m.Nzeta, 
                                                      nnz,nnz,nnz)   
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        #print(np.max(indexI), np.max(indexJ), ndof)
        H = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return H  

        
    def RHS(self,g,m,nbippe,U, parallel=True):
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        # r = ar.DVC_RHS_Structured_TrilinearInterp(g.pix, g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, g.dx, g.dy, g.dz,
        #                                           self.fip, self.dfdxip , self.dfdyip , self.dfdzip ,
        #                                           m.pp[0], m.pp[1], m.pp[2],
        #                                           m.xxsi[0], m.xxsi[1], m.xxsi[2],
        #                                           int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
        #                                           U, ndof) 
        
        if parallel==True: 
            if self.pixEvalMethod=='trilinear':
                r = ar.DVC_RHS_Structured_TrilinearInterp_Parallel(g.pix, 
                                                          g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                          self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                          m.pp[0], m.pp[1], m.pp[2],
                                                          m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                          m.noelem, 
                                                          int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                          m.pix1d, m.piy1d, m.piz1d,
                                                          m.Nxi, m.Neta, m.Nzeta, U, ndof  )
            else: 
                raise ValueError('Check input ')
        else : 
            if self.pixEvalMethod=='trilinear':
                r = ar.DVC_RHS_Structured_TrilinearInterp(g.pix, 
                                                          g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                          self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                          m.pp[0], m.pp[1], m.pp[2],
                                                          m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                          m.noelem, 
                                                          int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                          m.pix1d, m.piy1d, m.piz1d,
                                                          m.Nxi, m.Neta, m.Nzeta, U, ndof  )
            elif self.pixEvalMethod=='cBspline2':
                r = ar.DVC_RHS_Structured_CBspline2(g.pix, 
                                                    g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                                    g.dx, g.dy, g.dz, g.pxc, g.pyc, g.pzc, 
                                                    self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                    m.pp[0], m.pp[1], m.pp[2],
                                                    m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                    m.noelem, 
                                                    int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                    m.pix1d, m.piy1d, m.piz1d,
                                                    m.Nxi, m.Neta, m.Nzeta, U, ndof  )
            elif self.pixEvalMethod=='cBspline3':
                r = ar.DVC_RHS_Structured_CBspline3(g.pix, 
                                                    g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                                    g.dx, g.dy, g.dz, g.pxc, g.pyc, g.pzc, 
                                                    self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                    m.pp[0], m.pp[1], m.pp[2],
                                                    m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                    m.noelem, 
                                                    int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                    m.pix1d, m.piy1d, m.piz1d,
                                                    m.Nxi, m.Neta, m.Nzeta, U, ndof  )   
            else: 
                raise ValueError('Sub-voxel evaluation method not defined') 
            
        ssd = r[0]
        rhs = r[1] 
        return rhs, ssd  
    
    def RHS_ZN(self,g,m,nbippe,U):
        ne = int(np.prod(m.n_elems)) 
        nbf = m.Get_nbf()
        ndof = 3*nbf 
        if self.pixEvalMethod=='trilinear':
            r = ar.DVC_RHS_ZN_Structured_TrilinearInterp(g.pix, 
                                                         g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                         self.fip, self.dfdxip, self.dfdyip, self.dfdzip,
                                                         self.feMean, self.feStd, self.feDyn, 
                                                         m.pp[0], m.pp[1], m.pp[2],
                                                         m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                         m.noelem, 
                                                         int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                         m.pix1d, m.piy1d, m.piz1d,
                                                         m.Nxi, m.Neta, m.Nzeta, U, ndof, ne  ) 
            rhs        = r[0]
            elementRes = r[1]
        return rhs, elementRes 
    
    def GLR_Structured(self,g,m,nbippe,U,parallel=True):
        nvoxels = int(np.product(m.n_elems)*np.product(nbippe))
        
        if parallel==True: 
            if self.pixEvalMethod=='trilinear':
                res = ar.GLR_Structured_TrilinearInterp_Parallel(self.fip,g.pix, 
                                                        g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                        m.pp[0], m.pp[1], m.pp[2], 
                                                        m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                        int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                        m.pix1d, m.piy1d, m.piz1d,
                                                        m.Nxi, m.Neta, m.Nzeta, 
                                                        U, nvoxels) 
            else: 
                raise ValueError('Check input ')
        
        else: 
            if self.pixEvalMethod=='trilinear':
                res = ar.GLR_Structured_TrilinearInterp(self.fip,g.pix, 
                                                        g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                        m.pp[0], m.pp[1], m.pp[2], 
                                                        m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                        int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                        m.pix1d, m.piy1d, m.piz1d,
                                                        m.Nxi, m.Neta, m.Nzeta, 
                                                        U, nvoxels) 
            elif self.pixEvalMethod=='cBspline2':
                res = ar.GLR_Structured_CBspline2(self.fip, g.pix, 
                                                  g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                                  g.dx, g.dy, g.dz, 
                                                  g.pxc, g.pyc, g.pzc, 
                                                  m.pp[0], m.pp[1], m.pp[2], 
                                                  m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                  int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                  m.pix1d, m.piy1d, m.piz1d,
                                                  m.Nxi, m.Neta, m.Nzeta, 
                                                  U, nvoxels) 
            elif self.pixEvalMethod=='cBspline3':
                res = ar.GLR_Structured_CBspline3(self.fip, g.pix, 
                                                  g.xmin, g.xmax, g.ymin, g.ymax, g.zmin, g.zmax, 
                                                  g.dx, g.dy, g.dz, 
                                                  g.pxc, g.pyc, g.pzc, 
                                                  m.pp[0], m.pp[1], m.pp[2], 
                                                  m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                  int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                  m.pix1d, m.piy1d, m.piz1d,
                                                  m.Nxi, m.Neta, m.Nzeta, 
                                                  U, nvoxels) 
        return res 
    
    
    def GLR_ZN_Structured(self,g,m,nbippe,U):
        nvoxels = int(np.product(m.n_elems)*np.product(nbippe))
        if self.pixEvalMethod=='trilinear':
            res = ar.GLR_ZN_Structured_TrilinearInterp(self.fip, self.feMean, self.feStd,
                                                       g.pix, g.knotXi1, g.knotEta1, g.knotZeta1, 
                                                       m.pp[0], m.pp[1], m.pp[2], 
                                                       m.xxsi[0], m.xxsi[1], m.xxsi[2],
                                                       int(nbippe[0]) , int(nbippe[1]) ,  int(nbippe[2]),
                                                       m.pix1d, m.piy1d, m.piz1d,
                                                       m.Nxi, m.Neta, m.Nzeta, 
                                                       U, nvoxels) 
        return res 
    
    
 
        
        
    def ComputeLHS(self,f,m):
        if self.pixEvalMethod=='trilinear-structured':
            self.f,dfx,dfy,dfz = f.EvaluateTrilinearInterpolationAndGradientStructured(m.pix1d,m.piy1d,m.piz1d)
        elif self.pixEvalMethod=='trilinear-unstructured':
            self.f, dfx,dfy,dfz = f.EvaluateTrilinearInterpolationAndGradient(m.pix,m.piy,m.piz) 
        elif self.pixEvalMethod=='cBspline2-structured':
            self.f,dfx,dfy,dfz = f.EvaluateCardinalBsplineAndGradientStructured(m.pix1d,m.piy1d,m.piz1d,degree=2)
        elif self.pixEvalMethod=='cBspline2-unstructured':
            self.f,dfx,dfy,dfz = f.EvaluateCardinalBsplineAndGradient(m.pix,m.piy,m.piz,degree=2) 
        elif self.pixEvalMethod=='cBspline3-structured':
            self.f,dfx,dfy,dfz = f.EvaluateCardinalBsplineAndGradientStructured(m.pix1d,m.piy1d,m.piz1d,degree=3)
        elif self.pixEvalMethod=='cBspline3-unstructured':
            self.f,dfx,dfy,dfz = f.EvaluateCardinalBsplineAndGradient(m.pix,m.piy,m.piz,degree=3) 
        elif self.pixEvalMethod=='lvlstBspline-structured':
            self.f,dfx,dfy,dfz = f.EvaluateBsplineLevelSetAndGradientStructured(m.pix1d,m.piy1d,m.piz1d) 
        elif self.pixEvalMethod=='lvlstBspline-unstructured':
            self.f,dfx,dfy,dfz = f.EvaluateBsplineLevelSetAndGradient(m.pix,m.piy,m.piz) 
        else: 
            raise ValueError('Wrong image evaluation method')
            
        phidf = sps.diags(dfx).dot(m.phix)+sps.diags(dfy).dot(m.phiy)+sps.diags(dfz).dot(m.phiz) 
        self.dyn   = np.max(self.f)-np.min(self.f) 
        self.mean0 = np.mean(self.f)
        self.std0  = np.std(self.f)
        #self.f    -= self.mean0
        
        if type(m.wg) == np.float64 : 
            self.wphidf = m.wg*phidf # Here the integration weight is constant (one voxel volume) 
        else: 
            self.wphidf = m.wg.dot(phidf) 
            
        return phidf.T.dot(self.wphidf) 


    
    def ComputeRHS(self,g,m,U):
        
        x =  m.pix + m.phix.dot(U)
        y =  m.piy + m.phiy.dot(U)
        z =  m.piz + m.phiz.dot(U)
        
        if self.pixEvalMethod=='trilinear':
            res = g.EvaluateTrilinearInterpolation(x,y,z) 
        elif self.pixEvalMethod=='cBspline2': 
            res = g.EvaluateCardinalBspline(x,y,z,degree=2) 
        elif self.pixEvalMethod=='cBspline3':
            res = g.EvaluateCardinalBspline(x,y,z,degree=3) 
        elif self.pixEvalMethod=='lvlstBspline':
            res = g.EvaluateBsplineLevelSet(x,y,z)
        else: 
            raise ValueError('Wrong image evaluation method')
        #res -=  np.mean(res) 
        #std1 =  np.std(res) 
        #res  =  self.f-self.std0/std1*res
        res  = self.f - res 
        b = self.wphidf.T.dot(res)
        return b,res  
    
    def GnInit(self,f,m):
        self.f= f.EvaluateTrilinearInterpolation(m.pix,m.piy,m.piz)
        self.dyn=np.max(self.f)-np.min(self.f)
        self.mean0 = np.mean(self.f)
        self.std0  = np.std(self.f)
        self.f    -= self.mean0
        
    def ComputeGnMembers(self,f,g,m,U):
        x =  m.pix + m.phix.dot(U)
        y =  m.piy + m.phiy.dot(U)
        z =  m.piz + m.phiz.dot(U)
        
        if self.pixEvalMethod=='trilinear':
            res, dgx, dgy, dgz = g.EvaluateTrilinearInterpolationAndGradient(x,y,z) 
        elif self.pixEvalMethod=='cBspline2': 
            res, dgx, dgy, dgz = g.EvaluateCardinalBsplineAndGradient(x,y,z,degree=2) 
        elif self.pixEvalMethod=='cBspline3':
            res, dgx, dgy, dgz = g.EvaluateCardinalBsplineAndGradient(x,y,z,degree=3) 
        elif self.pixEvalMethod=='lvlstBspline':
            res, dgx, dgy, dgz = g.EvaluateBsplineLevelSetAndGradient(x,y,z)
        else: 
            raise ValueError('Wrong image evaluation method')
 
        res -=  np.mean(res) 
        std1 =  np.std(res) 
        res  =  self.f-self.std0/std1*res
        
        phidg  = sps.diags(dgx).dot(m.phix)+sps.diags(dgy).dot(m.phiy)+sps.diags(dgz).dot(m.phiz) 
        if type(m.wg) == np.float64 : 
            wphidg = m.wg*phidg # Here the integration weight is constant (one voxel volume) 
        else: 
            wphidg = m.wg.dot(phidg)
        d = wphidg.T.dot(res)
        H = phidg.T.dot(wphidg)  
        return H,d,res
    
    
    
    
#%% Tetra class
class TetraMesh1:
    def __init__(self,elem=None,nodes=None):
        if elem is not None and nodes is not None:
            self.n = nodes
            self.e  = elem.astype(np.int32) 
        else:
            self.n = None 
            self.e = None 
        
        self.pix = 0   # x coordinates of the integration points 
        self.piy = 0   # y coordinates of the integration points  
        self.piz = 0   # z coordinates of the integration points 
        self.wg =  0   # Integration weightsxabs(detJ) diagonal matrix  
        
        self.phi = 0          
        self.dphidx = 0      
        self.dphidy = 0       
        self.dphidz = 0      
 
    
        self.phix    = 0
        self.phiy    = 0
        self.phiz    = 0
        
        self.dphixdx = 0
        self.dphixdy = 0
        self.dphixdz = 0
        
        self.dphiydx = 0
        self.dphiydy = 0
        self.dphiydz = 0
        
        self.dphizdx = 0
        self.dphizdy = 0
        self.dphizdz = 0
    
    
    def Load(self,filename,removeLowerCells=True):
        mesh  = meshio.read(filename)
        if removeLowerCells==True:  
            print('Meshio: removing lower dimensional cells')
            mesh.remove_lower_dimensional_cells() 
        self.n = np.copy(mesh.points).astype(np.float64)  # I MUST PAY ATTENTION WITH DTYPES !!! ;(
        self.e = np.copy(mesh.cells_dict['tetra']).astype(np.int32) 
    
        
    def CleanMesh(self,removeIsolatedComponents=True, epsVolume=0):
        volElements = self.ComputeVolumeOfElements()
        ieToremove = np.where(volElements<epsVolume)[0]
        self.e   = np.delete(self.e, ieToremove, axis=0)
        print(str(len(ieToremove))+' elements with a volume smaller than '+str(epsVolume)+' were removed')
        
        # jacobE_terms = self.GetJacobianOfElements(np.arange(self.e.shape[0]))
        # detJelements = jacobE_terms[-1]
        # ieToremove = np.where(np.abs(detJelements)>epsJacobian)[0]
        # self.e   = np.delete(self.e, ieToremove, axis=0)
        # print(str(len(ieToremove))+' elements with a |Jac| greater than '+str(epsJacobian)+' were removed')
        
        nold = self.n.shape[0]
        mesh     = pymesh.form_mesh(self.n, faces = None,  voxels=self.e)
        mesh,_     = pymesh.remove_isolated_vertices(mesh)
        self.n   = mesh.vertices 
        self.e   = mesh.voxels 
        if removeIsolatedComponents==True: 
            print(str(nold-self.n.shape[0])+' isolated nodes were removed')
            list_meshes     = pymesh.separate_mesh(mesh, connectivity_type='voxel')
            ne = np.array([l.voxels.shape[0] for l in list_meshes ])
            iMesh = int(np.argmax(ne))
            print(str(len(ne)-1)+' mesh components were removed')
            self.n   = list_meshes[iMesh].vertices
            self.e   = list_meshes[iMesh].voxels  
      
            
 
    def RemoveElementsBasedOnGrayLevelThreshold(self, f, pixEvalMethod, thrsh_gl, thrsh_vf):
        """
        Removes the elements that have a volume fraction (based on thrsh_gl)
        smaller than the value thrsh_vf
        thrsh_gl: grey-level threshold 
        thrsh_vf: volume fraction threshold 
        Returns a new finite element mesh 
        """
        self.SetVoxelTetraIntegrationRule(f.dx)
        dvc = DVCEngine(pixEvalMethod)
        dvc.GetImage_Mean_Std_FE(f, self)
        # Get the volume fraction for each element 
        nvoxels = self.N.shape[0] # Constant number of voxels per element 
        voxels_thrsh = ((dvc.fip>thrsh_gl)*1).reshape((self.e.shape[0],nvoxels))
        volume_fraction = np.sum(voxels_thrsh,axis=1)/nvoxels
        print('Smallest volume_fraction value: '+str(np.min(volume_fraction)))
        print('Greatest volume_fraction value: '+str(np.max(volume_fraction)))
        ik = np.where(volume_fraction>thrsh_vf)[0] # indices of elements to be keeped  
        print(str(self.e.shape[0]-len(ik))+' elements were removed based on the volume fraction')
        m = TetraMesh1(self.e[ik],self.n)
        m.CleanMesh(removeIsolatedComponents=False)
        return m 
        
        
        
 
    def Connectivity(self):
        """ Compute connectivity """
        self.conn=-np.ones(self.n.shape[0],dtype=np.int)
        c=0
        self.nvm=0
        for je in range(self.e.shape[0]):
            newnodes,=np.where(self.conn[self.e[je,:]]<0)
            self.conn[self.e[je,newnodes]]=c+np.arange(len(newnodes))
            c+=newnodes.shape[0]               
            self.nvm+=4*4**2
        self.conn=np.c_[self.conn,self.conn+c*(self.conn>=0),self.conn+2*c*(self.conn>=0)]
        self.ndof=c*3   
        self.conn = self.conn.astype(np.int32)
        
    def SetFaceTetrahedronConnectivity(self):
        
        # Old method 
        # ip1   =  self.e[:,0] ; ip2 =  self.e[:,1] ; ip3 =  self.e[:,2] ; ip4 =  self.e[:,3] 
        # face1 =  np.c_[ip1,ip2,ip3]
        # face2 =  np.c_[ip1,ip2,ip4]
        # face3 =  np.c_[ip1,ip3,ip4]
        # face4 =  np.c_[ip2,ip3,ip4]     
        # faces = np.vstack((face1,face2,face3,face4))
        # faceElem = np.kron(np.ones(4), np.arange(self.e.shape[0])).astype('int')
        # faces = np.sort(faces,axis=1) 
        # unique_faces = np.unique(faces,axis=0) 
        # conn_list = [np.where((faces==r).all(axis=1))[0] for r in unique_faces]
        # conn = -1*np.ones((len(conn_list),2))
        # for i,c in enumerate(conn_list):
        #     conn[i,:len(c)] = faceElem[c]
        # self.tetFaces      = unique_faces 
        # self.connFaces  = conn.astype(np.int32)
        # self.faceMesh   = pymesh.form_mesh(self.n, self.tetFaces )
        
        
        ip1   =  self.e[:,0] ; ip2 =  self.e[:,1] ; ip3 =  self.e[:,2] ; ip4 =  self.e[:,3] 
        face1 =  np.c_[ip1,ip2,ip3]
        face2 =  np.c_[ip1,ip2,ip4]
        face3 =  np.c_[ip1,ip3,ip4]
        face4 =  np.c_[ip2,ip3,ip4]     
        faces = np.vstack((face1,face2,face3,face4))
        faceElem = np.kron(np.ones(4), np.arange(self.e.shape[0])).astype('int')
        faces = np.sort(faces,axis=1) 
        face_groups = tools.group_duplicate_index_v2(faces) 
        ngroups = len(face_groups)
        nfaces = faces.shape[0] - ngroups
        iDup    = np.array([item for sublist in face_groups for item in sublist])
        iNonDup = np.setdiff1d(np.arange(faces.shape[0]),iDup )
        eDup    = faceElem[iDup].reshape((ngroups,2))
        eNonDup = np.vstack((faceElem[iNonDup],-1*np.ones(nfaces-ngroups))).T 
        conn = np.vstack((eDup,eNonDup)) 
        self.tetFaces  = np.vstack((faces[iDup[::2],:],faces[iNonDup])) 
        self.connFaces = conn.astype(np.int32)
        self.faceMesh  = pymesh.form_mesh(self.n, self.tetFaces )
 
        
        
    
    def SetGaussIntegrationRule(self):
        
        self.N       = np.array([[1./4, 1./4, 1./4, 1./4]])
        self.dNdxi   = np.array([[-1,1,0,0]])
        self.dNdeta  = np.array([[-1,0,1,0]]) 
        self.dNdzeta = np.array([[-1,0,0,1]])
        self.iw = np.array([1./6])
        self.ip = np.array([[1./4, 1./4, 1./4, 1./4]])

    def ComputeVolume(self):
        x1 = self.n[self.e[:,0],0] ; y1 = self.n[self.e[:,0],1] ; z1 = self.n[self.e[:,0],2]
        x2 = self.n[self.e[:,1],0] ; y2 = self.n[self.e[:,1],1] ; z2 = self.n[self.e[:,1],2]
        x3 = self.n[self.e[:,2],0] ; y3 = self.n[self.e[:,2],1] ; z3 = self.n[self.e[:,2],2] 
        x4 = self.n[self.e[:,3],0] ; y4 = self.n[self.e[:,3],1] ; z4 = self.n[self.e[:,3],2]
        volE = (1./6.)*( np.abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                   (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                   (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) )
        return np.sum(volE) 
    def ComputeVolumeOfElements(self):
        x1 = self.n[self.e[:,0],0] ; y1 = self.n[self.e[:,0],1] ; z1 = self.n[self.e[:,0],2]
        x2 = self.n[self.e[:,1],0] ; y2 = self.n[self.e[:,1],1] ; z2 = self.n[self.e[:,1],2]
        x3 = self.n[self.e[:,2],0] ; y3 = self.n[self.e[:,2],1] ; z3 = self.n[self.e[:,2],2] 
        x4 = self.n[self.e[:,3],0] ; y4 = self.n[self.e[:,3],1] ; z4 = self.n[self.e[:,3],2]
        volE = (1./6.)*( np.abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                   (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                   (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) )
        return volE 
 
                
    def SetVoxelTetraIntegrationRule(self, voxel_size, confirm=False):
        meanVol = self.ComputeVolume()/self.e.shape[0] # Mean element volume 
        N = int(np.ceil( meanVol**(1/3)/(voxel_size) )) 
        
        dx = 1/N
        x = dx/2  + np.arange(N)*dx 
        y = dx/2  + np.arange(N)*dx 
        z = dx/2  + np.arange(N)*dx 
        X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
        iTetra = (Z <= 1-X-Y) 
        
        xi  =  X[iTetra]
        eta = Y[iTetra]
        zeta = Z[iTetra]

        self.iw = np.ones(xi.size)/(6*xi.size)
        self.ip = np.concatenate((xi,eta,zeta)).reshape((3,len(xi))).T 

        one  = np.ones(self.ip.shape[0])
        zero = np.zeros(self.ip.shape[0])
        self.N       = np.c_[ xi, eta, zeta, 1-xi-eta-zeta ] 
        self.dNdxi   = np.c_[ one, zero, zero, -one ]
        self.dNdeta  = np.c_[ zero, one, zero, -one ]
        self.dNdzeta = np.c_[ zero, zero, one, -one ]
        
        self.nodesRef  = None # The nodes of the subdivided reference triangle 
        self.elemsRef  = None   # The elements of the subdivided reference triangle 
        self.NnodesRef = None 
        
        meanIt = np.mean(self.iw)*meanVol # Mean integration tetrahedral volume 
        print('Mean volume of integration tetra ', meanIt )
        
        
    
    def SetVoxelTetraIntegrationRuleOld(self, N, confirm=False):
        
        meanVol = self.ComputeVolume()/self.e.shape[0] # Mean element volume 
        # N = meanVol**(1/3)/(voxelsize)
        print('N new')
        print(N) # Virer N
        
        lc1 = 1/N   
        lc2 = 1/N  
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add("subdivisionTetra")
        gmsh.model.geo.addPoint( 1, 0, 0, lc2, 1)
        gmsh.model.geo.addPoint( 0, 1, 0, lc2, 2)
        gmsh.model.geo.addPoint( 0, 0 ,1, lc2, 3)
        gmsh.model.geo.addPoint( 0, 0 ,0, lc1, 4)    
        gmsh.model.geo.addLine(4,3,1)
        gmsh.model.geo.addLine(4,2,2)
        gmsh.model.geo.addLine(2,3,3)
        gmsh.model.geo.addLine(3,1,4)  
        gmsh.model.geo.addLine(4,1,5)
        gmsh.model.geo.addLine(1,2,6)        
        gmsh.model.geo.addCurveLoop([2,-6,-5],1)
        gmsh.model.geo.addCurveLoop([1,4,-5],2)
        gmsh.model.geo.addCurveLoop([3,4,6],3)
        gmsh.model.geo.addCurveLoop([2,3,-1],4)
        gmsh.model.geo.addPlaneSurface([1],1)    
        gmsh.model.geo.addPlaneSurface([2],2)    
        gmsh.model.geo.addPlaneSurface([3],3)    
        gmsh.model.geo.addPlaneSurface([4],4)    
        gmsh.model.geo.addSurfaceLoop([2, 4, 1, 3],1)
        gmsh.model.geo.addVolume([1],1)
         
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
 
        nums,nodes,e=gmsh.model.mesh.getNodes()
        nodes=nodes.reshape((len(nums),3))

        nums,els = gmsh.model.mesh.getElementsByType(4)
        nnd=len(els)//len(nums)
        els=els.reshape((len(nums),nnd))-1
        
        # gmsh.finalize()
        
        print(els.shape)
        print(nodes.shape)
        
        ip = np.sum( nodes[els], axis=1)/4
        # Volume of each tetrahedron 
        x1 = nodes[els[:,0],0] ; y1 = nodes[els[:,0],1] ; z1 = nodes[els[:,0],2]
        x2 = nodes[els[:,1],0] ; y2 = nodes[els[:,1],1] ; z2 = nodes[els[:,1],2]
        x3 = nodes[els[:,2],0] ; y3 = nodes[els[:,2],1] ; z3 = nodes[els[:,2],2] 
        x4 = nodes[els[:,3],0] ; y4 = nodes[els[:,3],1] ; z4 = nodes[els[:,3],2]
        iw = (1./6.)*( np.abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                   (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                   (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) )  
        self.ip = ip  # list of integration poitns 
        self.iw = iw  # list of integration weights 
        xi   = ip[:,0] 
        eta  = ip[:,1]
        zeta = ip[:,2]
        one  = np.ones(ip.shape[0])
        zero = np.zeros(ip.shape[0])
        self.N       = np.c_[ xi, eta, zeta, 1-xi-eta-zeta ] 
        self.dNdxi   = np.c_[ one, zero, zero, -one ]
        self.dNdeta  = np.c_[ zero, one, zero, -one ]
        self.dNdzeta = np.c_[ zero, zero, one, -one ]
        
        self.nodesRef  = nodes # The nodes of the subdivided reference triangle 
        self.elemsRef  = els   # The elements of the subdivided reference triangle 
        self.NnodesRef = np.c_[ nodes[:,0], nodes[:,1], nodes[:,2], 1-nodes[:,0]-nodes[:,1]-nodes[:,2]  ] 
        
        
        meanVol = self.ComputeVolume()/self.e.shape[0] # Mean element volume 
        meanIt = np.mean(iw)*meanVol # Mean integration tetrahedral volume 
        print('Mean volume of integration tetra ', meanIt )
        if confirm==True: 
            cont = input('are you sure to choose this discretization ? (yes or no): ')
            if cont=='no':
                sys.exit(1) 
                
  
        
 
        fig = plt.figure() 
        ax = fig.gca( projection ='3d')
        ax.plot(nodes[:,0],nodes[:,1],nodes[:,2],'o', color='red')
        ax.plot([0,1],[0,0],[0,0],'k-') 
        ax.plot([0,0],[0,1],[0,0],'k-') 
        ax.plot([0,0],[0,0],[0,1],'k-') 
        ax.plot([1,0],[0,0],[0,1],'k-') 
        ax.plot([0,0],[0,1],[1,0],'k-') 
        ax.plot([1,0],[0,1],[0,0],'k-') 
        ax.plot(ip[:,0],ip[:,1],ip[:,2], 'o', markersize =2, color='blue')
        plt.title('Reference tetrahedron integration rule')
 
 
             
        
    
    def UniformIntegration(self):
        ne = 3 
        x = np.linspace(0,1,ne+1)
        X,Y = np.meshgrid(x,x)
        Xr = X.ravel() 
        Yr = Y.ravel()

        i  = np.arange(ne)
        I,J = np.meshgrid(i,i,indexing='ij')
        
        i = np.kron(np.ones(ne), np.arange(ne) ) + np.kron( np.linspace(0,ne**2-1,ne), np.ones(ne))
        i = i.astype('int')
        
        xp1 = (Xr[i]+Xr[i+1]+Xr[i+ne+1])/3
        yp1 = (Yr[i]+Yr[i+1]+Yr[i+ne+1])/3
        
        xp2 = (Xr[i+1]+Xr[i+ne+1]+Xr[i+ne+2])/3
        yp2 = (Yr[i+1]+Yr[i+ne+1]+Yr[i+ne+2])/3
        
        iT1 = yp1<=1-xp1
        iT2 = yp2<=1-xp2
        
        xp1 = xp1[iT1]
        yp1 = yp1[iT1]
        xp2 = xp2[iT2]
        yp2 = yp2[iT2]
        
        plt.figure() 
        plt.plot(Xr,Yr,'o')
        plt.plot(X,Y,color='black')
        plt.plot(Y,X,color='black')
        plt.plot(xp1,yp1,'.', color='green', markersize=2)
        plt.plot(xp2,yp2,'.', color='red', markersize=2)
        
        Nx = 10  
        Ny = 10
        Nz = 10
        dx = 1/Nx 
        dy = 1/Ny 
        dz = 1/Nz 
        w = dx*dy*dz
        x = dx/2  + np.arange(Nx)*dx 
        y = dy/2  + np.arange(Ny)*dx 
        z = dz/2  + np.arange(Nz)*dx 
        X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
        iTetra = (Z <= 1-X-Y) 
        
        X = X[iTetra]
        Y = Y[iTetra]
        Z = Z[iTetra]
        
        n = len(X)
        print('Reference tetrahedron area approximation')
        print(n*w)
 

        fig = plt.figure() 
        ax = fig.gca( projection ='3d')
        ax.plot(X,Y,Z,'o', color='red')
        ax.plot([0,1],[0,0],[0,0],'k-')
        ax.plot([0,0],[0,1],[0,0],'k-') 
        ax.plot([0,0],[0,0],[0,1],'k-') 
        ax.plot([1,0],[0,0],[0,1],'k-') 
        ax.plot([0,0],[0,1],[1,0],'k-') 
        ax.plot([1,0],[0,1],[0,0],'k-') 
 
 
    def LaplacianEA(self):
        n_elems  = self.e.shape[0] 
        nbf_elem = 4 
        nnz = int( 3*(nbf_elem)**2*n_elems ) 
        
        r = ar.Laplacian_FE(self.e, self.n, self.conn, self.N, self.dNdxi, self.dNdeta, self.dNdzeta, self.iw, nnz, nnz, nnz)
 
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        L = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (self.ndof,self.ndof ))
        return L 


    
    def StiffnessEA(self, E, nu): 
        n_elems  = self.e.shape[0] 
        nbf_elem = 4         
        nnz = int( 9*(nbf_elem)**2*n_elems ) 

        r = ar.Stiffness_FE(E, nu, self.e, self.n, self.conn, self.N, self.dNdxi, self.dNdeta, self.dNdzeta, self.iw, nnz, nnz, nnz)
        
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        K = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (self.ndof,self.ndof )) 
        
        return K 
        
    

    
    def StressFromStrain(self,hooke, exx,eyy,ezz,exy,exz,eyz): 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(ezz) + hooke[0,3]*(2*exy) + hooke[0,4]*(2*exz) + hooke[0,5]*(2*eyz)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(ezz) + hooke[1,3]*(2*exy) + hooke[1,4]*(2*exz) + hooke[1,5]*(2*eyz)
        szz = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(ezz) + hooke[2,3]*(2*exy) + hooke[2,4]*(2*exz) + hooke[2,5]*(2*eyz)
        sxy = hooke[3,0]*(exx) + hooke[3,1]*(eyy) + hooke[3,2]*(ezz) + hooke[3,3]*(2*exy) + hooke[3,4]*(2*exz) + hooke[3,5]*(2*eyz)
        sxz = hooke[4,0]*(exx) + hooke[4,1]*(eyy) + hooke[4,2]*(ezz) + hooke[4,3]*(2*exy) + hooke[4,4]*(2*exz) + hooke[4,5]*(2*eyz)
        syz = hooke[5,0]*(exx) + hooke[5,1]*(eyy) + hooke[5,2]*(ezz) + hooke[5,3]*(2*exy) + hooke[5,4]*(2*exz) + hooke[5,5]*(2*eyz)
        return sxx,syy,szz,sxy,sxz,syz
    
    def ExportMeshToVtk(self, filename, output_fields_p=None, output_fields_c=None, file_format='vtk'):
        """ 
        output_fields_p : dictionary of point data 
        output_fields_c : dictionary of cell data 
        """
        points = self.n
        cells  =  [("tetra",self.e)]
        mesh = meshio.Mesh( points, cells ) 
        if output_fields_p is not None: 
            mesh.point_data = output_fields_p 
        if output_fields_c is not None:
            mesh.cell_data = output_fields_c 
        mesh.write(filename+'.'+file_format) 
        print(filename+'.'+file_format + ' written')
       
        
    
    def ExportGLRField(self, filename, dvc, res):
        """ 
        This method allows to create a new mesh 
        by merging all the mapped reference tetrahedrons 
        """
        # Taking one element 
        # Setting the nodes of the sub-triangles of the element in the physical space 
        pnodes = self.NnodesRef.dot( np.column_stack(self.n[self.e]) ) # Nodes of the reference tetrahedron in the physical space 
        list_meshes = [] 
        for i in range(self.e.shape[0]):
            e_pnodes = pnodes[:,3*i:3*(i+1)]  # projected nodes of the element 
            # faces = np.array([[0,1,3],[0,2,3],[0,1,2],[1,2,3]]) 
            local_mesh     = pymesh.form_mesh(e_pnodes, faces = None,  voxels=self.elemsRef)
            list_meshes.append(local_mesh)
        merged_mesh = pymesh.merge_meshes(list_meshes)
        # Export merged mesh into vtk file 
        points = merged_mesh.vertices
        cells  =  {"tetra": merged_mesh.voxels }
        mesh = meshio.Mesh( points, cells )
        mesh.cell_data = { "f":  dvc.fip, "f-go(id+u)": res   }
        # mesh.remove_lower_dimensional_cells()
        file_format = 'vtk'
        mesh.write(filename+'.'+file_format) 
        print(filename+'.'+file_format + ' written')   
        
    
    def LocatePoints(self,x,y,z,eps=1.e-15):
        """
        Finds to which element belongs the points (x,y,z)
        Method1: 
        loop over point, loop over element 
        compute barycentric coordinates if all >0 the point is in the tetrahedron
        
        """
        ie = tr.LocatePointsInTetraFE_Mesh1(x, y, z, self.e, self.n, len(x),eps )
        return ie 
    def LocatePoints2(self,x,y,z,eps=1.e-15):
        """
        Finds to which element belongs the points (x,y,z)
        Method2: 
        loop over point, loop over element 
        computes the tetrahedron faces level-set if all <0 the point is in the tetrahedron 
        """
        ie = tr.LocatePointsInTetraFE_Mesh2(x, y, z, self.e, self.n, len(x), eps )
        return ie 
    def LocatePoints3(self,x,y,z,eps=1.e-15):
        """
        Finds the nearest face of the tetrahedron 
        Then peforms element test location 
        """
        pts   = np.c_[x,y,z]
        _, fi, _ = pymesh.distance_to_mesh(self.faceMesh, pts)
        fi = fi.reshape(-1)
        # Now that we have the closest faces to the point
        # the query for each point is performed only for at most 2 tetrahedrons that 
        # share the closest face of each point 
        ie = tr.LocatePointsInTetraFE_Mesh3(x,y,z,fi,self.e,self.n,self.connFaces,len(x),eps)
        return ie 
    
    
    def GetJacobianOfElements(self,ie):
        """ 
        Returns the terms of the Jacobian of each fe element 
        and the determinant 
        """        
        x1 = self.n[self.e[ie,0],0]   ; y1 = self.n[self.e[ie,0],1] ; z1 = self.n[self.e[ie,0],2]
        x2 = self.n[self.e[ie,1],0]   ; y2 = self.n[self.e[ie,1],1] ; z2 = self.n[self.e[ie,1],2]
        x3 = self.n[self.e[ie,2],0]   ; y3 = self.n[self.e[ie,2],1] ; z3 = self.n[self.e[ie,2],2]
        x4 = self.n[self.e[ie,3],0]   ; y4 = self.n[self.e[ie,3],1] ; z4 = self.n[self.e[ie,3],2] 
        
        
        dxdxi   =  x1-x4
        dxdeta  =  x2-x4 
        dxdzeta =  x3-x4
        
        dydxi   =  y1-y4
        dydeta  =  y2-y4
        dydzeta =  y3-y4
        
        dzdxi   =  z1-z4
        dzdeta  =  z2-z4
        dzdzeta =  z3-z4
        

        det = dxdxi*dydeta*dzdzeta + \
              dxdeta*dydzeta*dzdxi + \
              dxdzeta*dydxi*dzdeta - \
              dxdzeta*dydeta*dzdxi - \
              dxdeta*dydxi*dzdzeta - \
              dxdxi*dydzeta*dzdeta
        
        return dxdxi, dxdeta, dxdzeta,  dydxi, dydeta, dydzeta,  dzdxi, dzdeta, dzdzeta, det 
        
        
    
    def InverseIsoparamMapping(self, x, y, z, ie, returnJacobian=False):
        """ Inversion of the isoparametric transformation 
            between the reference tetrahedron and the physical tetrahedron 
        """
        # ie[i] is the index of the element containing (x[i],y[i],z[i]) 
        
        dxdxi, dxdeta, dxdzeta,  dydxi, dydeta, dydzeta,  dzdxi, dzdeta, dzdzeta, det = self.GetJacobianOfElements(ie)
        
        
        detInv = 1/det 
        
        x4 = self.n[self.e[ie,3],0]   ; y4 = self.n[self.e[ie,3],1] ; z4 = self.n[self.e[ie,3],2] 

        x_x4 = x-x4
        y_y4 = y-y4 
        z_z4 = z-z4
        
        xi   =  detInv*(  (dydeta*dzdzeta-dydzeta*dzdeta)*x_x4  - (dxdeta*dzdzeta-dxdzeta*dzdeta)*y_y4 + (dxdeta*dydzeta-dxdzeta*dydeta)*z_z4 )
        eta  =  detInv*(  -(dydxi*dzdzeta-dydzeta*dzdxi)*x_x4  + (dxdxi*dzdzeta-dxdzeta*dzdxi)*y_y4 - (dxdxi*dydzeta-dxdzeta*dydxi)*z_z4 )
        zeta =  detInv*(  (dydxi*dzdeta-dydeta*dzdxi)*x_x4  -(dxdxi*dzdeta-dxdeta*dzdxi)*y_y4 + (dxdxi*dydeta-dxdeta*dydxi)*z_z4 )

        if returnJacobian == False :
            return xi,eta,zeta
        else: 
            return xi,eta,zeta,  dxdxi, dxdeta, dxdzeta,  dydxi, dydeta, dydzeta,  dzdxi, dzdeta, dzdzeta, det 

                            
        
    def EvalUOnPoints(self, x,y,z,ie, U):
        xi,eta,zeta = self.InverseIsoparamMapping(x, y, z, ie)
        N1 = xi ; N2 = eta ; N3 = zeta ; N4 = 1-xi-eta-zeta ; 
 
        ux = N1*U[self.conn[self.e[ie,0],0]] + \
             N2*U[self.conn[self.e[ie,1],0]] + \
             N3*U[self.conn[self.e[ie,2],0]] + \
             N4*U[self.conn[self.e[ie,3],0]] 
             
        uy = N1*U[self.conn[self.e[ie,0],1]] + \
             N2*U[self.conn[self.e[ie,1],1]] + \
             N3*U[self.conn[self.e[ie,2],1]] + \
             N4*U[self.conn[self.e[ie,3],1]]   
 
        uz = N1*U[self.conn[self.e[ie,0],2]] + \
             N2*U[self.conn[self.e[ie,1],2]] + \
             N3*U[self.conn[self.e[ie,2],2]] + \
             N4*U[self.conn[self.e[ie,3],2]] 

        return ux,uy,uz 

    def EvalUGradUOnPoints(self, x,y,z,ie, U, returnIsoParamCoord=False):
        xi,eta,zeta,  dxdxi, dxdeta, dxdzeta,  dydxi, dydeta, dydzeta,  dzdxi, dzdeta, dzdzeta, detJ \
        = self.InverseIsoparamMapping(x, y, z, ie, returnJacobian=True)
        
        N1 = xi ; N2 = eta ; N3 = zeta ; N4 = 1-xi-eta-zeta ; 
        
        
        dN1dxi   =  1 *np.ones(len(xi)) 
        # dN1deta  =  0
        # dN1dzeta =  0
        
        # dN2dxi   =  0
        dN2deta  =  1 *np.ones(len(xi))
        # dN2dzeta =  0 
        
        # dN3dxi   =  0 
        # dN3deta  =  0
        dN3dzeta =  1 *np.ones(len(xi))
        
        dN4dxi   =  -1 *np.ones(len(xi))
        dN4deta  =  -1 *np.ones(len(xi)) 
        dN4dzeta =  -1 *np.ones(len(xi))
        
        
        dN1dx = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dN1dxi  
        dN1dy = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dN1dxi 
        dN1dz = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dN1dxi 
                                 
                                 
        dN2dx =  (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dN2deta
        dN2dy =  ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dN2deta 
        dN2dz =  (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dN2deta
 
                             
        dN3dx = ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dN3dzeta 
        dN3dy = (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dN3dzeta 
        dN3dz = ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dN3dzeta      
                                 
 
        dN4dx = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dN4dxi + \
                (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dN4deta + \
                ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dN4dzeta  
                                 
        dN4dy = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dN4dxi + \
                ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dN4deta + \
                (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dN4dzeta  
                                 
        dN4dz = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dN4dxi + \
                (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dN4deta + \
                ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dN4dzeta      
                                 
 

        ux = N1*U[self.conn[self.e[ie,0],0]] + \
             N2*U[self.conn[self.e[ie,1],0]] + \
             N3*U[self.conn[self.e[ie,2],0]] + \
             N4*U[self.conn[self.e[ie,3],0]] 
             
        duxdx = dN1dx*U[self.conn[self.e[ie,0],0]] + \
                dN2dx*U[self.conn[self.e[ie,1],0]] + \
                dN3dx*U[self.conn[self.e[ie,2],0]] + \
                dN4dx*U[self.conn[self.e[ie,3],0]]     

        duxdy = dN1dy*U[self.conn[self.e[ie,0],0]] + \
                dN2dy*U[self.conn[self.e[ie,1],0]] + \
                dN3dy*U[self.conn[self.e[ie,2],0]] + \
                dN4dy*U[self.conn[self.e[ie,3],0]]  

        duxdz = dN1dz*U[self.conn[self.e[ie,0],0]] + \
                dN2dz*U[self.conn[self.e[ie,1],0]] + \
                dN3dz*U[self.conn[self.e[ie,2],0]] + \
                dN4dz*U[self.conn[self.e[ie,3],0]]                   
             
             

        uy = N1*U[self.conn[self.e[ie,0],1]] + \
             N2*U[self.conn[self.e[ie,1],1]] + \
             N3*U[self.conn[self.e[ie,2],1]] + \
             N4*U[self.conn[self.e[ie,3],1]] 
             
        duydx = dN1dx*U[self.conn[self.e[ie,0],1]] + \
                dN2dx*U[self.conn[self.e[ie,1],1]] + \
                dN3dx*U[self.conn[self.e[ie,2],1]] + \
                dN4dx*U[self.conn[self.e[ie,3],1]]     

        duydy = dN1dy*U[self.conn[self.e[ie,0],1]] + \
                dN2dy*U[self.conn[self.e[ie,1],1]] + \
                dN3dy*U[self.conn[self.e[ie,2],1]] + \
                dN4dy*U[self.conn[self.e[ie,3],1]]  

        duydz = dN1dz*U[self.conn[self.e[ie,0],1]] + \
                dN2dz*U[self.conn[self.e[ie,1],1]] + \
                dN3dz*U[self.conn[self.e[ie,2],1]] + \
                dN4dz*U[self.conn[self.e[ie,3],1]]               

                

        uz = N1*U[self.conn[self.e[ie,0],2]] + \
             N2*U[self.conn[self.e[ie,1],2]] + \
             N3*U[self.conn[self.e[ie,2],2]] + \
             N4*U[self.conn[self.e[ie,3],2]] 
             
        duzdx = dN1dx*U[self.conn[self.e[ie,0],2]] + \
                dN2dx*U[self.conn[self.e[ie,1],2]] + \
                dN3dx*U[self.conn[self.e[ie,2],2]] + \
                dN4dx*U[self.conn[self.e[ie,3],2]]     

        duzdy = dN1dy*U[self.conn[self.e[ie,0],2]] + \
                dN2dy*U[self.conn[self.e[ie,1],2]] + \
                dN3dy*U[self.conn[self.e[ie,2],2]] + \
                dN4dy*U[self.conn[self.e[ie,3],2]]  

        duzdz = dN1dz*U[self.conn[self.e[ie,0],2]] + \
                dN2dz*U[self.conn[self.e[ie,1],2]] + \
                dN3dz*U[self.conn[self.e[ie,2],2]] + \
                dN4dz*U[self.conn[self.e[ie,3],2]]  
                
        if returnIsoParamCoord==False: 
            return ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz  
        else: 
            return xi,eta,zeta, ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz  


    
    def EvalUGradUOnGaussPoints(self,U):
        """
        Strain is constant per tetrahedral element, 
        therefore we compute it at the barycenter of each element 
        """
        npg = self.e.shape[0] # number of elements 
        ie = np.arange(npg)
        dxdxi, dxdeta, dxdzeta,  dydxi, dydeta, dydzeta,  dzdxi, dzdeta, dzdzeta, detJ  = self.GetJacobianOfElements(ie) 
        
        N1 = 1/4 
        N2 = 1/4 
        N3 = 1/4 
        N4 = 1/4
        
        dN1dxi   =  1 *np.ones(npg) 
        
        dN2deta  =  1 *np.ones(npg)
        
        dN3dzeta =  1 *np.ones(npg)
        
        dN4dxi   =  -1 *np.ones(npg)
        dN4deta  =  -1 *np.ones(npg) 
        dN4dzeta =  -1 *np.ones(npg)


        dN1dx = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dN1dxi  
        dN1dy = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dN1dxi 
        dN1dz = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dN1dxi 
                                 
                                 
        dN2dx =  (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dN2deta
        dN2dy =  ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dN2deta 
        dN2dz =  (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dN2deta
 
                             
        dN3dx = ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dN3dzeta 
        dN3dy = (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dN3dzeta 
        dN3dz = ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dN3dzeta      
                                 
 
        dN4dx = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dN4dxi + \
                (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dN4deta + \
                ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dN4dzeta  
                                 
        dN4dy = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dN4dxi + \
                ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dN4deta + \
                (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dN4dzeta  
                                 
        dN4dz = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dN4dxi + \
                (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dN4deta + \
                ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dN4dzeta          
        

        ux = N1*U[self.conn[self.e[ie,0],0]] + \
             N2*U[self.conn[self.e[ie,1],0]] + \
             N3*U[self.conn[self.e[ie,2],0]] + \
             N4*U[self.conn[self.e[ie,3],0]] 
             
        duxdx = dN1dx*U[self.conn[self.e[ie,0],0]] + \
                dN2dx*U[self.conn[self.e[ie,1],0]] + \
                dN3dx*U[self.conn[self.e[ie,2],0]] + \
                dN4dx*U[self.conn[self.e[ie,3],0]]     

        duxdy = dN1dy*U[self.conn[self.e[ie,0],0]] + \
                dN2dy*U[self.conn[self.e[ie,1],0]] + \
                dN3dy*U[self.conn[self.e[ie,2],0]] + \
                dN4dy*U[self.conn[self.e[ie,3],0]]  

        duxdz = dN1dz*U[self.conn[self.e[ie,0],0]] + \
                dN2dz*U[self.conn[self.e[ie,1],0]] + \
                dN3dz*U[self.conn[self.e[ie,2],0]] + \
                dN4dz*U[self.conn[self.e[ie,3],0]]                   
             
             

        uy = N1*U[self.conn[self.e[ie,0],1]] + \
             N2*U[self.conn[self.e[ie,1],1]] + \
             N3*U[self.conn[self.e[ie,2],1]] + \
             N4*U[self.conn[self.e[ie,3],1]] 
             
        duydx = dN1dx*U[self.conn[self.e[ie,0],1]] + \
                dN2dx*U[self.conn[self.e[ie,1],1]] + \
                dN3dx*U[self.conn[self.e[ie,2],1]] + \
                dN4dx*U[self.conn[self.e[ie,3],1]]     

        duydy = dN1dy*U[self.conn[self.e[ie,0],1]] + \
                dN2dy*U[self.conn[self.e[ie,1],1]] + \
                dN3dy*U[self.conn[self.e[ie,2],1]] + \
                dN4dy*U[self.conn[self.e[ie,3],1]]  

        duydz = dN1dz*U[self.conn[self.e[ie,0],1]] + \
                dN2dz*U[self.conn[self.e[ie,1],1]] + \
                dN3dz*U[self.conn[self.e[ie,2],1]] + \
                dN4dz*U[self.conn[self.e[ie,3],1]]               

                

        uz = N1*U[self.conn[self.e[ie,0],2]] + \
             N2*U[self.conn[self.e[ie,1],2]] + \
             N3*U[self.conn[self.e[ie,2],2]] + \
             N4*U[self.conn[self.e[ie,3],2]] 
             
        duzdx = dN1dx*U[self.conn[self.e[ie,0],2]] + \
                dN2dx*U[self.conn[self.e[ie,1],2]] + \
                dN3dx*U[self.conn[self.e[ie,2],2]] + \
                dN4dx*U[self.conn[self.e[ie,3],2]]     

        duzdy = dN1dy*U[self.conn[self.e[ie,0],2]] + \
                dN2dy*U[self.conn[self.e[ie,1],2]] + \
                dN3dy*U[self.conn[self.e[ie,2],2]] + \
                dN4dy*U[self.conn[self.e[ie,3],2]]  

        duzdz = dN1dz*U[self.conn[self.e[ie,0],2]] + \
                dN2dz*U[self.conn[self.e[ie,1],2]] + \
                dN3dz*U[self.conn[self.e[ie,2],2]] + \
                dN4dz*U[self.conn[self.e[ie,3],2]]  
                

        return ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz  
 
    def StrainAtGaussPoints2(self,U)  :
        ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz = self.EvalUGradUOnGaussPoints(U)
        exx = duxdx 
        eyy = duydy 
        ezz = duzdz 
        exy = 0.5*(duydx+duxdy)
        exz = 0.5*(duzdx+duxdz)
        eyz = 0.5*(duzdy+duydz)
        return exx,eyy,ezz,exy,exz,eyz    
    
    def DispAndStrainAtPoints(self,x,y,z,ie,U,returnIsoParamCoord=False):
        if returnIsoParamCoord == False:
            ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz = self.EvalUGradUOnPoints(x, y, z, ie, U,returnIsoParamCoord) 
            exx = duxdx 
            eyy = duydy 
            ezz = duzdz 
            exy = 0.5*(duydx+duxdy)
            exz = 0.5*(duzdx+duxdz)
            eyz = 0.5*(duzdy+duydz)
            return ux,uy,uz,exx,eyy,ezz,exy,exz,eyz    
        else: 
            xi, eta, zeta, ux,uy,uz, duxdx, duxdy, duxdz, duydx, duydy, duydz, duzdx, duzdy, duzdz = self.EvalUGradUOnPoints(x, y, z, ie, U,returnIsoParamCoord) 
            exx = duxdx 
            eyy = duydy 
            ezz = duzdz 
            exy = 0.5*(duydx+duxdy)
            exz = 0.5*(duzdx+duxdz)
            eyz = 0.5*(duzdy+duydz)            
            return xi,eta,zeta,ux,uy,uz,exx,eyy,ezz,exy,exz,eyz   

#%%% Vectorized rootines 
    def LocatePointsInMesh(self, x,y,z):
        
        pts   = np.c_[x,y,z]
 
# =============================================================================
#         Fail attempt to locate a 3d point in tetrahedral finite element mesh 
#         ip1   =  self.e[:,0] ; ip2 =  self.e[:,1] ; ip3 =  self.e[:,2] ; ip4 =  self.e[:,3] 
#         face1 =  np.c_[ip1,ip2,ip3]
#         face2 =  np.c_[ip1,ip2,ip4]
#         face3 =  np.c_[ip1,ip3,ip4]
#         face4 =  np.c_[ip2,ip3,ip4]
#         faces = np.vstack((face1,face2,face3,face4))
#         faceElem = np.kron(np.ones(4), np.arange(self.e.shape[0])).astype('int') # Knowing to which tetrahedron belongs a triangular face  
#         _,faces, info = pymesh.remove_duplicated_faces_raw(self.n, faces)
#         mesh = pymesh.form_mesh(self.n, faces , voxels=self.e)
#         faceElem = faceElem[info['ori_face_index']]
#         _, fi, _ = pymesh.distance_to_mesh(mesh, pts)
#         # We look for the elements defined by the closest faces 
# =============================================================================

        
        # Example locate a point in a finite element mesh 
        
        ie = 2135
        x1 = self.n[self.e[ie,0],0] ; y1 = self.n[self.e[ie,0],1] ; z1 = self.n[self.e[ie,0],2]
        x2 = self.n[self.e[ie,1],0] ; y2 = self.n[self.e[ie,1],1] ; z2 = self.n[self.e[ie,1],2]
        x3 = self.n[self.e[ie,2],0] ; y3 = self.n[self.e[ie,2],1] ; z3 = self.n[self.e[ie,2],2] 
        x4 = self.n[self.e[ie,3],0] ; y4 = self.n[self.e[ie,3],1] ; z4 = self.n[self.e[ie,3],2]
        
        
        x = 2 + 1.e-15
        y = 0 + 1.e-15
        z = 0 + 1.e-15
        
        x1 = 0 ; y1 = 0 ; z1 = 0
        x2 = 1 ; y2 = 0 ; z2 = 0
        x3 = 0 ; y3 = 1 ; z3 = 0 
        x4 = 0 ; y4 = 0 ; z4 = 1 
        
        A = np.array([[x1,x2,x3,x4],
                      [y1,y2,y3,y4],
                      [z1,z2,z3,z4],
                      [1, 1, 1, 1]])
        b = np.array([x,y,z,1])
        l = np.linalg.solve(A,b)
        if (l> 0).all() :
            print('Point is inside element ')
            
 
    
    def GaussIntegration(self):
        """ Gauss integration: One point per element""" 
        # The integration points and the elements are aranged in the same order 
        x1 = self.n[self.e[:,0],0]   ; y1 = self.n[self.e[:,0],1] ; z1 = self.n[self.e[:,0],2]
        x2 = self.n[self.e[:,1],0]   ; y2 = self.n[self.e[:,1],1] ; z2 = self.n[self.e[:,1],2]
        x3 = self.n[self.e[:,2],0]   ; y3 = self.n[self.e[:,2],1] ; z3 = self.n[self.e[:,2],2]
        x4 = self.n[self.e[:,3],0]   ; y4 = self.n[self.e[:,3],1] ; z4 = self.n[self.e[:,3],2] 
        
        self.pix =  (x1+x2+x3+x4)/4
        self.piy =  (y1+y2+y3+y4)/4
        self.piz =  (z1+z2+z3+z4)/4
        
        self.npg = len(self.pix)
        
        # Basis function and derivatives evaluated at the integration point (reference element) 
        N_ref       = np.array([1/4,1/4,1/4,1/4])
        dNdxi_ref   = np.array([-1,1,0,0])
        dNdeta_ref  = np.array([-1,0,1,0])
        dNdzeta_ref = np.array([-1,0,0,1]) 
        
        indexI = np.kron( np.arange(self.npg),np.ones(4))
        indexJ = self.conn[self.e.ravel(),0]
        
        values_phi       = np.kron( np.ones(self.npg),  N_ref )
        
        values_dphidxi   = np.kron( np.ones(self.npg),  dNdxi_ref )
        values_dphideta  = np.kron( np.ones(self.npg),  dNdeta_ref )
        values_dphidzeta =  np.kron( np.ones(self.npg), dNdzeta_ref )     
        
        phi       = sps.csc_matrix(( values_phi, (indexI,indexJ)), shape = (self.npg,self.ndof//3) )
        dphidxi   = sps.csc_matrix(( values_dphidxi, (indexI,indexJ)), shape = (self.npg,self.ndof//3) )
        dphideta  = sps.csc_matrix(( values_dphideta, (indexI,indexJ)), shape = (self.npg,self.ndof//3) )
        dphidzeta = sps.csc_matrix(( values_dphidzeta, (indexI,indexJ)), shape = (self.npg,self.ndof//3) )
        

        dxdxi   = x2-x1
        dxdeta  = x3-x1 
        dxdzeta = x4-x1
        
        dydxi   = y2-y1 
        dydeta  = y3-y1 
        dydzeta = y4-y1
        
        dzdxi   = z2-z1 
        dzdeta  = z3-z1 
        dzdzeta = z4-z1 
 
 
        detJ   = dxdxi*dydeta*dzdzeta + \
                 dydxi*dzdeta*dxdzeta + \
                 dzdxi*dxdeta*dydzeta - \
                 dzdxi*dydeta*dxdzeta - \
                 dydxi*dxdeta*dzdzeta - \
                 dxdxi*dzdeta*dydzeta
 
           
        dphidx = sps.diags((dydeta*dzdzeta-dzdeta*dydzeta)/detJ).dot(dphidxi) +\
                 sps.diags(-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ).dot(dphideta) +\
                 sps.diags((dydxi*dzdeta-dzdxi*dydeta)/detJ).dot(dphidzeta)
            
    
        dphidy = sps.diags((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ).dot(dphidxi) +\
                 sps.diags((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ).dot(dphideta) +\
                 sps.diags(-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ).dot(dphidzeta)
                 
        dphidz = sps.diags((dxdeta*dydzeta-dydeta*dxdzeta)/detJ).dot(dphidxi) +\
                 sps.diags(-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ).dot(dphideta) +\
                 sps.diags((dxdxi*dydeta-dydxi*dxdeta)/detJ).dot(dphidzeta)
                 
        wg = (1./6)*np.abs(detJ) 
                     
                     
        self.wg = sps.diags(wg)        
        zero    = sps.csr_matrix((self.npg,self.ndof//3))
        
        self.phix    = sps.hstack((phi,zero,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi,zero)   ,  'csc')
        self.phiz    = sps.hstack((zero,zero,phi)  ,  'csc')
        
        self.dphixdx = sps.hstack((dphidx,zero,zero),  'csc')
        self.dphixdy = sps.hstack((dphidy,zero,zero),  'csc')
        self.dphixdz = sps.hstack((dphidz,zero,zero),  'csc')
        
        self.dphiydx = sps.hstack((zero,dphidx,zero),  'csc')
        self.dphiydy = sps.hstack((zero,dphidy,zero),  'csc')
        self.dphiydz = sps.hstack((zero,dphidz,zero),  'csc')
        
        self.dphizdx = sps.hstack((zero,zero,dphidx),  'csc')
        self.dphizdy = sps.hstack((zero,zero,dphidy),  'csc')
        self.dphizdz = sps.hstack((zero,zero,dphidz),  'csc')                      
                     
    def Laplacian(self):
        L = self.dphixdx.T.dot(self.wg.dot(self.dphixdx)) + \
            self.dphixdy.T.dot(self.wg.dot(self.dphixdy)) + \
            self.dphixdz.T.dot(self.wg.dot(self.dphixdz)) + \
            self.dphiydx.T.dot(self.wg.dot(self.dphiydx)) + \
            self.dphiydy.T.dot(self.wg.dot(self.dphiydy)) + \
            self.dphiydz.T.dot(self.wg.dot(self.dphiydz)) + \
            self.dphizdx.T.dot(self.wg.dot(self.dphizdx)) + \
            self.dphizdy.T.dot(self.wg.dot(self.dphizdy)) + \
            self.dphizdz.T.dot(self.wg.dot(self.dphizdz))      
        return L  

    def Stiffness(self,hooke):
       
        Bxy = self.dphixdy + self.dphiydx
        Bxz = self.dphixdz + self.dphizdx
        Byz = self.dphiydz + self.dphizdy 
        
        K =  hooke[0,0]*self.dphixdx.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[0,1]*self.dphiydy.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[0,2]*self.dphizdz.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[1,0]*self.dphixdx.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[1,1]*self.dphiydy.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[1,2]*self.dphizdz.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[2,0]*self.dphixdx.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[2,1]*self.dphiydy.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[2,2]*self.dphizdz.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[3,3]*Bxy.T.dot(self.wg.dot(Bxy)) + \
             hooke[4,4]*Bxz.T.dot(self.wg.dot(Bxz)) + \
             hooke[5,5]*Byz.T.dot(self.wg.dot(Byz))
        return K 

    def StrainAtGaussPoints(self,U):
        # n = self.ndof//3 
        # Ux = U[:n]; Uy =U[n:2*n]; Uz=U[2*n:]
        exx = self.dphixdx.dot(U)
        eyy = self.dphiydy.dot(U) 
        ezz = self.dphizdz.dot(U)
        exy = 0.5*(self.dphiydx.dot(U)+self.dphixdy.dot(U))
        exz = 0.5*(self.dphizdx.dot(U)+self.dphixdz.dot(U))
        eyz = 0.5*(self.dphizdy.dot(U)+self.dphiydz.dot(U))
        return exx,eyy,ezz,exy,exz,eyz 
#%% Structured Mesh 
class StructuredMesh3d: 
    def __init__(self,degree):
        self.pp    = degree 
        self.n_elems   = None 
        self.xxsi      = None    
        self.noelem    = None 
        self.pix = 0 # x coordinates of the integration points 
        self.piy = 0 # y coordinates of the integration points  
        self.piz = 0 # z coordinates of the integration points 
        self.wg = 0        # Integration weights diagonal matrix 
 
 
        """ fixed parameters for integration """
        self.nbg_xi  =  0
        self.nbg_eta  = 0
        self.nbg_zeta  = 0
        
        self.oneXi   = 0  
        self.oneEta  = 0
        self.oneZeta = 0  
        
        self.Gauss_xi  =  0
        self.Gauss_eta =  0
        self.Gauss_zeta =  0
        self.wgRef =  0
        self.GaussTetra =  0
        
        # Saved univariate basis functions and derivatives 
        self.Nxi = 0 
        self.dNxidxi = 0 
        self.Neta = 0 
        self.dNetadeta = 0 
        self.Nzeta = 0 
        self.dNzetadzeta = 0 
        
        
        self.integrationCellsCoord = 0 #list of integration cells: rectangular cuboids and tetrahedrons for plotting and verification  
        

    def Get_nbf(self):
        return self.Btot.shape[1]  
    def Get_Btot(self):
        return self.Btot 
    
    
    def GetBoundaryIndices(self,eps=1.e-12):
        P = self.Get_Btot() 
        b1 = np.where( np.abs ( P[0,:]  - np.min(P[0,:])  )  < eps  )[0] # x=xmin plane 1
        b2 = np.where( np.abs ( P[0,:]  - np.max(P[0,:])  )  < eps  )[0] # x=xmax plane 2
        b3 = np.where( np.abs ( P[1,:]  - np.min(P[1,:])  )  < eps  )[0] # y=ymin plane 3
        b4 = np.where( np.abs ( P[1,:]  - np.max(P[1,:])  )  < eps  )[0] # y=ymax plane 4
        b5 = np.where( np.abs ( P[2,:]  - np.min(P[2,:])  )  < eps  )[0] # z=zmin plane 5
        b6 = np.where( np.abs ( P[2,:]  - np.max(P[2,:])  )  < eps  )[0] # z=zmax plane 6
        return b1,b2,b3,b4,b5,b6 
    
    def GetBoundaryOperator(self,planes=[1,1,1,1,1,1]):
        b = self.GetBoundaryIndices()
        ib = np.array([])
        nbf = self.Get_nbf() ; ndof = 3*nbf 
        for i in range(6):
            if planes[i]==1:
                ib = np.r_[ib,b[i]]
        ib = np.unique(ib) 
        index = np.r_[ib, ib+nbf, ib+2*nbf]
        D = tools.GetBinaryMatrix(index, ndof)
        return D
            
        
    
    
    def SetContinuousKnotVector(self,Xbounds,Ybounds,Zbounds,n_elems):
        self.n_elems = n_elems  
        xmin = Xbounds[0]  
        xmax = Xbounds[1]  
        ymin = Ybounds[0]  
        ymax = Ybounds[1]  
        zmin = Zbounds[0]  
        zmax = Zbounds[1]  
        self.xxsi = dict()
        self.xxsi[0] = np.array([xmin,xmin,xmax,xmax])    # Xi knot vector ( direction x)
        self.xxsi[1] = np.array([ymin,ymin,ymax,ymax])    # Eta knot vector ( direction y)   
        self.xxsi[2] = np.array([zmin,zmin,zmax,zmax])    # Zeta knot vector (direction z) 
        #Initial B-spline mesh 
        x = np.array([xmin,xmax]) ; x = x.reshape((1,-1)) 
        y = np.array([ymin,ymax]) ; y = y.reshape((1,-1))           
        z = np.array([zmin,zmax]) ; z = z.reshape((1,-1)) 
        """ Degree elevation """ 
        # Elevalting from degree 1 to the desired degree 
        x, self.xxsi[0] = nb.bspdegelev(1,x,self.xxsi[0],self.pp[0]-1)
        y, self.xxsi[1] = nb.bspdegelev(1,y,self.xxsi[1],self.pp[1]-1)  
        z, self.xxsi[2] = nb.bspdegelev(1,z,self.xxsi[2],self.pp[2]-1)          
        """ Knot refinement """ 
        ubar_xi   = np.linspace(xmin,xmax,self.n_elems[0]+1)[1:-1]
        ubar_eta  = np.linspace(ymin,ymax,self.n_elems[1]+1)[1:-1] 
        ubar_zeta = np.linspace(zmin,zmax,self.n_elems[2]+1)[1:-1]
        if (len(ubar_xi)!=0):
            x, self.xxsi[0] = nb.bspkntins(self.pp[0], x, self.xxsi[0], ubar_xi )
        if (len(ubar_eta)!=0):
            y, self.xxsi[1] = nb.bspkntins(self.pp[1], y, self.xxsi[1], ubar_eta )
        if (len(ubar_zeta)!=0):
            z, self.xxsi[2] = nb.bspkntins(self.pp[2], z, self.xxsi[2], ubar_zeta )
 
        x = x[0]; y = y[0]; z=z[0]
        nbf_xi = x.shape[0]; nbf_eta = y.shape[0]; nbf_zeta = z.shape[0]
 
        Xtot = np.kron(np.ones(nbf_eta*nbf_zeta), x)
        Ytot = np.kron(np.ones(nbf_zeta), np.kron(y,np.ones(nbf_xi)))
        Ztot = np.kron(z,np.ones(nbf_eta*nbf_xi))
        self.Btot = np.c_[Xtot,Ytot,Ztot].T
        self.noelem = nb.computeNOELEM3D(self.n_elems[0], self.n_elems[1], self.n_elems[2], self.pp[0], self.pp[1], self.pp[2])
        # self.X, self.Y, self.Z  = np.meshgrid(x,y,z)
        # self.Btot = np.c_[self.X.ravel(),self.Y.ravel(),self.Z.ravel()].T   
 
    
    def BezierUniformIntegrationRule(self,nx,ny,nz):
        """ Computes the Bezier basis functions 
        on the referece [-1,1]^3 cube 
        and the Bezier extraction operator
        maps the bezier control points as b-spline control points 
        """
        # Getting the quadrature rule 
        x1d = np.linspace(-1,1,2*nx+1)[1::2]
        y1d = np.linspace(-1,1,2*ny+1)[1::2]
        z1d = np.linspace(-1,1,2*nz+1)[1::2]
        # Getting the trivariate Bernstein basis functions 
        Bx, _ = nb.BernsteinRef(x1d,self.pp[0])
        By, _ = nb.BernsteinRef(y1d,self.pp[1])
        Bz, _ = nb.BernsteinRef(z1d,self.pp[1])
        B = np.kron(Bz,np.kron(By,Bx))
        # Getting the Bezier extraction operator 
        # for each element         
        Cx = nb.bezier_extraction_nurbs_1d( self.xxsi[0], len(self.xxsi[0]), self.pp[0] )
        Cy = nb.bezier_extraction_nurbs_1d( self.xxsi[1], len(self.xxsi[1]), self.pp[1] )
        Cz = nb.bezier_extraction_nurbs_1d( self.xxsi[2], len(self.xxsi[2]), self.pp[2] )
        ne = np.product(self.n_elems)
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1)
        C = np.zeros((ne,nbf_elem,nbf_elem))
        ie=0
        for k in range(self.n_elems[2]):
            for j in range(self.n_elems[1]):
                for i in range(self.n_elems[0]):
                    C[ie, :,:] = np.kron(Cz[k,:,:],np.kron(Cx[j,:,:],Cy[i,:,:]))
                    ie+=1 
        self.C = C 
        self.B = B 

                
    def SetUnivariateBasisFunctionsAndDerivatives(self,nbippe):

        xi   = np.linspace(self.xxsi[0][0], self.xxsi[0][-1], 2*self.n_elems[0]*nbippe[0]+1)[1::2] 
        eta  = np.linspace(self.xxsi[1][0], self.xxsi[1][-1], 2*self.n_elems[1]*nbippe[1]+1)[1::2] 
        zeta = np.linspace(self.xxsi[2][0], self.xxsi[2][-1], 2*self.n_elems[2]*nbippe[2]+1)[1::2] 
        
        self.pix1d = xi 
        self.piy1d = eta 
        self.piz1d = zeta 
        
        nbf_xi   = len(self.xxsi[0]) -1 - self.pp[0] ; 
        nbf_eta  = len(self.xxsi[1]) -1 - self.pp[1] ;   
        nbf_zeta = len(self.xxsi[2]) -1 - self.pp[2] ;
        
        self.Nxi         = np.zeros((len(xi),self.pp[0]+1))
        self.dNxidxi     = np.zeros((len(xi),self.pp[0]+1)) 
        self.Neta        = np.zeros((len(eta),self.pp[1]+1))
        self.dNetadeta   = np.zeros((len(eta),self.pp[1]+1))
        self.Nzeta       = np.zeros((len(zeta),self.pp[2]+1)) 
        self.dNzetadzeta = np.zeros((len(zeta),self.pp[2]+1)) 
        
        for i in range(len(xi)):
            span_xi  = nb.findspan(nbf_xi,  self.pp[0], xi[i] , self.xxsi[0])
            self.Nxi[i,:], self.dNxidxi[i,:] = nb.derbasisfuns(span_xi,self.pp[0],self.xxsi[0],1,xi[i])
        for i in range(len(eta)):
            span_eta = nb.findspan(nbf_eta, self.pp[1], eta[i], self.xxsi[1])
            self.Neta[i,:], self.dNetadeta[i,:] = nb.derbasisfuns(span_eta,self.pp[1],self.xxsi[1],1,eta[i])
        for i in range(len(zeta)):
            span_zeta = nb.findspan(nbf_zeta, self.pp[2], zeta[i], self.xxsi[2])
            self.Nzeta[i,:], self.dNzetadzeta[i,:] = nb.derbasisfuns(span_zeta,self.pp[2],self.xxsi[2],1,zeta[i])            
        
        nip_xi   =  self.n_elems[0]*nbippe[0] 
        nip_eta  =  self.n_elems[1]*nbippe[1]  
        nip_zeta =  self.n_elems[2]*nbippe[2]    

        self.pix = np.kron(np.ones(nip_zeta*nip_eta),xi)
        self.piy = np.kron(np.ones(nip_zeta),np.kron(eta, np.ones(nip_xi)))
        self.piz = np.kron(zeta,np.ones(nip_eta*nip_xi))
        self.npg = self.pix.shape[0]         

    def ComputeUnivariateIntegratedDiffTensors(self, nip, ne):
      
        """ to be translated to C++ """
        # nip : [nipx,nipy] number of integration points per integration element 
        # ne : [nex,ney] number of integration elements per knot element (knot span) 
        
        e_xi   = np.unique(self.xxsi[0]) 
        e_eta  = np.unique(self.xxsi[1])   
        e_zeta = np.unique(self.xxsi[2])  
 
        nbf_elem_xi    =  self.pp[0]+1
        nbf_elem_eta   =  self.pp[1]+1
        nbf_elem_zeta  =  self.pp[2]+1
        
 
        self.dxi_dxi = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.dxi_xi    = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.xi_dxi    = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.xi_xi       = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))

        self.deta_deta = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.deta_eta     = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.eta_deta      = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.eta_eta           = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        
        self.dzeta_dzeta = np.zeros((self.n_elems[2]*ne[2],nbf_elem_zeta,nbf_elem_zeta))
        self.dzeta_zeta       = np.zeros((self.n_elems[2]*ne[2],nbf_elem_zeta,nbf_elem_zeta))
        self.zeta_dzeta      = np.zeros((self.n_elems[2]*ne[2],nbf_elem_zeta,nbf_elem_zeta))
        self.zeta_zeta             = np.zeros((self.n_elems[2]*ne[2],nbf_elem_zeta,nbf_elem_zeta))        
        
        Gauss_xi   =  nb.GaussLegendre(nip[0])
        Gauss_eta  =  nb.GaussLegendre(nip[1])
        Gauss_zeta =  nb.GaussLegendre(nip[2])
        
        xig    = Gauss_xi[0]
        wgxi   = Gauss_xi[1]
        etag   = Gauss_eta[0]
        wgeta  = Gauss_eta[1]
        zetag  = Gauss_zeta[0]
        wgzeta = Gauss_zeta[1]        
        
        
        mes_xi   =  (e_xi[1]  - e_xi[0] )   / ne[0]
        mes_eta  =  (e_eta[1] - e_eta[0])   / ne[1]
        mes_zeta =  (e_zeta[1] - e_zeta[0]) / ne[2]
 
        
        # Loop over basis elements (xi direction)
        for i in range(self.n_elems[0]):
            # loop over integration elements of the current element 
            for j in range(ne[0]):
                xi_min = e_xi[i] + j*mes_xi  
                xi_p   = xi_min  + 0.5*(xig+1)*mes_xi
                for ip in range(nip[0]):
                     N, dN = nb.derbasisfuns(i+self.pp[0], self.pp[0], self.xxsi[0], 1, xi_p[ip])
                     self.dxi_dxi[j+i*ne[0],:,:] += wgxi[ip]*np.outer(dN,dN)*mes_xi/2
                     self.dxi_xi[j+i*ne[0],:,:]  += wgxi[ip]*np.outer(dN,N)*mes_xi/2
                     self.xi_dxi[j+i*ne[0],:,:]  += wgxi[ip]*np.outer(N,dN)*mes_xi/2  
                     self.xi_xi[j+i*ne[0],:,:]   += wgxi[ip]*np.outer(N,N)*mes_xi/2  
                     
        # Loop over basis elements (eta direction)
        for i in range(self.n_elems[1]): 
            # loop over integration elements of the current element 
            for j in range(ne[1]):
                eta_min = e_eta[i] + j*mes_eta 
                eta_p   = eta_min  + 0.5*(etag+1)*mes_eta
                for ip in range(nip[1]):
                     N, dN = nb.derbasisfuns(i+self.pp[1], self.pp[1], self.xxsi[1], 1, eta_p[ip])
                     self.deta_deta[j+i*ne[1],:,:] += wgeta[ip]*np.outer(dN,dN)*mes_eta/2
                     self.deta_eta[j+i*ne[1],:,:]  += wgeta[ip]*np.outer(dN,N)*mes_eta/2
                     self.eta_deta[j+i*ne[1],:,:]  += wgeta[ip]*np.outer(N,dN)*mes_eta/2                      
                     self.eta_eta[j+i*ne[1],:,:]   += wgeta[ip]*np.outer(N,N)*mes_eta/2  
                     

        # Loop over basis elements (zeta direction)
        for i in range(self.n_elems[2]): 
            # loop over integration elements of the current element 
            for j in range(ne[2]):
                zeta_min = e_zeta[i] + j*mes_zeta 
                zeta_p   = zeta_min  + 0.5*(zetag+1)*mes_zeta
                for ip in range(nip[2]):
                     N, dN = nb.derbasisfuns(i+self.pp[2], self.pp[2], self.xxsi[2], 1, zeta_p[ip])
                     self.dzeta_dzeta[j+i*ne[2],:,:] += wgzeta[ip]*np.outer(dN,dN)*mes_zeta/2
                     self.dzeta_zeta[j+i*ne[2],:,:]  += wgzeta[ip]*np.outer(dN,N)*mes_zeta/2
                     self.zeta_dzeta[j+i*ne[2],:,:]  += wgzeta[ip]*np.outer(N,dN)*mes_zeta/2                      
                     self.zeta_zeta[j+i*ne[2],:,:]   += wgzeta[ip]*np.outer(N,N)*mes_zeta/2 
                     
                     

    def ComputeBivariateIntegratedDiffTensors(self):
        """ to be translated to C++ """
        nbf_elem_xi    =  self.pp[0]+1
        nbf_elem_eta   =  self.pp[1]+1
        nbf_elem2d     =  nbf_elem_eta*nbf_elem_xi
 
        nei_xi  = self.dxi_dxi.shape[0]//self.n_elems[0] # Number of integration elements per xi element 
        nei_eta = self.deta_deta.shape[0]//self.n_elems[1] # Number of integration elements per eta element 
        
        nei_xi_tot = self.dxi_dxi.shape[0]
        
        # Computes for each 2d element the products of the 1d integrals 
        self.eta_eta_dxi_dxi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.eta_deta_dxi_xi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.deta_eta_xi_dxi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.deta_deta_xi_xi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        
        
        # Loop over basis elements  
        for j in range(self.n_elems[1]):
            for i in range(self.n_elems[0]):
               # For each basis element 
               # Loop over integration elements 
               for ji in range(nei_eta):
                   ej = ji + j*nei_eta
                   for ii in range(nei_xi):
                       ei = ii + i*nei_xi 
                       iie = ei + ej*nei_xi_tot
                       self.eta_eta_dxi_dxi[iie,:,:] = np.kron(self.eta_eta[ej,:,:], self.dxi_dxi[ei,:,:])
                       self.eta_deta_dxi_xi[iie,:,:] = np.kron(self.eta_deta[ej,:,:], self.dxi_xi[ei,:,:]) 
                       self.deta_eta_xi_dxi[iie,:,:] = np.kron(self.deta_eta[ej,:,:], self.xi_dxi[ei,:,:]) 
                       self.deta_deta_xi_xi[iie,:,:] = np.kron(self.deta_deta[ej,:,:], self.xi_xi[ei,:,:])
                    
                      
    def ComputeIntegratedDiffTensors(self, nip, ne):
        
        nbf_elem_xi    =  self.pp[0]+1
        nbf_elem_eta   =  self.pp[1]+1
        nbf_elem_zeta  =  self.pp[2]+1
        
        nbf_elem2d     =  nbf_elem_eta*nbf_elem_xi
    
        s_xi     = self.n_elems[0]*ne[0]*nbf_elem_xi*nbf_elem_xi 
        s_eta    = self.n_elems[1]*ne[1]*nbf_elem_eta*nbf_elem_eta 
        s_zeta   = self.n_elems[2]*ne[2]*nbf_elem_zeta*nbf_elem_zeta 
        s_eta_xi = self.n_elems[0]*self.n_elems[1]*ne[0]*ne[1]*nbf_elem2d*nbf_elem2d 
  
        r = ar.ComputeIntegratedDiffTensors(self.xxsi[0], self.xxsi[1], self.xxsi[2], 
                                             self.pp[0], self.pp[1], self.pp[2],
                                             nip[0], nip[1], nip[2], ne[0], ne[1], ne[2],
                                             s_xi,s_xi,s_xi,s_xi,
                                             s_eta,s_eta,s_eta,s_eta,
                                             s_zeta, s_zeta, s_zeta, s_zeta, 
                                             s_eta_xi, s_eta_xi, s_eta_xi, s_eta_xi )  
 
        # dxi_dxi =  r[0]
        # dxi_xi  =  r[1]
        # xi_dxi  =  r[2]
        # xi_xi   =  r[3]    
 
        # deta_deta =  r[4]
        # deta_eta  =  r[5]
        # eta_deta  =  r[6]                 
        # eta_eta   =  r[7] 


        # dzeta_dzeta = r[8]
        # dzeta_zeta  = r[9]
        # zeta_dzeta  = r[10]               
        # zeta_zeta   = r[11]

        # eta_eta_dxi_dxi = r[12]
        # eta_deta_dxi_xi = r[13] 
        # deta_eta_xi_dxi = r[14] 
        # deta_deta_xi_xi = r[15]   
        
        # dxi_dxi, dxi_xi, xi_dxi, xi_xi, deta_deta, deta_eta, eta_deta, eta_eta, \   
        #        dzeta_dzeta, zeta_zeta, zeta_dzeta, zeta_zeta, \   
        #        eta_eta_dxi_dxi, eta_deta_dxi_xi, deta_eta_xi_dxi, deta_deta_xi_xi     
        
 
                     
        return  r  
 
        
                     
        
    
        
    def L2Projection(self,m,U):
        n_elems  = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1)
        nnz = int( 3*(nbf_elem)**2*n_elems )
        nbf = self.Get_nbf()
        ndof = 3*nbf 
        
        
        r = ar.L2Projection(self.pp[0], self.pp[1], self.pp[2], 
                        self.xxsi[0], self.xxsi[1], self.xxsi[2],
                        m.pp[0], m.pp[1], m.pp[2], 
                        m.xxsi[0], m.xxsi[1], m.xxsi[2],
                        U, nnz, nnz, nnz, ndof )
        
        indexI     = r[0]
        indexJ     = r[1]
        nnz_values = r[2] 
        rhs        = r[3]
        M = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        print('Solving the L2 projection linear system')
        # M_LU  = splalg.splu(M) 
        # Up    = M_LU.solve(rhs) # L2 projection solution  
        M_LLt = cholesky(M)
        Up    = M_LLt.solve_A(rhs)
        print('L2 projection finished')
        return Up 
    
    def LinearProjection(self,m,U):
        """ Works only for linear functions """ 
        isLinear = all( p==1 for p in self.pp)
        
        if not isLinear:
            raise ValueError('I should be a linear mesh !')
        
        xminM  = m.xxsi[0][0] ; xmaxM = m.xxsi[0][-1]; 
        yminM  = m.xxsi[1][0] ; ymaxM = m.xxsi[1][-1]; 
        zminM  = m.xxsi[2][0] ; zmaxM = m.xxsi[2][-1];  

        xmin  = self.xxsi[0][0] ; xmax = self.xxsi[0][-1]; 
        ymin  = self.xxsi[1][0] ; ymax = self.xxsi[1][-1]; 
        zmin  = self.xxsi[2][0] ; zmax = self.xxsi[2][-1];
        
        errorDomain = (xmin<xminM) or (ymin<yminM) or (zmin<zminM) or (xmax>xmaxM) or (ymax>ymaxM) or (zmax>zmaxM) 
        if errorDomain:
            raise ValueError('My domain sould be included in the mesh domain that you want to project !')

        nbf = m.Get_nbf();  
        cpx1d = np.linspace(xmin, xmax, self.n_elems[0]+1)
        cpy1d = np.linspace(ymin, ymax, self.n_elems[1]+1)
        cpz1d = np.linspace(zmin, zmax, self.n_elems[2]+1)
        
        Ux = m.EvaluateBsplineScalarStructured(cpx1d, cpy1d, cpz1d, U[:nbf])
        Uy = m.EvaluateBsplineScalarStructured(cpx1d, cpy1d, cpz1d, U[nbf:2*nbf])
        Uz = m.EvaluateBsplineScalarStructured(cpx1d, cpy1d, cpz1d, U[2*nbf:])

        Us = np.r_[Ux,Uy,Uz]
        return Us 
    
    def ProjectDispOnMesh(self,m,U):
        im, = np.where(m.conn[:,0]!=-1)
        x = m.n[im,0]
        y = m.n[im,1]
        z = m.n[im,2]
        
        xmin  = np.min(x) ; xmax = np.max(x); 
        ymin  = np.min(y) ; ymax = np.max(y); 
        zmin  = np.min(z) ; zmax = np.max(z);  

        xminM  = self.xxsi[0][0] ; xmaxM = self.xxsi[0][-1]; 
        yminM  = self.xxsi[1][0] ; ymaxM = self.xxsi[1][-1]; 
        zminM  = self.xxsi[2][0] ; zmaxM = self.xxsi[2][-1];
        
        errorDomain = (xmin<xminM) or (ymin<yminM) or (zmin<zminM) or (xmax>xmaxM) or (ymax>ymaxM) or (zmax>zmaxM) 
        if errorDomain:
            raise ValueError('My domain sould be included in the mesh domain that you want to project !')
        
        ndof = int(len(U)/3) 
        Ux = self.EvaluateBsplineScalar(x,y,z,U[:ndof])
        Uy = self.EvaluateBsplineScalar(x,y,z,U[ndof:2*ndof])
        Uz = self.EvaluateBsplineScalar(x,y,z,U[2*ndof:])
        
        Up = np.zeros(m.ndof)
        Up[m.conn[im,0]] = Ux
        Up[m.conn[im,1]] = Uy 
        Up[m.conn[im,2]] = Uz 
        
        return Up 
 


 

       
    def GetFcmLevelSetCellTesselationCoord(self,Phi,lvlmax):
        self.integrationCellsCoord = [] 
        #self.pix = []
        #self.piy = []
        #self.piz = [] 
        #self.wg  = []  
        e_xi      =  self.xxsi[0][self.pp[0]:-self.pp[0]]  
        e_eta     =  self.xxsi[1][self.pp[1]:-self.pp[1]]  
        e_zeta    =  self.xxsi[2][self.pp[2]:-self.pp[2]]
        nx = self.n_elems[0]
        ny = self.n_elems[1]
        nz = self.n_elems[2]
        
        for ne in range(nx*ny*nz):
            #print('---Element ', k)
            k = int(np.floor(ne/(nx*ny)))
            j = int(np.floor((ne-k*nx*ny)/nx))
            i = ne - j*nx -k*nx*ny 
            xmin   =  e_xi[i]
            xmax   =  e_xi[i+1]
            ymin   =  e_eta[j]
            ymax   =  e_eta[j+1]
            zmin   =  e_zeta[k]
            zmax   =  e_zeta[k+1]
            C = Cell3d(xmin,xmax,ymin,ymax,zmin,zmax)
            C.DecomposeLevelSetCoord(self, Phi, lvlmax )
            
    def FcmLevelSetIntegrationTesselation(self,Phi,lvlmax):
        self.nbg_xi   = self.pp[0]+1
        self.nbg_eta  = self.pp[1]+1
        self.nbg_zeta = self.pp[2]+1 
        
        self.oneXi   = np.ones(self.nbg_xi)  
        self.oneEta  = np.ones(self.nbg_eta)
        self.oneZeta = np.ones(self.nbg_zeta)
        
        self.Gauss_xi   =  nb.GaussLegendre(self.nbg_xi)
        self.Gauss_eta  =  nb.GaussLegendre(self.nbg_eta)
        self.Gauss_zeta =  nb.GaussLegendre(self.nbg_zeta)
        
 
        self.wgRef = np.kron( self.Gauss_zeta[1],np.kron(self.Gauss_eta[1], self.Gauss_xi[1]) )

        self.GaussTetra = nb.GaussTetrahedron((max(self.pp[0],self.pp[1],self.pp[2])))

        self.pix = []
        self.piy = []
        self.piz = [] 
        self.wg  = []  
        e_xi      =  self.xxsi[0][self.pp[0]:-self.pp[0]]  
        e_eta     =  self.xxsi[1][self.pp[1]:-self.pp[1]]  
        e_zeta    =  self.xxsi[2][self.pp[2]:-self.pp[2]]
        nx = self.n_elems[0]
        ny = self.n_elems[1]
        nz = self.n_elems[2]
        
        for ne in range(nx*ny*nz):
            k = int(np.floor(ne/(nx*ny)))
            j = int(np.floor((ne-k*nx*ny)/nx))
            i = ne - j*nx -k*nx*ny 
            xmin    =  e_xi[i]
            xmax    =  e_xi[i+1]
            ymin   =  e_eta[j]
            ymax   =  e_eta[j+1]
            zmin   =  e_zeta[k]
            zmax   =  e_zeta[k+1]
            C = Cell3d(xmin,xmax,ymin,ymax,zmin,zmax)
            C.DecomposeLevelSetIntegration(self, Phi, lvlmax )
            
        self.pix = np.array(self.pix)
        self.piy = np.array(self.piy)
        self.piz = np.array(self.piz)
        self.wg  = sps.diags ( np.array(self.wg) ) 
        self.npg = len(self.pix)
    
    
    def FcmLevelSetIntegrationTesselationImage(self, im, thrsh, lvlmax):
        r = ar.FcmIntegrationTrilinearInterp(im.pix, im.knotXi1, im.knotEta1, im.knotZeta1, thrsh, 
                                         self.pp[0], self.pp[1], self.pp[2], self.xxsi[0], self.xxsi[1], self.xxsi[2], lvlmax) 
        rv = np.array(r)
        rv = rv.reshape((-1,4)) # Each point is defined by (x,y,z,w) 
        self.pix = rv[:,0]
        self.piy = rv[:,1]
        self.piz = rv[:,2]
        self.wg  = sps.diags ( rv[:,3] )
        self.npg = self.pix.shape[0]                                        

     

    
    def GetGlobalBasisFunctionsMatrixAtStructuredPoints(self,x,y,z):
        """ 
        Returns the B-spline shape function matrix 
        Each row corresponds to a point of the structured grid (tensor product z-y-x)
        """
  
        nbf = self.Get_nbf()  
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1) 
        # n_elems = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        # nnz = int(n_elems*(3*nbf_elem)**3)
        npoints = len(x)*len(y)*len(z)
        nnz = int(nbf_elem*npoints)
        r = ar.GetBsplineFunctionsMatrixStructured(x,y,z,self.xxsi[0],self.xxsi[1],self.xxsi[2],
                                                    self.pp[0],self.pp[1],self.pp[2],
                                                    nnz,nnz,nnz)
        phiValues = r[0]
        indexI    = r[1]
        indexJ    = r[2] 
        phi    = sps.csc_matrix(( phiValues,    (indexI,indexJ)), shape = (npoints,nbf))
        return phi 
    def GetGlobalBasisFunctionsAndDerivMatrixAtStructuredPoints(self,x,y,z):
        """ 
        Returns the B-spline shape function matrices and derivatives  
        Each row corresponds to a point of the structured grid (tensor product z-y-x)
        """
        
        nbf = self.Get_nbf() 
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1) 
        # n_elems = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        # nnz = int(n_elems*(3*nbf_elem)**3)
        npoints = len(x)*len(y)*len(z)
        nnz = int(nbf_elem*npoints)
        r = ar.GetBsplineFunctionsAndDerivativesMatrixStructured(x,y,z,self.xxsi[0],self.xxsi[1],self.xxsi[2],
                                                    self.pp[0],self.pp[1],self.pp[2],
                                                    nnz,nnz,nnz,nnz,nnz,nnz)
        phiValues    = r[0]
        dphidxValues = r[1]
        dphidyValues = r[2]
        dphidzValues = r[3]
        indexI    = r[4]
        indexJ    = r[5] 
        phi       = sps.csc_matrix(( phiValues,    (indexI,indexJ)), shape = (npoints,nbf))
        dphidx    = sps.csc_matrix(( dphidxValues,    (indexI,indexJ)), shape = (npoints,nbf))
        dphidy    = sps.csc_matrix(( dphidyValues,    (indexI,indexJ)), shape = (npoints,nbf))
        dphidz    = sps.csc_matrix(( dphidzValues,    (indexI,indexJ)), shape = (npoints,nbf))
        return phi,dphidx,dphidy,dphidz  
    
    
    def EvaluateBsplineScalar(self,x,y,z,U1):
        u = tr.EvaluateBspline3D(x,y,z, self.xxsi[0],self.xxsi[1], self.xxsi[2],
                                  self.pp[0], self.pp[1], self.pp[2], U1, len(x))
        return u 
    def EvaluateBsplineAndDerivativesScalar(self,x,y,z,U1):
        r = tr.EvaluateBsplineAndGradient3D(x,y,z, self.xxsi[0],self.xxsi[1], self.xxsi[2],
                                  self.pp[0], self.pp[1], self.pp[2], U1, len(x),len(x),len(x),len(x))
        return r[0],r[1],r[2],r[3] 
    
    def EvaluateBsplineScalarStructured(self,x,y,z,U1):
        noutput = len(x)*len(y)*len(z)
        u = tr.EvaluateBsplineStructured3D(x,y,z,self.xxsi[0], self.xxsi[1], self.xxsi[2],
                                           self.pp[0], self.pp[1], self.pp[2], U1, noutput)
        
        return u 
    def EvaluateBsplineAndDerivativesScalarStructured(self,x,y,z,U1):
        noutput = len(x)*len(y)*len(z)
        r = tr.EvaluateBsplineAndGradientStructured3D(x,y,z,self.xxsi[0], self.xxsi[1], self.xxsi[2],
                                           self.pp[0], self.pp[1], self.pp[2], U1, noutput, noutput, noutput, noutput)
        return r[0],r[1],r[2],r[3] 
    
          
 
    def SetBasisFunctionsAtIntegrationPoints(self):
        phi, dphidx, dphidy, dphidz = nb.Get3dBasisFunctionsAtPts(self.pix, self.piy, self.piz,
                                                                  self.xxsi[0], self.xxsi[1], self.xxsi[2],
                                                                  self.pp[0], self.pp[1], self.pp[2])
        nbf  =  self.Get_nbf() 
        zero    = sps.csr_matrix((self.npg,nbf))
        self.phiMatrix = phi      

        self.phix    = sps.hstack((phi,zero,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi,zero)   ,  'csc')
        self.phiz    = sps.hstack((zero,zero,phi)  ,  'csc')
        
        self.dphixdx = sps.hstack((dphidx,zero,zero),  'csc')
        self.dphixdy = sps.hstack((dphidy,zero,zero),  'csc')
        self.dphixdz = sps.hstack((dphidz,zero,zero),  'csc')
        
        self.dphiydx = sps.hstack((zero,dphidx,zero),  'csc')
        self.dphiydy = sps.hstack((zero,dphidy,zero),  'csc')
        self.dphiydz = sps.hstack((zero,dphidz,zero),  'csc')
        
        self.dphizdx = sps.hstack((zero,zero,dphidx),  'csc')
        self.dphizdy = sps.hstack((zero,zero,dphidy),  'csc')
        self.dphizdz = sps.hstack((zero,zero,dphidz),  'csc')
    
    def SetUnivariateBasisFunctionsAndDerivatives(self,nbippe):

        xi   = np.linspace(self.xxsi[0][0], self.xxsi[0][-1], 2*self.n_elems[0]*nbippe[0]+1)[1::2] 
        eta  = np.linspace(self.xxsi[1][0], self.xxsi[1][-1], 2*self.n_elems[1]*nbippe[1]+1)[1::2] 
        zeta = np.linspace(self.xxsi[2][0], self.xxsi[2][-1], 2*self.n_elems[2]*nbippe[2]+1)[1::2] 
        
        self.pix1d = xi 
        self.piy1d = eta 
        self.piz1d = zeta 
        
        nbf_xi   = len(self.xxsi[0]) -1 - self.pp[0] ; 
        nbf_eta  = len(self.xxsi[1]) -1 - self.pp[1] ;   
        nbf_zeta = len(self.xxsi[2]) -1 - self.pp[2] ;
        
        self.Nxi         = np.zeros((len(xi),self.pp[0]+1))
        self.dNxidxi     = np.zeros((len(xi),self.pp[0]+1)) 
        self.Neta        = np.zeros((len(eta),self.pp[1]+1))
        self.dNetadeta   = np.zeros((len(eta),self.pp[1]+1))
        self.Nzeta       = np.zeros((len(zeta),self.pp[2]+1)) 
        self.dNzetadzeta = np.zeros((len(zeta),self.pp[2]+1)) 
        
        for i in range(len(xi)):
            span_xi  = nb.findspan(nbf_xi,  self.pp[0], xi[i] , self.xxsi[0])
            self.Nxi[i,:], self.dNxidxi[i,:] = nb.derbasisfuns(span_xi,self.pp[0],self.xxsi[0],1,xi[i])
        for i in range(len(eta)):
            span_eta = nb.findspan(nbf_eta, self.pp[1], eta[i], self.xxsi[1])
            self.Neta[i,:], self.dNetadeta[i,:] = nb.derbasisfuns(span_eta,self.pp[1],self.xxsi[1],1,eta[i])
        for i in range(len(zeta)):
            span_zeta = nb.findspan(nbf_zeta, self.pp[2], zeta[i], self.xxsi[2])
            self.Nzeta[i,:], self.dNzetadzeta[i,:] = nb.derbasisfuns(span_zeta,self.pp[2],self.xxsi[2],1,zeta[i])            
        
        nip_xi   =  self.n_elems[0]*nbippe[0] 
        nip_eta  =  self.n_elems[1]*nbippe[1]  
        nip_zeta =  self.n_elems[2]*nbippe[2]    

        self.pix = np.kron(np.ones(nip_zeta*nip_eta),xi)
        self.piy = np.kron(np.ones(nip_zeta),np.kron(eta, np.ones(nip_xi)))
        self.piz = np.kron(zeta,np.ones(nip_eta*nip_xi))
        self.npg = self.pix.shape[0] 
        
    

        
   
    def VoxelIntegration(self,nbippe):
        """
        Integration rule: Riemann sum (rectangle rule)
        Performs uniform sub-division of the parametric space 
        nbipe : number of integration points per element in each direction 
        """
        
        nip_xi   =  self.n_elems[0]*nbippe[0] 
        nip_eta  =  self.n_elems[1]*nbippe[1]  
        nip_zeta =  self.n_elems[2]*nbippe[2]    

        xi   = np.linspace(self.xxsi[0][0], self.xxsi[0][-1], 2*self.n_elems[0]*nbippe[0]+1)[1::2] 
        eta  = np.linspace(self.xxsi[1][0], self.xxsi[1][-1], 2*self.n_elems[1]*nbippe[1]+1)[1::2] 
        zeta = np.linspace(self.xxsi[2][0], self.xxsi[2][-1], 2*self.n_elems[2]*nbippe[2]+1)[1::2] 

        mes_xi   =  (self.xxsi[0][self.pp[0]+1]-self.xxsi[0][self.pp[0]])/nbippe[0] 
        mes_eta  =  (self.xxsi[1][self.pp[1]+1]-self.xxsi[1][self.pp[1]])/nbippe[1]  
        mes_zeta =  (self.xxsi[2][self.pp[2]+1]-self.xxsi[2][self.pp[2]])/nbippe[2]  
        
        # Constant measure of an integration cuboid 
        self.wg = mes_xi*mes_eta*mes_zeta # Here wg is a scalar 
        
        """ Integration points and weights """
        self.pix = np.kron(np.ones(nip_zeta*nip_eta),xi)
        self.piy = np.kron(np.ones(nip_zeta),np.kron(eta, np.ones(nip_xi)))
        self.piz = np.kron(zeta,np.ones(nip_eta*nip_xi))
        self.npg = self.pix.shape[0]
        
        self.pix1d = xi 
        self.piy1d = eta 
        self.piz1d = zeta 
        
        
        phi = self.GetGlobalBasisFunctionsMatrixAtStructuredPoints(xi,eta,zeta)
        
        

        nbf  =  self.Get_nbf()
        zero    = sps.csr_matrix((self.npg,nbf))
        
        self.phiMatrix = phi      
 
        self.phix    = sps.hstack((phi,zero,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi,zero)   ,  'csc')
        self.phiz    = sps.hstack((zero,zero,phi)   ,  'csc')   
 

    
    def Stiffness(self,E,nu):
        """ 
        Stiffness Matrix """
        hooke = E/((1+nu)*(1-2*nu))*np.array([[1-nu, nu, nu,     0          , 0 ,    0],
                                              [nu, 1-nu, nu,     0          , 0 ,    0],
                                              [nu, nu, 1-nu,     0          , 0 ,    0],
                                              [0,  0,   0,   0.5*(1-2*nu)   , 0 ,    0],
                                              [0,  0,   0,       0  , 0.5*(1-2*nu),  0],
                                              [0,  0,   0,       0,           0,  0.5*(1-2*nu)]]) 
        
        Bxy = self.dphixdy + self.dphiydx
        Bxz = self.dphixdz + self.dphizdx
        Byz = self.dphiydz + self.dphizdy 
        
        K =  hooke[0,0]*self.dphixdx.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[0,1]*self.dphiydy.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[0,2]*self.dphizdz.T.dot(self.wg.dot(self.dphixdx)) + \
             hooke[1,0]*self.dphixdx.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[1,1]*self.dphiydy.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[1,2]*self.dphizdz.T.dot(self.wg.dot(self.dphiydy)) + \
             hooke[2,0]*self.dphixdx.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[2,1]*self.dphiydy.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[2,2]*self.dphizdz.T.dot(self.wg.dot(self.dphizdz)) + \
             hooke[3,3]*Bxy.T.dot(self.wg.dot(Bxy)) + \
             hooke[4,4]*Bxz.T.dot(self.wg.dot(Bxz)) + \
             hooke[5,5]*Byz.T.dot(self.wg.dot(Byz))  
        
        return K
    
    def Laplacian(self):
        """ Homogeneous Laplacian definied on all structured mesh domain """ 
        n_elems  = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1)
        nnz = int( 3*(nbf_elem)**2*n_elems )
        nbf = self.Get_nbf()
        ndof = 3*nbf 
        
        r = ar.Laplacian_Structured(self.pp[0], self.pp[1], self.pp[2], self.xxsi[0], self.xxsi[1], self.xxsi[2],nnz,nnz,nnz) 
 
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        L = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return L 
    
    def HomogeneousStiffnessEA(self,E,nu):
        """ Homogeneous stiffness matrix """
        n_elems  = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1)
        nnz = int( 9*(nbf_elem)**2*n_elems )
        nbf = self.Get_nbf()
        ndof = 3*nbf  
        
        r = ar.HomogeneousStiffness(E,nu,self.pp[0], self.pp[1], self.pp[2], self.xxsi[0], self.xxsi[1], self.xxsi[2],nnz,nnz,nnz)
        
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2] 
        K = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        
        return K
    
    def FcmStiffnessAssemblyAndBfCriteria(self, im, thrsh, lvlmax, E, nu, Parallel=True):
        """ Stiffness matrix using the FCM scheme
        Integration only on the domain defined by the level-set
        This matrix is in general singular
        Returns also the integral of each B-spline basis function over
        the level-set domain domain
        """
        n_elems  = self.n_elems[0]*self.n_elems[1]*self.n_elems[2]
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)*(self.pp[2]+1)
        nnz = int( 9*(nbf_elem)**2*n_elems )
        nbf = self.Get_nbf()
        ndof = 3*nbf
        
        if Parallel==True:
            r = ar.FcmTrilinearInterpStiffnessAndBfIntegralParallel(im.pix, im.knotXi1, im.knotEta1, im.knotZeta1,
                                                        thrsh, self.pp[0], self.pp[1], self.pp[2],
                                                        self.xxsi[0], self.xxsi[1], self.xxsi[2],
                                                        lvlmax,E,nu, nnz,nnz,nnz, nbf)    
        else:
            
            r = ar.FcmTrilinearInterpStiffnessAndBfIntegral(im.pix, im.knotXi1, im.knotEta1, im.knotZeta1,
                                                        thrsh, self.pp[0], self.pp[1], self.pp[2],
                                                        self.xxsi[0], self.xxsi[1], self.xxsi[2],
                                                        lvlmax,E,nu, nnz,nnz,nnz, nbf)
        
        indexI = r[0]
        indexJ = r[1]
        nnz_values = r[2]      
        intBf = r[3]
        
        K = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        
        return K,intBf
    
    
    def HomogeneousGaussLaplacian(self):
        """
        Gauss Integration 
        """         
        nbg_xi    = self.pp[0]+1 
        nbg_eta   = self.pp[1]+1  
        nbg_zeta  = self.pp[2]+1 
        
        Gauss_xi   =  nb.GaussLegendre(nbg_xi)
        Gauss_eta  =  nb.GaussLegendre(nbg_eta)
        Gauss_zeta =  nb.GaussLegendre(nbg_zeta)
        
        nbf = self.Get_nbf() 
        
        e_xi   = np.unique(self.xxsi[0]) ; ne_xi   = e_xi.shape[0]-1
        e_eta  = np.unique(self.xxsi[1]) ; ne_eta  = e_eta.shape[0]-1
        e_zeta = np.unique(self.xxsi[2]) ; ne_zeta = e_zeta.shape[0]-1
        
        xi_min    =  np.kron(e_xi[:-1],np.ones(nbg_xi))    
        eta_min   =  np.kron(e_eta[:-1],np.ones(nbg_eta))  
        zeta_min  =  np.kron(e_zeta[:-1],np.ones(nbg_zeta))  
        
        
        xi_g      =  np.kron(np.ones(ne_xi)  , Gauss_xi[0])
        eta_g     =  np.kron(np.ones(ne_eta) , Gauss_eta[0])
        zeta_g    =  np.kron(np.ones(ne_zeta) , Gauss_zeta[0])
        
        """ Measures of elements """
        mes_xi   = e_xi[1:]  - e_xi[:-1]
        mes_eta  = e_eta[1:] - e_eta[:-1]
        mes_zeta = e_zeta[1:] - e_zeta[:-1]
        
        mes_xi   = np.kron(mes_xi,np.ones(nbg_xi))   
        mes_eta  = np.kron(mes_eta,np.ones(nbg_eta)) 
        mes_zeta = np.kron(mes_zeta,np.ones(nbg_zeta)) 
        
         
        """ Going from the reference element to the parametric space  """ 
        xi         = xi_min   + 0.5*(xi_g+1)*mes_xi       # Aranged gauss points in  xi direction  
        eta        = eta_min  + 0.5*(eta_g+1)*mes_eta     # Aranged gauss points in  eta direction 
        zeta       = zeta_min + 0.5*(zeta_g+1)*mes_zeta   # Aranged gauss points in  zeta direction 
        
        phi, dphidx, dphidy, dphidz = self.GetGlobalBasisFunctionsAndDerivMatrixAtStructuredPoints(xi,eta,zeta)
        
        npg        = dphidx.shape[0]
        
        wg_xi        =  np.kron(np.ones(ne_xi) , Gauss_xi[1])
        wg_eta       =  np.kron(np.ones(ne_eta), Gauss_eta[1])
        wg_zeta      =  np.kron(np.ones(ne_zeta), Gauss_zeta[1])
        
        oneXi   = np.ones(xi.shape[0])
        oneEta  = np.ones(eta.shape[0])
        oneZeta = np.ones(zeta.shape[0])
        
        
        
        mes_xi   = np.kron(oneZeta, np.kron(oneEta,mes_xi)) 
        mes_eta  = np.kron(oneZeta, np.kron(mes_eta,oneXi)) 
        mes_zeta = np.kron(mes_zeta, np.kron(oneEta,oneXi))  
          
        wg     =  np.kron(wg_zeta,np.kron(wg_eta, wg_xi))*mes_xi*mes_eta*mes_zeta/8
        

        wg = sps.diags(wg)        
        zero    = sps.csr_matrix((npg,nbf))

  
        dphixdx = sps.hstack((dphidx,zero,zero),  'csc')
        dphixdy = sps.hstack((dphidy,zero,zero),  'csc')
        dphixdz = sps.hstack((dphidz,zero,zero),  'csc')
        
        dphiydx = sps.hstack((zero,dphidx,zero),  'csc')
        dphiydy = sps.hstack((zero,dphidy,zero),  'csc')
        dphiydz = sps.hstack((zero,dphidz,zero),  'csc')
        
        dphizdx = sps.hstack((zero,zero,dphidx),  'csc')
        dphizdy = sps.hstack((zero,zero,dphidy),  'csc')
        dphizdz = sps.hstack((zero,zero,dphidz),  'csc')
        

        """ 
        Assembly 
        """
        
        L = dphixdx.T.dot(wg.dot(dphixdx)) + \
            dphixdy.T.dot(wg.dot(dphixdy)) + \
            dphixdz.T.dot(wg.dot(dphixdz)) + \
            dphiydx.T.dot(wg.dot(dphiydx)) + \
            dphiydy.T.dot(wg.dot(dphiydy)) + \
            dphiydz.T.dot(wg.dot(dphiydz)) + \
            dphizdx.T.dot(wg.dot(dphizdx)) + \
            dphizdy.T.dot(wg.dot(dphizdy)) + \
            dphizdz.T.dot(wg.dot(dphizdz))
            
        return L    
 
    def HomogeneousGaussStiffness(self,E,nu):
        
        """
        Gauss Integration 
        """ 
        nbg_xi    = self.pp[0]+1 
        nbg_eta   = self.pp[1]+1  
        nbg_zeta  = self.pp[2]+1 

        Gauss_xi   =  nb.GaussLegendre(nbg_xi)
        Gauss_eta  =  nb.GaussLegendre(nbg_eta)
        Gauss_zeta =  nb.GaussLegendre(nbg_zeta)
        
 
        nbf = self.Get_nbf() 
 

        e_xi   = np.unique(self.xxsi[0]) ; ne_xi   = e_xi.shape[0]-1
        e_eta  = np.unique(self.xxsi[1]) ; ne_eta  = e_eta.shape[0]-1
        e_zeta = np.unique(self.xxsi[2]) ; ne_zeta = e_zeta.shape[0]-1
        
        xi_min    =  np.kron(e_xi[:-1],np.ones(nbg_xi))    
        eta_min   =  np.kron(e_eta[:-1],np.ones(nbg_eta))  
        zeta_min  =  np.kron(e_zeta[:-1],np.ones(nbg_zeta))  
        
        
        xi_g      =  np.kron(np.ones(ne_xi)  , Gauss_xi[0])
        eta_g     =  np.kron(np.ones(ne_eta) , Gauss_eta[0])
        zeta_g    =  np.kron(np.ones(ne_zeta) , Gauss_zeta[0])
        
        """ Measures of elements """
        mes_xi   = e_xi[1:]  - e_xi[:-1]
        mes_eta  = e_eta[1:] - e_eta[:-1]
        mes_zeta = e_zeta[1:] - e_zeta[:-1]
        
        mes_xi   = np.kron(mes_xi,np.ones(nbg_xi))   
        mes_eta  = np.kron(mes_eta,np.ones(nbg_eta)) 
        mes_zeta = np.kron(mes_zeta,np.ones(nbg_zeta)) 
        
         
        """ Going from the reference element to the parametric space  """ 
        xi         = xi_min   + 0.5*(xi_g+1)*mes_xi       # Aranged gauss points in  xi direction  
        eta        = eta_min  + 0.5*(eta_g+1)*mes_eta     # Aranged gauss points in  eta direction 
        zeta       = zeta_min + 0.5*(zeta_g+1)*mes_zeta   # Aranged gauss points in  zeta direction 
        
        # phi_xi   , dphi_xi    = nb.global_basisfuns(self.pp[0],self.xxsi[0],xi)
        # phi_eta  , dphi_eta   = nb.global_basisfuns(self.pp[1],self.xxsi[1],eta)
        # phi_zeta , dphi_zeta  = nb.global_basisfuns(self.pp[2],self.xxsi[2],zeta)
        
        
        # phiEta_phiXi    =  sps.kron(phi_eta, phi_xi,  'csc') 
        # phiEta_dphidXi  =  sps.kron(phi_eta,dphi_xi,  'csc') 
        # dphidx         =  sps.kron(phi_zeta, phiEta_dphidXi,  'csc') 
        # dphiEta_phiXi   =  sps.kron(dphi_eta, phi_xi,  'csc') 
        # dphidy        =  sps.kron(phi_zeta, dphiEta_phiXi,  'csc') 
        # dphidz       =  sps.kron(dphi_zeta, phiEta_phiXi,  'csc') 
        
        phi, dphidx, dphidy, dphidz = self.GetGlobalBasisFunctionsAndDerivMatrixAtStructuredPoints(xi,eta,zeta)

        
        npg        = dphidx.shape[0]
        
        wg_xi        =  np.kron(np.ones(ne_xi) , Gauss_xi[1])
        wg_eta       =  np.kron(np.ones(ne_eta), Gauss_eta[1])
        wg_zeta      =  np.kron(np.ones(ne_zeta), Gauss_zeta[1])
        
        oneXi   = np.ones(xi.shape[0])
        oneEta  = np.ones(eta.shape[0])
        oneZeta = np.ones(zeta.shape[0])
        
        
        
        mes_xi   = np.kron(oneZeta, np.kron(oneEta,mes_xi)) 
        mes_eta  = np.kron(oneZeta, np.kron(mes_eta,oneXi)) 
        mes_zeta = np.kron(mes_zeta, np.kron(oneEta,oneXi))  
          
        wg     =  np.kron(wg_zeta,np.kron(wg_eta, wg_xi))*mes_xi*mes_eta*mes_zeta/8
        

        wg = sps.diags(wg)        
        zero    = sps.csr_matrix((npg,nbf))

  
        dphixdx = sps.hstack((dphidx,zero,zero),  'csc')
        dphixdy = sps.hstack((dphidy,zero,zero),  'csc')
        dphixdz = sps.hstack((dphidz,zero,zero),  'csc')
        
        dphiydx = sps.hstack((zero,dphidx,zero),  'csc')
        dphiydy = sps.hstack((zero,dphidy,zero),  'csc')
        dphiydz = sps.hstack((zero,dphidz,zero),  'csc')
        
        dphizdx = sps.hstack((zero,zero,dphidx),  'csc')
        dphizdy = sps.hstack((zero,zero,dphidy),  'csc')
        dphizdz = sps.hstack((zero,zero,dphidz),  'csc')
        

        """ 
        Assembly 
        """
        hooke = E/((1+nu)*(1-2*nu))*np.array([[1-nu, nu, nu,     0          , 0 ,    0],
                                              [nu, 1-nu, nu,     0          , 0 ,    0],
                                              [nu, nu, 1-nu,     0          , 0 ,    0],
                                              [0,  0,   0,   0.5*(1-2*nu)   , 0 ,    0],
                                              [0,  0,   0,       0  , 0.5*(1-2*nu),  0],
                                              [0,  0,   0,       0,           0,  0.5*(1-2*nu)]]) 
                
        Bxy = dphixdy + dphiydx
        Bxz = dphixdz + dphizdx
        Byz = dphiydz + dphizdy 
        
        K =  hooke[0,0]*dphixdx.T.dot(wg.dot(dphixdx)) + \
             hooke[0,1]*dphiydy.T.dot(wg.dot(dphixdx)) + \
             hooke[0,2]*dphizdz.T.dot(wg.dot(dphixdx)) + \
             hooke[1,0]*dphixdx.T.dot(wg.dot(dphiydy)) + \
             hooke[1,1]*dphiydy.T.dot(wg.dot(dphiydy)) + \
             hooke[1,2]*dphizdz.T.dot(wg.dot(dphiydy)) + \
             hooke[2,0]*dphixdx.T.dot(wg.dot(dphizdz)) + \
             hooke[2,1]*dphiydy.T.dot(wg.dot(dphizdz)) + \
             hooke[2,2]*dphizdz.T.dot(wg.dot(dphizdz)) + \
             hooke[3,3]*Bxy.T.dot(wg.dot(Bxy)) + \
             hooke[4,4]*Bxz.T.dot(wg.dot(Bxz)) + \
             hooke[5,5]*Byz.T.dot(wg.dot(Byz)) 
            
        return K    
    

        
         
            

    def VtkEmbeddingElements(self,path,U=None): 
        x       =  self.xxsi[0][self.pp[0]:-self.pp[0]]  
        y       =  self.xxsi[1][self.pp[1]:-self.pp[1]]  
        z       =  self.xxsi[2][self.pp[2]:-self.pp[2]]
        
        xd1 = np.kron( np.ones(len(z)),  np.kron( np.ones(len(y)), np.array([x[0], x[-1]])) ) 
        yd1 = np.kron( np.ones(len(z)), np.kron(y,np.ones(2)))
        zd1 = np.kron( z, np.ones(2*len(y)) )  

        yd2 = np.kron( np.ones(len(z)),  np.kron( np.ones(len(x)), np.array([y[0], y[-1]])) ) 
        xd2 = np.kron( np.ones(len(z)), np.kron(x,np.ones(2)))
        zd2 = np.kron( z, np.ones(2*len(x)) )  
        
                
        zd3 = np.kron( np.ones(len(x)),  np.kron( np.ones(len(y)), np.array([z[0], z[-1]])) ) 
        yd3 = np.kron( np.ones(len(x)), np.kron(y,np.ones(2)))
        xd3 = np.kron( x, np.ones(2*len(y)) ) 
        
        if U is not None :
            nbf = self.Get_nbf() 
            Ux = U[:nbf]
            Uy = U[nbf:2*nbf]
            Uz = U[2*nbf:]
            
            uxd1 = self.EvaluateBsplineScalar(xd1,yd1,zd1,Ux)
            uyd1 = self.EvaluateBsplineScalar(xd1,yd1,zd1,Uy)
            uzd1 = self.EvaluateBsplineScalar(xd1,yd1,zd1,Uz)
            
            uxd2 = self.EvaluateBsplineScalar(xd2,yd2,zd2,Ux)
            uyd2 = self.EvaluateBsplineScalar(xd2,yd2,zd2,Uy)
            uzd2 = self.EvaluateBsplineScalar(xd2,yd2,zd2,Uz)
            
            uxd3 = self.EvaluateBsplineScalar(xd3,yd3,zd3,Ux)
            uyd3 = self.EvaluateBsplineScalar(xd3,yd3,zd3,Uy)
            uzd3 = self.EvaluateBsplineScalar(xd3,yd3,zd3,Uz)
            
          
            x  = np.r_[xd1,xd2,xd3]
            y  = np.r_[yd1,yd2,yd3]
            z  = np.r_[zd1,zd2,zd3]
            
            ux = np.r_[uxd1,uxd2,uxd3]
            uy = np.r_[uyd1,uyd2,uyd3]
            uz = np.r_[uzd1,uzd2,uzd3] 
 
            
            npoints = len(x)
            ncells = int(len(x) / 2.0)
 
            # create some temporary arrays to write grid topology
            offsets = np.arange(start = 2, step = 2, stop = npoints + 1, dtype = 'int32')   # index of last node in each cell
            connectivity = np.arange(npoints, dtype = 'int32')                              # each point is only connected to itself
            cell_types = np.empty(npoints, dtype = 'uint8') 
           
            cell_types[:] = VtkLine.tid
        
            w = VtkFile(path, VtkUnstructuredGrid)

            w.openGrid()
            w.openPiece(ncells = ncells, npoints = npoints)
            
            w.openElement("Points")
            w.addData("points", (x,y,z))
            w.closeElement("Points")
            w.openElement("Cells")
            w.addData("connectivity", connectivity)
            w.addData("offsets", offsets)
            w.addData("types", cell_types)
            w.closeElement("Cells")
            
            w.openData("Point")
            w.addData("disp",(ux,uy,uz))
            w.closeData("Point")
        
            w.closePiece()
            w.closeGrid()
            w.appendData( (x,y,z) )
            w.appendData(connectivity).appendData(offsets).appendData(cell_types)
            w.appendData( (ux,uy,uz))
 
            w.save()
            print(w.getFileName()+' saved') 
        else: 
            linesToVTK(path, np.r_[xd1,xd2,xd3], np.r_[yd1,yd2,yd3], np.r_[zd1,zd2,zd3])  
        
        

            
 
    def VtkIntegrationCells(self,path):
        """ Getting the lines of each integration cell"""
        lx = np.array([])
        ly = np.array([])
        lz = np.array([])
        for c in self.integrationCellsCoord : 
            if c[0]=='t':
                for i in range(c[2].shape[0]):
                    # Loop over simplices 
                    pts = c[1][c[2][i,:]]   
                    lx = np.r_[lx, pts[[0,1],0], pts[[0,2],0], pts[[0,3],0], pts[[2,3],0], pts[[1,3],0], pts[[1,2],0] ]
                    ly = np.r_[ly, pts[[0,1],1], pts[[0,2],1], pts[[0,3],1], pts[[2,3],1], pts[[1,3],1], pts[[1,2],1] ]
                    lz = np.r_[lz, pts[[0,1],2], pts[[0,2],2], pts[[0,3],2], pts[[2,3],2], pts[[1,3],2], pts[[1,2],2] ]
            if c[0]=='c':
                xmin = c[1] ; xmax = c[2] ; ymin = c[3] ; ymax = c[4]; zmin = c[5]; zmax = c[6]
                lx = np.r_[lx, np.array([xmin,xmax,
                                          xmax,xmax,
                                          xmax,xmin,
                                          xmin,xmin,
                                          xmin,xmax,
                                          xmax,xmax,
                                          xmax,xmin,
                                          xmin,xmin,
                                          xmax,xmax,
                                          xmax,xmax,
                                          xmin,xmin,
                                          xmin,xmin]) ]
                ly = np.r_[ly, np.array([ymin,ymin,
                                          ymin,ymax,
                                          ymax,ymax,
                                          ymin,ymax,
                                          ymin,ymin,
                                          ymin,ymax,
                                          ymax,ymax,
                                          ymin,ymax,
                                          ymin,ymin,
                                          ymax,ymax,
                                          ymin,ymin,
                                          ymax,ymax]) ]
                lz = np.r_[lz, np.array([zmin,zmin,
                                          zmin,zmin,
                                          zmin,zmin,
                                          zmin,zmin,
                                          zmax,zmax,
                                          zmax,zmax,
                                          zmax,zmax,
                                          zmax,zmax,
                                          zmin,zmax,
                                          zmin,zmax,
                                          zmin,zmax,
                                          zmin,zmax]) ]   
                
        """ Exporting to VTK using EVTK library """ 
        linesToVTK(path, lx,ly,lz   )  
 
    def get_U_Eps_Sigma_OnPointCloud_Seq(self, x,y,z, U, hooke):
        
        nbf = self.Get_nbf() 
        dim = 3 
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:] # displacement at the control points 
        
        """ Displacement """ 
        ux,duxdx, duxdy, duxdz =  self.EvaluateBsplineAndDerivativesScalar(x,y,z,Ux) 
        uy,duydx, duydy, duydz =  self.EvaluateBsplineAndDerivativesScalar(x,y,z,Uy) 
        uz,duzdx, duzdy, duzdz =  self.EvaluateBsplineAndDerivativesScalar(x,y,z,Uz)
        
 
        """ Strain """ 
        exx = duxdx 
        eyy = duydy 
        ezz = duzdz 
        exy = 0.5*(duydx+duxdy)
        exz = 0.5*(duzdx+duxdz)
        eyz = 0.5*(duzdy+duzdz)   
        
        """ Stress """ 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(ezz) + hooke[0,3]*(2*exy) + hooke[0,4]*(2*exz) + hooke[0,5]*(2*eyz)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(ezz) + hooke[1,3]*(2*exy) + hooke[1,4]*(2*exz) + hooke[1,5]*(2*eyz)
        szz = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(ezz) + hooke[2,3]*(2*exy) + hooke[2,4]*(2*exz) + hooke[2,5]*(2*eyz)
        sxy = hooke[3,0]*(exx) + hooke[3,1]*(eyy) + hooke[3,2]*(ezz) + hooke[3,3]*(2*exy) + hooke[3,4]*(2*exz) + hooke[3,5]*(2*eyz)
        sxz = hooke[4,0]*(exx) + hooke[4,1]*(eyy) + hooke[4,2]*(ezz) + hooke[4,3]*(2*exy) + hooke[4,4]*(2*exz) + hooke[4,5]*(2*eyz)
        syz = hooke[5,0]*(exx) + hooke[5,1]*(eyy) + hooke[5,2]*(ezz) + hooke[5,3]*(2*exy) + hooke[5,4]*(2*exz) + hooke[5,5]*(2*eyz)
        
        return ux,uy,uz,exx,eyy,ezz,exy,exz,eyz,sxx,syy,szz,sxy,sxz,syz  

    def get_U_Eps_Sigma_OnPointCloud(self,x,y,z,U,hooke):
 
        nbf = self.Get_nbf() 
        ndof =3*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
            
        phi, dphidx, dphidy, dphidz = nb.Get3dBasisFunctionsAtPts(x,y,z, self.xxsi[0], self.xxsi[1], self.xxsi[2], self.pp[0], self.pp[1], self.pp[2])

        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:]  # displacement at the control points 
   
        """ Displacement """
        ux = phi.dot(Ux)
        uy = phi.dot(Uy)
        uz = phi.dot(Uz)
        """ Strain """ 
        exx = dphidx.dot(Ux)
        eyy = dphidy.dot(Uy) 
        ezz = dphidz.dot(Uz)
        exy = 0.5*(dphidx.dot(Uy)+dphidy.dot(Ux))
        exz = 0.5*(dphidx.dot(Uz)+dphidz.dot(Ux))
        eyz = 0.5*(dphidy.dot(Uz)+dphidz.dot(Uy))
        """ Stress """ 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(ezz) + hooke[0,3]*(2*exy) + hooke[0,4]*(2*exz) + hooke[0,5]*(2*eyz)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(ezz) + hooke[1,3]*(2*exy) + hooke[1,4]*(2*exz) + hooke[1,5]*(2*eyz)
        szz = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(ezz) + hooke[2,3]*(2*exy) + hooke[2,4]*(2*exz) + hooke[2,5]*(2*eyz)
        sxy = hooke[3,0]*(exx) + hooke[3,1]*(eyy) + hooke[3,2]*(ezz) + hooke[3,3]*(2*exy) + hooke[3,4]*(2*exz) + hooke[3,5]*(2*eyz)
        sxz = hooke[4,0]*(exx) + hooke[4,1]*(eyy) + hooke[4,2]*(ezz) + hooke[4,3]*(2*exy) + hooke[4,4]*(2*exz) + hooke[4,5]*(2*eyz)
        syz = hooke[5,0]*(exx) + hooke[5,1]*(eyy) + hooke[5,2]*(ezz) + hooke[5,3]*(2*exy) + hooke[5,4]*(2*exz) + hooke[5,5]*(2*eyz)
        
        return ux,uy,uz,exx,eyy,ezz,exy,exz,eyz,sxx,syy,szz,sxy,sxz,syz  
 

    def VtkPlot(self, path , U, neval):
        eps = 0
        xi   = np.linspace(self.xxsi[0][self.pp[0]]+eps  , self.xxsi[0][-self.pp[0]]-eps , neval[0])
        eta  = np.linspace(self.xxsi[1][self.pp[1]]+eps  , self.xxsi[1][-self.pp[1]]-eps , neval[1])
        zeta = np.linspace(self.xxsi[2][self.pp[2]]+eps  , self.xxsi[2][-self.pp[2]]-eps , neval[2])
 
        nzeta = len(zeta)
        neta  = len(eta)
        nxi   = len(xi)
 
        # phi, dphidx, dphidy, dphidz = self.GetGlobalBasisFunctionsAndDerivMatrixAtStructuredPoints(xi, eta, zeta)
        nbf = self.Get_nbf() 
        dim = 3 
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:] # displacement at the control points 
        P = self.Get_Btot() # control points 
        
        # """ Displacement """
        # ux = phi.dot(Ux)
        # uy = phi.dot(Uy)
        # uz = phi.dot(Uz)
        # """ Strain """ 
        # exx = dphidx.dot(Ux)
        # eyy = dphidy.dot(Uy) 
        # ezz = dphidz.dot(Uz)
        # exy = 0.5*(dphidx.dot(Uy)+dphidy.dot(Ux))
        # exz = 0.5*(dphidx.dot(Uz)+dphidz.dot(Ux))
        # eyz = 0.5*(dphidy.dot(Uz)+dphidz.dot(Uy))
        
        # """ Solid points """ 
        # x  = phi.dot(P[0,:]) 
        # y  = phi.dot(P[1,:])
        # z  = phi.dot(P[2,:])
        
        """ Displacement """ 
        ux,duxdx, duxdy, duxdz =  self.EvaluateBsplineAndDerivativesScalarStructured(xi,eta,zeta,Ux) 
        uy,duydx, duydy, duydz =  self.EvaluateBsplineAndDerivativesScalarStructured(xi,eta,zeta,Uy) 
        uz,duzdx, duzdy, duzdz =  self.EvaluateBsplineAndDerivativesScalarStructured(xi,eta,zeta,Uz)
        """ Strain """ 
        exx = duxdx 
        eyy = duydy 
        ezz = duzdz 
        exy = 0.5*(duydx+duxdy)
        exz = 0.5*(duzdx+duxdz)
        eyz = 0.5*(duzdy+duydz)        
        """ Solid points """
        x = self.EvaluateBsplineScalarStructured(xi, eta, zeta, P[0,:])
        y = self.EvaluateBsplineScalarStructured(xi, eta, zeta, P[1,:])
        z = self.EvaluateBsplineScalarStructured(xi, eta, zeta, P[2,:])
  
        # Reshaping the fields  
        ux = np.ascontiguousarray( ux.reshape((nzeta,neta,nxi)))
        uy = np.ascontiguousarray( uy.reshape((nzeta,neta,nxi)))
        uz = np.ascontiguousarray( uz.reshape((nzeta,neta,nxi)))
        
        exx = np.ascontiguousarray( exx.reshape((nzeta,neta,nxi)))
        eyy = np.ascontiguousarray( eyy.reshape((nzeta,neta,nxi)))
        ezz = np.ascontiguousarray( ezz.reshape((nzeta,neta,nxi)))
        exy = np.ascontiguousarray( exy.reshape((nzeta,neta,nxi)))
        exz = np.ascontiguousarray( exz.reshape((nzeta,neta,nxi)))
        eyz = np.ascontiguousarray( eyz.reshape((nzeta,neta,nxi)))
 
        x  =  np.ascontiguousarray( x.reshape((nzeta,neta,nxi)))
        y  =  np.ascontiguousarray( y.reshape((nzeta,neta,nxi)))
        z  =  np.ascontiguousarray( z.reshape((nzeta,neta,nxi)))

 
        
        """ Exporting to VTK using EVTK library """ 
 
        start = (0,0,0)
        end = (nzeta-1,neta-1,nxi-1)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("U",(ux,uy,uz))
        w.addData("Exx", exx) 
        w.addData("Eyy", eyy)
        w.addData("Ezz", ezz) 
        w.addData("Exy", exy)
        w.addData("Exz", exz)
        w.addData("Eyz", eyz) 
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz)) 
        w.appendData(exx) 
        w.appendData(eyy)
        w.appendData(ezz) 
        w.appendData(exy)
        w.appendData(exz)
        w.appendData(eyz) 
 

        w.save()   
    
    def DVC_VtkPlot_StructuredFromFEMesh(self, path, dvcSt, dvcFe, g,  mFe, Ufe,mask=None):
        """
        This function performs the processing of a fe calculation on a structured grid 
        that contains the finite element mesh
        
        dvcSt is the structured dvc 
        dvcFe is the tetra finite element dvc 
        """
        # ie = mFe.LocatePoints(self.pix, self.piy, self.piz
        # in1 = self.pix > np.min(mFe.n[:,0]) 
        # in2 = self.pix < np.max(mFe.n[:,0]) 
        # in3 = self.piy > np.min(mFe.n[:,1]) 
        # in4 = self.piy < np.max(mFe.n[:,1]) 
        # in5 = self.piz > np.min(mFe.n[:,2])  
        # in6 = self.piz < np.max(mFe.n[:,2]) 
        # iIn = np.where(in1*in2*in3*in4*in5*in6)[0]

        if mask =='Disk':
            r  = (np.max(mFe.n[:,0])-np.min(mFe.n[:,0]))/2
            x0 = np.max(mFe.n[:,0]) -r 
            z0 = np.max(mFe.n[:,2]) -r 
            print('Disk estimation')
            print('x0='+str(x0)+', z0='+str(z0)+', r='+str(r))
            d = np.sqrt((self.pix-x0)**2+(self.piz-z0)**2)
            iInDisk = np.where(d<r)[0]
            ie = mFe.LocatePoints3(self.pix[iInDisk], self.piy[iInDisk], self.piz[iInDisk])
            iInE  = np.where(ie!=-1)[0]  # indices of the points of the disk that are in the fe mesh  
            iIn  = iInDisk[iInE] # indices of the points of the visualization grid that are in the fe mesh
        else:
            ie = mFe.LocatePoints3(self.pix, self.piy, self.piz)
            iIn  = np.where(ie!=-1)[0]
            iInE = iIn 
        
        print('Finished location of points')
        
     
        xi,eta,zeta,uxIn,uyIn,uzIn,exxIn,eyyIn,ezzIn,exyIn,exzIn,eyzIn = mFe.DispAndStrainAtPoints(self.pix[iIn],self.piy[iIn],self.piz[iIn],ie[iInE],Ufe,returnIsoParamCoord=True)
        goPhiIn = dvcFe.GoIdu(g, mFe , Ufe, xi,eta,zeta,ie[iInE])

        print('Writing file...')

        nzeta   = len(self.piz1d) 
        neta    = len(self.piy1d)
        nxi     = len(self.pix1d)   
        npoints = len(self.pix)
        """ Displacement and strain """ 
        ux    = np.zeros(npoints)
        uy    = np.zeros(npoints)
        uz    = np.zeros(npoints)
        exx   = np.zeros(npoints)
        eyy   = np.zeros(npoints)
        ezz   = np.zeros(npoints)
        exy   = np.zeros(npoints)
        exz   = np.zeros(npoints)
        eyz   = np.zeros(npoints)
        
        goPhi = np.zeros(npoints)
        
        # meshMask = np.zeros(npoints)
        # meshMask[iIn] = 1 
 
 
        ux[iIn]  = uxIn
        uy[iIn]  = uyIn 
        uz[iIn]  = uzIn 
        exx[iIn] = exxIn 
        eyy[iIn] = eyyIn 
        ezz[iIn] = ezzIn 
        exy[iIn] = exyIn 
        exz[iIn] = exzIn
        eyz[iIn] = eyzIn 
        goPhi[iIn] = goPhiIn
        


        # Fields  
        ux = np.ascontiguousarray( ux.reshape((nzeta,neta,nxi)))
        uy = np.ascontiguousarray( uy.reshape((nzeta,neta,nxi)))
        uz = np.ascontiguousarray( uz.reshape((nzeta,neta,nxi)))
        
        exx = np.ascontiguousarray( exx.reshape((nzeta,neta,nxi)))
        eyy = np.ascontiguousarray( eyy.reshape((nzeta,neta,nxi)))
        ezz = np.ascontiguousarray( ezz.reshape((nzeta,neta,nxi)))
        exy = np.ascontiguousarray( exy.reshape((nzeta,neta,nxi)))
        exz = np.ascontiguousarray( exz.reshape((nzeta,neta,nxi)))
        eyz = np.ascontiguousarray( eyz.reshape((nzeta,neta,nxi)))
     
        x  =  np.ascontiguousarray( self.pix.reshape((nzeta,neta,nxi)))
        y  =  np.ascontiguousarray( self.piy.reshape((nzeta,neta,nxi)))
        z  =  np.ascontiguousarray( self.piz.reshape((nzeta,neta,nxi)))
        
        goPhi = np.ascontiguousarray( goPhi.reshape((nzeta,neta,nxi)))
        
        if mask=='Disk':
            distCenter = np.ascontiguousarray( d.reshape((nzeta,neta,nxi)))

        
        # Rearanging the voxels in the tensor (x,y,z) order 
        nx = self.n_elems[0]
        ny = self.n_elems[1]
        nz = self.n_elems[2]
        nixe = int ( len(self.pix1d)/nx ) 
        niye = int ( len(self.piy1d)/ny )
        nize = int ( len(self.piz1d)/nz )
        
        
        px  = np.kron( np.ones(nx), np.arange(nixe) ) + np.kron( np.arange(nx)*(nixe*niye*nize), np.ones(nixe))   
        pex = np.kron( np.ones(niye), px ) + np.kron( np.arange(niye)*nixe, np.ones(nx*nixe) )
        p  = np.kron( np.ones(ny), pex ) + np.kron( np.arange(ny)*(nixe*niye*nize*nx), np.ones(nixe*niye*nx) )
        ph = np.kron(np.ones(nize), p )
        t  = np.kron( np.arange(nize)*(nixe*niye), np.ones(nx*ny*nixe*niye) )
        i = ph + t    
        I = np.kron( np.ones(nz), i ) + np.kron(  np.arange(nz)*(nixe*niye*nize*nx*ny), np.ones(nx*ny*nixe*niye*nize) )
        I = I.astype('int')
        f     = np.ascontiguousarray( dvcSt.fip[I].reshape((nzeta,neta,nxi)) ) # image f 
        
        res = f - goPhi 
             
        """ Exporting to VTK using EVTK library """ 
     
        start = (0,0,0)
        end = (nzeta-1,neta-1,nxi-1)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)
    
        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("U",(ux,uy,uz))
        w.addData("Exx", exx) 
        w.addData("Eyy", eyy)
        w.addData("Ezz", ezz) 
        w.addData("Exy", exy)
        w.addData("Exz", exz)
        w.addData("Eyz", eyz) 
        w.addData("Ref image", f)
        w.addData("residual", res)
        # w.addData("mesh mask", meshMask)
        if mask=='Disk':
            w.addData("distance Center", distCenter)
        w.closeData("Point")
     
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz)) 
        w.appendData(exx) 
        w.appendData(eyy)
        w.appendData(ezz) 
        w.appendData(exy)
        w.appendData(exz)
        w.appendData(eyz) 
        w.appendData(f)
        w.appendData(res)
        # w.appendData(meshMask)
        if mask=='Disk':
            w.appendData(distCenter)
        
        w.save()   
        print(w.getFileName()+' saved') 

        
        
        
        
        
        
        

    def DVC_VtkPlot_Structured(self, path, dvc, U, res):
        nbf = self.Get_nbf() 
        ndof = 3*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:] # displacement at the control points 
     
        nzeta = len(self.piz1d) 
        neta  = len(self.piy1d)
        nxi   = len(self.pix1d)
        """ Displacement """ 
        ux,duxdx, duxdy, duxdz =  self.EvaluateBsplineAndDerivativesScalarStructured(self.pix1d,self.piy1d,self.piz1d,Ux) 
        uy,duydx, duydy, duydz =  self.EvaluateBsplineAndDerivativesScalarStructured(self.pix1d,self.piy1d,self.piz1d,Uy) 
        uz,duzdx, duzdy, duzdz =  self.EvaluateBsplineAndDerivativesScalarStructured(self.pix1d,self.piy1d,self.piz1d,Uz)
        """ Strain """ 
        exx = duxdx 
        eyy = duydy 
        ezz = duzdz 
        exy = 0.5*(duydx+duxdy)
        exz = 0.5*(duzdx+duxdz)
        eyz = 0.5*(duzdy+duydz)    
        
        
        # Fields  
        ux = np.ascontiguousarray( ux.reshape((nzeta,neta,nxi)))
        uy = np.ascontiguousarray( uy.reshape((nzeta,neta,nxi)))
        uz = np.ascontiguousarray( uz.reshape((nzeta,neta,nxi)))
        
        exx = np.ascontiguousarray( exx.reshape((nzeta,neta,nxi)))
        eyy = np.ascontiguousarray( eyy.reshape((nzeta,neta,nxi)))
        ezz = np.ascontiguousarray( ezz.reshape((nzeta,neta,nxi)))
        exy = np.ascontiguousarray( exy.reshape((nzeta,neta,nxi)))
        exz = np.ascontiguousarray( exz.reshape((nzeta,neta,nxi)))
        eyz = np.ascontiguousarray( eyz.reshape((nzeta,neta,nxi)))
     
        x  =  np.ascontiguousarray( self.pix.reshape((nzeta,neta,nxi)))
        y  =  np.ascontiguousarray( self.piy.reshape((nzeta,neta,nxi)))
        z  =  np.ascontiguousarray( self.piz.reshape((nzeta,neta,nxi)))
        
        res3d  = np.ascontiguousarray( res.reshape((nzeta,neta,nxi)) )
        
        # Rearanging the voxels in the tensor (x,y,z) order 
        nx = self.n_elems[0]
        ny = self.n_elems[1]
        nz = self.n_elems[2]
        nixe = int ( len(self.pix1d)/nx ) 
        niye = int ( len(self.piy1d)/ny )
        nize = int ( len(self.piz1d)/nz )
        
        
        px  = np.kron( np.ones(nx), np.arange(nixe) ) + np.kron( np.arange(nx)*(nixe*niye*nize), np.ones(nixe))   
        pex = np.kron( np.ones(niye), px ) + np.kron( np.arange(niye)*nixe, np.ones(nx*nixe) )
        p  = np.kron( np.ones(ny), pex ) + np.kron( np.arange(ny)*(nixe*niye*nize*nx), np.ones(nixe*niye*nx) )
        ph = np.kron(np.ones(nize), p )
        t  = np.kron( np.arange(nize)*(nixe*niye), np.ones(nx*ny*nixe*niye) )
        i = ph + t    
        I = np.kron( np.ones(nz), i ) + np.kron(  np.arange(nz)*(nixe*niye*nize*nx*ny), np.ones(nx*ny*nixe*niye*nize) )
        I = I.astype('int')
        gl     = np.ascontiguousarray( dvc.fip[I].reshape((nzeta,neta,nxi)) )
             
        """ Exporting to VTK using EVTK library """ 
     
        start = (0,0,0)
        end = (nzeta-1,neta-1,nxi-1)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)
    
        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("U",(ux,uy,uz))
        w.addData("Exx", exx) 
        w.addData("Eyy", eyy)
        w.addData("Ezz", ezz) 
        w.addData("Exy", exy)
        w.addData("Exz", exz)
        w.addData("Eyz", eyz) 
        w.addData("Ref image", gl)
        w.addData("residual", res3d)
        w.closeData("Point")
     
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz)) 
        w.appendData(exx) 
        w.appendData(eyy)
        w.appendData(ezz) 
        w.appendData(exy)
        w.appendData(exz)
        w.appendData(eyz) 
        w.appendData(gl)
        w.appendData(res3d)
        
        w.save()   
        print(w.getFileName()+' saved') 

#%% Finite cell
class Cell3d : 
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
        self.xmin = xmin 
        self.xmax = xmax 
        self.ymin = ymin 
        self.ymax = ymax 
        self.zmin = zmin 
        self.zmax = zmax 
        self.lvl = 0 
        self.bottomBackLeft   = None
        self.bottomBackRight  = None
        self.bottomFrontLeft  = None
        self.bottomFrontRight = None
        self.topBackLeft      = None
        self.topBackRight     = None
        self.topFrontLeft     = None
        self.topFrontRight    = None     
        
    def DecomposeLevelSetIntegration(self,mesh,Phi,lvlmax):
        if self.lvl == lvlmax : 
            # if the cell is homogeneous 
            # take all the cell 
            eps = 1.e-8 
            xc =  np.array([self.xmin+eps, self.xmax-eps, self.xmax-eps, self.xmin+eps, self.xmin+eps, self.xmax-eps, self.xmax-eps, self.xmin+eps])
            yc =  np.array([self.ymin+eps, self.ymin+eps, self.ymax-eps, self.ymax-eps, self.ymin+eps, self.ymin+eps, self.ymax-eps, self.ymax-eps])
            zc =  np.array([self.zmin+eps, self.zmin+eps, self.zmin+eps, self.zmin+eps, self.zmax-eps, self.zmax-eps, self.zmax-eps, self.zmax-eps])
            ls = Phi(xc,yc,zc)  
            if (ls>0).all():
                xim   = np.kron(self.xmin,np.ones(mesh.nbg_xi))
                etam  = np.kron(self.ymin,np.ones(mesh.nbg_eta))
                zetam = np.kron(self.zmin,np.ones(mesh.nbg_zeta))
                
                mesx = self.xmax-self.xmin
                mesy = self.ymax-self.ymin
                mesz = self.zmax-self.zmin
                
                xi         = xim  + 0.5*(mesh.Gauss_xi[0]+1)*mesx
                eta        = etam + 0.5*(mesh.Gauss_eta[0]+1)*mesy
                zeta       = zetam + 0.5*(mesh.Gauss_zeta[0]+1)*mesz 
                                
                xi        = np.kron(mesh.oneZeta, np.kron(mesh.oneEta, xi))
                eta       = np.kron(mesh.oneZeta, np.kron(eta, mesh.oneXi))
                zeta      = np.kron(zeta, np.kron(mesh.oneEta, mesh.oneXi))
                
                
                
                mesh.pix.extend(list(xi))
                mesh.piy.extend(list(eta))
                mesh.piz.extend(list(zeta))
                mesh.wg.extend(mesh.wgRef*mesx*mesy*mesz/8)    
 
                
            elif  (ls<0).all() ==False :
 
                # perform the tesselation procedure using delaunay algorithm 
                
                TriCorners = [] 
                for i in range(8):
                    if ls[i] > 0 :
                        TriCorners.append(np.array([ xc[i],yc[i],zc[i] ]))
                        
                # Distance linearization of the interface 
             
                # Edge 1  0-1 - x
                if ( ls[0] > 0 and ls[1] < 0 ) or ( ls[0] < 0 and ls[1] > 0) :
                   a,b = nb.interpolateLinearly(xc[0],xc[1],ls[0],ls[1])
                   xb =  -b/a 
                   yb =  self.ymin 
                   zb =  self.zmin 
                   TriCorners.append( np.array([xb,yb,zb]))    
                   
                   
                # Edge 2  1-2 - y
                if ( ls[1] > 0 and ls[2] < 0 ) or ( ls[1] < 0 and ls[2] > 0) :
                   a,b = nb.interpolateLinearly(yc[1],yc[2],ls[1],ls[2])
                   xb =  self.xmax
                   yb =  -b/a  
                   zb =  self.zmin
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 3  2-3 - x
                if ( ls[2] > 0 and ls[3] < 0 ) or ( ls[2] < 0 and ls[3] > 0) :
                   a,b = nb.interpolateLinearly(xc[2],xc[3],ls[2],ls[3])
                   xb =  -b/a
                   yb =  self.ymax  
                   zb =  self.zmin 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                   
                # Edge 4  0-3 - y
                if ( ls[0] > 0 and ls[3] < 0 ) or ( ls[0] < 0 and ls[3] > 0) :
                   a,b = nb.interpolateLinearly(yc[0],yc[3],ls[0],ls[3])
                   xb =  self.xmin
                   yb =  -b/a 
                   zb =  self.zmin  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 5  4-5 - x
                if ( ls[4] > 0 and ls[5] < 0 ) or ( ls[4] < 0 and ls[5] > 0) :
                   a,b = nb.interpolateLinearly(xc[4],xc[5],ls[4],ls[5])
                   xb =  -b/a 
                   yb =  self.ymin 
                   zb =  self.zmax 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 6  5-6 - y
                if ( ls[5] > 0 and ls[6] < 0 ) or ( ls[5] < 0 and ls[6] > 0) :
                   a,b = nb.interpolateLinearly(yc[5],yc[6],ls[5],ls[6])
                   xb =  self.xmax
                   yb =  -b/a  
                   zb =  self.zmax   
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 7  6-7 - x
                if ( ls[6] > 0 and ls[7] < 0 ) or ( ls[6] < 0 and ls[7] > 0) :
                   a,b = nb.interpolateLinearly(xc[6],xc[7],ls[6],ls[7])
                   xb =  -b/a
                   yb =  self.ymax  
                   zb =  self.zmax  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 8  7-4 - y
                if ( ls[7] > 0 and ls[4] < 0 ) or ( ls[7] < 0 and ls[4] > 0) :
                   a,b = nb.interpolateLinearly(yc[7],yc[4],ls[7],ls[4])
                   xb =  self.xmin
                   yb =  -b/a 
                   zb =  self.zmax 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 9  0-4 - z
                if ( ls[0] > 0 and ls[4] < 0 ) or ( ls[0] < 0 and ls[4] > 0) :
                   a,b = nb.interpolateLinearly(zc[0],zc[4],ls[0],ls[4])
                   xb =  self.xmin
                   yb =  self.ymin 
                   zb =  -b/a 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 10 1-5 - z
                if ( ls[1] > 0 and ls[5] < 0 ) or ( ls[1] < 0 and ls[5] > 0) :
                   a,b = nb.interpolateLinearly(zc[1],zc[5],ls[1],ls[5])
                   xb =  self.xmax
                   yb =  self.ymin  
                   zb =  -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 11  2-6 -z
                if ( ls[2] > 0 and ls[6] < 0 ) or ( ls[2] < 0 and ls[6] > 0) :
                   a,b = nb.interpolateLinearly(zc[2],zc[6],ls[2],ls[6])
                   xb = self.xmax
                   yb = self.ymax  
                   zb = -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
           
                # Edge 12  3-7 -z 
                if ( ls[3] > 0 and ls[7] < 0 ) or ( ls[3] < 0 and ls[7] > 0) :
                   a,b = nb.interpolateLinearly(zc[3],zc[7],ls[3],ls[7])
                   xb =  self.xmin
                   yb =  self.ymax 
                   zb =  -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
            
                nodes = np.array([TriCorners])[0] 
                tri = Delaunay(nodes)  
                # mesh.integrationCellsCoord.append(['t',tri.points,tri.simplices])   
                for i in range(tri.simplices.shape[0]):
                    n1 = tri.simplices[i,0] ; n2 = tri.simplices[i,1]  ; n3=  tri.simplices[i,2]; n4 = tri.simplices[i,3]
                    x1 = tri.points[n1,0] ; y1 = tri.points[n1,1]; z1 = tri.points[n1,2]  
                    x2 = tri.points[n2,0] ; y2 = tri.points[n2,1]; z2 = tri.points[n2,2]  
                    x3 = tri.points[n3,0] ; y3 = tri.points[n3,1]; z3 = tri.points[n3,2]                     
                    x4 = tri.points[n4,0] ; y4 = tri.points[n4,1]; z4 = tri.points[n4,2]   
                    
 
                    pg = mesh.GaussTetra[0].dot(np.c_[np.array([x1,x2,x3,x4]),
                                                 np.array([y1,y2,y3,y4]),
                                                 np.array([z1,z2,z3,z4])])
                    
                    V = 1/6*( np.abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                               (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                               (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) ) 
                              
                    wg = mesh.GaussTetra[1]*V
                    mesh.pix.extend(list(pg[:,0]))
                    mesh.piy.extend(list(pg[:,1]))
                    mesh.piz.extend(list(pg[:,2]))
                    mesh.wg.extend(wg) 
 
    

        else : 
            nbp = 2**(lvlmax - self.lvl)+1 
            xi   = np.linspace(self.xmin +1.e-8, self.xmax -1.e-8, nbp)
            eta  = np.linspace(self.ymin +1.e-8, self.ymax -1.e-8, nbp)
            zeta = np.linspace(self.zmin +1.e-8, self.zmax -1.e-8, nbp)
            
            # Uniform discretization of the 6 faces 
            fxy, fxz = np.meshgrid(eta,zeta)
            fyx, fyz = np.meshgrid(xi,zeta)
            fzx, fzy = np.meshgrid(xi,eta)
            
            fxyr = fxy.ravel()
            fxzr = fxz.ravel()
            fyxr = fyx.ravel()
            fyzr = fyz.ravel()
            fzxr = fzx.ravel() 
            fzyr = fzy.ravel() 
            
            one = np.ones(nbp**2)
            # Setting all the face points 
            x = np.r_[fyxr, fyxr, one*self.xmin, one*self.xmax, fzxr, fzxr ]
            y = np.r_[one*self.ymin, one*self.ymax, fxyr, fxyr, fzyr, fzyr ]
            z = np.r_[fyzr, fyzr, fxzr, fxzr, one*self.zmin, one*self.zmin]

            ls = Phi(x,y,z)  
            
            if (ls>0).all() :
                # integrate the cell 
                xim   = np.kron(self.xmin,np.ones(mesh.nbg_xi))
                etam  = np.kron(self.ymin,np.ones(mesh.nbg_eta))
                zetam = np.kron(self.zmin,np.ones(mesh.nbg_zeta))
                
                mesx = self.xmax-self.xmin
                mesy = self.ymax-self.ymin
                mesz = self.zmax-self.zmin
                
                xi         = xim  + 0.5*(mesh.Gauss_xi[0]+1)*mesx
                eta        = etam + 0.5*(mesh.Gauss_eta[0]+1)*mesy
                zeta       = zetam + 0.5*(mesh.Gauss_zeta[0]+1)*mesz 
                
             
                xi        = np.kron(mesh.oneZeta, np.kron(mesh.oneEta, xi))
                eta       = np.kron(mesh.oneZeta, np.kron(eta, mesh.oneXi))
                zeta      = np.kron(zeta, np.kron(mesh.oneEta, mesh.oneXi))
 
                
                mesh.pix.extend(list(xi))
                mesh.piy.extend(list(eta))
                mesh.piz.extend(list(zeta))
                mesh.wg.extend(mesh.wgRef*mesx*mesy*mesz/8)   
                
 
            
            elif (ls<0).all() == False : 
                # cut the cell : recursive call
                self.bottomBackLeft   = Cell3d(self.xmin, (self.xmin+self.xmax)/2, self.ymin, (self.ymin+self.ymax)/2, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomBackRight  = Cell3d((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomFrontLeft  = Cell3d(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomFrontRight = Cell3d((self.xmin+self.xmax)/2, self.xmax, (self.ymin+self.ymax)/2, self.ymax, self.zmin, (self.zmin+self.zmax)/2)
                self.topBackLeft      = Cell3d(self.xmin, (self.xmin+self.xmax)/2, self.ymin, (self.ymin+self.ymax)/2, (self.zmin+self.zmax)/2, self.zmax)
                self.topBackRight     = Cell3d((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2, (self.zmin+self.zmax)/2, self.zmax)
                self.topFrontLeft     = Cell3d(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax, (self.zmin+self.zmax)/2, self.zmax)
                self.topFrontRight    = Cell3d((self.xmin+self.xmax)/2, self.xmax, (self.ymin+self.ymax)/2, self.ymax, (self.zmin+self.zmax)/2, self.zmax)
 
                
                self.bottomBackLeft.lvl = self.lvl + 1 
                self.bottomBackRight.lvl = self.lvl + 1 
                self.bottomFrontLeft.lvl = self.lvl + 1 
                self.bottomFrontRight.lvl = self.lvl + 1 
                self.topBackLeft.lvl = self.lvl + 1 
                self.topBackRight.lvl = self.lvl + 1 
                self.topFrontLeft.lvl = self.lvl + 1 
                self.topFrontRight.lvl = self.lvl + 1 
  
                self.bottomBackLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.bottomBackRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.bottomFrontLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.bottomFrontRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.topBackLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.topBackRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.topFrontLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.topFrontRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
 
 

    def DecomposeLevelSetCoord(self,mesh,Phi,lvlmax):
        if self.lvl == lvlmax : 
            # if the cell is homogeneous 
            # take all the cell 
            eps = 1.e-8 
            xc =  np.array([self.xmin+eps, self.xmax-eps, self.xmax-eps, self.xmin+eps, self.xmin+eps, self.xmax-eps, self.xmax-eps, self.xmin+eps])
            yc =  np.array([self.ymin+eps, self.ymin+eps, self.ymax-eps, self.ymax-eps, self.ymin+eps, self.ymin+eps, self.ymax-eps, self.ymax-eps])
            zc =  np.array([self.zmin+eps, self.zmin+eps, self.zmin+eps, self.zmin+eps, self.zmax-eps, self.zmax-eps, self.zmax-eps, self.zmax-eps])
            ls = Phi(xc,yc,zc)  
            if (ls>0).all():
                mesh.integrationCellsCoord.append([ 'c', self.xmin,self.xmax, self.ymin, self.ymax, self.zmin, self.zmax ]) 
            elif  (ls<0).all() ==False :
                # perform the tesselation procedure using delaunay algorithm 
                
                TriCorners = [] 
                for i in range(8):
                    if ls[i] > 0 :
                        TriCorners.append(np.array([ xc[i],yc[i],zc[i] ]))
                        
                # Distance linearization of the interface 
             
                # Edge 1  0-1 - x
                if ( ls[0] > 0 and ls[1] < 0 ) or ( ls[0] < 0 and ls[1] > 0) :
                   a,b = nb.interpolateLinearly(xc[0],xc[1],ls[0],ls[1])
                   xb =  -b/a 
                   yb =  self.ymin 
                   zb =  self.zmin 
                   TriCorners.append( np.array([xb,yb,zb]))    
                   
                   
                # Edge 2  1-2 - y
                if ( ls[1] > 0 and ls[2] < 0 ) or ( ls[1] < 0 and ls[2] > 0) :
                   a,b = nb.interpolateLinearly(yc[1],yc[2],ls[1],ls[2])
                   xb =  self.xmax
                   yb =  -b/a  
                   zb =  self.zmin
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 3  2-3 - x
                if ( ls[2] > 0 and ls[3] < 0 ) or ( ls[2] < 0 and ls[3] > 0) :
                   a,b = nb.interpolateLinearly(xc[2],xc[3],ls[2],ls[3])
                   xb =  -b/a
                   yb =  self.ymax  
                   zb =  self.zmin 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                   
                # Edge 4  0-3 - y
                if ( ls[0] > 0 and ls[3] < 0 ) or ( ls[0] < 0 and ls[3] > 0) :
                   a,b = nb.interpolateLinearly(yc[0],yc[3],ls[0],ls[3])
                   xb =  self.xmin
                   yb =  -b/a 
                   zb =  self.zmin  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 5  4-5 - x
                if ( ls[4] > 0 and ls[5] < 0 ) or ( ls[4] < 0 and ls[5] > 0) :
                   a,b = nb.interpolateLinearly(xc[4],xc[5],ls[4],ls[5])
                   xb =  -b/a 
                   yb =  self.ymin 
                   zb =  self.zmax 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 6  5-6 - y
                if ( ls[5] > 0 and ls[6] < 0 ) or ( ls[5] < 0 and ls[6] > 0) :
                   a,b = nb.interpolateLinearly(yc[5],yc[6],ls[5],ls[6])
                   xb =  self.xmax
                   yb =  -b/a  
                   zb =  self.zmax   
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 7  6-7 - x
                if ( ls[6] > 0 and ls[7] < 0 ) or ( ls[6] < 0 and ls[7] > 0) :
                   a,b = nb.interpolateLinearly(xc[6],xc[7],ls[6],ls[7])
                   xb =  -b/a
                   yb =  self.ymax  
                   zb =  self.zmax  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 8  7-4 - y
                if ( ls[7] > 0 and ls[4] < 0 ) or ( ls[7] < 0 and ls[4] > 0) :
                   a,b = nb.interpolateLinearly(yc[7],yc[4],ls[7],ls[4])
                   xb =  self.xmin
                   yb =  -b/a 
                   zb =  self.zmax 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 9  0-4 - z
                if ( ls[0] > 0 and ls[4] < 0 ) or ( ls[0] < 0 and ls[4] > 0) :
                   a,b = nb.interpolateLinearly(zc[0],zc[4],ls[0],ls[4])
                   xb =  self.xmin
                   yb =  self.ymin 
                   zb =  -b/a 
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 10 1-5 - z
                if ( ls[1] > 0 and ls[5] < 0 ) or ( ls[1] < 0 and ls[5] > 0) :
                   a,b = nb.interpolateLinearly(zc[1],zc[5],ls[1],ls[5])
                   xb =  self.xmax
                   yb =  self.ymin  
                   zb =  -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
                   
                # Edge 11  2-6 -z
                if ( ls[2] > 0 and ls[6] < 0 ) or ( ls[2] < 0 and ls[6] > 0) :
                   a,b = nb.interpolateLinearly(zc[2],zc[6],ls[2],ls[6])
                   xb = self.xmax
                   yb = self.ymax  
                   zb = -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
                   
           
                # Edge 12  3-7 -z 
                if ( ls[3] > 0 and ls[7] < 0 ) or ( ls[3] < 0 and ls[7] > 0) :
                   a,b = nb.interpolateLinearly(zc[3],zc[7],ls[3],ls[7])
                   xb =  self.xmin
                   yb =  self.ymax 
                   zb =  -b/a  
                   TriCorners.append( np.array([xb,yb,zb]))  
            
                nodes = np.array([TriCorners])[0] 
                tri = Delaunay(nodes)  
                mesh.integrationCellsCoord.append(['t',tri.points,tri.simplices])    

        else : 
            nbp = 2**(lvlmax - self.lvl)+1 
            xi   = np.linspace(self.xmin +1.e-8, self.xmax -1.e-8, nbp)
            eta  = np.linspace(self.ymin +1.e-8, self.ymax -1.e-8, nbp)
            zeta = np.linspace(self.zmin +1.e-8, self.zmax -1.e-8, nbp)
            
            # Uniform discretization of the 6 faces 
            fxy, fxz = np.meshgrid(eta,zeta)
            fyx, fyz = np.meshgrid(xi,zeta)
            fzx, fzy = np.meshgrid(xi,eta)
            
            fxyr = fxy.ravel()
            fxzr = fxz.ravel()
            fyxr = fyx.ravel()
            fyzr = fyz.ravel()
            fzxr = fzx.ravel() 
            fzyr = fzy.ravel() 
            
            one = np.ones(nbp**2)
            # Setting all the face points 
            x = np.r_[fyxr, fyxr, one*self.xmin, one*self.xmax, fzxr, fzxr ]
            y = np.r_[one*self.ymin, one*self.ymax, fxyr, fxyr, fzyr, fzyr ]
            z = np.r_[fyzr, fyzr, fxzr, fxzr, one*self.zmin, one*self.zmin]

            ls = Phi(x,y,z)  
            
            if (ls>0).all() :
                # integrate the cell 
                mesh.integrationCellsCoord.append([ 'c', self.xmin,self.xmax, self.ymin, self.ymax, self.zmin, self.zmax ])  
            elif (ls<0).all() == False : 
                # cut the cell : recursive call
                self.bottomBackLeft   = Cell3d(self.xmin, (self.xmin+self.xmax)/2, self.ymin, (self.ymin+self.ymax)/2, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomBackRight  = Cell3d((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomFrontLeft  = Cell3d(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax, self.zmin, (self.zmin+self.zmax)/2)
                self.bottomFrontRight = Cell3d((self.xmin+self.xmax)/2, self.xmax, (self.ymin+self.ymax)/2, self.ymax, self.zmin, (self.zmin+self.zmax)/2)
                self.topBackLeft      = Cell3d(self.xmin, (self.xmin+self.xmax)/2, self.ymin, (self.ymin+self.ymax)/2, (self.zmin+self.zmax)/2, self.zmax)
                self.topBackRight     = Cell3d((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2, (self.zmin+self.zmax)/2, self.zmax)
                self.topFrontLeft     = Cell3d(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax, (self.zmin+self.zmax)/2, self.zmax)
                self.topFrontRight    = Cell3d((self.xmin+self.xmax)/2, self.xmax, (self.ymin+self.ymax)/2, self.ymax, (self.zmin+self.zmax)/2, self.zmax)
 
                
                self.bottomBackLeft.lvl = self.lvl + 1 
                self.bottomBackRight.lvl = self.lvl + 1 
                self.bottomFrontLeft.lvl = self.lvl + 1 
                self.bottomFrontRight.lvl = self.lvl + 1 
                self.topBackLeft.lvl = self.lvl + 1 
                self.topBackRight.lvl = self.lvl + 1 
                self.topFrontLeft.lvl = self.lvl + 1 
                self.topFrontRight.lvl = self.lvl + 1 
  
                self.bottomBackLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.bottomBackRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.bottomFrontLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.bottomFrontRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.topBackLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.topBackRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.topFrontLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.topFrontRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
#%% 3D IGA 
class MeshVolume : 
    """ 
    Mesh class volume : IGA
    Can be ameliorated or integrated in the general class Mesh (to do in the future)
    """                 
    def __init__(self,X,Y,Z,W,pp,xxsi):
        self.pp = pp 
        self.xxsi = xxsi 
        self.X = X # x coordinates control point matrix  
        self.Y = Y # y coordinates control points matrix 
        self.Z = Z # z coordinates control points matrix 
        self.W = W # Nurb weight matrix 
        
        
        self.phi = 0          
        self.dphidx = 0      
        self.dphidy = 0       
        self.dphidz = 0      
        
        
        self.pix = 0  # Integration points in xi   = x direction 
        self.piy = 0  # Integration points in eta  = y direction   
        self.piz = 0  # Integration points in zeta = z direction  
        self.wg  = 0  # Integration weights diagonal matrix 
        
        
        self.phix    = 0
        self.phiy    = 0
        self.phiz    = 0
        
        self.dphixdx = 0
        self.dphixdy = 0
        self.dphixdz = 0
        
        self.dphiydx = 0
        self.dphiydy = 0
        self.dphiydz = 0
        
        self.dphizdx = 0
        self.dphizdy = 0
        self.dphizdz = 0
        

    def DegElevation(self,ppnew):
        ''' INPUT target degree, example: ppnew=[2,3]
        '''
        t = ppnew - self.pp
        # to homogeneous coordinates
        self.X = self.X*self.W
        self.Y = self.Y*self.W 
        self.Z = self.Z*self.W 
                
        if t[0] !=0 :
            # degree elevation in xi direction 
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xr = np.zeros((nbf_zeta,nbf_xi+t[0],nbf_eta))
            Yr = np.zeros((nbf_zeta,nbf_xi+t[0],nbf_eta))
            Zr = np.zeros((nbf_zeta,nbf_xi+t[0],nbf_eta))
            Wr = np.zeros((nbf_zeta,nbf_xi+t[0],nbf_eta))
        
            
            # loop over eta perform degree elevation on each iso-surface 
            for i in range(nbf_zeta):
                Xe, newKnotXi = nb.bspdegelev(self.pp[0],self.X[i].T,self.xxsi[0],t[0])
                Ye, newKnotXi = nb.bspdegelev(self.pp[0],self.Y[i].T,self.xxsi[0],t[0])
                Ze, newKnotXi = nb.bspdegelev(self.pp[0],self.Z[i].T,self.xxsi[0],t[0])
                We, newKnotXi = nb.bspdegelev(self.pp[0],self.W[i].T,self.xxsi[0],t[0])
                Xr[i] = Xe.T 
                Yr[i] = Ye.T 
                Zr[i] = Ze.T
                Wr[i] = We.T  
            self.X = Xr 
            self.Y = Yr 
            self.Z = Zr 
            self.W = Wr 
            self.xxsi[0] = newKnotXi
            
        if t[1] !=0 :
            # degree elevation in eta direction 
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+t[1]))
            Yr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+t[1]))
            Zr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+t[1]))
            Wr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+t[1]))
        
            
            # loop over eta perform degree elevation on each iso-surface 
            for i in range(nbf_zeta):
                Xe, newKnotEta = nb.bspdegelev(self.pp[1],self.X[i],self.xxsi[1],t[1])
                Ye, newKnotEta = nb.bspdegelev(self.pp[1],self.Y[i],self.xxsi[1],t[1])
                Ze, newKnotEta = nb.bspdegelev(self.pp[1],self.Z[i],self.xxsi[1],t[1])
                We, newKnotEta = nb.bspdegelev(self.pp[1],self.W[i],self.xxsi[1],t[1])
                Xr[i] = Xe 
                Yr[i] = Ye 
                Zr[i] = Ze 
                Wr[i] = We  
            self.X = Xr 
            self.Y = Yr 
            self.Z = Zr 
            self.W = Wr 
            self.xxsi[1] = newKnotEta
            
        if t[2] !=0 :
            # degree elevation in zeta direction 
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xrm = self.X.ravel().reshape((self.X.shape[0],self.X.shape[1]*self.X.shape[2])).T
            Yrm = self.Y.ravel().reshape((self.Y.shape[0],self.Y.shape[1]*self.Y.shape[2])).T
            Zrm = self.Z.ravel().reshape((self.Z.shape[0],self.Z.shape[1]*self.Z.shape[2])).T
            Wrm = self.W.ravel().reshape((self.W.shape[0],self.W.shape[1]*self.W.shape[2])).T
            
            Xe, newKnotZeta = nb.bspdegelev(self.pp[2],Xrm,self.xxsi[2],t[2])
            Ye, newKnotZeta = nb.bspdegelev(self.pp[2],Yrm,self.xxsi[2],t[2])
            Ze, newKnotZeta = nb.bspdegelev(self.pp[2],Zrm,self.xxsi[2],t[2])
            We, newKnotZeta = nb.bspdegelev(self.pp[2],Wrm,self.xxsi[2],t[2])
            
          
            self.X = Xe.T.reshape((Xe.shape[1],nbf_xi,nbf_eta))
            self.Y = Ye.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.Z = Ze.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.W = We.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.xxsi[2] = newKnotZeta 
            
            
        self.pp=ppnew
        # back to real coordinates
        self.X = self.X/self.W
        self.Y = self.Y/self.W
        self.Z = self.Z/self.W
    
    def KnotInsertion(self,u_refinement):
        '''
        INPUT
        ubar=dict()
        ubar[0] = ubarXi
        ubar[1] = ubarEta
        ubar[2] = ubarZeta 
        '''
 
        # to homogeneous coordinates
        self.X = self.X*self.W
        self.Y = self.Y*self.W 
        self.Z = self.Z*self.W 
        
        nXi   = np.size(u_refinement[0])
        nEta  = np.size(u_refinement[1])
        nZeta = np.size(u_refinement[2])
        
        if nXi !=0: 
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xr = np.zeros((nbf_zeta,nbf_xi+nXi,nbf_eta))
            Yr = np.zeros((nbf_zeta,nbf_xi+nXi,nbf_eta))
            Zr = np.zeros((nbf_zeta,nbf_xi+nXi,nbf_eta))
            Wr = np.zeros((nbf_zeta,nbf_xi+nXi,nbf_eta))
        
            
        
            for i in range(nbf_zeta):
                Xe, newKnotXi = nb.bspkntins(self.pp[0],self.X[i].T,self.xxsi[0],u_refinement[0])
                Ye, newKnotXi = nb.bspkntins(self.pp[0],self.Y[i].T,self.xxsi[0],u_refinement[0])
                Ze, newKnotXi = nb.bspkntins(self.pp[0],self.Z[i].T,self.xxsi[0],u_refinement[0])
                We, newKnotXi = nb.bspkntins(self.pp[0],self.W[i].T,self.xxsi[0],u_refinement[0])
                Xr[i] = Xe.T 
                Yr[i] = Ye.T 
                Zr[i] = Ze.T
                Wr[i] = We.T  
            self.X = Xr 
            self.Y = Yr 
            self.Z = Zr 
            self.W = Wr 
            self.xxsi[0] = newKnotXi
            
            
        if nEta != 0:
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+nEta))
            Yr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+nEta))
            Zr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+nEta))
            Wr = np.zeros((nbf_zeta ,nbf_xi  ,nbf_eta+nEta))
        
            
   
            for i in range(nbf_zeta):
                Xe, newKnotEta = nb.bspkntins(self.pp[1],self.X[i],self.xxsi[1],u_refinement[1])
                Ye, newKnotEta = nb.bspkntins(self.pp[1],self.Y[i],self.xxsi[1],u_refinement[1])
                Ze, newKnotEta = nb.bspkntins(self.pp[1],self.Z[i],self.xxsi[1],u_refinement[1])
                We, newKnotEta = nb.bspkntins(self.pp[1],self.W[i],self.xxsi[1],u_refinement[1])
                Xr[i] = Xe 
                Yr[i] = Ye 
                Zr[i] = Ze 
                Wr[i] = We  
            self.X = Xr 
            self.Y = Yr 
            self.Z = Zr 
            self.W = Wr 
            self.xxsi[1] = newKnotEta
            
            
        if nZeta !=0: 
            nbf_xi   = self.X.shape[1]
            nbf_eta  = self.X.shape[2]
            nbf_zeta = self.X.shape[0]
            
            Xrm = self.X.ravel().reshape((self.X.shape[0],self.X.shape[1]*self.X.shape[2])).T
            Yrm = self.Y.ravel().reshape((self.Y.shape[0],self.Y.shape[1]*self.Y.shape[2])).T
            Zrm = self.Z.ravel().reshape((self.Z.shape[0],self.Z.shape[1]*self.Z.shape[2])).T
            Wrm = self.W.ravel().reshape((self.W.shape[0],self.W.shape[1]*self.W.shape[2])).T
 
    
            Xe, newKnotZeta = nb.bspkntins(self.pp[2],Xrm,self.xxsi[2],u_refinement[2])
            Ye, newKnotZeta = nb.bspkntins(self.pp[2],Yrm,self.xxsi[2],u_refinement[2])
            Ze, newKnotZeta = nb.bspkntins(self.pp[2],Zrm,self.xxsi[2],u_refinement[2])
            We, newKnotZeta = nb.bspkntins(self.pp[2],Wrm,self.xxsi[2],u_refinement[2])
            
          
            self.X = Xe.T.reshape((Xe.shape[1],nbf_xi,nbf_eta))
            self.Y = Ye.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.Z = Ze.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.W = We.T.reshape((Xe.shape[1],nbf_xi,nbf_eta)) 
            self.xxsi[2] = newKnotZeta 
 
 
        # back to real coordinates
        self.X = self.X/self.W
        self.Y = self.Y/self.W
        self.Z = self.Z/self.W

        
    def Get_NumberOfunibf(self):
        """ Number of control points in u and v direction  """
        return self.X.shape[1], self.X.shape[2], self.X.shape[0]
    def Get_Btot(self):
        """List of control points """
        P=np.c_[(np.transpose(self.X,(0,2,1)).reshape(self.X.shape[2]*self.X.shape[0],-1).T ).ravel(order='F'), 
                (np.transpose(self.Y,(0,2,1)).reshape(self.Y.shape[2]*self.Y.shape[0],-1).T ).ravel(order='F'),
                (np.transpose(self.Z,(0,2,1)).reshape(self.Z.shape[2]*self.Z.shape[0],-1).T ).ravel(order='F'),
                (np.transpose(self.W,(0,2,1)).reshape(self.W.shape[2]*self.W.shape[0],-1).T ).ravel(order='F')]   
 
        return P.T
        
    def IsRational(self):
        """ Returns True if the B-spline is rational"""
        return (self.W!=1).any()

    def GaussIntegration(self):
        """ Gauss integration: build of the global differential operators """
        nbg_xi    = self.pp[0]+1 
        nbg_eta   = self.pp[1]+1  
        nbg_zeta  = self.pp[2]+1 

        Gauss_xi   =  nb.GaussLegendre(nbg_xi)
        Gauss_eta  =  nb.GaussLegendre(nbg_eta)
        Gauss_zeta =  nb.GaussLegendre(nbg_zeta)
        
        nbf_xi , nbf_eta, nbf_zeta = self.Get_NumberOfunibf() 
        nbf = nbf_xi * nbf_eta * nbf_zeta 
 

        e_xi   = np.unique(self.xxsi[0]) ; ne_xi   = e_xi.shape[0]-1
        e_eta  = np.unique(self.xxsi[1]) ; ne_eta  = e_eta.shape[0]-1
        e_zeta = np.unique(self.xxsi[2]) ; ne_zeta = e_zeta.shape[0]-1
        
        xi_min    =  np.kron(e_xi[:-1],np.ones(nbg_xi))    
        eta_min   =  np.kron(e_eta[:-1],np.ones(nbg_eta))  
        zeta_min  =  np.kron(e_zeta[:-1],np.ones(nbg_zeta))  
        
        
        xi_g      =  np.kron(np.ones(ne_xi)  , Gauss_xi[0])
        eta_g     =  np.kron(np.ones(ne_eta) , Gauss_eta[0])
        zeta_g    =  np.kron(np.ones(ne_zeta) , Gauss_zeta[0])
        
        """ Measures of elements """
        mes_xi   = e_xi[1:]  - e_xi[:-1]
        mes_eta  = e_eta[1:] - e_eta[:-1]
        mes_zeta = e_zeta[1:] - e_zeta[:-1]
        
        mes_xi   = np.kron(mes_xi,np.ones(nbg_xi))   
        mes_eta  = np.kron(mes_eta,np.ones(nbg_eta)) 
        mes_zeta = np.kron(mes_zeta,np.ones(nbg_zeta)) 
        
         
        """ Going from the reference element to the parametric space  """ 
        xi         = xi_min   + 0.5*(xi_g+1)*mes_xi       # Aranged gauss points in  xi direction  
        eta        = eta_min  + 0.5*(eta_g+1)*mes_eta     # Aranged gauss points in  eta direction 
        zeta       = zeta_min + 0.5*(zeta_g+1)*mes_zeta   # Aranged gauss points in  zeta direction 
        
        phi_xi   , dphi_xi    = nb.global_basisfuns(self.pp[0],self.xxsi[0],xi)
        phi_eta  , dphi_eta   = nb.global_basisfuns(self.pp[1],self.xxsi[1],eta)
        phi_zeta , dphi_zeta  = nb.global_basisfuns(self.pp[2],self.xxsi[2],zeta)
        
        
        phiEta_phiXi    =  sps.kron(phi_eta, phi_xi,  'csc') 
        phi             =  sps.kron(phi_zeta, phiEta_phiXi,  'csc') 
        
        phiEta_dphidXi  =  sps.kron(phi_eta,dphi_xi,  'csc') 
        dphidxi         =  sps.kron(phi_zeta, phiEta_dphidXi,  'csc') 
        
        dphiEta_phiXi   =  sps.kron(dphi_eta, phi_xi,  'csc') 
        dphideta        =  sps.kron(phi_zeta, dphiEta_phiXi,  'csc') 
   
        dphidzeta       =  sps.kron(dphi_zeta, phiEta_phiXi,  'csc') 
        
        self.npg        = phi.shape[0]
        
        wg_xi        =  np.kron(np.ones(ne_xi) , Gauss_xi[1])
        wg_eta       =  np.kron(np.ones(ne_eta), Gauss_eta[1])
        wg_zeta      =  np.kron(np.ones(ne_zeta), Gauss_zeta[1])
        
        oneXi   = np.ones(xi.shape[0])
        oneEta  = np.ones(eta.shape[0])
        oneZeta = np.ones(zeta.shape[0])
        
        
        
        mes_xi   = np.kron(oneZeta, np.kron(oneEta,mes_xi)) 
        mes_eta  = np.kron(oneZeta, np.kron(mes_eta,oneXi)) 
        mes_zeta = np.kron(mes_zeta, np.kron(oneEta,oneXi))  
          
        P = self.Get_Btot()
         
  
        isRational = self.IsRational()
 
        if isRational == False:
            
            """ Jacobian elements"""
            dxdxi   = dphidxi.dot(P[0,:])      
            dxdeta  = dphideta.dot(P[0,:])     
            dxdzeta = dphidzeta.dot(P[0,:])    
            
            dydxi   = dphidxi.dot(P[1,:])     
            dydeta  = dphideta.dot(P[1,:])     
            dydzeta = dphidzeta.dot(P[1,:])    
            
            dzdxi   = dphidxi.dot(P[2,:])      
            dzdeta  = dphideta.dot(P[2,:])     
            dzdzeta = dphidzeta.dot(P[2,:])    
             
 
            detJ   = dxdxi*dydeta*dzdzeta + \
                     dydxi*dzdeta*dxdzeta + \
                     dzdxi*dxdeta*dydzeta - \
                     dzdxi*dydeta*dxdzeta - \
                     dydxi*dxdeta*dzdzeta - \
                     dxdxi*dzdeta*dydzeta
 
            dphidx = sps.diags((dydeta*dzdzeta-dzdeta*dydzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ).dot(dphideta) +\
                     sps.diags((dydxi*dzdeta-dzdxi*dydeta)/detJ).dot(dphidzeta)
                
        
            dphidy = sps.diags((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ).dot(dphidxi) +\
                     sps.diags((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags(-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ).dot(dphidzeta)
                     
            dphidz = sps.diags((dxdeta*dydzeta-dydeta*dxdzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags((dxdxi*dydeta-dydxi*dxdeta)/detJ).dot(dphidzeta)
                     
            """ Integration weights x isoparametric transformation x physical transformation """ 
            wg     =  np.kron(wg_zeta,np.kron(wg_eta, wg_xi))*np.abs(detJ)*mes_xi*mes_eta*mes_zeta/8
            
        if isRational == True: 

             
            denom      = phi.dot(       P[3,:]   )
            denomdxi   = dphidxi.dot(   P[3,:]   )
            denomdeta  = dphideta.dot(  P[3,:]   )
            denomdzeta = dphidzeta.dot( P[3,:]   )

            Ddenom      = sps.diags(denom)
            Ddenom_squ  = sps.diags(1/denom**2)
            Ddenomdxi   = sps.diags(denomdxi)
            Ddenomdeta  = sps.diags(denomdeta)
            Ddenomdzeta = sps.diags(denomdzeta) 
            
            
            Wphi       = phi.multiply(P[3,:])
            Wdphidxi   = dphidxi.multiply(P[3,:])
            Wdphideta  = dphideta.multiply(P[3,:])
            Wdphidzeta = dphidzeta.multiply(P[3,:])
            
            dNdxi   = Ddenom_squ.dot (  Ddenom.dot(Wdphidxi)   - Ddenomdxi.dot(Wphi)   ) 
            dNdeta  = Ddenom_squ.dot (  Ddenom.dot(Wdphideta)  - Ddenomdeta.dot(Wphi)  )
            dNdzeta = Ddenom_squ.dot (  Ddenom.dot(Wdphidzeta) - Ddenomdzeta.dot(Wphi) )
            
            
            dxdxi   = dNdxi.dot(P[0,:])
            dxdeta  = dNdeta.dot(P[0,:])
            dxdzeta = dNdzeta.dot(P[0,:])
            
            dydxi   = dNdxi.dot(P[1,:])
            dydeta  = dNdeta.dot(P[1,:])
            dydzeta = dNdzeta.dot(P[1,:])
            
            dzdxi   = dNdxi.dot(P[2,:])
            dzdeta  = dNdeta.dot(P[2,:])
            dzdzeta = dNdzeta.dot(P[2,:])
            
            detJ   = dxdxi*dydeta*dzdzeta + \
                     dydxi*dzdeta*dxdzeta + \
                     dzdxi*dxdeta*dydzeta - \
                     dzdxi*dydeta*dxdzeta - \
                     dydxi*dxdeta*dzdzeta - \
                     dxdxi*dzdeta*dydzeta
 
            
            dphidx = sps.diags((dydeta*dzdzeta-dzdeta*dydzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ).dot(dphideta) +\
                     sps.diags((dydxi*dzdeta-dzdxi*dydeta)/detJ).dot(dphidzeta)
                
        
            dphidy = sps.diags((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ).dot(dphidxi) +\
                     sps.diags((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags(-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ).dot(dphidzeta)
                     
            dphidz = sps.diags((dxdeta*dydzeta-dydeta*dxdzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags((dxdxi*dydeta-dydxi*dxdeta)/detJ).dot(dphidzeta)
                     
            """ Integration weights x isoparametric transformation x physical transformation """ 
            wg     =  np.kron(wg_zeta,np.kron(wg_eta, wg_xi))*np.abs(detJ)*mes_xi*mes_eta*mes_zeta/8
        
        
        self.wg = sps.diags(wg)        
        zero    = sps.csr_matrix((self.npg,nbf))

        self.phix    = sps.hstack((phi,zero,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi,zero)   ,  'csc')
        self.phiz    = sps.hstack((zero,zero,phi)  ,  'csc')
        
        self.dphixdx = sps.hstack((dphidx,zero,zero),  'csc')
        self.dphixdy = sps.hstack((dphidy,zero,zero),  'csc')
        self.dphixdz = sps.hstack((dphidz,zero,zero),  'csc')
        
        self.dphiydx = sps.hstack((zero,dphidx,zero),  'csc')
        self.dphiydy = sps.hstack((zero,dphidy,zero),  'csc')
        self.dphiydz = sps.hstack((zero,dphidz,zero),  'csc')
        
        self.dphizdx = sps.hstack((zero,zero,dphidx),  'csc')
        self.dphizdy = sps.hstack((zero,zero,dphidy),  'csc')
        self.dphizdz = sps.hstack((zero,zero,dphidz),  'csc')
        
        
        
    

    def GetBfForGridPoints(self,xi,eta,zeta):
        """ xi,eta,zeta are the 1d points 
        This method computes the basis functions on the mesh-grid point 
        obtained from the 1d vector points xi, eta and zeta  
        """
 
        phi_xi   , dphi_xi    = nb.global_basisfuns(self.pp[0],self.xxsi[0],xi)
        phi_eta  , dphi_eta   = nb.global_basisfuns(self.pp[1],self.xxsi[1],eta)
        phi_zeta , dphi_zeta  = nb.global_basisfuns(self.pp[2],self.xxsi[2],zeta)
        
        
        phiEta_phiXi    =  sps.kron(phi_eta, phi_xi,  'csc') 
        phi             =  sps.kron(phi_zeta, phiEta_phiXi,  'csc') 
        
        phiEta_dphidXi  =  sps.kron(phi_eta,dphi_xi,  'csc') 
        dphidxi         =  sps.kron(phi_zeta, phiEta_dphidXi,  'csc') 
        
        dphiEta_phiXi   =  sps.kron(dphi_eta, phi_xi,  'csc') 
        dphideta        =  sps.kron(phi_zeta, dphiEta_phiXi,  'csc') 
   
        dphidzeta       =  sps.kron(dphi_zeta, phiEta_phiXi,  'csc') 
 
        
        P = self.Get_Btot()
        
        isRational = self.IsRational()
        if isRational == False:
            
            """ Jacobian elements"""
            dxdxi   = dphidxi.dot(P[0,:])      
            dxdeta  = dphideta.dot(P[0,:])     
            dxdzeta = dphidzeta.dot(P[0,:])    
            
            dydxi   = dphidxi.dot(P[1,:])     
            dydeta  = dphideta.dot(P[1,:])     
            dydzeta = dphidzeta.dot(P[1,:])    
            
            dzdxi   = dphidxi.dot(P[2,:])      
            dzdeta  = dphideta.dot(P[2,:])     
            dzdzeta = dphidzeta.dot(P[2,:])    
             
            
 
            detJ   = dxdxi*dydeta*dzdzeta + \
                     dydxi*dzdeta*dxdzeta + \
                     dzdxi*dxdeta*dydzeta - \
                     dzdxi*dydeta*dxdzeta - \
                     dydxi*dxdeta*dzdzeta - \
                     dxdxi*dzdeta*dydzeta
 
           
            dphidx = sps.diags((dydeta*dzdzeta-dzdeta*dydzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ).dot(dphideta) +\
                     sps.diags((dydxi*dzdeta-dzdxi*dydeta)/detJ).dot(dphidzeta)
                
        
            dphidy = sps.diags((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ).dot(dphidxi) +\
                     sps.diags((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags(-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ).dot(dphidzeta)
                     
            dphidz = sps.diags((dxdeta*dydzeta-dydeta*dxdzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags((dxdxi*dydeta-dydxi*dxdeta)/detJ).dot(dphidzeta)
 
            
            # Univariate basis functions if needed 
            
#            Nxi  = phi_xi  
#            Neta = phi_eta 
            N = phi 
            
        if isRational == True: 

             
            denom      = phi.dot(       P[3,:]   )
            denomdxi   = dphidxi.dot(   P[3,:]   )
            denomdeta  = dphideta.dot(  P[3,:]   )
            denomdzeta = dphidzeta.dot( P[3,:]   )

            Ddenom      = sps.diags(denom)
            Ddenom_squ  = sps.diags(1/denom**2)
            Ddenomdxi   = sps.diags(denomdxi)
            Ddenomdeta  = sps.diags(denomdeta)
            Ddenomdzeta = sps.diags(denomdzeta) 
            
            
            Wphi       = phi.multiply(P[3,:])
            Wdphidxi   = dphidxi.multiply(P[3,:])
            Wdphideta  = dphideta.multiply(P[3,:])
            Wdphidzeta = dphidzeta.multiply(P[3,:])
            
            dNdxi   = Ddenom_squ.dot (  Ddenom.dot(Wdphidxi)   - Ddenomdxi.dot(Wphi)   ) 
            dNdeta  = Ddenom_squ.dot (  Ddenom.dot(Wdphideta)  - Ddenomdeta.dot(Wphi)  )
            dNdzeta = Ddenom_squ.dot (  Ddenom.dot(Wdphidzeta) - Ddenomdzeta.dot(Wphi) )
            
            
            dxdxi   = dNdxi.dot(P[0,:])
            dxdeta  = dNdeta.dot(P[0,:])
            dxdzeta = dNdzeta.dot(P[0,:])
            
            dydxi   = dNdxi.dot(P[1,:])
            dydeta  = dNdeta.dot(P[1,:])
            dydzeta = dNdzeta.dot(P[1,:])
            
            dzdxi   = dNdxi.dot(P[2,:])
            dzdeta  = dNdeta.dot(P[2,:])
            dzdzeta = dNdzeta.dot(P[2,:])
            
            detJ   = dxdxi*dydeta*dzdzeta + \
                     dydxi*dzdeta*dxdzeta + \
                     dzdxi*dxdeta*dydzeta - \
                     dzdxi*dydeta*dxdzeta - \
                     dydxi*dxdeta*dzdzeta - \
                     dxdxi*dzdeta*dydzeta
    
            
            dphidx = sps.diags((dydeta*dzdzeta-dzdeta*dydzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ).dot(dphideta) +\
                     sps.diags((dydxi*dzdeta-dzdxi*dydeta)/detJ).dot(dphidzeta)
                
        
            dphidy = sps.diags((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ).dot(dphidxi) +\
                     sps.diags((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags(-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ).dot(dphidzeta)
                     
            dphidz = sps.diags((dxdeta*dydzeta-dydeta*dxdzeta)/detJ).dot(dphidxi) +\
                     sps.diags(-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ).dot(dphideta) +\
                     sps.diags((dxdxi*dydeta-dydxi*dxdeta)/detJ).dot(dphidzeta)
                     
 
            N = sps.diags(1/denom).dot(Wphi)
            
            # Univariate basis functions if needed 
            
#            w_xi  = self.ctrlPts[2][:,0] 
#            w_eta = self.ctrlPts[2][0,:]
#            denom_xi  =  phi_xi.dot(w_xi) 
#            denom_eta =  phi_eta.dot(w_eta) 
#            
#            Nxi  =   sps.diags(1/denom_xi).dot(phi_xi).multiply(w_xi) 
#            Neta =   sps.diags(1/denom_eta).dot(phi_eta).multiply(w_eta) 
              
        return N, dphidx, dphidy, dphidz  
    
 
    
    def VtkIsoLines(self, path, U, neval): 
        
        nbf_xi , nbf_eta, nbf_zeta = self.Get_NumberOfunibf() 
        nbf = nbf_xi * nbf_eta * nbf_zeta 
        dim = 3 
        ndof = dim*nbf    

        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:] # displacement at the control points 
        P = self.Get_Btot() # control points 
        
        xi   = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta  = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
        zeta = np.linspace(self.xxsi[2][self.pp[2]] , self.xxsi[2][-self.pp[2]] , neval[2])
        
          
        # Iso parameters for the univariate elemnts 
        xiu   = np.unique(self.xxsi[0])   
        etau  = np.unique(self.xxsi[1])   
        zetau = np.unique(self.xxsi[2])
        
        # Univariate basis functions 
        # (1,2,3)=(xi,eta,zeta)
        phi_xi1   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xi)
        phi_eta1  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],etau)
        phi_zeta1 = nb.global_basisfunsWd(self.pp[2],self.xxsi[2],zetau)
        
        phi_xi2   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xiu)
        phi_eta2  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],eta)
        phi_zeta2 = nb.global_basisfunsWd(self.pp[2],self.xxsi[2],zetau)
        
        phi_xi3   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xiu)
        phi_eta3  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],etau)
        phi_zeta3 = nb.global_basisfunsWd(self.pp[2],self.xxsi[2],zeta)
        
        
        isRational = self.IsRational()
        if isRational == False :
            phi1      = sps.kron( phi_zeta1, sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') )  
            phi2      = sps.kron( phi_zeta2, sps.kron( phi_eta2  ,  phi_xi2   ,  'csc') )  
            phi3      = sps.kron( phi_zeta3, sps.kron( phi_eta3  ,  phi_xi3   ,  'csc') )
        if isRational == True :
            phi1      = sps.kron( phi_zeta1, sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') )  
            phi2      = sps.kron( phi_zeta2, sps.kron( phi_eta2  ,  phi_xi2   ,  'csc') )  
            phi3      = sps.kron( phi_zeta3, sps.kron( phi_eta3  ,  phi_xi3   ,  'csc') )
            Wphi1     = phi1.multiply(P[3,:])
            Wphi2     = phi2.multiply(P[3,:])
            Wphi3     = phi3.multiply(P[3,:])
            denom1    = phi1.dot( P[3,:])
            denom2    = phi2.dot( P[3,:])
            denom3    = phi3.dot( P[3,:])
            phi1 = sps.diags(1/denom1).dot(Wphi1)
            phi2 = sps.diags(1/denom2).dot(Wphi2)
            phi3 = sps.diags(1/denom3).dot(Wphi3) 
        
        xe1  = phi1.dot(P[0,:]) 
        ye1  = phi1.dot(P[1,:])
        ze1  = phi1.dot(P[2,:])
        
        xe2  = phi2.dot(P[0,:]) 
        ye2  = phi2.dot(P[1,:])
        ze2  = phi2.dot(P[2,:])
        
        xe3  = phi3.dot(P[0,:]) 
        ye3  = phi3.dot(P[1,:])
        ze3  = phi3.dot(P[2,:])
        
 
   
        xe1 = xe1.reshape((zetau.size,etau.size,neval[0]))
        ye1 = ye1.reshape((zetau.size,etau.size,neval[0]))
        ze1 = ze1.reshape((zetau.size,etau.size,neval[0]))
        
        xe2 = xe2.reshape((zetau.size,neval[1],xiu.size))
        ye2 = ye2.reshape((zetau.size,neval[1],xiu.size))
        ze2 = ze2.reshape((zetau.size,neval[1],xiu.size))
        
        xe3 = xe3.reshape((neval[2],etau.size,xiu.size))
        ye3 = ye3.reshape((neval[2],etau.size,xiu.size))
        ze3 = ze3.reshape((neval[2],etau.size,xiu.size))
        
        
        uxe1 = phi1.dot(Ux)
        uye1 = phi1.dot(Uy)
        uze1 = phi1.dot(Uz)
        
        uxe2 = phi2.dot(Ux)
        uye2 = phi2.dot(Uy)
        uze2 = phi2.dot(Uz)
        
        uxe3 = phi3.dot(Ux)
        uye3 = phi3.dot(Uy)
        uze3 = phi3.dot(Uz)
        
        
        uxe1 = uxe1.reshape((zetau.size,etau.size,neval[0]))
        uye1 = uye1.reshape((zetau.size,etau.size,neval[0]))
        uze1 = uze1.reshape((zetau.size,etau.size,neval[0]))
        
        uxe2 = uxe2.reshape((zetau.size,neval[1],xiu.size)) 
        uye2 = uye2.reshape((zetau.size,neval[1],xiu.size)) 
        uze2 = uze2.reshape((zetau.size,neval[1],xiu.size)) 
        
        uxe3 = uxe3.reshape((neval[2],etau.size,xiu.size)) 
        uye3 = uye3.reshape((neval[2],etau.size,xiu.size)) 
        uze3 = uze3.reshape((neval[2],etau.size,xiu.size)) 
        
        
        """ Getting the Mesh Lines"""
        
        Ix_xi   =  xe1.ravel()
        Iy_xi   =  ye1.ravel()
        Iz_xi   =  ze1.ravel()
        
        ux_xi =  uxe1.ravel()
        uy_xi =  uye1.ravel()
        uz_xi =  uze1.ravel()
        
        Ix_eta   =  xe2.transpose(0,2,1).ravel()
        Iy_eta   =  ye2.transpose(0,2,1).ravel()
        Iz_eta   =  ze2.transpose(0,2,1).ravel()
        
        ux_eta =  uxe2.transpose(0,2,1).ravel()
        uy_eta =  uye2.transpose(0,2,1).ravel()
        uz_eta =  uze2.transpose(0,2,1).ravel()
        
        Ix_zeta   =  xe3.transpose(1,2,0).ravel()
        Iy_zeta   =  ye3.transpose(1,2,0).ravel()
        Iz_zeta   =  ze3.transpose(1,2,0).ravel()
        
        ux_zeta =  uxe3.transpose(1,2,0).ravel()
        uy_zeta =  uye3.transpose(1,2,0).ravel()
        uz_zeta =  uze3.transpose(1,2,0).ravel()        
        
        
        """ Exporting to VTK using EVTK library """ 
        Ix  = np.r_[Ix_xi, Ix_eta, Ix_zeta]
        Iy  = np.r_[Iy_xi, Iy_eta, Iy_zeta]
        Iz  = np.r_[Iz_xi, Iz_eta, Iz_zeta]
        Iux = np.r_[ux_xi, ux_eta, ux_zeta]
        Iuy = np.r_[uy_xi, uy_eta, uy_zeta]
        Iuz = np.r_[uz_xi, uz_eta, uz_zeta]
        pointsPerLine = np.r_[np.repeat(neval[0],etau.size*zetau.size), 
                              np.repeat(neval[1],xiu.size*zetau.size),
                              np.repeat(neval[2],xiu.size*etau.size)]
                
        npoints = Ix.size
        ncells = pointsPerLine.size
        # create some temporary arrays to write grid topology
        offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells):
            ii += pointsPerLine[i]
            offsets[i] = ii
        connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive
        cell_types = np.empty(npoints, dtype = 'uint8') 
        cell_types[:] = VtkPolyLine.tid
        
        w = VtkFile(path, VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells = ncells, npoints = npoints)
        
        w.openElement("Points")
        w.addData("points", (Ix,Iy,Iz))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        
        w.openData("Point")
        w.addData("disp",(Iux,Iuy,Iuz))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (Ix,Iy,Iz) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)
        w.appendData( (Iux,Iuy,Iuz))

        w.save()
        
        
        
        
        
    
    def VtkPlot(self, path , hooke, U, neval):
        # eps = 1.e-8
        # xipRef   =  np.linspace(0+eps,1-eps,nevalperElem[0])
        # etapRef  =  np.linspace(0+eps,1-eps,nevalperElem[1])
        # zetapRef =  np.linspace(0+eps,1-eps,nevalperElem[2])
 
        # xipRefMin   = xipRef[0]
        # xipRefMax   = xipRef[-1]
        # etapRefMin  = etapRef[0]
        # etapRefMax  = etapRef[-1]
        # zetapRefMin = zetapRef[0]
        # zetapRefMax = zetapRef[-1]
          
        # e_xi   = np.unique(self.xxsi[0]) ; ne_xi   = e_xi.shape[0]-1
        # e_eta  = np.unique(self.xxsi[1]) ; ne_eta  = e_eta.shape[0]-1
        # e_zeta = np.unique(self.xxsi[2]) ; ne_zeta = e_zeta.shape[0]-1
        
        # print('Total number of visualization elements')
        # print('xi:',   ne_xi*nevalperElem[0])
        # print('eta:',  ne_eta*nevalperElem[1])
        # print('zeta:', ne_zeta*nevalperElem[2])
        
        # xi_min    =  np.kron(e_xi[:-1],np.ones(nevalperElem[0]))    
        # eta_min   =  np.kron(e_eta[:-1],np.ones(nevalperElem[1]))  
        # zeta_min  =  np.kron(e_zeta[:-1],np.ones(nevalperElem[2]))  
        
        # xi_max    =  np.kron(e_xi[1:],np.ones(nevalperElem[0]))    
        # eta_max   =  np.kron(e_eta[1:],np.ones(nevalperElem[1]))  
        # zeta_max  =  np.kron(e_zeta[1:],np.ones(nevalperElem[2]))  
        
        
        # xi_p      =  np.kron(np.ones(ne_xi)   , xipRef)
        # eta_p     =  np.kron(np.ones(ne_eta)  , etapRef)
        # zeta_p    =  np.kron(np.ones(ne_zeta) , zetapRef)

        # # Mapped visualization points in the xi parametric space 
        # xi         = (xi_min-xi_max)/(xipRefMin-xipRefMax)* xi_p + (xipRefMin*xi_max-xipRefMax*xi_min)/(xipRefMin-xipRefMax)    
        # eta        = (eta_min-eta_max)/(etapRefMin-etapRefMax)* eta_p + (etapRefMin*eta_max-etapRefMax*eta_min)/(etapRefMin-etapRefMax)    
        # zeta       = (zeta_min-zeta_max)/(zetapRefMin-zetapRefMax)* zeta_p + (zetapRefMin*zeta_max-zetapRefMax*zeta_min)/(zetapRefMin-zetapRefMax)    

 
        eps = 0
        xi   = np.linspace(self.xxsi[0][self.pp[0]]+eps  , self.xxsi[0][-self.pp[0]]-eps , neval[0])
        eta  = np.linspace(self.xxsi[1][self.pp[1]]+eps  , self.xxsi[1][-self.pp[1]]-eps , neval[1])
        zeta = np.linspace(self.xxsi[2][self.pp[2]]+eps  , self.xxsi[2][-self.pp[2]]-eps , neval[2])
        
    
        xi   = np.setdiff1d(xi,self.xxsi[0])        
        eta  = np.setdiff1d(eta,self.xxsi[1])    
        zeta = np.setdiff1d(zeta,self.xxsi[2])    
 
 
        
        nzeta = len(zeta)
        neta  = len(eta)
        nxi   = len(xi)
        print(nxi,neta,nzeta)
        
        phi, dphidx, dphidy, dphidz  = self.GetBfForGridPoints(xi, eta, zeta) 
        nbf_xi , nbf_eta, nbf_zeta = self.Get_NumberOfunibf() 
        nbf = nbf_xi * nbf_eta * nbf_zeta 
        dim = 3 
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] ; Uz = U[2*nbf:] # displacement at the control points 
        P = self.Get_Btot() # control points 
        
        """ Displacement """
        ux = phi.dot(Ux)
        uy = phi.dot(Uy)
        uz = phi.dot(Uz)
        """ Strain """ 
        exx = dphidx.dot(Ux)
        eyy = dphidy.dot(Uy) 
        ezz = dphidz.dot(Uz)
        exy = 0.5*(dphidx.dot(Uy)+dphidy.dot(Ux))
        exz = 0.5*(dphidx.dot(Uz)+dphidz.dot(Ux))
        eyz = 0.5*(dphidy.dot(Uz)+dphidz.dot(Uy))
        """ Stress """ 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(ezz) + hooke[0,3]*(2*exy) + hooke[0,4]*(2*exz) + hooke[0,5]*(2*eyz)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(ezz) + hooke[1,3]*(2*exy) + hooke[1,4]*(2*exz) + hooke[1,5]*(2*eyz)
        szz = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(ezz) + hooke[2,3]*(2*exy) + hooke[2,4]*(2*exz) + hooke[2,5]*(2*eyz)
        sxy = hooke[3,0]*(exx) + hooke[3,1]*(eyy) + hooke[3,2]*(ezz) + hooke[3,3]*(2*exy) + hooke[3,4]*(2*exz) + hooke[3,5]*(2*eyz)
        sxz = hooke[4,0]*(exx) + hooke[4,1]*(eyy) + hooke[4,2]*(ezz) + hooke[4,3]*(2*exy) + hooke[4,4]*(2*exz) + hooke[4,5]*(2*eyz)
        syz = hooke[5,0]*(exx) + hooke[5,1]*(eyy) + hooke[5,2]*(ezz) + hooke[5,3]*(2*exy) + hooke[5,4]*(2*exz) + hooke[5,5]*(2*eyz)
        
        """ Solid points """ 
        x  = phi.dot(P[0,:]) 
        y  = phi.dot(P[1,:])
        z  = phi.dot(P[2,:])
 
  
        # Reshaping the fields  
        ux = np.ascontiguousarray( ux.reshape((nzeta,neta,nxi)))
        uy = np.ascontiguousarray( uy.reshape((nzeta,neta,nxi)))
        uz = np.ascontiguousarray( uz.reshape((nzeta,neta,nxi)))
        
        exx = np.ascontiguousarray( exx.reshape((nzeta,neta,nxi)))
        eyy = np.ascontiguousarray( eyy.reshape((nzeta,neta,nxi)))
        ezz = np.ascontiguousarray( ezz.reshape((nzeta,neta,nxi)))
        exy = np.ascontiguousarray( exy.reshape((nzeta,neta,nxi)))
        exz = np.ascontiguousarray( exz.reshape((nzeta,neta,nxi)))
        eyz = np.ascontiguousarray( eyz.reshape((nzeta,neta,nxi)))
        
        sxx = np.ascontiguousarray( sxx.reshape((nzeta,neta,nxi)))
        syy = np.ascontiguousarray( syy.reshape((nzeta,neta,nxi)))
        szz = np.ascontiguousarray( szz.reshape((nzeta,neta,nxi)))
        sxy = np.ascontiguousarray( sxy.reshape((nzeta,neta,nxi)))
        sxz = np.ascontiguousarray( sxz.reshape((nzeta,neta,nxi)))
        syz = np.ascontiguousarray( syz.reshape((nzeta,neta,nxi)))
        
        x  =  np.ascontiguousarray( x.reshape((nzeta,neta,nxi)))
        y  =  np.ascontiguousarray( y.reshape((nzeta,neta,nxi)))
        z  =  np.ascontiguousarray( z.reshape((nzeta,neta,nxi)))

 
        
        """ Exporting to VTK using EVTK library """ 
 
        start = (0,0,0)
        end = (nzeta-1,neta-1,nxi-1)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("U",(ux,uy,uz))
        w.addData("Sxx", sxx) 
        w.addData("Syy", syy)
        w.addData("Szz", szz) 
        w.addData("Sxy", sxy)
        w.addData("Sxz", sxz)
        w.addData("Syz", syz) 
        w.addData("Exx", exx) 
        w.addData("Eyy", eyy)
        w.addData("Ezz", ezz) 
        w.addData("Exy", exy)
        w.addData("Exz", exz)
        w.addData("Eyz", eyz) 
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz))
        w.appendData(sxx ) 
        w.appendData(syy)
        w.appendData(szz) 
        w.appendData(sxy)
        w.appendData(sxz)
        w.appendData(syz) 
        w.appendData(exx) 
        w.appendData(eyy)
        w.appendData(ezz) 
        w.appendData(exy)
        w.appendData(exz)
        w.appendData(eyz) 
 

        w.save()   
        
   
 
 
 
    def VtkVolumeAndControlGrid(self,path,neval):
        P = self.Get_Btot() # control points 
        nbf_xi , nbf_eta, nbf_zeta = self.Get_NumberOfunibf() 
        nbf = nbf_xi * nbf_eta * nbf_zeta                 
        
        """ Exporting the control points cloud """ 
        xp = np.ascontiguousarray(P[0,:])
        yp = np.ascontiguousarray(P[1,:])
        zp = np.ascontiguousarray(P[2,:])
        tools.tvtk.pointsToVTK(path+'-control-points', xp, yp, zp, data= {"point": np.zeros(nbf)} ) 
        
        """ Exporting the control grid """
        # Defining lines for each iso-grid 
 
        # Zeta iso-grids 
        Ix_zeta  = np.array([])
        Iy_zeta  = np.array([])
        Iz_zeta  = np.array([])        
 
        Xt = self.X 
        Yt = self.Y 
        Zt = self.Z 
        for i in range(nbf_zeta):
            Ixs = np.r_[Xt[i,:,:].ravel(order='F'),Xt[i,:,:].ravel(order='C')]  
            Iys = np.r_[Yt[i,:,:].ravel(order='F'),Yt[i,:,:].ravel(order='C')]  
            Izs = np.r_[Zt[i,:,:].ravel(order='F'),Zt[i,:,:].ravel(order='C')]  
            Ix_zeta = np.concatenate((Ix_zeta,Ixs))
            Iy_zeta = np.concatenate((Iy_zeta,Iys))
            Iz_zeta = np.concatenate((Iz_zeta,Izs)) 

        # Eta iso-grids 
        Ix_eta  = np.array([])
        Iy_eta  = np.array([])
        Iz_eta  = np.array([])        
 
        Xt = self.X 
        Yt = self.Y 
        Zt = self.Z 
        for i in range(nbf_eta):
            Ixs = np.r_[Xt[:,:,i].ravel(order='F'),Xt[:,:,i].ravel(order='C')]  
            Iys = np.r_[Yt[:,:,i].ravel(order='F'),Yt[:,:,i].ravel(order='C')]  
            Izs = np.r_[Zt[:,:,i].ravel(order='F'),Zt[:,:,i].ravel(order='C')]  
            Ix_eta = np.concatenate((Ix_eta,Ixs))
            Iy_eta = np.concatenate((Iy_eta,Iys))
            Iz_eta = np.concatenate((Iz_eta,Izs)) 
            
 
        pointsPerLine_zeta = np.kron(np.ones(nbf_zeta), np.r_[np.repeat(nbf_xi,nbf_eta), np.repeat(nbf_eta,nbf_xi)]  )
        pointsPerLine_eta  = np.kron(np.ones(nbf_eta),  np.r_[np.repeat(nbf_zeta,nbf_xi), np.repeat(nbf_xi,nbf_zeta)] )  
        
 
          
        npoints_zeta = Ix_zeta.size
        ncells_zeta = pointsPerLine_zeta.size
        # create some temporary arrays to write grid topology
        offsets_zeta = np.zeros(ncells_zeta, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells_zeta):
            ii += pointsPerLine_zeta[i]
            offsets_zeta[i] = ii
        connectivity_zeta = np.arange(npoints_zeta, dtype = 'int32')      # each line connects points that are consecutive
        cell_types_zeta = np.empty(npoints_zeta, dtype = 'uint8') 
        cell_types_zeta[:] = VtkPolyLine.tid
        
        
        npoints_eta = Ix_eta.size
        ncells_eta = pointsPerLine_eta.size
        # create some temporary arrays to write grid topology
        offsets_eta = np.zeros(ncells_eta, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells_eta):
            ii += pointsPerLine_eta[i]
            offsets_eta[i] = ii
        connectivity_eta = np.arange(npoints_eta, dtype = 'int32')      # each line connects points that are consecutive
        cell_types_eta = np.empty(npoints_eta, dtype = 'uint8') 
        cell_types_eta[:] = VtkPolyLine.tid
        
        
        
        w = VtkFile(path+'-control-grid', VtkUnstructuredGrid)
        
        w.openGrid()
        
        w.openPiece(ncells = ncells_zeta, npoints = npoints_zeta)
        w.openElement("Points")
        w.addData("points", (Ix_zeta,Iy_zeta,Iz_zeta))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity_zeta)
        w.addData("offsets", offsets_zeta)
        w.addData("types", cell_types_zeta)
        w.closeElement("Cells")
        w.closePiece()
        
        
        w.openPiece(ncells = ncells_eta, npoints = npoints_eta)
        w.openElement("Points")
        w.addData("points", (Ix_eta,Iy_eta,Iz_eta))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity_eta)
        w.addData("offsets", offsets_eta)
        w.addData("types", cell_types_eta)
        w.closeElement("Cells")
        w.closePiece()
        
        w.closeGrid()
        
 
        w.appendData( (Ix_zeta,Iy_zeta,Iz_zeta) )
        w.appendData(connectivity_zeta).appendData(offsets_zeta).appendData(cell_types_zeta)        

        w.appendData( (Ix_eta,Iy_eta,Iz_eta) )
        w.appendData(connectivity_eta).appendData(offsets_eta).appendData(cell_types_eta)      
 
        w.save()
        
        
        
        """Exporting the volume """  
        eps = 0
        xi   = np.linspace(self.xxsi[0][self.pp[0]]+eps  , self.xxsi[0][-self.pp[0]]-eps , neval[0])
        eta  = np.linspace(self.xxsi[1][self.pp[1]]+eps  , self.xxsi[1][-self.pp[1]]-eps , neval[1])
        zeta = np.linspace(self.xxsi[2][self.pp[2]]+eps  , self.xxsi[2][-self.pp[2]]-eps , neval[2])
        
 
        
        nzeta = len(zeta)
        neta  = len(eta)
        nxi   = len(xi)
 
        
        phi, dphidx, dphidy, dphidz  = self.GetBfForGridPoints(xi, eta, zeta) 
        nbf_xi , nbf_eta, nbf_zeta = self.Get_NumberOfunibf() 
        nbf = nbf_xi * nbf_eta * nbf_zeta 
 
        P = self.Get_Btot() # control points 
 
        x  = phi.dot(P[0,:]) 
        y  = phi.dot(P[1,:])
        z  = phi.dot(P[2,:])
       
        x  =  np.ascontiguousarray( x.reshape((nzeta,neta,nxi)))
        y  =  np.ascontiguousarray( y.reshape((nzeta,neta,nxi)))
        z  =  np.ascontiguousarray( z.reshape((nzeta,neta,nxi)))
 
 
        start = (0,0,0)
        end = (nzeta-1,neta-1,nxi-1)
        
        w =  VtkFile(path+'-volume', VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
 
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
 
        w.save()           

 
        
    
    
    def Stiffness(self,hooke):
        """ 
        Stiffness Matrix """
        
        if isinstance(hooke[0,0],sps.dia.dia_matrix) :
            """
            If Hook's law depends on the spatial position of each integration point 
            """
            Bxy = self.dphixdy + self.dphiydx
            Bxz = self.dphixdz + self.dphizdx
            Byz = self.dphiydz + self.dphizdy 
        
            K =  self.dphixdx.T.dot(hooke[0,0].dot(self.wg.dot(self.dphixdx))) + \
                 self.dphiydy.T.dot(hooke[0,1].dot(self.wg.dot(self.dphixdx))) + \
                 self.dphizdz.T.dot(hooke[0,2].dot(self.wg.dot(self.dphixdx))) + \
                 Bxy.T.dot(hooke[0,3].dot(self.wg.dot(self.dphixdx)))          + \
                 Bxz.T.dot(hooke[0,4].dot(self.wg.dot(self.dphixdx)))          + \
                 Byz.T.dot(hooke[0,5].dot(self.wg.dot(self.dphixdx)))          + \
                 self.dphixdx.T.dot(hooke[1,0].dot(self.wg.dot(self.dphiydy))) + \
                 self.dphiydy.T.dot(hooke[1,1].dot(self.wg.dot(self.dphiydy))) + \
                 self.dphizdz.T.dot(hooke[1,2].dot(self.wg.dot(self.dphiydy))) + \
                 Bxy.T.dot(hooke[1,3].dot(self.wg.dot(self.dphiydy)))          + \
                 Bxz.T.dot(hooke[1,4].dot(self.wg.dot(self.dphiydy)))          + \
                 Byz.T.dot(hooke[1,5].dot(self.wg.dot(self.dphiydy)))          + \
                 self.dphixdx.T.dot(hooke[2,0].dot(self.wg.dot(self.dphizdz))) + \
                 self.dphiydy.T.dot(hooke[2,1].dot(self.wg.dot(self.dphizdz))) + \
                 self.dphizdz.T.dot(hooke[2,2].dot(self.wg.dot(self.dphizdz))) + \
                 Bxy.T.dot(hooke[2,3].dot(self.wg.dot(self.dphizdz)))          + \
                 Bxz.T.dot(hooke[2,4].dot(self.wg.dot(self.dphizdz)))          + \
                 Byz.T.dot(hooke[2,5].dot(self.wg.dot(self.dphizdz)))          + \
                 self.dphixdx.T.dot(hooke[3,0].dot(self.wg.dot(Bxy)))          + \
                 self.dphiydy.T.dot(hooke[3,1].dot(self.wg.dot(Bxy)))          + \
                 self.dphizdz.T.dot(hooke[3,2].dot(self.wg.dot(Bxy)))          + \
                 Bxy.T.dot(hooke[3,3].dot(self.wg.dot(Bxy)))                   + \
                 Bxz.T.dot(hooke[3,4].dot(self.wg.dot(Bxy)))                   + \
                 Byz.T.dot(hooke[3,5].dot(self.wg.dot(Bxy)))                   + \
                 self.dphixdx.T.dot(hooke[4,0].dot(self.wg.dot(Bxz)))          + \
                 self.dphiydy.T.dot(hooke[4,1].dot(self.wg.dot(Bxz)))          + \
                 self.dphizdz.T.dot(hooke[4,2].dot(self.wg.dot(Bxz)))          + \
                 Bxy.T.dot(hooke[4,3].dot(self.wg.dot(Bxz)))                   + \
                 Bxz.T.dot(hooke[4,4].dot(self.wg.dot(Bxz)))                   + \
                 Byz.T.dot(hooke[4,5].dot(self.wg.dot(Bxz)))                   + \
                 self.dphixdx.T.dot(hooke[5,0].dot(self.wg.dot(Byz)))          + \
                 self.dphiydy.T.dot(hooke[5,1].dot(self.wg.dot(Byz)))          + \
                 self.dphizdz.T.dot(hooke[5,2].dot(self.wg.dot(Byz)))          + \
                 Bxy.T.dot(hooke[5,3].dot(self.wg.dot(Byz)))                   + \
                 Bxz.T.dot(hooke[5,4].dot(self.wg.dot(Byz)))                   + \
                 Byz.T.dot(hooke[5,5].dot(self.wg.dot(Byz))) 
                 
                 
          
        else: 
            
            Bxy = self.dphixdy + self.dphiydx
            Bxz = self.dphixdz + self.dphizdx
            Byz = self.dphiydz + self.dphizdy 
        
            K =  hooke[0,0]*self.dphixdx.T.dot(self.wg.dot(self.dphixdx)) + \
                 hooke[0,1]*self.dphiydy.T.dot(self.wg.dot(self.dphixdx)) + \
                 hooke[0,2]*self.dphizdz.T.dot(self.wg.dot(self.dphixdx)) + \
                 hooke[0,3]*Bxy.T.dot(self.wg.dot(self.dphixdx))          + \
                 hooke[0,4]*Bxz.T.dot(self.wg.dot(self.dphixdx))          + \
                 hooke[0,5]*Byz.T.dot(self.wg.dot(self.dphixdx))          + \
                 hooke[1,0]*self.dphixdx.T.dot(self.wg.dot(self.dphiydy)) + \
                 hooke[1,1]*self.dphiydy.T.dot(self.wg.dot(self.dphiydy)) + \
                 hooke[1,2]*self.dphizdz.T.dot(self.wg.dot(self.dphiydy)) + \
                 hooke[1,3]*Bxy.T.dot(self.wg.dot(self.dphiydy))          + \
                 hooke[1,4]*Bxz.T.dot(self.wg.dot(self.dphiydy))          + \
                 hooke[1,5]*Byz.T.dot(self.wg.dot(self.dphiydy))          + \
                 hooke[2,0]*self.dphixdx.T.dot(self.wg.dot(self.dphizdz)) + \
                 hooke[2,1]*self.dphiydy.T.dot(self.wg.dot(self.dphizdz)) + \
                 hooke[2,2]*self.dphizdz.T.dot(self.wg.dot(self.dphizdz)) + \
                 hooke[2,3]*Bxy.T.dot(self.wg.dot(self.dphizdz))          + \
                 hooke[2,4]*Bxz.T.dot(self.wg.dot(self.dphizdz))          + \
                 hooke[2,5]*Byz.T.dot(self.wg.dot(self.dphizdz))          + \
                 hooke[3,0]*self.dphixdx.T.dot(self.wg.dot(Bxy))          + \
                 hooke[3,1]*self.dphiydy.T.dot(self.wg.dot(Bxy))          + \
                 hooke[3,2]*self.dphizdz.T.dot(self.wg.dot(Bxy))          + \
                 hooke[3,3]*Bxy.T.dot(self.wg.dot(Bxy))                   + \
                 hooke[3,4]*Bxz.T.dot(self.wg.dot(Bxy))                   + \
                 hooke[3,5]*Byz.T.dot(self.wg.dot(Bxy))                   + \
                 hooke[4,0]*self.dphixdx.T.dot(self.wg.dot(Bxz))          + \
                 hooke[4,1]*self.dphiydy.T.dot(self.wg.dot(Bxz))          + \
                 hooke[4,2]*self.dphizdz.T.dot(self.wg.dot(Bxz))          + \
                 hooke[4,3]*Bxy.T.dot(self.wg.dot(Bxz))                   + \
                 hooke[4,4]*Bxz.T.dot(self.wg.dot(Bxz))                   + \
                 hooke[4,5]*Byz.T.dot(self.wg.dot(Bxz))                   + \
                 hooke[5,0]*self.dphixdx.T.dot(self.wg.dot(Byz))          + \
                 hooke[5,1]*self.dphiydy.T.dot(self.wg.dot(Byz))          + \
                 hooke[5,2]*self.dphizdz.T.dot(self.wg.dot(Byz))          + \
                 hooke[5,3]*Bxy.T.dot(self.wg.dot(Byz))                   + \
                 hooke[5,4]*Bxz.T.dot(self.wg.dot(Byz))                   + \
                 hooke[5,5]*Byz.T.dot(self.wg.dot(Byz)) 
            
        return K     
    
#%% 2D CODE 
#%%% Finite cell 
class CellTesselation :
    def __init__(self,xmin,xmax,ymin,ymax):
        self.xmin = xmin 
        self.xmax = xmax
        self.ymin = ymin 
        self.ymax = ymax 
        self.mesx = xmax - xmin 
        self.mesy = ymax - ymin 
        self.lvl = 0 
        self.topRight = None  
        self.topLeft = None 
        self.bottomRight = None 
        self.bottomLeft = None 

    def DecomposeLevelSetIntegration(self,mesh,Phi,lvlmax):
        """ Tesselation based FCM integration """
        """ Getting the different cells and triangles  """ 
        if self.lvl == lvlmax : 
            # First checking if the four corners are homogeneous 
            # in this case, we perform an integration on a the square cell 
            xc =  np.array([self.xmin+1.e-8 , self.xmax-1.e-8, self.xmax-1.e-8, self.xmin+1.e-8 ])
            yc = np.array([self.ymin +1.e-8, self.ymin+1.e-8 , self.ymax-1.e-8, self.ymax-1.e-8 ])
            ls = Phi(xc,yc)  
            if ls[0] >= 0  and  ls[1] >= 0 and ls[2] >= 0 and ls[3] >= 0 :
                xim  = np.kron(self.xmin,np.ones(mesh.nbg_xi))
                etam = np.kron(self.ymin,np.ones(mesh.nbg_eta))
                xi        = xim  + 0.5*(mesh.Gauss_xi[0]+1)*self.mesx        
                eta       = etam + 0.5*(mesh.Gauss_eta[0]+1)*self.mesy
                xi        = np.kron(np.ones(mesh.nbg_eta), xi)
                eta       = np.kron(eta, np.ones(mesh.nbg_xi))                
                mesh.pix.extend(list(xi))
                mesh.piy.extend(list(eta))
                mesh.wg.extend(mesh.wgRef*self.mesx*self.mesy/4 )                
                # mesh.integrationCellsCoord.append([ 'r', self.xmin,self.xmax, self.ymin, self.ymax ])
                
            elif (ls[0] < 0 and ls[1] < 0 and ls[2] < 0 and ls[3] < 0) == False :
                # perform the tesselation procedure using delaunay algorithm 
                TriCorners = [] 
                if ls[0] > 0 :
                     TriCorners.append( np.array([ xc[0],yc[0] ]) )
                if ls[1] > 0 :
                    TriCorners.append( np.array([ xc[1],yc[1] ]) )
                if ls[2] > 0 :
                    TriCorners.append( np.array([ xc[2],yc[2] ]) )
                if ls[3] > 0 :
                    TriCorners.append( np.array([ xc[3],yc[3] ]) )
        
                # bottom segment  
                if ( ls[0] > 0 and ls[1] < 0 ) or ( ls[0] < 0 and ls[1] > 0) :
                    a,b = nb.interpolateLinearly(xc[0],xc[1],ls[0],ls[1])
                    xb = -b/a 
                    yb =  self.ymin 
                    TriCorners.append( np.array([xb,yb]))
                # right segment  
                if ( ls[1] > 0 and ls[2] < 0 ) or ( ls[1] < 0 and ls[2] > 0) :
                    a,b = nb.interpolateLinearly(yc[1],yc[2],ls[1],ls[2])
                    xr =  self.xmax
                    yr =  -b/a 
                    TriCorners.append( np.array([xr,yr]))
                # top segment 
                if ( ls[2] > 0 and ls[3] < 0 ) or ( ls[2] < 0 and ls[3] > 0) :
                    a,b = nb.interpolateLinearly(xc[2],xc[3],ls[2],ls[3])
                    xt =  -b/a 
                    yt =  self.ymax
                    TriCorners.append( np.array([xt,yt]) )
                # left segment 
                if ( ls[3] > 0 and ls[0] < 0 ) or ( ls[3] < 0 and ls[0] > 0) :
                    a,b = nb.interpolateLinearly(yc[3],yc[0],ls[3],ls[0])
                    xl =  self.xmin  
                    yl =  -b/a 
                    TriCorners.append( np.array([xl,yl]))
                nodes = np.array([TriCorners])[0] 
#                print (TriCorners)
#                print(self.xmin,self.xmax,self.ymin,self.ymax)
                # Now we give the tesselation nodes to the delaunay algorithm in order 
                # to obtain a tesselation at the boundary 
                tri = Delaunay(nodes)  
                # mesh.integrationCellsCoord.append(['t',tri.points,tri.simplices]) 
                for i in range(tri.simplices.shape[0]):
                    n1 = tri.simplices[i,0] ; n2 = tri.simplices[i,1]  ; n3=  tri.simplices[i,2]
                    x1 = tri.points[n1,0] ; y1 = tri.points[n1,1] 
                    x2 = tri.points[n2,0] ; y2 = tri.points[n2,1] 
                    x3 = tri.points[n3,0] ; y3 = tri.points[n3,1]
                    
                    xgTref = mesh.refGaussTriangle[0][:,0]
                    ygTref = mesh.refGaussTriangle[0][:,1]
                    wgTref = mesh.refGaussTriangle[1]
                    # Mapping from the reference triangle to the physical triangle 
                    xi  = (x2-x1)*xgTref + (x3-x1)*ygTref + x1 
                    eta = (y2-y1)*xgTref + (y3-y1)*ygTref + y1
                    # Getting the area of the triangle 
                    areaT = 0.5*np.abs( (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
                    mesh.pix.extend(list(xi))
                    mesh.piy.extend(list(eta))
                    mesh.wg.extend(wgTref*2*areaT ) 
        else : 
            nbp = 2**(lvlmax - self.lvl)+1 
            xi  = np.linspace(self.xmin +1.e-8, self.xmax -1.e-8, nbp)
            eta = np.linspace(self.ymin +1.e-8, self.ymax -1.e-8, nbp)
            x = np.r_[xi,xi,np.kron(xi[0],np.ones(nbp)),np.kron(xi[-1],np.ones(nbp))]
            y = np.r_[np.kron(eta[0],np.ones(nbp)),np.kron(eta[-1],np.ones(nbp)),eta,eta]
            ls = Phi(x,y)  
            if (ls>0).all():
                # integrate the cell 
                # mesh.integrationCellsCoord.append([ 'r', self.xmin,self.xmax, self.ymin, self.ymax ])
                xim  = np.kron(self.xmin,np.ones(mesh.nbg_xi))
                etam = np.kron(self.ymin,np.ones(mesh.nbg_eta))
                xi        = xim  + 0.5*(mesh.Gauss_xi[0]+1)*self.mesx        
                eta       = etam + 0.5*(mesh.Gauss_eta[0]+1)*self.mesy
                xi        = np.kron(np.ones(mesh.nbg_eta), xi)
                eta       = np.kron(eta, np.ones(mesh.nbg_xi))                
                mesh.pix.extend(list(xi))
                mesh.piy.extend(list(eta))
                mesh.wg.extend(mesh.wgRef*self.mesx*self.mesy/4 )    
            elif (ls<0).all() == False : 
                # cut the cell : recursive call
                self.topRight    = CellTesselation((self.xmin+self.xmax)/2, self.xmax,(self.ymin+self.ymax)/2, self.ymax)
                self.topLeft     = CellTesselation(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax)
                self.bottomRight = CellTesselation((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2)
                self.bottomLeft  = CellTesselation(self.xmin, (self.xmin+self.xmax)/2,self.ymin,(self.ymin+self.ymax)/2)
                
                self.topRight.lvl = self.lvl + 1 
                self.topLeft.lvl  = self.lvl + 1 
                self.bottomRight.lvl  = self.lvl + 1 
                self.bottomLeft.lvl   = self.lvl + 1 
                
                self.topRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.topLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.bottomRight.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)
                self.bottomLeft.DecomposeLevelSetIntegration(mesh,Phi,lvlmax)     
        
    def DecomposeLevelSetCoord(self,mesh,Phi,lvlmax):
        """ Getting the different cells and triangles  """ 
        if self.lvl == lvlmax : 
            # First checking if the four corners are homogeneous 
            # in this case, we perform an integration on a the square cell 
            xc =  np.array([self.xmin+1.e-8 , self.xmax-1.e-8, self.xmax-1.e-8, self.xmin+1.e-8 ])
            yc = np.array([self.ymin +1.e-8, self.ymin+1.e-8 , self.ymax-1.e-8, self.ymax-1.e-8 ])
            ls = Phi(xc,yc)  
            if ls[0] > 0  and  ls[1] > 0 and ls[2] > 0 and ls[3] > 0 :
                mesh.integrationCellsCoord.append([ 'r', self.xmin,self.xmax, self.ymin, self.ymax ]) 
            elif (ls[0] < 0 and ls[1] < 0 and ls[2] < 0 and ls[3] < 0) == False :
                # perform the tesselation procedure using delaunay algorithm 
                TriCorners = [] 
                if ls[0] > 0 :
                     TriCorners.append( np.array([ xc[0],yc[0] ]) )
                if ls[1] > 0 :
                    TriCorners.append( np.array([ xc[1],yc[1] ]) )
                if ls[2] > 0 :
                    TriCorners.append( np.array([ xc[2],yc[2] ]) )
                if ls[3] > 0 :
                    TriCorners.append( np.array([ xc[3],yc[3] ]) )
        
                # bottom segment  
                if ( ls[0] > 0 and ls[1] < 0 ) or ( ls[0] < 0 and ls[1] > 0) :
                    a,b = nb.interpolateLinearly(xc[0],xc[1],ls[0],ls[1])
                    xb =  -b/a 
                    yb =  self.ymin 
                    TriCorners.append( np.array([xb,yb]))
                # right segment  
                if ( ls[1] > 0 and ls[2] < 0 ) or ( ls[1] < 0 and ls[2] > 0) :
                    a,b = nb.interpolateLinearly(yc[1],yc[2],ls[1],ls[2])
                    xr =  self.xmax
                    yr =   -b/a 
                    TriCorners.append( np.array([xr,yr]))
                # top segment 
                if ( ls[2] > 0 and ls[3] < 0 ) or ( ls[2] < 0 and ls[3] > 0) :
                    a,b = nb.interpolateLinearly(xc[2],xc[3],ls[2],ls[3])
                    xt =  -b/a 
                    yt =  self.ymax
                    TriCorners.append( np.array([xt,yt]) )
                # left segment 
                if ( ls[3] > 0 and ls[0] < 0 ) or ( ls[3] < 0 and ls[0] > 0) :
                    a,b = nb.interpolateLinearly(yc[3],yc[0],ls[3],ls[0])
                    xl =  self.xmin  
                    yl =   -b/a  
                    TriCorners.append( np.array([xl,yl]))
                nodes = np.array([TriCorners])[0] 
#                print (TriCorners)
#                print(self.xmin,self.xmax,self.ymin,self.ymax)
                # Now we give the tesselation nodes to the delaunay algorithm in order 
                # to obtain a tesselation at the boundary   
                tri = Delaunay(nodes)  
                mesh.integrationCellsCoord.append(['t',tri.points,tri.simplices])    
        else : 
            nbp = 2**(lvlmax - self.lvl)+1 
            xi  = np.linspace(self.xmin +1.e-8, self.xmax -1.e-8, nbp)
            eta = np.linspace(self.ymin +1.e-8, self.ymax -1.e-8, nbp)
            x = np.r_[xi,xi,np.kron(xi[0],np.ones(nbp)),np.kron(xi[-1],np.ones(nbp))]
            y = np.r_[np.kron(eta[0],np.ones(nbp)),np.kron(eta[-1],np.ones(nbp)),eta,eta]
            ls = Phi(x,y)  
            if (ls>=0).all():
                # integrate the cell 
                mesh.integrationCellsCoord.append([ 'r', self.xmin,self.xmax, self.ymin, self.ymax ]) 
            elif (ls<0).all() == False : 
                # cut the cell : recursive call
                self.topRight    = CellTesselation((self.xmin+self.xmax)/2, self.xmax,(self.ymin+self.ymax)/2, self.ymax)
                self.topLeft     = CellTesselation(self.xmin, (self.xmin+self.xmax)/2, (self.ymin+self.ymax)/2, self.ymax)
                self.bottomRight = CellTesselation((self.xmin+self.xmax)/2, self.xmax, self.ymin, (self.ymin+self.ymax)/2)
                self.bottomLeft  = CellTesselation(self.xmin, (self.xmin+self.xmax)/2,self.ymin,(self.ymin+self.ymax)/2)
                
                self.topRight.lvl = self.lvl + 1 
                self.topLeft.lvl  = self.lvl + 1 
                self.bottomRight.lvl  = self.lvl + 1 
                self.bottomLeft.lvl   = self.lvl + 1 
                
                self.topRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.topLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.bottomRight.DecomposeLevelSetCoord(mesh,Phi,lvlmax)
                self.bottomLeft.DecomposeLevelSetCoord(mesh,Phi,lvlmax)    
                
#%%% IGA  
class Mesh:
    """
    Mesh class : IGA  
    """
    def __init__(self,ctrlPts,pp,xxsi):
        self.ctrlPts=ctrlPts    # Inhomogeneous Control Points Coordinates and weights
        self.pp=pp              # degree [ppu,ppv]
        self.xxsi=xxsi          # knot vector p-uplet (XXsiu,XXsiv)
        self.ien=0              # NURBS Connectivity p-uplet (IENu,IENv)
        self.noelem=0           # Connectivity: Control Point number (column) of each element (line)
        self.tripleien=0        # Correspondance between 2D elements and 1D elements
        self.iperM=0            # Sparse matrix for coincident control points
        
        self.wdetJmes = 0     # Dictionary containing the value det(Jacobian)*Gauss_Weights*ElementMeasure/4 
        self.mes = 0         # Dictionary containing the measure of the elements in the ordre of listeltot
        
        self.phi = 0         # Dictionary of basis functions evaluated at gauss points 
        self.dphidx = 0      # Dictionary containing the derivative of Basis function in x direction
        self.dphidy = 0      # Dictionary containing the derivative of Basis function in y direction
        
        
        """ Attributes when using vectorization  """ 
        """ In this case, the implicit connectivity of the structured B-spline parametric space is used """ 
        self.npg = 0 
        self.phix = 0      # Matrix (N,0)
        self.phiy = 0      # Matrix (0,N)
        self.dphixdx =  0  # Matrix (dNdx,0)
        self.dphixdy =  0  # Matrix (dNdy,0)
        self.dphiydx =  0  # Matrix (0,dNdx)
        self.dphiydy =  0  # Matrix (0,dNdy) 
        self.wg = 0        # Integration weights diagonal matrix 
        
         
        self.phiMatrix = 0 
        self.n_elems = 0 
        self.pix = 0
        self.piy = 0  
        self.piz = 0 
        self.integrationCellsCoord = 0 
        """ fixed parameters for integration """
        self.nbg_xi  =  0
        self.nbg_eta  = 0
        self.Gauss_xi  =  0
        self.Gauss_eta =  0
        self.wgRef =  0
        self.refGaussTriangle =  0
       
       
        

    def Get_nunv(self):
        """ Number of control points in u and v direction  """
        return self.ctrlPts.shape[1:]
    def Get_Btot(self):
        """List of control points """
        P=np.c_[self.ctrlPts[0].ravel(order='F'),self.ctrlPts[1].ravel(order='F'),self.ctrlPts[2].ravel(order='F')]            
        if self.ctrlPts.shape[0]==4:
            P=np.c_[P,self.ctrlPts[3].ravel(order='F')]
        return P.T
    def Get_listeltot(self):
        """ Indices of elements """
        return np.arange(self.ien[0].shape[1]*self.ien[1].shape[1])
    def Get_listBtot(self):
        """ Indices of control points """
        return np.arange(self.Get_Btot().shape[1])
    def Get_nenunv(self):
        """ Number of basis functions for one element in u and v directions """
        return self.pp+1
    def Get_dim(self):
        """ Dimension of the control points """
        return self.ctrlPts.shape[0] - 1
    def Get_ndof(self):
        """ Number of degrees of freedom"""
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        return self.Get_dim()*nbf 
    def Get_nbf(self):
        return self.Get_Btot().shape[1]
    def IsRational(self):
        """ Returns True if the B-spline is rational"""
        if self.ctrlPts.shape[0] == 3:
            return (self.ctrlPts[2]!=1).any()
        if self.ctrlPts.shape[0] == 4:
            return (self.ctrlPts[3]!=1).any()
        
        
    def DegElevation(self,ppnew):
        ''' INPUT target degree, example: ppnew=[2,3]
        '''
        t = ppnew - self.pp
        # to homogeneous coordinates
        self.ctrlPts[:-1,:,:] = self.ctrlPts[:-1,:,:]*self.ctrlPts[-1,:,:]
        # NURBS represents a surface
        dim = self.Get_dim()        # dimension de l'espace physique (3 si on utilise x,y)
 
        num1 = self.ctrlPts.shape[1]       # nombre de points de controle selon u (nu)
        num2 = self.ctrlPts.shape[2]       # nombre de points de controle selon v (nv)
        # Degree elevate along the v direction
        if t[1] != 0:
            coefs = self.ctrlPts.reshape(((dim+1)*num1,num2),order="F")  # fait un tableau 2D pour pouvoir utiliser bspdegelev
            self.ctrlPts,self.xxsi[1] = nb.bspdegelev(self.pp[1],coefs,self.xxsi[1],t[1])
            num2 = self.ctrlPts.shape[1]                             # nouvelle valeur de nv (on a ajouté des pts de controle)
            self.ctrlPts = self.ctrlPts.reshape(((dim+1),num1,num2),order="F") # on remet sous la forme de matrice 3D
        # Degree elevate along the u direction
        if t[0] != 0:
            coefs = np.transpose(self.ctrlPts,(0,2,1))   # pour pouvoir utiliser bspdegelev, on met la dimension sur laquelle on travaille (u) au bon endroit
            coefs = coefs.reshape(((dim+1)*num2,num1),order="F")
            self.ctrlPts,self.xxsi[0] = nb.bspdegelev(self.pp[0],coefs,self.xxsi[0],t[0])
            num1 = self.ctrlPts.shape[1]
            self.ctrlPts = self.ctrlPts.reshape(((dim+1),num2,num1),order="F")
            self.ctrlPts = np.transpose(self.ctrlPts,(0,2,1)) # on remet les dimensions dans le bon ordre
 
        self.pp=ppnew
        # back to real coordinates
        self.ctrlPts[:-1,:,:] = self.ctrlPts[:-1,:,:]/self.ctrlPts[-1,:,:]
        
    def KnotInsertion(self,u_refinement):
        '''
        INPUT
        ubar=dict()
        ubar[0] = ubaru
        ubar[1] = ubarv
        '''
        # to homogeneous coordinates
        self.ctrlPts[:-1,:,:] = self.ctrlPts[:-1,:,:]*self.ctrlPts[-1,:,:]
        # NURBS represents a surface
        dim = self.Get_dim()
        num1 = self.ctrlPts.shape[1] # nombre de points de controle selon u (nu)
        num2 = self.ctrlPts.shape[2] # nombre de points de controle selon v (nv)
        nu = np.size(u_refinement[0])
        nv = np.size(u_refinement[1])
        # Degree elevate along the v direction
        if nv != 0:
            coefs = self.ctrlPts.reshape(((dim+1)*num1,num2),order="F") # fait un tableau 2D pour pouvoir utiliser bspdegelev
            self.ctrlPts,self.xxsi[1] = nb.bspkntins(self.pp[1],coefs,self.xxsi[1],u_refinement[1])
            num2 = self.ctrlPts.shape[1] # nouvelle valeur de nv (on a ajouté des pts de controle)
            self.ctrlPts = self.ctrlPts.reshape(((dim+1),num1,num2),order="F") # on remet sous la forme de matrice 3D
        # Degree elevate along the u direction
        if nu != 0:
            coefs = np.transpose(self.ctrlPts,(0,2,1)) # pour pouvoir utiliser bspdegelev, on met la dimension sur laquelle on travaille (u) au bon endroit
            coefs = coefs.reshape(((dim+1)*num2,num1),order="F")
            self.ctrlPts, self.xxsi[0] = nb.bspkntins(self.pp[0], coefs, self.xxsi[0], u_refinement[0])
            num1 = self.ctrlPts.shape[1]
            self.ctrlPts = self.ctrlPts.reshape(((dim+1),num2,num1),order="F")
            self.ctrlPts = np.transpose(self.ctrlPts,(0,2,1)) # on remets les dimensions dans le bon ordre
            
        # back to real coordinates
        self.ctrlPts[:-1,:,:] = self.ctrlPts[:-1,:,:]/self.ctrlPts[-1,:,:]

    def Connectivity(self):
        dim=self.Get_dim()
        self.ien=dict()
        for i in range(dim):
            self.ien[i] = nb.nubsconnect(self.pp[i],self.Get_nunv()[i])
 
        ''' connect_2D(IENu,IENv,nu) '''
        nu=self.Get_nunv()[0]
        # 4 boucles for : marche avec input : nelx,nely,nenx,neny,nx   ### Attention : ancienne version pas avec indices Python qui partent de 0
        #    NOELEM = np.zeros((nely*nelx,neny*nenx),dtype='int')
        #    for j in range(nely):   #1:nely
        #        for i in range(nelx):   #1:nelx
        #            for jj in range(neny):   #1:neny
        #                for ii in range (nenx):  #1:nenx
        #                    NOELEM[i+j*nelx,ii+jj*nenx]=ii+1+nx*jj+i+j*nx
        decal = np.fliplr(self.ien[1].T)*nu
        comptage = np.fliplr(self.ien[0].T)
        tramev = np.ones_like(decal)
        trameu = np.ones_like(comptage)
        self.noelem = np.kron(tramev,comptage)+np.kron(decal,trameu)
        ''' IEN_2D(nelx,nely,listelx,listely) '''
        #    TRIPLEIEN = np.zeros((nelx*nely,2))    ### Attention : ancienne version pas avec indices Python qui partent de 0
        #    for j in range(nely):   #1:nely
        #        for i in range(nelx):   #1:nelx
        #            TRIPLEIEN[i+nelx*j,0] = listelx[i]
        #            TRIPLEIEN[i+nelx*j,1] = listely[j]    
        #    return TRIPLEIEN
        listelx=np.arange(self.ien[0].shape[1])
        listely=np.arange(self.ien[1].shape[1])
        xx,yy = np.meshgrid(listelx,listely)
        self.tripleien = np.array([xx.ravel(),yy.ravel()]).T
        '''IPERm'''
        Btot=self.Get_Btot()
        MasterPts,indices = np.unique(Btot[:-1,:].T,axis=0,return_inverse = True)
        IPERm  = sps.csc_matrix((np.ones(len(indices)), (np.arange(len(indices)),indices)), shape = (len(indices),MasterPts.shape[0] ))
        self.iperM = sps.kron(np.eye(2),IPERm) # 2 ddl par point de contrôle
 
            
    def GetBoundaryIndices(self):
        P = self.Get_Btot() 
        left   =   np.where( np.abs ( P[0,:]  - np.min(P[0,:])  )  <= 1.e-12  )[0]
        right  =   np.where( np.abs ( P[0,:]  - np.max(P[0,:])  )  <= 1.e-12  )[0]
        bottom =   np.where( np.abs ( P[1,:]  - np.min(P[1,:])  )  <= 1.e-12  )[0]
        top    =   np.where( np.abs ( P[1,:]  - np.max(P[1,:])  )  <= 1.e-12  )[0]
        sleft = np.size(left)
        stop = np.size(top)
        sbottom = np.size(bottom)
        sright = np.size(right)
        return bottom,top,right,left,sbottom,stop,sright,sleft 
 
    def ComputeUnivariateIntegratedDiffTensors(self, nip, ne):
        # nip : [nipx,nipy] number of integration points per integration element 
        # ne : [nex,ney] number of integration elements per knot element (knot span) 
        
        e_xi  = np.unique(self.xxsi[0])  
        e_eta = np.unique(self.xxsi[1])  
        
        self.n_elems = [ e_xi.shape[0]-1, e_eta.shape[0]-1 ]
 
        nbf_elem_xi  =  self.pp[0]+1
        nbf_elem_eta =  self.pp[1]+1
 
        self.dxi_dxi = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.dxi_xi  = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.xi_dxi  = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))
        self.xi_xi   = np.zeros((self.n_elems[0]*ne[0],nbf_elem_xi,nbf_elem_xi))

        self.deta_deta = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.deta_eta  = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.eta_deta  = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        self.eta_eta   = np.zeros((self.n_elems[1]*ne[1],nbf_elem_eta,nbf_elem_eta))
        
        Gauss_xi  =  nb.GaussLegendre(nip[0])
        Gauss_eta =  nb.GaussLegendre(nip[1])
        
        xig   = Gauss_xi[0]
        wgxi  = Gauss_xi[1]
        etag  = Gauss_eta[0]
        wgeta = Gauss_eta[1]
        
        mes_xi  =  (e_xi[1]  - e_xi[0] ) / ne[0]
        mes_eta =  (e_eta[1] - e_eta[0]) / ne[1]
 
        
        # Loop over basis elements (xi direction)
        for i in range(self.n_elems[0]):
            # loop over integration elements of the current element 
            for j in range(ne[0]):
                xi_min = e_xi[i] + j*mes_xi  
                xi_p   = xi_min  + 0.5*(xig+1)*mes_xi
                for ip in range(nip[0]):
                     N, dN = nb.derbasisfuns(i+self.pp[0], self.pp[0], self.xxsi[0], 1, xi_p[ip])
                     self.dxi_dxi[j+i*ne[0],:,:] += wgxi[ip]*np.outer(dN,dN)*mes_xi/2
                     self.dxi_xi[j+i*ne[0],:,:]  += wgxi[ip]*np.outer(dN,N)*mes_xi/2
                     self.xi_dxi[j+i*ne[0],:,:]  += wgxi[ip]*np.outer(N,dN)*mes_xi/2  
                     self.xi_xi[j+i*ne[0],:,:]   += wgxi[ip]*np.outer(N,N)*mes_xi/2  
                     
        # Loop over basis elements (eta direction)
        for i in range(self.n_elems[1]): 
            # loop over integration elements of the current element 
            for j in range(ne[1]):
                eta_min = e_eta[i] + j*mes_eta 
                eta_p   = eta_min  + 0.5*(etag+1)*mes_eta
                for ip in range(nip[1]):
                     N, dN = nb.derbasisfuns(i+self.pp[1], self.pp[1], self.xxsi[1], 1, eta_p[ip])
                     self.deta_deta[j+i*ne[1],:,:] += wgeta[ip]*np.outer(dN,dN)*mes_eta/2
                     self.deta_eta[j+i*ne[1],:,:]  += wgeta[ip]*np.outer(dN,N)*mes_eta/2
                     self.eta_deta[j+i*ne[1],:,:]  += wgeta[ip]*np.outer(N,dN)*mes_eta/2                      
                     self.eta_eta[j+i*ne[1],:,:]   += wgeta[ip]*np.outer(N,N)*mes_eta/2  
        

    def ComputeBivariateIntegratedDiffTensors(self):
        nbf_elem_xi    =  self.pp[0]+1
        nbf_elem_eta   =  self.pp[1]+1
        nbf_elem2d     =  nbf_elem_eta*nbf_elem_xi
 
        nei_xi  = self.dxi_dxi.shape[0]//self.n_elems[0] # Number of integration elements per xi element 
        nei_eta = self.deta_deta.shape[0]//self.n_elems[1] # Number of integration elements per eta element 
        
        nei_xi_tot = self.dxi_dxi.shape[0]
        
        # Computes for each 2d element the products of the 1d integrals 
        self.eta_eta_dxi_dxi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.eta_deta_dxi_xi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.deta_eta_xi_dxi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        self.deta_deta_xi_xi = np.zeros(( self.n_elems[0]*self.n_elems[1]*nei_xi*nei_eta, nbf_elem2d, nbf_elem2d ))
        
        
        # Loop over basis elements  
        for j in range(self.n_elems[1]):
            for i in range(self.n_elems[0]):
               # For each basis element 
               # Loop over integration elements 
               for ji in range(nei_eta):
                   ej = ji + j*nei_eta
                   for ii in range(nei_xi):
                       ei = ii + i*nei_xi 
                       iie = ei + ej*nei_xi_tot
                       self.eta_eta_dxi_dxi[iie,:,:] = np.kron(self.eta_eta[ej,:,:], self.dxi_dxi[ei,:,:])
                       self.eta_deta_dxi_xi[iie,:,:] = np.kron(self.eta_deta[ej,:,:], self.dxi_xi[ei,:,:]) 
                       self.deta_eta_xi_dxi[iie,:,:] = np.kron(self.deta_eta[ej,:,:], self.xi_dxi[ei,:,:]) 
                       self.deta_deta_xi_xi[iie,:,:] = np.kron(self.deta_deta[ej,:,:], self.xi_xi[ei,:,:])
                           
    
            
    def GaussIntegrationVectorized(self):
        """ Gauss integration: build of the global differential operators """
        nbg_xi    = self.pp[0]+1 
        nbg_eta   = self.pp[1]+1  

        Gauss_xi  =  nb.GaussLegendre(nbg_xi)
        Gauss_eta =  nb.GaussLegendre(nbg_eta)
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]

        e_xi  = np.unique(self.xxsi[0]) ; ne_xi  = e_xi.shape[0]-1
        e_eta = np.unique(self.xxsi[1]) ; ne_eta = e_eta.shape[0]-1
        xi_min    =  np.kron(e_xi[:-1],np.ones(nbg_xi))    
        eta_min   =  np.kron(e_eta[:-1],np.ones(nbg_eta))  
        xi_g      =  np.kron(np.ones(ne_xi)  , Gauss_xi[0])
        eta_g     =  np.kron(np.ones(ne_eta) , Gauss_eta[0])
        
        """ Measures of elements """
        mes_xi  = e_xi[1:]  - e_xi[:-1]
        mes_eta = e_eta[1:] - e_eta[:-1]
        
        mes_xi  = np.kron(mes_xi,np.ones(nbg_xi))   
        mes_eta = np.kron(mes_eta,np.ones(nbg_eta)) 
        
         
        """ Going from the reference element to the parametric space  """ 
        xi        = xi_min  + 0.5*(xi_g+1)*mes_xi     # Aranged gauss points in  xi direction  
        eta       = eta_min + 0.5*(eta_g+1)*mes_eta    # Aranged gauss points in  eta direction 
        
        phi_xi  , dphi_xi   = nb.global_basisfuns(self.pp[0],self.xxsi[0],xi)
        phi_eta , dphi_eta  = nb.global_basisfuns(self.pp[1],self.xxsi[1],eta)
        
        phi        = sps.kron( phi_eta  ,  phi_xi   ,  'csc') 
        dphidxi    = sps.kron( phi_eta  ,  dphi_xi  ,  'csc')
        dphideta   = sps.kron( dphi_eta ,  phi_xi   ,  'csc')
        self.npg       = phi.shape[0]
        
        wg_xi       =  np.kron(np.ones(ne_xi) , Gauss_xi[1])
        wg_eta      =  np.kron(np.ones(ne_eta), Gauss_eta[1])
        
        mes_xi  =  np.kron(np.ones(eta.shape[0]), mes_xi)
        mes_eta =  np.kron(mes_eta, np.ones(xi.shape[0])) 
          
        P = self.Get_Btot()
         
  
        isRational = self.IsRational()
 
        if isRational == False:
            
            """ Jacobian elements"""
            dxdxi  = dphidxi.dot(P[0,:])
            dxdeta = dphideta.dot(P[0,:])
            dydxi  = dphidxi.dot(P[1,:])
            dydeta = dphideta.dot(P[1,:])
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            """ Spatial derivatives """
            dphidx = sps.diags(dydeta/detJ).dot(dphidxi) + sps.diags(-dydxi/detJ).dot(dphideta)
            dphidy = sps.diags(-dxdeta/detJ).dot(dphidxi) + sps.diags(dxdxi/detJ).dot(dphideta)
            """ Integration weights + measures + Jacobian of the transformation """ 
            wg     =  np.kron(wg_eta, wg_xi)*np.abs(detJ)*mes_xi*mes_eta/4
            
        if isRational == True: 

             
            denom      = phi.dot( P[2,:]   )
            denomdxi   = dphidxi.dot( P[2,:]  )
            denomdeta  = dphideta.dot( P[2,:]   )

            Ddenom     = sps.diags(denom)
            Ddenom_squ = sps.diags(1/denom**2)
            Ddenomdxi  = sps.diags(denomdxi)
            Ddenomdeta = sps.diags(denomdeta)
            
            Wphi      = phi.multiply(P[2,:])
            Wdphidxi  = dphidxi.multiply(P[2,:])
            Wdphideta = dphideta.multiply(P[2,:])
            
            dNdxi  = Ddenom_squ.dot (  Ddenom.dot(Wdphidxi)  - Ddenomdxi.dot(Wphi) ) 
            dNdeta = Ddenom_squ.dot (  Ddenom.dot(Wdphideta) - Ddenomdeta.dot(Wphi) )
            
            
            dxdxi  = dNdxi.dot(P[0,:])
            dxdeta = dNdeta.dot(P[0,:])
            dydxi  = dNdxi.dot(P[1,:])
            dydeta = dNdeta.dot(P[1,:])
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            
            
            """ Spatial derivatives """
            dphidx = sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)
            dphidy = sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)
            """ Integration weights + measures + Jacobian of the transformation """ 
            wg     =  np.kron(wg_eta, wg_xi)*np.abs(detJ)*mes_xi*mes_eta/4
        
        self.wg = sps.diags(wg)
        zero    = sps.csr_matrix((self.npg,nbf))
        self.phi     = phi 
        self.phix    = sps.hstack((phi,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi)   ,  'csc')
        self.dphixdx = sps.hstack((dphidx,zero),  'csc') 
        self.dphixdy = sps.hstack((dphidy,zero),  'csc')
        self.dphiydx = sps.hstack((zero,dphidx),  'csc')
        self.dphiydy = sps.hstack((zero,dphidy),  'csc')       
        
    def GaussIntegration1(self):
        """ Build the integration scheme for the construction of the global differential operators """
        # Gauss quadrature points and weights  for domain elements 
        nbg_xi  = self.pp[0]+1 
        nbg_eta = self.pp[1]+1  
        Gauss_xi  =  nb.GaussLegendre(nbg_xi)
        Gauss_eta = nb.GaussLegendre(nbg_eta)
        w_g = np.kron(Gauss_eta[1], Gauss_xi[1])
        
        # ones_xi  = np.ones((nbg_xi,1)) # Useful for kronecker product ( similar to tile and repeat) 
        # ones_eta = np.ones((nbg_eta,1))
        
        nbg_elem = nbg_xi*nbg_eta
        
        # Number of elements 
        n_elems = self.ien[0].shape[1]*self.ien[1].shape[1]
        # Total number of Gauss integration points 
        self.npg = n_elems*nbg_elem
        nbvelem = nbg_elem**2
        
        indexI = np.zeros(self.npg*nbg_elem)
        indexJx = indexI.copy()
        indexJy = indexI.copy()
        v_phi = indexI.copy()
        v_dphidx = indexI.copy()
        v_dphidy = indexI.copy() 
        wg        = np.zeros(self.npg)
 
        
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = 2 
        ndof = dim*nbf 
        
        
 
 
        
        self.phi = dict()
        self.dphidx = dict()
        self.dphidy = dict()
        self.wdetJmes = dict() 
        self.mes = dict()
        
        listel_tot = self.Get_listeltot()
        P = self.Get_Btot() 
        
 
        isRational = self.IsRational()
        if isRational == True :

            for e in listel_tot: 
                ne_xi  = self.tripleien[e,0]
                ni_xi  = self.ien[0][0,ne_xi]
                ne_eta = self.tripleien[e,1]
                ni_eta = self.ien[1][0,ne_eta]
                
                xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
                eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
                # Treating elements of non zero measures only 
                mes_xi = xi_max-xi_min
                mes_eta = eta_max-eta_min
                self.mes[e] = mes_xi*mes_eta
                if self.mes[e] != 0:
                    # Mapping to the knot span 
                    xi_p   = xi_min  + 0.5*(Gauss_xi[0]+1)*mes_xi
                    eta_p  = eta_min + 0.5*(Gauss_eta[0]+1)*mes_eta
                    # Basis functions and dertivatives evaluated on Gauss points 
                    Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],xi_p,nbg_xi,ni_xi,1)
                    Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],eta_p,nbg_eta,ni_eta,1)

                    
                    Neta_Dot_Nxi      = np.kron(Neta,Nxi)
                    Neta_Dot_dNxidxi  = np.kron(Neta,dNxidxi)
                    dNetadeta_Dot_Nxi = np.kron(dNetadeta,Nxi)
                    
                    
                    num_vector     = Neta_Dot_Nxi*P[2,self.noelem[e,:]]
                    num_vector_xi  = Neta_Dot_dNxidxi*P[2,self.noelem[e,:]]
                    num_vector_eta = dNetadeta_Dot_Nxi*P[2,self.noelem[e,:]]
                    
                    denom     = np.sum(num_vector, axis = -1 )
                    denom_xi  = np.sum(num_vector_xi, axis = -1)
                    denom_eta = np.sum( num_vector_eta, axis = -1)
            
                    denom_g           =  sps.diags(denom)
                    inv_squ_denom_g   =  sps.diags(1./denom**2)
                    
                    dNdxi  = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_xi)  - sps.diags(denom_xi).dot(num_vector)  ) 
                    dNdeta = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_eta) - sps.diags(denom_eta).dot(num_vector) ) 
                    
                    # Jacobian elements 
                    dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]])
                    dxdeta = dNdeta.dot(P[0,self.noelem[e,:]])
                    dydxi  = dNdxi.dot(P[1,self.noelem[e,:]])
                    dydeta = dNdeta.dot(P[1,self.noelem[e,:]])
                    # Deterinant of Jacobian 
                    detJ   = dxdxi*dydeta - dydxi*dxdeta 
                  
                    indexI[e*nbvelem  : (e+1)*nbvelem]  = np.kron(np.arange(e*nbg_elem,(e+1)*nbg_elem), np.ones(nbg_elem))
                    indexJx[e*nbvelem  : (e+1)*nbvelem] = np.kron( np.ones(nbg_elem), self.noelem[e,:] )  
                    indexJy[e*nbvelem  : (e+1)*nbvelem] = np.kron( np.ones(nbg_elem), self.noelem[e,:]+nbf )
                    
                    v_phi[e*nbvelem  : (e+1)*nbvelem]    =  ( num_vector/(np.array([denom]).T) ).ravel()
                    v_dphidx[e*nbvelem  : (e+1)*nbvelem] =  ( sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)).ravel()
                    v_dphidy[e*nbvelem  : (e+1)*nbvelem] =  ( sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)).ravel() 
                    wg[e*nbg_elem:(e+1)*nbg_elem]  = w_g*np.abs(detJ)*self.mes[e]/4
     
                    
        if isRational == False : 
            
            for e in listel_tot: 
                ne_xi  = self.tripleien[e,0]
                ni_xi  = self.ien[0][0,ne_xi]
                ne_eta = self.tripleien[e,1]
                ni_eta = self.ien[1][0,ne_eta]
                
                xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
                eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
                # Treating elements of non zero measures only
                mes_xi = xi_max-xi_min
                mes_eta = eta_max-eta_min
                self.mes[e] = mes_xi*mes_eta
                if self.mes[e] != 0:
                    # Mapping to the knot span 
                    xi_p   = xi_min  + 0.5*(Gauss_xi[0]+1)*mes_xi
                    eta_p  = eta_min + 0.5*(Gauss_eta[0]+1)*mes_eta
                    # Basis functions and dertivatives evaluated on Gauss points 
                    Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],xi_p,nbg_xi,ni_xi,1)
                    Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],eta_p,nbg_eta,ni_eta,1)
              
 
                    Neta_Dot_Nxi      = np.kron(Neta,Nxi)
                    dNdxi             = np.kron(Neta,dNxidxi)
                    dNdeta            = np.kron(dNetadeta,Nxi)
        
                    
                    # Jacobian elements 
                    dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]])
                    dxdeta = dNdeta.dot(P[0,self.noelem[e,:]])
                    dydxi  = dNdxi.dot(P[1,self.noelem[e,:]])
                    dydeta = dNdeta.dot(P[1,self.noelem[e,:]])
                    # Determinant of Jacobian 
                    detJ   = dxdxi*dydeta - dydxi*dxdeta 
                  
                    indexI[e*nbvelem  : (e+1)*nbvelem] = np.kron(np.arange(e*nbg_elem,(e+1)*nbg_elem), np.ones(nbg_elem))
                    indexJx[e*nbvelem  : (e+1)*nbvelem] = np.kron( np.ones(nbg_elem), self.noelem[e,:] )  
                    indexJy[e*nbvelem  : (e+1)*nbvelem] = np.kron( np.ones(nbg_elem), self.noelem[e,:]+nbf )
                    
                    v_phi[e*nbvelem  : (e+1)*nbvelem]    =  Neta_Dot_Nxi.ravel()
                    v_dphidx[e*nbvelem  : (e+1)*nbvelem] =  ( sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)).ravel()
                    v_dphidy[e*nbvelem  : (e+1)*nbvelem] =  ( sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)).ravel() 
                    wg[e*nbg_elem:(e+1)*nbg_elem]  = w_g*np.abs(detJ)*self.mes[e]/4
        
        self.wg      =  sps.diags(wg)
        self.phix    =  sps.csc_matrix(( v_phi, (indexI,indexJx)), shape = (self.npg,ndof) )
        self.phiy    =  sps.csc_matrix(( v_phi, (indexI,indexJy)), shape = (self.npg,ndof) )
        self.dphixdx =  sps.csc_matrix(( v_dphidx, (indexI,indexJx)), shape = (self.npg,ndof) )
        self.dphixdy =  sps.csc_matrix(( v_dphidy, (indexI,indexJx)), shape = (self.npg,ndof) )
        self.dphiydx =  sps.csc_matrix(( v_dphidx, (indexI,indexJy)), shape = (self.npg,ndof) )
        self.dphiydy =  sps.csc_matrix(( v_dphidy, (indexI,indexJy)), shape = (self.npg,ndof) )          


    def GetGaussIntegrationPoints(self,nbg):
            """ Build the integration scheme """
            # Gauss quadrature points and weights  for domain elements 
            nbg_xi  = nbg[0]
            nbg_eta = nbg[1]  
            Gauss_xi  =  nb.GaussLegendre(nbg_xi)
            Gauss_eta = nb.GaussLegendre(nbg_eta)
            w_g = np.kron(Gauss_eta[1], Gauss_xi[1])
            ones_xi  = np.ones((nbg_xi,1)) # Useful for kronecker product ( similar to tile and repeat) 
            ones_eta = np.ones((nbg_eta,1))
            
            pix = np.array([])
            piy = np.array([]) 
            wgdetJ  = np.array([]) 
            
            self.mes = dict()
            
            listel_tot = self.Get_listeltot()
            P = self.Get_Btot() 
            isRational = self.IsRational()
             
            if isRational == True :
    
                for e in listel_tot: 
                    ne_xi  = self.tripleien[e,0]
                    ni_xi  = self.ien[0][0,ne_xi]
                    ne_eta = self.tripleien[e,1]
                    ni_eta = self.ien[1][0,ne_eta]
                    
                    xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
                    eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
                    # Treating elements of non zero measures only 
                    mes_xi = xi_max-xi_min
                    mes_eta = eta_max-eta_min
                    self.mes[e] = mes_xi*mes_eta
                    if self.mes[e] != 0:
                        # Mapping to the knot span 
                        xi_p   = xi_min  + 0.5*(Gauss_xi[0]+1)*mes_xi
                        eta_p  = eta_min + 0.5*(Gauss_eta[0]+1)*mes_eta
                        # Basis functions and dertivatives evaluated on Gauss points 
                        Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],xi_p,nbg_xi,ni_xi,1)
                        Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],eta_p,nbg_eta,ni_eta,1)
                        
                        r_Nxi       = np.kron(ones_eta,Nxi)
                        r_dNxidxi   = np.kron(ones_eta,dNxidxi)
                        r_Neta      = np.kron(Neta, ones_xi )
                        r_dNetadeta = np.kron(dNetadeta, ones_xi)
                        
                        Neta_Dot_Nxi      = np.einsum('nk,nl->nkl',r_Neta,r_Nxi).reshape(r_Neta.shape[0],-1)
                        Neta_Dot_dNxidxi  = np.einsum('nk,nl->nkl',r_Neta,r_dNxidxi).reshape(r_Neta.shape[0],-1)
                        dNetadeta_Dot_Nxi = np.einsum('nk,nl->nkl',r_dNetadeta,r_Nxi).reshape(r_dNetadeta.shape[0],-1) 
          
                        
                        num_vector     = Neta_Dot_Nxi*P[2,self.noelem[e,:]]
                        num_vector_xi  = Neta_Dot_dNxidxi*P[2,self.noelem[e,:]]
                        num_vector_eta = dNetadeta_Dot_Nxi*P[2,self.noelem[e,:]]
                        
                        denom     = np.sum(num_vector, axis = -1 )
                        denom_xi  = np.sum(num_vector_xi, axis = -1)
                        denom_eta = np.sum( num_vector_eta, axis = -1)
                
                        denom_g           =  sps.diags(denom)
                        inv_squ_denom_g   =  sps.diags(1./denom**2)
                        
                        dNdxi  = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_xi)  - sps.diags(denom_xi).dot(num_vector)  ) 
                        dNdeta = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_eta) - sps.diags(denom_eta).dot(num_vector) ) 
                        
                        # Jacobian elements 
                        dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]])
                        dxdeta = dNdeta.dot(P[0,self.noelem[e,:]])
                        dydxi  = dNdxi.dot(P[1,self.noelem[e,:]])
                        dydeta = dNdeta.dot(P[1,self.noelem[e,:]])
                        # Deterinant of Jacobian 
                        detJ   = dxdxi*dydeta - dydxi*dxdeta 
                        # Physical coordinates of the integration points
                        pix_coord = np.sum(num_vector*P[0,self.noelem[e,:]], axis =-1)/denom
                        piy_coord = np.sum(num_vector*P[1,self.noelem[e,:]], axis =-1)/denom
                        pix = np.r_[pix,pix_coord]
                        piy = np.r_[piy,piy_coord]
                        wgdetJ = np.r_[wgdetJ, w_g*np.abs(detJ)*self.mes[e]/4]
                        
            return pix,piy,wgdetJ   
        
        
    def SetBasisFunctionsAtIntegrationPoints(self):
        phi, dphidx, dphidy = nb.Get2dBasisFunctionsAtPts(self.pix, self.piy, self.xxsi[0], self.xxsi[1], self.pp[0], self.pp[1])
        nbf  =  self.Get_ndof()//2
        self.wg = sps.diags(self.wg)
        zero    = sps.csr_matrix((self.npg,nbf))
        self.phiMatrix = phi      
        self.phix    = sps.hstack((phi,zero)   ,  'csc')
        self.phiy    = sps.hstack((zero,phi)   ,  'csc')
        self.dphixdx = sps.hstack((dphidx,zero),  'csc') 
        self.dphixdy = sps.hstack((dphidy,zero),  'csc')
        self.dphiydx = sps.hstack((zero,dphidx),  'csc')
        self.dphiydy = sps.hstack((zero,dphidy),  'csc')   
    
    def get_U_Eps_Sigma_OnPointCloud(self,x,y,U,hooke):
 
        ndof = self.Get_ndof() 
        nbf = ndof//2 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
            
        phi, dphidx, dphidy = nb.Get2dBasisFunctionsAtPts(x,y, self.xxsi[0], self.xxsi[1], self.pp[0], self.pp[1])

        Ux = U[:nbf]; Uy = U[nbf:] # displacement at the control points 
   
        """ Displacement """
        ux = phi.dot(Ux)
        uy = phi.dot(Uy)
        """ Strain """ 
        exx = dphidx.dot(Ux)
        eyy = dphidy.dot(Uy)
        exy = 0.5*( dphidy.dot(Ux) + dphidx.dot(Uy) )
        """ Stress """ 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(2*exy)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(2*exy)
        sxy = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(2*exy)
        
        return ux,uy,exx,eyy,exy,sxx,syy,sxy 
    
    def StressFromStrain(self,U, hooke, epsxx, epsyy, epsxy):
        sigxx = hooke[0,0]*(epsxx) + hooke[0,1]*(epsyy) + hooke[0,2]*(2*epsxy)
        sigyy = hooke[1,0]*(epsxx) + hooke[1,1]*(epsyy) + hooke[1,2]*(2*epsxy)
        sigxy = hooke[2,0]*(epsxx) + hooke[2,1]*(epsyy) + hooke[2,2]*(2*epsxy)
        return sigxx,sigyy,sigxy
    

    def GetFcmLevelSetCellTesselationCoord(self,Phi,lvlmax):
        self.integrationCellsCoord = [] 
        self.pix = []
        self.piy = []
        self.wg  = []  
        e_xi      =  self.xxsi[0][self.pp[0]:-self.pp[0]]  
        e_eta     =  self.xxsi[1][self.pp[1]:-self.pp[1]]  
        for k in range(self.n_elems[0]*self.n_elems[1]):
            #print('---Element ', k)
            j = int(np.floor(k/self.n_elems[0]))
            i = k - j*self.n_elems[0]
            xmin    =  e_xi[i]
            xmax    =  e_xi[i+1]
            ymin   =  e_eta[j]
            ymax   =  e_eta[j+1]
            C = CellTesselation(xmin,xmax,ymin,ymax)
            C.DecomposeLevelSetCoord(self, Phi, lvlmax )
            
    def FcmLevelSetIntegrationTesselation(self,Phi,lvlmax):
        self.nbg_xi  = self.pp[0]+1
        self.nbg_eta  =self.pp[1]+1
        self.Gauss_xi  =  nb.GaussLegendre(self.nbg_xi)
        self.Gauss_eta =  nb.GaussLegendre(self.nbg_eta)
        self.wgRef = np.kron(self.Gauss_eta[1], self.Gauss_xi[1])
        self.refGaussTriangle = nb.GaussTriangle((max(self.pp[0],self.pp[1])))
        
        self.integrationCellsCoord = [] 
        self.pix = []
        self.piy = []
        self.wg  = []  
        e_xi      =  self.xxsi[0][self.pp[0]:-self.pp[0]]  
        e_eta     =  self.xxsi[1][self.pp[1]:-self.pp[1]]  
        for k in range(self.n_elems[0]*self.n_elems[1]):
#            print('---Element ', k)
            j = int(np.floor(k/self.n_elems[0]))
            i = k - j*self.n_elems[0]
            xmin    =  e_xi[i]
            xmax    =  e_xi[i+1]
            ymin   =  e_eta[j]
            ymax   =  e_eta[j+1]
            C = CellTesselation(xmin,xmax,ymin,ymax)
            C.DecomposeLevelSetIntegration(self, Phi, lvlmax )
        self.pix = np.array(self.pix)
        self.piy = np.array(self.piy)
        self.wg  = np.array(self.wg)
        self.npg = len(self.pix)
 
        
    def FcmLevelSetIntegrationTesselation2(self,Phi,lvlmax):
        self.nbg_xi  = self.pp[0]+1
        self.nbg_eta  =self.pp[1]+1
        self.Gauss_xi  =  nb.GaussLegendre(self.nbg_xi)
        self.Gauss_eta =  nb.GaussLegendre(self.nbg_eta)
        self.wgRef = np.kron(self.Gauss_eta[1], self.Gauss_xi[1])
        self.refGaussTriangle = nb.GaussTriangle((max(self.pp[0],self.pp[1])))

        nbf_elem = np.prod(self.pp+1)
       
        self.phi = dict()
        self.dphidx = dict()
        self.dphidy = dict()
        self.wdetJmes = dict() 
        self.mes = dict()
        
        self.pixe = dict() 
        self.piye = dict() 
        
        listel_tot = self.Get_listeltot()
 
        for e in listel_tot: 
            ne_xi  = self.tripleien[e,0]
            ni_xi  = self.ien[0][0,ne_xi]
            ne_eta = self.tripleien[e,1]
            ni_eta = self.ien[1][0,ne_eta]
            
            xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
            eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
            
            # Treating elements of non zero measures only 
            mes_xi = xi_max-xi_min
            mes_eta = eta_max-eta_min
            self.mes[e] = mes_xi*mes_eta
            
            if self.mes[e] != 0:
                # Here we should clear in the integraion points list that we filled for the previous element  
                self.pix = []
                self.piy = [] 
                self.wg  = []
                
                C = CellTesselation(xi_min, xi_max, eta_min, eta_max)
                C.DecomposeLevelSetIntegration(self, Phi, lvlmax )
                
                # Ater decomposing this cell  we get the integration points and weights 
                
                # Saving the coordinates of the integration points for each element if needed 
                self.pixe[e] = self.pix 
                self.piye[e] = self.piy 
                
                self.wdetJmes[e] = sps.diags(self.wg )  # Not necessary : just so that Stiffness2 works without changing anything 
                # Computing the basis functions at the integration points of the element 
                Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],self.pixe[e],len(self.pixe[e]),ni_xi,1)
                Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],self.piye[e],len(self.pixe[e]),ni_eta,1)
                
                self.phi[e]    = np.zeros( (len(self.pixe[e]),nbf_elem) )
                self.dphidx[e] = np.zeros( (len(self.pixe[e]),nbf_elem) )
                self.dphidy[e] = np.zeros( (len(self.pixe[e]),nbf_elem) )
                               
                # Performing row wise kronecker product  
                # Loop over integration points 
                for i in range(len(self.pixe[e])):
                    self.phi[e][i,:]    =  np.kron(Neta[i,:],Nxi[i,:])
                    self.dphidx[e][i,:] =  np.kron(Neta[i,:],dNxidxi[i,:])
                    self.dphidy[e][i,:] =  np.kron(dNetadeta[i,:],Nxi[i,:]) 

                    
        self.pix = []
        self.piy = []
        self.wg  = []            

                    
 
   
    def PlotFcmIntegrationElements(self,ax): 
        lines = [] 
        for c in self.integrationCellsCoord : 
            if c[0]=='r':
                xmin = c[1] ; xmax = c[2] ; ymin = c[3] ; ymax = c[4]
                lines.append (  [ (xmin,ymin),(xmax,ymin) ] )
                lines.append (  [ (xmin,ymin),(xmin,ymax) ] )  
                lines.append (  [ (xmax,ymin),(xmax,ymax) ] ) 
                lines.append (  [(xmin,ymax),(xmax,ymax)  ] ) 
            if c[0] =='t':
                # loop over element of the tesselation 
                for i in range(c[2].shape[0]):
                    n1 = c[2][i,0] ; n2 = c[2][i,1]  ; n3=  c[2][i,2]
                    x1 = c[1][n1,0] ; y1 = c[1][n1,1] 
                    x2 = c[1][n2,0] ; y2 = c[1][n2,1] 
                    x3 = c[1][n3,0] ; y3 = c[1][n3,1]
                    lines.append(  [(x1,y1),(x2,y2)])
                    lines.append(  [(x1,y1),(x3,y3)])
                    lines.append(  [(x2,y2),(x3,y3)])
        lc = cols.LineCollection(lines, linewidths=2) 
        ax.add_collection(lc)  
        
    def GetVisualizationPointsFromIntegrationElements(self,neval):
        # neval = [number of visu pts in x direction, number of visualization pts in y direction]
        lpx =  self.xxsi[0][-1]-self.xxsi[0][0]  # length of the FCM mesh in x and y directions 
        lpy =  self.xxsi[1][-1]-self.xxsi[1][0]
        area = lpx*lpy 
        nevalTotal = neval[0]*neval[1]
        xv = np.array([])
        yv = np.array([]) 
        for c in self.integrationCellsCoord:
            if c[0]=='r':
                xmin = c[1] ; xmax = c[2] ; ymin = c[3] ; ymax = c[4]
                lx = xmax -xmin
                ly = ymax -ymin 
                nvx = int(np.ceil(lx/lpx*neval[0])) 
                nvy = int(np.ceil(ly/lpy*neval[1])) 
                xp = np.linspace(xmin+1.e-12, xmax-1.e-12, nvx+1)
                yp = np.linspace(ymin+1.e-12, ymax-1.e-12 ,nvy+1)
                Xp,Yp = np.meshgrid(xp,yp)
                xv = np.r_[xv,Xp.ravel()]
                yv = np.r_[yv,Yp.ravel()] 
            if c[0] == 't':
                # loop over the elements of the tesselation 
                for i in range(c[2].shape[0]):
                    n1 = c[2][i,0]  ; n2 = c[2][i,1]  ; n3=  c[2][i,2]
                    x1 = c[1][n1,0] ; y1 = c[1][n1,1] 
                    x2 = c[1][n2,0] ; y2 = c[1][n2,1] 
                    x3 = c[1][n3,0] ; y3 = c[1][n3,1]
 
                    A = 0.5* ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)) # Triangle area 
                    nv = int(np.ceil(np.sqrt(nevalTotal*A/area)))*2 
                    
                    xp = np.linspace(0+1.e-12, 1-1.e-12, nv+1)
                    yp = np.linspace(0+1.e-12, 1-1.e-12 ,nv+1)
                    Xp,Yp = np.meshgrid(xp,yp)
                    
                    xp = Xp[Yp<1-Xp].ravel()
                    yp = Yp[Yp<1-Xp].ravel() 

                    # Mapping of the visualization points from the reference triangle to the physical triangle
                    xve  = (x2-x1)*xp + (x3-x1)*yp + x1 ;
                    yve  = (y2-y1)*xp + (y3-y1)*yp + y1 ;
 
                    xv = np.r_[xv, xve ]
                    yv = np.r_[yv, yve ] 
        return xv,yv 
                    
            
        
    
    
    def StiffnessVectorized(self,hooke):
        """ 
        Stiffness Matrix """
        
        if isinstance(hooke[0,0],sps.dia.dia_matrix) :
            """
            If Hook's law depends on the spatial position of each integration point 
            """
            Bxy=self.dphixdy+self.dphiydx
            K =  self.dphixdx.T.dot(hooke[0,0].dot(self.wg.dot(self.dphixdx))) +     \
                 self.dphiydy.T.dot(hooke[1,1].dot(self.wg.dot(self.dphiydy))) +   \
                 Bxy.T.dot(hooke[2,2].dot(self.wg).dot(Bxy)) + \
                 self.dphixdx.T.dot(hooke[0,1].dot(self.wg.dot(self.dphiydy))) +   \
                 self.dphixdx.T.dot(hooke[0,2].dot(self.wg.dot(Bxy))) +  \
                 self.dphiydy.T.dot(hooke[1,2].dot(self.wg.dot(Bxy))) +  \
                 self.dphiydy.T.dot(hooke[1,0].dot(self.wg.dot(self.dphixdx))) +   \
                 Bxy.T.dot(hooke[2,0].dot(self.wg.dot(self.dphixdx)))+  \
                 Bxy.T.dot(hooke[2,1].dot(self.wg.dot(self.dphiydy)))
                 
        else: 
            Bxy=self.dphixdy+self.dphiydx
            K =  hooke[0,0]*self.dphixdx.T.dot(self.wg.dot(self.dphixdx)) +   \
                 hooke[1,1]*self.dphiydy.T.dot(self.wg.dot(self.dphiydy)) +   \
                 hooke[2,2]*Bxy.T.dot(self.wg.dot(Bxy)) + \
                 hooke[0,1]*self.dphixdx.T.dot(self.wg.dot(self.dphiydy)) +   \
                 hooke[0,2]*self.dphixdx.T.dot(self.wg.dot(Bxy)) +  \
                 hooke[1,2]*self.dphiydy.T.dot(self.wg.dot(Bxy)) +  \
                 hooke[1,0]*self.dphiydy.T.dot(self.wg.dot(self.dphixdx)) +   \
                 hooke[2,0]*Bxy.T.dot(self.wg.dot(self.dphixdx)) +  \
                 hooke[2,1]*Bxy.T.dot(self.wg.dot(self.dphiydy))
        return K   
    
 
    def Stiffness1(self, hooke):
        """
        Returns stiffness matrix 
        Int ( Epsilon* Sigma) = Sum_(On elements of ) Int ( B'HB)
        """
        listel_tot = self.Get_listeltot()
        nbf_elem_uv = self.Get_nenunv()
        nbf_elem = nbf_elem_uv[0]*nbf_elem_uv[1]
        
        n_elems = listel_tot[-1] + 1
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = 2 
        ndof = dim*nbf 
 
        
        
        # Indices and values for the sparse stiffness matrix K  
        nnz  = (dim*nbf_elem)**2*n_elems
        indexI = np.zeros(nnz)
        indexJ = np.zeros(nnz)
        nnz_values = np.zeros(nnz)
        Ke = np.zeros((dim*nbf_elem,dim*nbf_elem)) # Elementary stiffness matrix 
        
        sp_count = 0 # Index couter for the sparse values

        for e in listel_tot :
             if self.mes[e] !=0 :
    
                 # Dot product of derivatives of basis functions 
                         
                 dNdx_dNdx = self.dphidx[e].T.dot(self.wdetJmes[e].dot(self.dphidx[e]))
                 dNdx_dNdy = self.dphidx[e].T.dot(self.wdetJmes[e].dot(self.dphidy[e]))
                 dNdy_dNdx = self.dphidy[e].T.dot(self.wdetJmes[e].dot(self.dphidx[e]))
                 dNdy_dNdy = self.dphidy[e].T.dot(self.wdetJmes[e].dot(self.dphidy[e]))
                
                 # 4 blocs of the elementary stiffness matrix 
                
                 # Bloc 0,0
                 Ke[:nbf_elem,:nbf_elem] =  hooke[0,0]*dNdx_dNdx + \
                                           hooke[2,0]*dNdy_dNdx + \
                                           hooke[0,2]*dNdx_dNdy + \
                                           hooke[2,2]*dNdy_dNdy 
                 # Bloc 0,1
                 Ke[:nbf_elem,nbf_elem:] =  hooke[0,1]*dNdx_dNdy +\
                                           hooke[2,1]*dNdy_dNdy +\
                                           hooke[0,2]*dNdx_dNdx +\
                                           hooke[2,2]*dNdy_dNdx 
                 # Bloc 1,0
                 Ke[nbf_elem:,:nbf_elem] =  hooke[1,0]*dNdy_dNdx +\
                                           hooke[2,0]*dNdx_dNdx +\
                                           hooke[1,2]*dNdy_dNdy +\
                                           hooke[2,2]*dNdx_dNdy 
                 # Bloc 1,1 
                 Ke[nbf_elem:,nbf_elem:] =  hooke[1,1]*dNdy_dNdy +\
                                           hooke[2,1]*dNdx_dNdy +\
                                           hooke[1,2]*dNdy_dNdx +\
                                           hooke[2,2]*dNdx_dNdx 
                                   
 
                                   
                 rep=np.r_[self.noelem[e,:],self.noelem[e,:]+nbf]
                 [repi,repj]=np.meshgrid(rep,rep)  
                 repi = repi.ravel()
                 repj = repj.ravel()
                 indexI[sp_count+np.arange(len(repi))] = repi
                 indexJ[sp_count+np.arange(len(repj))] = repj
                 nnz_values[sp_count+np.arange(len(repj))] = Ke.ravel()        
                 sp_count+=len(repj)
                 
        # Sparse stiffness matrix 
        K  = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return  K

    def Stiffness3(self,hooke):

        dim = 2    
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        ndof = dim*nbf 
 
        
        nei_xi  = self.dxi_dxi.shape[0]//self.n_elems[0]
        nei_eta = self.deta_deta.shape[0]//self.n_elems[1]
        
        nbf_elem = (self.pp[0]+1)*(self.pp[1]+1)
        
        
        # Indices and values for the sparse stiffness matrix K  
        n_elems = self.n_elems[0] * self.n_elems[1]
        
        nei_xi_tot = self.dxi_dxi.shape[0]
        
        nnz  = (dim*nbf_elem)**2*n_elems
        indexI = np.zeros(nnz)
        indexJ = np.zeros(nnz)
        nnz_values = np.zeros(nnz)
        Ke = np.zeros((dim*nbf_elem,dim*nbf_elem)) # Elementary stiffness matrix 
        
        sp_count = 0 # Index couter for the sparse values
        
        
        # Loop over basis elements  
        for j in range(self.n_elems[1]):
            for i in range(self.n_elems[0]):
               Ke = np.zeros((dim*nbf_elem,dim*nbf_elem)) # Elementary stiffness matrix 
               e = i + j*self.n_elems[0] 
               # For each basis element 
               # Loop over integration elements 
               for ji in range(nei_eta):
                   ej = ji + j*nei_eta
                   for ii in range(nei_xi):
                       ei = ii + i*nei_xi 
                       iie = ei + ej*nei_xi_tot
                       
                       # 4 blocs of the elementary stiffness matrix 
                       dNdx_dNdx = self.eta_eta_dxi_dxi[iie,:,:]     
                       dNdx_dNdy = self.eta_deta_dxi_xi[iie,:,:]       
                       dNdy_dNdx = self.deta_eta_xi_dxi[iie,:,:]                        
                       dNdy_dNdy = self.deta_deta_xi_xi[iie,:,:]                      
                                             
                       # Bloc 0,0
                       Ke[:nbf_elem,:nbf_elem] += hooke[0,0]*dNdx_dNdx + \
                                                  hooke[2,0]*dNdy_dNdx + \
                                                  hooke[0,2]*dNdx_dNdy + \
                                                  hooke[2,2]*dNdy_dNdy 
                       # Bloc 0,1
                       Ke[:nbf_elem,nbf_elem:] += hooke[0,1]*dNdx_dNdy +\
                                                  hooke[2,1]*dNdy_dNdy +\
                                                  hooke[0,2]*dNdx_dNdx +\
                                                  hooke[2,2]*dNdy_dNdx 
                       # Bloc 1,0
                       Ke[nbf_elem:,:nbf_elem] += hooke[1,0]*dNdy_dNdx +\
                                                  hooke[2,0]*dNdx_dNdx +\
                                                  hooke[1,2]*dNdy_dNdy +\
                                                  hooke[2,2]*dNdx_dNdy 
                       # Bloc 1,1 
                       Ke[nbf_elem:,nbf_elem:] += hooke[1,1]*dNdy_dNdy +\
                                                  hooke[2,1]*dNdx_dNdy +\
                                                  hooke[1,2]*dNdy_dNdx +\
                                                  hooke[2,2]*dNdx_dNdx 
                                           
               rep=np.r_[self.noelem[e,:],self.noelem[e,:]+nbf]
               [repi,repj]=np.meshgrid(rep,rep)  
               repi = repi.ravel()
               repj = repj.ravel()
               indexI[sp_count+np.arange(len(repi))] = repi
               indexJ[sp_count+np.arange(len(repj))] = repj
               nnz_values[sp_count+np.arange(len(repj))] = Ke.ravel()        
               sp_count+=len(repj)
                 
        # Sparse stiffness matrix 
        K  = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return  K
        
        
        
    
    def DispAtIp(self,U):
        ux = self.phix.dot(U)
        uy = self.phiy.dot(U)
        return ux, uy  
    def StrainAtIP(self,U):
        epsx=self.dphixdx.dot(U)
        epsy=self.dphiydy.dot(U)
        epsxy=0.5*self.dphixdy.dot(U)+0.5*self.dphiydx.dot(U)
        return epsx,epsy,epsxy  

    def GetBfForGridPoints(self,xi,eta):
        """ xi and eta are the 1d points 
        This method computes the basis functions on the mesh-grid point 
        obtained from the 1d vector points xi and eta 
        """

        phi_xi  , dphi_xi   = nb.global_basisfuns(self.pp[0],self.xxsi[0],xi)
        phi_eta , dphi_eta  = nb.global_basisfuns(self.pp[1],self.xxsi[1],eta)
        
        phi        = sps.kron( phi_eta  ,  phi_xi   ,  'csc') 
        dphidxi    = sps.kron( phi_eta  ,  dphi_xi  ,  'csc')
        dphideta   = sps.kron( dphi_eta ,  phi_xi   ,  'csc')
        
        P = self.Get_Btot()
        
        isRational = self.IsRational()
        if isRational == False:
            
            """ Jacobian elements"""
            dxdxi  = dphidxi.dot(P[0,:])
            dxdeta = dphideta.dot(P[0,:])
            dydxi  = dphidxi.dot(P[1,:])
            dydeta = dphideta.dot(P[1,:])
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            """ Spatial derivatives """
            dphidx = sps.diags(dydeta/detJ).dot(dphidxi) + sps.diags(-dydxi/detJ).dot(dphideta)
            dphidy = sps.diags(-dxdeta/detJ).dot(dphidxi) + sps.diags(dxdxi/detJ).dot(dphideta)
            
            # Univariate basis functions if needed 
            
#            Nxi  = phi_xi  
#            Neta = phi_eta 
            N = phi 
            
        if isRational == True: 

             
            denom      = phi.dot( P[2,:]   )
            denomdxi   = dphidxi.dot( P[2,:]  )
            denomdeta  = dphideta.dot( P[2,:]   )

            Ddenom     = sps.diags(denom)
            Ddenom_squ = sps.diags(1/denom**2)
            Ddenomdxi  = sps.diags(denomdxi)
            Ddenomdeta = sps.diags(denomdeta)
            
            Wphi      = phi.multiply(P[2,:])
            Wdphidxi  = dphidxi.multiply(P[2,:])
            Wdphideta = dphideta.multiply(P[2,:])
            
            dNdxi  = Ddenom_squ.dot (  Ddenom.dot(Wdphidxi)  - Ddenomdxi.dot(Wphi) ) 
            dNdeta = Ddenom_squ.dot (  Ddenom.dot(Wdphideta) - Ddenomdeta.dot(Wphi) )
            
            
            dxdxi  = dNdxi.dot(P[0,:])
            dxdeta = dNdeta.dot(P[0,:])
            dydxi  = dNdxi.dot(P[1,:])
            dydeta = dNdeta.dot(P[1,:])
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            
            """ Spatial derivatives """
            dphidx = sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)
            dphidy = sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)
            
            N = sps.diags(1/denom).dot(Wphi)
            
            # Univariate basis functions if needed 
            
#            w_xi  = self.ctrlPts[2][:,0] 
#            w_eta = self.ctrlPts[2][0,:]
#            denom_xi  =  phi_xi.dot(w_xi) 
#            denom_eta =  phi_eta.dot(w_eta) 
#            
#            Nxi  =   sps.diags(1/denom_xi).dot(phi_xi).multiply(w_xi) 
#            Neta =   sps.diags(1/denom_eta).dot(phi_eta).multiply(w_eta) 
              
        return N, dphidx, dphidy 
  
 
    def MatplotlibMeshVectorized(self, U, neval, color ):
        """ Physical elements = Image of the parametric elements on Python """
        P = self.Get_Btot() # control points 
        nbf = len(P[0,:]) #Number of control points (=ndof/2)
        Pxm = P[0,:] + U[:nbf]
        Pym = P[1,:] + U[nbf:]
        Pw = P[2,:]
        xi  = np.linspace( self.xxsi[0][self.pp[0]], self.xxsi[0][-self.pp[0]], neval[0] )
        eta = np.linspace( self.xxsi[1][self.pp[1]], self.xxsi[1][-self.pp[1]], neval[1] ) 
        # Iso parameters for the elemnts 
        xiu  = np.unique(self.xxsi[0]);   
        etau = np.unique(self.xxsi[1]);  
    
        # Basis functions 
        phi_xi1   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xiu)
        phi_eta1  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],eta)
        phi_xi2   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xi)
        phi_eta2  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],etau)    
        
        isRational = self.IsRational()
        if isRational == False :
            phi1      = sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') 
            phi2      = sps.kron( phi_eta2  ,  phi_xi2   ,  'csc')
        if isRational == True :
            phi1      = sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') 
            phi2      = sps.kron( phi_eta2  ,  phi_xi2   ,  'csc')
            Wphi1     = phi1.multiply(Pw)
            Wphi2     = phi2.multiply(Pw)
            denom1    = phi1.dot( Pw   )
            denom2    = phi2.dot( Pw   )
            phi1 = sps.diags(1/denom1).dot(Wphi1)
            phi2 = sps.diags(1/denom2).dot(Wphi2)
            
        xe1  = phi1.dot(Pxm) 
        ye1  = phi1.dot(Pym)
        xe2  = phi2.dot(Pxm) 
        ye2  = phi2.dot(Pym)  
        
        xe1 = xe1.reshape((xiu.size,neval[1]),order='F')
        ye1 = ye1.reshape((xiu.size,neval[1]),order='F')
        xe2 = xe2.reshape((neval[0],etau.size),order='F')
        ye2 = ye2.reshape((neval[0],etau.size),order='F')
        
        for i in range(xiu.size):
            # loop on xi 
            # Getting one eta iso-curve
            plt.plot( xe1[i,:],ye1[i,:], color=color   )
 
        for i in range(etau.size):
            # loop on eta 
            # Getting one xi iso-curve      
            plt.plot( xe2[:,i], ye2[:,i], color=color )  
            
    def VtkPlotVectorized(self, path , hooke, U, neval ):
        
        """ Interpolating fields """ 
        
        xi  = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
        
        phi, dphidx, dphidy  = self.GetBfForGridPoints(xi, eta) 
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = 2 
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        Ux = U[:nbf]; Uy = U[nbf:] # displacement at the control points 
        P = self.Get_Btot() # control points 
        
        """ Displacement """
        ux = phi.dot(Ux)
        uy = phi.dot(Uy)
        """ Strain """ 
        exx = dphidx.dot(Ux)
        eyy = dphidy.dot(Uy)
        exy = 0.5*( dphidy.dot(Ux) + dphidx.dot(Uy) )
        """ Stress """ 
        sxx = hooke[0,0]*(exx) + hooke[0,1]*(eyy) + hooke[0,2]*(2*exy)
        syy = hooke[1,0]*(exx) + hooke[1,1]*(eyy) + hooke[1,2]*(2*exy)
        sxy = hooke[2,0]*(exx) + hooke[2,1]*(eyy) + hooke[2,2]*(2*exy)
        """ Surface points """ 
        x  = phi.dot(P[0,:]) 
        y  = phi.dot(P[1,:])
 
  
        # Reshaping the fields  
        ux = ux.reshape((neval[0],neval[1],1),order='F')
        uy = uy.reshape((neval[0],neval[1],1),order='F')
        uz = ux *0 
        
        exx = exx.reshape((neval[0],neval[1],1),order='F')
        eyy = eyy.reshape((neval[0],neval[1],1),order='F')
        exy = exy.reshape((neval[0],neval[1],1),order='F')
        
        sxx = sxx.reshape((neval[0],neval[1],1),order='F')
        syy = syy.reshape((neval[0],neval[1],1),order='F')
        sxy = sxy.reshape((neval[0],neval[1],1),order='F')
        
        x = x.reshape((neval[0],neval[1],1),order='F')
        y = y.reshape((neval[0],neval[1],1),order='F')
        z = 0*x
        
        """ Exporting to VTK using EVTK library """ 
        start = (0,0,0)
        end = (neval[0]-1, neval[1]-1, 0)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("disp",(ux,uy,uz))
        w.addData("stress",(sxx,syy,sxy))
        w.addData("strain",(exx,eyy,exy))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz))
        w.appendData((sxx,syy,sxy))
        w.appendData((exx,eyy,exy))

        w.save()   
    
 
        

    def VtkSurfaceAndControlGrid(self,path,neval): 
        P = self.Get_Btot() # control points 
        dim = self.Get_dim()
        nbf_uv =  self.Get_nunv()
        nbf_xi  = nbf_uv[0] 
        nbf_eta = nbf_uv[1]
        nbf =  nbf_xi * nbf_eta 
        """ Exporting the control points could """ 
        xp = np.ascontiguousarray(P[0,:])
        yp = np.ascontiguousarray(P[1,:])
        if dim ==2 :
            tools.tvtk.pointsToVTK(path+'control-points',xp , yp, np.zeros(nbf), data= {"point": np.zeros(nbf)} ) 
        elif dim==3 : 
            zp = np.ascontiguousarray(P[2,:])
            tools.tvtk.pointsToVTK(path+'control-points', xp, yp, zp, data= {"point": np.zeros(nbf)} ) 
            
        """ Exporting the control grid """ 
        if dim==2 :
            Ix  = np.r_[self.ctrlPts[0,:,:].ravel(order='F'),self.ctrlPts[0,:,:].ravel(order='C')]    
            Iy  = np.r_[self.ctrlPts[1,:,:].ravel(order='F'),self.ctrlPts[1,:,:].ravel(order='C')]    
            Iz = Ix*0 
        elif dim==3 : 
            Ix  = np.r_[self.ctrlPts[0,:,:].ravel(order='F'),self.ctrlPts[0,:,:].ravel(order='C')]    
            Iy  = np.r_[self.ctrlPts[1,:,:].ravel(order='F'),self.ctrlPts[1,:,:].ravel(order='C')]  
            Iz  = np.r_[self.ctrlPts[2,:,:].ravel(order='F'),self.ctrlPts[2,:,:].ravel(order='C')]  
        pointsPerLine = np.r_[np.repeat(nbf_xi,nbf_eta), np.repeat(nbf_eta,nbf_xi)]

        npoints = Ix.size
        ncells = pointsPerLine.size
        # create some temporary arrays to write grid topology
        offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells):
            ii += pointsPerLine[i]
            offsets[i] = ii
        connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive
        cell_types = np.empty(npoints, dtype = 'uint8') 
        cell_types[:] = VtkPolyLine.tid
        
        w = VtkFile(path+'control-grid', VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells = ncells, npoints = npoints)
        
        w.openElement("Points")
        w.addData("points", (Ix,Iy,Iz))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (Ix,Iy,Iz) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)

        w.save()
        
        """ Exporting the surface """
        xi  = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
        
        phi, _, _  = self.GetBfForGridPoints(xi, eta) 
        
        x  = phi.dot(P[0,:]) 
        y  = phi.dot(P[1,:])
        if dim==2:
            z = np.zeros(len(x))
        elif dim ==3: 
            z = phi.dot(P[2,:]) 
            
        x = x.reshape((neval[0],neval[1],1),order='F')
        y = y.reshape((neval[0],neval[1],1),order='F')
        z = z.reshape((neval[0],neval[1],1),order='F')
        
        """ Exporting to VTK using EVTK library """ 
        start = (0,0,0)
        end = (neval[0]-1, neval[1]-1, 0)
        
        w =  VtkFile(path+'surface', VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )

        w.save()               
        
        
        
        

        
            
        
        

    def VtkMeshVectorized(self, path, U, neval ):
        """ Physical elements = Image of the parametric elements """ 
        
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = self.Get_dim()
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        
        Ux = U[:nbf]; Uy = U[nbf:2*nbf] # displacement at the control points 
        if dim==3:
            Uz = U[2*nbf:]
            
        P = self.Get_Btot() # control points 
        
        xi  = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
        
        
        
        # Iso parameters for the elemnts 
        xiu  = np.unique(self.xxsi[0]);  
        etau = np.unique(self.xxsi[1]);  
        
        # Basis functions 
        phi_xi1   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xiu)
        phi_eta1  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],eta)
        phi_xi2   = nb.global_basisfunsWd(self.pp[0],self.xxsi[0],xi)
        phi_eta2  = nb.global_basisfunsWd(self.pp[1],self.xxsi[1],etau)
        
        isRational = self.IsRational()
        if isRational == False :
            phi1      = sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') 
            phi2      = sps.kron( phi_eta2  ,  phi_xi2   ,  'csc')
        if isRational == True :
            phi1      = sps.kron( phi_eta1  ,  phi_xi1   ,  'csc') 
            phi2      = sps.kron( phi_eta2  ,  phi_xi2   ,  'csc')
            if dim==2:
                Wphi1     = phi1.multiply(P[2,:])
                Wphi2     = phi2.multiply(P[2,:])
                denom1    = phi1.dot( P[2,:]   )
                denom2    = phi2.dot( P[2,:]   )
            elif dim ==3: 
                Wphi1     = phi1.multiply(P[3,:])
                Wphi2     = phi2.multiply(P[3,:])
                denom1    = phi1.dot( P[3,:]   )
                denom2    = phi2.dot( P[3,:]   )

            phi1 = sps.diags(1/denom1).dot(Wphi1)
            phi2 = sps.diags(1/denom2).dot(Wphi2)
            
        
        xe1  = phi1.dot(P[0,:]) 
        ye1  = phi1.dot(P[1,:])
        xe2  = phi2.dot(P[0,:]) 
        ye2  = phi2.dot(P[1,:])
        
        xe1 = xe1.reshape((xiu.size,neval[1]),order='F')
        ye1 = ye1.reshape((xiu.size,neval[1]),order='F')
        xe2 = xe2.reshape((neval[0],etau.size),order='F')
        ye2 = ye2.reshape((neval[0],etau.size),order='F')
        
        uxe1 = phi1.dot(Ux)
        uye1 = phi1.dot(Uy)
        uxe1 = uxe1.reshape((xiu.size,neval[1]),order='F')
        uye1 = uye1.reshape((xiu.size,neval[1]),order='F')
        
        uxe2 = phi2.dot(Ux)
        uye2 = phi2.dot(Uy)
        uxe2 = uxe2.reshape((neval[0],etau.size),order='F')
        uye2 = uye2.reshape((neval[0],etau.size),order='F')
        
        
        """ Getting the Mesh Lines""" 
        if dim==2:
            Ix_eta = np.zeros( xiu.size*neval[1] ) 
            Iy_eta = np.zeros( xiu.size*neval[1] ) 
            ux_eta = np.zeros( xiu.size*neval[1] ) 
            uy_eta = np.zeros( xiu.size*neval[1] ) 
            for i in range(xiu.size):
                # loop on xi 
                # Getting one eta iso-curve
                Ix_eta[i*neval[1]:(i+1)*neval[1]]  = xe1[i,:]
                Iy_eta[i*neval[1]:(i+1)*neval[1]]  = ye1[i,:]
                ux_eta[i*neval[1]:(i+1)*neval[1]]  = uxe1[i,:]
                uy_eta[i*neval[1]:(i+1)*neval[1]]  = uye1[i,:]
            
            Ix_xi = np.zeros( etau.size*neval[0] ) 
            Iy_xi = np.zeros( etau.size*neval[0] ) 
            ux_xi = np.zeros( etau.size*neval[0] ) 
            uy_xi = np.zeros( etau.size*neval[0] ) 
            for i in range(etau.size):
                # loop on eta 
                # Getting one xi iso-curve             
                Ix_xi[i*neval[0]:(i+1)*neval[0]] = xe2[:,i]
                Iy_xi[i*neval[0]:(i+1)*neval[0]] = ye2[:,i]
                ux_xi[i*neval[0]:(i+1)*neval[0]] = uxe2[:,i]
                uy_xi[i*neval[0]:(i+1)*neval[0]] = uye2[:,i]
                
            
            Ix  = np.r_[Ix_xi, Ix_eta]
            Iy  = np.r_[Iy_xi, Iy_eta]
            Iz  = Ix*0
     
            Iux = np.r_[ux_xi,ux_eta] 
            Iuy = np.r_[uy_xi,uy_eta]
            Iuz = Iux*0 
            
        elif dim==3:
            ze1  = phi1.dot(P[2,:]) 
            ze2  = phi2.dot(P[2,:])
            
            ze1 = ze1.reshape((xiu.size,neval[1]),order='F')
            ze2 = ze2.reshape((neval[0],etau.size),order='F')
            
            uze1 = phi1.dot(Uz)
            uze1 = uze1.reshape((xiu.size,neval[1]),order='F')
 
            uze2 = phi2.dot(Uz)
            uze2 = uze2.reshape((neval[0],etau.size),order='F')
            
            Ix_eta = np.zeros( xiu.size*neval[1] ) 
            Iy_eta = np.zeros( xiu.size*neval[1] ) 
            Iz_eta = np.zeros( xiu.size*neval[1] ) 
            
            ux_eta = np.zeros( xiu.size*neval[1] ) 
            uy_eta = np.zeros( xiu.size*neval[1] ) 
            uz_eta = np.zeros( xiu.size*neval[1] )
            
            for i in range(xiu.size):
                # loop on xi 
                # Getting one eta iso-curve
                Ix_eta[i*neval[1]:(i+1)*neval[1]]  = xe1[i,:]
                Iy_eta[i*neval[1]:(i+1)*neval[1]]  = ye1[i,:]
                Iz_eta[i*neval[1]:(i+1)*neval[1]]  = ze1[i,:]
                
                ux_eta[i*neval[1]:(i+1)*neval[1]]  = uxe1[i,:]
                uy_eta[i*neval[1]:(i+1)*neval[1]]  = uye1[i,:]
                uz_eta[i*neval[1]:(i+1)*neval[1]]  = uze1[i,:]
            
            Ix_xi = np.zeros( etau.size*neval[0] ) 
            Iy_xi = np.zeros( etau.size*neval[0] ) 
            Iz_xi = np.zeros( etau.size*neval[0] ) 
            
            ux_xi = np.zeros( etau.size*neval[0] ) 
            uy_xi = np.zeros( etau.size*neval[0] ) 
            uz_xi = np.zeros( etau.size*neval[0] ) 
            for i in range(etau.size):
                # loop on eta 
                # Getting one xi iso-curve             
                Ix_xi[i*neval[0]:(i+1)*neval[0]] = xe2[:,i]
                Iy_xi[i*neval[0]:(i+1)*neval[0]] = ye2[:,i]
                Iz_xi[i*neval[0]:(i+1)*neval[0]] = ze2[:,i]
                
                ux_xi[i*neval[0]:(i+1)*neval[0]] = uxe2[:,i]
                uy_xi[i*neval[0]:(i+1)*neval[0]] = uye2[:,i]
                uz_xi[i*neval[0]:(i+1)*neval[0]] = uze2[:,i]
                
            
            Ix  = np.r_[Ix_xi, Ix_eta]
            Iy  = np.r_[Iy_xi, Iy_eta]
            Iz  = np.r_[Iz_xi, Iz_eta]
     
            Iux = np.r_[ux_xi,ux_eta] 
            Iuy = np.r_[uy_xi,uy_eta]
            Iuz = np.r_[uz_xi,uz_eta]
        
 
        
        """ Exporting to VTK using EVTK library """ 
        pointsPerLine = np.r_[np.repeat(neval[0],etau.size), np.repeat(neval[1],xiu.size)]
                
        npoints = Ix.size
        ncells = pointsPerLine.size
        # create some temporary arrays to write grid topology
        offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells):
            ii += pointsPerLine[i]
            offsets[i] = ii
        connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive
        cell_types = np.empty(npoints, dtype = 'uint8') 
        cell_types[:] = VtkPolyLine.tid
        
        w = VtkFile(path, VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells = ncells, npoints = npoints)
        
        w.openElement("Points")
        w.addData("points", (Ix,Iy,Iz))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        
        w.openData("Point")
        w.addData("disp",(Iux,Iuy,Iuz))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (Ix,Iy,Iz) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)
        w.appendData( (Iux,Iuy,Iuz))

        w.save()
        
 
    
    def GaussIntegration2(self):
        """ Build the integration scheme using the loop over elements"""
        # Gauss quadrature points and weights  for domain elements 
        nbg_xi  = self.pp[0]+1 
        nbg_eta = self.pp[1]+1  
        Gauss_xi  =  nb.GaussLegendre(nbg_xi)
        Gauss_eta = nb.GaussLegendre(nbg_eta)
        w_g = np.kron(Gauss_eta[1], Gauss_xi[1])
        # ones_xi  = np.ones((nbg_xi,1)) # Useful for kronecker product ( similar to tile and repeat) 
        # ones_eta = np.ones((nbg_eta,1))
 
 
        
        self.phi = dict()
        self.dphidx = dict()
        self.dphidy = dict()
        self.wdetJmes = dict() 
        self.mes = dict()
        
        listel_tot = self.Get_listeltot()
        P = self.Get_Btot() 
        
        self.phi_xi = dict()
        self.phi_eta = dict()
        self.dphi_xi = dict()
        self.dphi_eta = dict()
        i_xi = 0 ; i_eta = 0
        
 
        isRational = self.IsRational()
  
        if isRational == True :

            for e in listel_tot: 
                ne_xi  = self.tripleien[e,0]
                ni_xi  = self.ien[0][0,ne_xi]
                ne_eta = self.tripleien[e,1]
                ni_eta = self.ien[1][0,ne_eta]
                
                xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
                eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
                # Treating elements of non zero measures only 
                mes_xi = xi_max-xi_min
                mes_eta = eta_max-eta_min
                self.mes[e] = mes_xi*mes_eta
                if self.mes[e] != 0:
                    # Mapping to the knot span 
                    xi_p   = xi_min  + 0.5*(Gauss_xi[0]+1)*mes_xi
                    eta_p  = eta_min + 0.5*(Gauss_eta[0]+1)*mes_eta
                    # Basis functions and dertivatives evaluated on Gauss points 
                    Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],xi_p,nbg_xi,ni_xi,1)
                    Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],eta_p,nbg_eta,ni_eta,1)
                    
                    
                    Neta_Dot_Nxi      = np.kron(Neta,Nxi)     # (Matrix kronecker product)
                    Neta_Dot_dNxidxi  = np.kron(Neta,dNxidxi)
                    dNetadeta_Dot_Nxi = np.kron(dNetadeta,Nxi)
      
                    
                    num_vector     = Neta_Dot_Nxi*P[2,self.noelem[e,:]]
                    num_vector_xi  = Neta_Dot_dNxidxi*P[2,self.noelem[e,:]]
                    num_vector_eta = dNetadeta_Dot_Nxi*P[2,self.noelem[e,:]]
                    
                    denom     = np.sum(num_vector, axis = -1 )
                    denom_xi  = np.sum(num_vector_xi, axis = -1)
                    denom_eta = np.sum( num_vector_eta, axis = -1)
            
                    denom_g           =  sps.diags(denom)
                    inv_squ_denom_g   =  sps.diags(1./denom**2)
                    
                    dNdxi  = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_xi)  - sps.diags(denom_xi).dot(num_vector)  ) 
                    dNdeta = inv_squ_denom_g.dot  ( denom_g.dot(num_vector_eta) - sps.diags(denom_eta).dot(num_vector) ) 
                    
                    # Jacobian elements 
                    dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]] )
                    dxdeta = dNdeta.dot(P[0,self.noelem[e,:]] )
                    dydxi  = dNdxi.dot(P[1,self.noelem[e,:]] )
                    dydeta = dNdeta.dot(P[1,self.noelem[e,:]] )
                    # Determinant of Jacobian 
                    detJ   = dxdxi*dydeta - dydxi*dxdeta 
                  
  
                    
                    self.phi[e]  = num_vector/(np.array([denom]).T) 
                    self.dphidx[e]   = sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)
                    self.dphidy[e]   = sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)
                    self.wdetJmes[e] = sps.diags(w_g*np.abs(detJ)*self.mes[e]/4)
                    
                    # Saving univariate basis functions  
                    if eta_min == self.xxsi[1][0]:
                        self.phi_xi[i_xi] = Nxi
                        self.dphi_xi[i_xi] = dNxidxi
                        i_xi += 1 
                    if xi_min ==  self.xxsi[0][0]:
                        self.phi_eta[i_eta] = Neta
                        self.dphi_eta[i_eta] = dNetadeta   
                        i_eta +=1 
                    
 
                        
        if isRational == False : 
            
            for e in listel_tot: 
                ne_xi  = self.tripleien[e,0]
                ni_xi  = self.ien[0][0,ne_xi]
                ne_eta = self.tripleien[e,1]
                ni_eta = self.ien[1][0,ne_eta]
                
                xi_min  = self.xxsi[0][ni_xi] ; xi_max  = self.xxsi[0][ni_xi+1]
                eta_min = self.xxsi[1][ni_eta]; eta_max = self.xxsi[1][ni_eta+1]
                # Treating elements of non zero measures only
                mes_xi = xi_max-xi_min
                mes_eta = eta_max-eta_min
                self.mes[e] = mes_xi*mes_eta
                if self.mes[e] != 0:
                    # Mapping to the knot span 
                    xi_p   = xi_min  + 0.5*(Gauss_xi[0]+1)*mes_xi
                    eta_p  = eta_min + 0.5*(Gauss_eta[0]+1)*mes_eta
                    # Basis functions and dertivatives evaluated on Gauss points 
                    Nxi  , dNxidxi    = nb.derbasisfuncVectorInput(self.pp[0],self.xxsi[0],xi_p,nbg_xi,ni_xi,1)
                    Neta , dNetadeta  = nb.derbasisfuncVectorInput(self.pp[1],self.xxsi[1],eta_p,nbg_eta,ni_eta,1)
              
                    
                    Neta_Dot_Nxi      = np.kron(Neta,Nxi)
                    dNdxi  = np.kron(Neta,dNxidxi)
                    dNdeta = np.kron(dNetadeta,Nxi)
        
                    
                    # Jacobian elements 
                    dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]])
                    dxdeta = dNdeta.dot(P[0,self.noelem[e,:]])
                    dydxi  = dNdxi.dot(P[1,self.noelem[e,:]])
                    dydeta = dNdeta.dot(P[1,self.noelem[e,:]])
                    # Determinant of Jacobian 
                    detJ   = dxdxi*dydeta - dydxi*dxdeta 
                  
                    self.phi[e]  = Neta_Dot_Nxi 
                    self.dphidx[e]   = sps.diags(dydeta/detJ).dot(dNdxi) + sps.diags(-dydxi/detJ).dot(dNdeta)
                    self.dphidy[e]   = sps.diags(-dxdeta/detJ).dot(dNdxi) + sps.diags(dxdxi/detJ).dot(dNdeta)
                    self.wdetJmes[e] = sps.diags(w_g*np.abs(detJ)*self.mes[e]/4)
                    
                    # Saving univariate basis functions 
                    if eta_min == self.xxsi[1][0]:
                        self.phi_xi[i_xi] = Nxi
                        self.dphi_xi[i_xi] = dNxidxi
                        i_xi += 1 
                    if xi_min ==  self.xxsi[0][0]:
                        self.phi_eta[i_eta] = Neta
                        self.dphi_eta[i_eta] = dNetadeta   
                        i_eta +=1         
                        
    
    def Stiffness2(self,hooke):
        
        listel_tot = self.Get_listeltot()
        nbf_elem_uv = self.Get_nenunv()
        nbf_elem = nbf_elem_uv[0]*nbf_elem_uv[1]
        
        n_elems = listel_tot[-1] + 1
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = 2 
        ndof = dim*nbf 
        
        
    
        # Indices and values for the sparse stiffness matrix K  
        nnz  = (dim*nbf_elem)**2*n_elems
        indexI = np.zeros(nnz)
        indexJ = np.zeros(nnz)
        nnz_values = np.zeros(nnz)
        Ke = np.zeros((dim*nbf_elem,dim*nbf_elem)) # Elementary stiffness matrix 
        Be = np.zeros((3,dim*nbf_elem))   # Elementary differential matrix 
        
        
        sp_count = 0 # Index couter for the sparse values
        for e in listel_tot :
             # Loop over elements 
             Ke = Ke*0 
             if self.mes[e] !=0 :
                 # Getting the interation weight with the jacobian of the transformation 
                 wdetJ = self.wdetJmes[e].data[0]
                 # Loop over integration points 
                 nbg_elem = self.phi[e].shape[0] # Number of integration points per element 
                 for ig in range(nbg_elem):
                     
                    # Filling the elementary differential matrix
                    Be[0,:nbf_elem] =  self.dphidx[e][ig] 
                    Be[1,nbf_elem:] =  self.dphidy[e][ig]   
                    Be[2,:nbf_elem] =  self.dphidy[e][ig]   
                    Be[2,nbf_elem:] =  self.dphidx[e][ig]  
                    
                    # Elemental stiffness matrix 
                    Ke += wdetJ[ig]*Be.T.dot(hooke.dot(Be))
                    
 
                 
                 # Elementary contributions to the global stiffness matrix 
                 for ibf in range(nbf_elem):
                     for jbf in range(nbf_elem):
                         
                         # output of NOELEM 
                         I = self.noelem[e,ibf]
                         J = self.noelem[e,jbf]
             
                         indexI[sp_count] = I
                         indexJ[sp_count] = J
                         nnz_values[sp_count] = Ke[ibf,jbf]
                         sp_count += 1
                        
                         indexI[sp_count] = I
                         indexJ[sp_count] = J+nbf
                         nnz_values[sp_count] = Ke[ibf,jbf+nbf_elem]
                         sp_count += 1
                        
                         indexI[sp_count] = I+nbf
                         indexJ[sp_count] = J
                         nnz_values[sp_count] = Ke[ibf+nbf_elem,jbf]
                         sp_count += 1
                        
                         indexI[sp_count] = I + nbf
                         indexJ[sp_count] = J + nbf
                         nnz_values[sp_count] = Ke[ibf+nbf_elem,jbf+nbf_elem]
                         sp_count += 1
# =============================================================================
#                # The previous loop for adding the elementary contributions 
                 # can be replaced by the following lines : 
                         
                 # rep=np.r_[self.noelem[e,:],self.noelem[e,:]+nbf]
                 # [repi,repj]=np.meshgrid(rep,rep)  
                 # repi = repi.ravel()
                 # repj = repj.ravel()
                 # indexI[sp_count+np.arange(len(repi))] = repi
                 # indexJ[sp_count+np.arange(len(repj))] = repj
                 # nnz_values[sp_count+np.arange(len(repj))] = Ke.ravel()        
                 # sp_count+=len(repj)  
# =============================================================================

                 
        # Sparse stiffness matrix 
        K  = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (ndof,ndof ))
        return  K    
 
 
    def DisplacementAtIP_OfOneElem(self,indexE, indexIp, U):
        # IndexE is the index of the element 
        # indexIp is the index of the integration point in the element 
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        
        # Displacement of control points corresponding to the basis functions 
        # that support the desired element 
        Uxe = U[self.noelem[indexE,:]]
        Uye = U[self.noelem[indexE,:]+nbf]
        
        ux = self.phi[indexE][indexIp,:].dot(Uxe)
        uy = self.phi[indexE][indexIp,:].dot(Uye)
                                                           
        return ux,uy 
 
    
    def StrainAtIP_OfOneElem(self,indexE, indexIp, U):
        # IndexE is the index of the element 
        # indexIp is the index of the integration point in the element 
        
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        
        # Displacement of control points corresponding to the basis functions 
        # that support the desired element 
        Uxe = U[self.noelem[indexE,:]]
        Uye = U[self.noelem[indexE,:]+nbf]

        epsx  = self.dphidx[indexE][indexIp,:].dot(Uxe)
        epsy  = self.dphidy[indexE][indexIp,:].dot(Uye)
        epsxy = 0.5*( self.dphidy[indexE][indexIp,:].dot(Uxe)  +  self.dphidx[indexE][indexIp,:].dot(Uye)) 
        
        return epsx, epsy, epsxy 
    
    
    def GetBfOnOnePoint(self,xi,eta):
        """ Returns the non-zero basis basis functions only 
        one parametric point (xi,eta)
        """
        nbf_xi, nbf_eta  = self.Get_nunv()
        span_xi  = nb.findspan(nbf_xi,  self.pp[0], xi , self.xxsi[0])
        span_eta = nb.findspan(nbf_eta, self.pp[1], eta, self.xxsi[1])
        
        # Getting the index of the parametric element containing the point (xi,eta)
        ne_xi  = np.where(self.ien[0][0,:]==span_xi)[0][-1]
        ne_eta = np.where(self.ien[1][0,:]==span_eta)[0][-1]
        e     = ne_xi+ne_eta*self.ien[0].shape[1]
        
        N_xi  =  nb.derbasisfuns(span_xi,self.pp[0],self.xxsi[0],0,xi)  
        N_eta =  nb.derbasisfuns(span_eta,self.pp[1],self.xxsi[1],0,eta) 
        
        phi      = np.kron( N_eta[0,:], N_xi[0,:] ) 
    
        P = self.Get_Btot() 
        isRational = self.IsRational()
        
        if isRational == False : 
            N    = phi 
        
        if isRational == True :
            num     =  P[2,self.noelem[e,:]]*phi
            denom      = np.sum(num)
            N   = num/denom 

        return N, e   
 
    
    def GetBfAndDerivOnOnePoint(self,xi,eta):
        """ Returns the non-zero basis basis function and their derivatives
        one parametric point (xi,eta)
        """
        nbf_xi, nbf_eta  = self.Get_nunv()
        span_xi  = nb.findspan(nbf_xi,  self.pp[0], xi , self.xxsi[0])
        span_eta = nb.findspan(nbf_eta, self.pp[1], eta, self.xxsi[1])
        
        # Getting the index of the parametric element containing the point (xi,eta)
        ne_xi  = np.where(self.ien[0][0,:]==span_xi)[0][-1]
        ne_eta = np.where(self.ien[1][0,:]==span_eta)[0][-1]
        e     = ne_xi+ne_eta*self.ien[0].shape[1]
        
        ders_xi  =  nb.derbasisfuns(span_xi,self.pp[0],self.xxsi[0],1,xi)  
        ders_eta =  nb.derbasisfuns(span_eta,self.pp[1],self.xxsi[1],1,eta) 
        
        phi      = np.kron( ders_eta[0,:], ders_xi[0,:] ) 
        dphidxi  = np.kron( ders_eta[0,:], ders_xi[1,:] )
        dphideta = np.kron( ders_eta[1,:], ders_xi[0,:] )
        
        P = self.Get_Btot() 
        isRational = self.IsRational()
        
        if isRational == False : 
            
            # Jacobian elements 
            dxdxi  = dphidxi.dot(P[0,self.noelem[e,:]] )
            dxdeta = dphideta.dot(P[0,self.noelem[e,:]] )
            dydxi  = dphidxi.dot(P[1,self.noelem[e,:]] )
            dydeta = dphideta.dot(P[1,self.noelem[e,:]] )
            # Determinant of Jacobian 
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            
            dNdx = (dydeta/detJ)*(dphidxi)  - (dydxi/detJ)*(dphideta)
            dNdy = (-dxdeta/detJ)*(dphidxi) + (dxdxi/detJ)*(dphideta)
            N    = phi 

 
        if isRational == True :
            
            num     =  P[2,self.noelem[e,:]]*phi
            num_xi  =  P[2,self.noelem[e,:]]*dphidxi
            num_eta =  P[2,self.noelem[e,:]]*dphideta
            
            denom      = np.sum(num) # denoted the denominator x2 or y2 in the manuscript 
            denom_dxi  = np.sum(num_xi)
            denom_deta = np.sum(num_eta)
            
            dNdxi  =    (num_xi*denom  -  num*denom_dxi) /denom**2 
            dNdeta =    (num_eta*denom -  num*denom_deta)/denom**2 
            
            # Jacobian elements 
            dxdxi  = dNdxi.dot(P[0,self.noelem[e,:]] )
            dxdeta = dNdeta.dot(P[0,self.noelem[e,:]] )
            dydxi  = dNdxi.dot(P[1,self.noelem[e,:]] )
            dydeta = dNdeta.dot(P[1,self.noelem[e,:]] )
            # Determinant of Jacobian 
            detJ   = dxdxi*dydeta - dydxi*dxdeta 
            
           
            dNdx = (dydeta/detJ)*(dNdxi)  - (dydxi/detJ)*(dNdeta)
            dNdy = (-dxdeta/detJ)*(dNdxi) + (dxdxi/detJ)*(dNdeta)
            N   = num/denom 
            
        return N, dNdx, dNdy, e
    
    
    def DispAndSurfaceAtOneParametricPoint(self,xi,eta,U):
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        
        P = self.Get_Btot() 
        
        N, e = self.GetBfOnOnePoint(xi, eta)
        Uxe = U[self.noelem[e,:]]
        Uye = U[self.noelem[e,:]+nbf]
        
        ux = N.dot(Uxe)
        uy = N.dot(Uye)
        
        x  = N.dot(P[0,self.noelem[e,:]])
        y  = N.dot(P[1,self.noelem[e,:]])
        
        return x, y, ux, uy
   
    
    def DispStrainAndSurfaceAtOneParametricPoint(self,xi,eta,U):
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        
        P = self.Get_Btot() 
        
        N, dNdx, dNdy, e = self.GetBfAndDerivOnOnePoint(xi, eta)
        Uxe = U[self.noelem[e,:]]
        Uye = U[self.noelem[e,:]+nbf]
        
        ux = N.dot(Uxe)
        uy = N.dot(Uye)
        
        x  = N.dot(P[0,self.noelem[e,:]])
        y  = N.dot(P[1,self.noelem[e,:]])
        
        epsx  = dNdx.dot(Uxe)
        epsy  = dNdy.dot(Uye)
        epsxy = 0.5*( dNdy.dot(Uxe)  +  dNdx.dot(Uye) ) 
        
        return x, y, ux, uy, epsx, epsy, epsxy 
        


     
        
    def MatplotlibMesh1(self, U, neval, color):
        
        # Visualization points 
        xi  = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
    
        # Iso parameters for the elemnts 
        xiu  = np.unique(self.xxsi[0]);  
        etau = np.unique(self.xxsi[1]);  
        

        xe1  = np.zeros((xiu.size,neval[1]))
        ye1  = np.zeros((xiu.size,neval[1]))
        uxe1 = np.zeros((xiu.size,neval[1]))
        uye1 = np.zeros((xiu.size,neval[1]))
        
        xe2  = np.zeros((etau.size,neval[0]))
        ye2  = np.zeros((etau.size,neval[0]))
        uxe2 = np.zeros((etau.size,neval[0]))
        uye2 = np.zeros((etau.size,neval[0]))
        
        for i in range(xiu.size):
            # loop on xi knots  
            # Getting one eta iso-curve by looping over eta visualization points
            for j in range(neval[1]):
                xe1[i,j], ye1[i,j], uxe1[i,j], uye1[i,j] = self.DispAndSurfaceAtOneParametricPoint(xiu[i], eta[j], U)
                
        for i in range(etau.size):
            # loop on eta knots
            # Getting one xi iso-curve over xi visualization points 
            for j in range(neval[0]):
                xe2[i,j], ye2[i,j], uxe2[i,j], uye2[i,j] = self.DispAndSurfaceAtOneParametricPoint(xi[j], etau[i], U)
        
        for i in range(xiu.size):
            # loop on xi 
            # Getting one eta iso-curve
            plt.plot( xe1[i,:]+uxe1[i,:],ye1[i,:]+uye1[i,:], color=color   )
 
        for i in range(etau.size):
            # loop on eta 
            # Getting one xi iso-curve      
            plt.plot( xe2[i,:]+uxe2[i,:], ye2[i,:]+uye2[i,:], color=color )  
                
        
    def VtkPlot1(self, path, hooke, U, neval):
        """ Interpolating fields """ 
        nbf_uv =  self.Get_nunv()
        nbf =  nbf_uv[0]*nbf_uv[1]
        dim = 2 
        ndof = dim*nbf 
        if ndof!=U.shape[0] : 
            raise ValueError('Verify the shape of U')
        
        xi  = np.linspace(self.xxsi[0][self.pp[0]], self.xxsi[0][-self.pp[0]], neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]], self.xxsi[1][-self.pp[1]], neval[1])
        Xi,Eta = np.meshgrid(xi,eta)
        Xir  = np.ravel(Xi)
        Etar = np.ravel(Eta) 
        
        NvisPoints = len(Xir)
        x   = np.zeros(NvisPoints)
        y   = np.zeros(NvisPoints)
        z   = np.zeros(NvisPoints) 
        ux  = np.zeros(NvisPoints)
        uy  = np.zeros(NvisPoints)
        uz  = np.zeros(NvisPoints) 
        exx = np.zeros(NvisPoints)
        eyy = np.zeros(NvisPoints)
        exy = np.zeros(NvisPoints) 
        sxx = np.zeros(NvisPoints)
        syy = np.zeros(NvisPoints)
        sxy = np.zeros(NvisPoints)   
        
        for i in range(NvisPoints): 
            x[i], y[i], ux[i], uy[i], exx[i], eyy[i], exy[i] = self.DispStrainAndSurfaceAtOneParametricPoint(Xir[i],Etar[i],U)
        
        sxx,syy,sxy = self.StressFromStrain(U, hooke, exx, eyy, exy) 
        
        """ Exporting to VTK using EVTK library """ 
        start = (0,0,0)
        end = (neval[0]-1, neval[1]-1, 0)
        
        w =  VtkFile(path, VtkStructuredGrid)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        
        w.openData("Point")
        w.addData("disp",(ux,uy,uz))
        w.addData("stress",(sxx,syy,sxy))
        w.addData("strain",(exx,eyy,exy))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData( (ux,uy,uz))
        w.appendData((sxx,syy,sxy))
        w.appendData((exx,eyy,exy))

        w.save()        
        
        

     
    
    def VtkMesh1(self,path,U,neval):
        
        
        # Visualization points 
        xi  = np.linspace(self.xxsi[0][self.pp[0]] , self.xxsi[0][-self.pp[0]] , neval[0])
        eta = np.linspace(self.xxsi[1][self.pp[1]] , self.xxsi[1][-self.pp[1]] , neval[1])
    
        # Iso parameters for the elemnts 
        xiu  = np.unique(self.xxsi[0]);  
        etau = np.unique(self.xxsi[1]);  
        

        xe1  = np.zeros((xiu.size,neval[1]))
        ye1  = np.zeros((xiu.size,neval[1]))
        uxe1 = np.zeros((xiu.size,neval[1]))
        uye1 = np.zeros((xiu.size,neval[1]))
        
        xe2  = np.zeros((etau.size,neval[0]))
        ye2  = np.zeros((etau.size,neval[0]))
        uxe2 = np.zeros((etau.size,neval[0]))
        uye2 = np.zeros((etau.size,neval[0]))
        
        for i in range(xiu.size):
            # loop on xi knots  
            # Getting one eta iso-curve by looping over eta visualization points
            for j in range(neval[1]):
                xe1[i,j], ye1[i,j], uxe1[i,j], uye1[i,j] = self.DispAndSurfaceAtOneParametricPoint(xiu[i], eta[j], U)
                
        for i in range(etau.size):
            # loop on eta knots
            # Getting one xi iso-curve over xi visualization points 
            for j in range(neval[0]):
                xe2[i,j], ye2[i,j], uxe2[i,j], uye2[i,j] = self.DispAndSurfaceAtOneParametricPoint(xi[j], etau[i], U)
                
 
        """ Getting the Mesh Lines""" 
        Ix_eta = np.zeros( xiu.size*neval[1] ) 
        Iy_eta = np.zeros( xiu.size*neval[1] ) 
        ux_eta = np.zeros( xiu.size*neval[1] ) 
        uy_eta = np.zeros( xiu.size*neval[1] ) 
        for i in range(xiu.size):
            # loop on xi 
            # Getting one eta iso-curve
            Ix_eta[i*neval[1]:(i+1)*neval[1]]  = xe1[i,:]
            Iy_eta[i*neval[1]:(i+1)*neval[1]]  = ye1[i,:]
            ux_eta[i*neval[1]:(i+1)*neval[1]]  = uxe1[i,:]
            uy_eta[i*neval[1]:(i+1)*neval[1]]  = uye1[i,:]
        
        Ix_xi = np.zeros( etau.size*neval[0] ) 
        Iy_xi = np.zeros( etau.size*neval[0] ) 
        ux_xi = np.zeros( etau.size*neval[0] ) 
        uy_xi = np.zeros( etau.size*neval[0] ) 
        for i in range(etau.size):
            # loop on eta 
            # Getting one xi iso-curve             
            Ix_xi[i*neval[0]:(i+1)*neval[0]] = xe2[i,:]
            Iy_xi[i*neval[0]:(i+1)*neval[0]] = ye2[i,:]
            ux_xi[i*neval[0]:(i+1)*neval[0]] = uxe2[i,:]
            uy_xi[i*neval[0]:(i+1)*neval[0]] = uye2[i,:]
            
        
        """ Exporting to VTK using EVTK library """ 
        Ix  = np.r_[Ix_xi, Ix_eta]
        Iy  = np.r_[Iy_xi, Iy_eta]
        Iz  = Ix*0
        Iux = np.r_[ux_xi,ux_eta] 
        Iuy = np.r_[uy_xi,uy_eta]
        Iuz = Iux*0 
        pointsPerLine = np.r_[np.repeat(neval[0],etau.size), np.repeat(neval[1],xiu.size)]
                
        npoints = Ix.size
        ncells = pointsPerLine.size
        # create some temporary arrays to write grid topology
        offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
        ii = 0
        for i in range(ncells):
            ii += pointsPerLine[i]
            offsets[i] = ii
        connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive
        cell_types = np.empty(npoints, dtype = 'uint8') 
        cell_types[:] = VtkPolyLine.tid
        
        w = VtkFile(path, VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells = ncells, npoints = npoints)
        
        w.openElement("Points")
        w.addData("points", (Ix,Iy,Iz))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        
        w.openData("Point")
        w.addData("disp",(Iux,Iuy,Iuz))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (Ix,Iy,Iz) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)
        w.appendData( (Iux,Iuy,Iuz))

        w.save()
        

#%%% STATIC METHODS     
def Vtk_Disp_Strain_AtPoints(x,y,ux,uy,exx,eyy,exy,sxx,syy,sxy,filename):
        npoints = len(x)
        offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32')   # index of last node in each cell
        connectivity = np.arange(npoints, dtype = 'int32')                    # each point is only connected to itself
        cell_types = np.empty(npoints, dtype = 'uint8') 
   
        cell_types[:] = VtkVertex.tid
    
        w = VtkFile(filename, VtkUnstructuredGrid)
        w.openGrid()
        w.openPiece(ncells = npoints, npoints = npoints)
        
        w.openElement("Points")
        w.addData("points", (x,y,0*x ))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")
        
        w.openData("Point")
        w.addData("disp",(ux ,uy ,0*ux ))
        w.addData("strain",(exx,eyy,exy))
        w.addData("stress",(sxx,syy,sxy ))
        w.closeData("Point")
 
        w.closePiece()
        w.closeGrid()
        w.appendData( (x ,y ,0*x ) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)
        w.appendData( (ux ,uy ,0*ux ))
        w.appendData((exx,eyy,exy))
        w.appendData((sxx,syy,sxy))

        w.save() 
        


 

        
        
        
        
        
        
   
