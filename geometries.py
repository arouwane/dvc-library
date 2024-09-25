#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 18:47:42 2018

@author: rouwane
"""
 
import numpy as np 
import classIgaMesh as cim 
import gmsh 


def RectangularPlate(xmin,xmax,ymin,ymax, ne_xi, ne_eta, degree_xi, degree_eta): 
    
    # Parametric space properties 
    p=1; q=1
    Xi =   np.concatenate((np.repeat(0,p+1),np.repeat(1,p+1)))*(xmax-xmin) + xmin
    Eta =  np.concatenate((np.repeat(0,q+1),np.repeat(1,q+1)))*(ymax-ymin) + ymin 
     
    
    # Control points for a recangular plate 
    x = np.array([[xmin,xmin],
                  [xmax,xmax]]) 
    y  = np.array([[ymin,ymax],
                   [ymin,ymax]]) 
    w =  np.ones((q+1,p+1))
 
    ctrlPts = np.array([x,y,w])
  
    
    # Dictionary for the knot vector 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
     
     
    m = cim.Mesh(ctrlPts,np.array([p,q]),knot_vector)
    
 
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta]))
 
    # Knot refinement  
    ubar    = dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)    *(xmax-xmin)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)  *(ymax-ymin)
 
    m.KnotInsertion(ubar) 
  
    # Building connectivity 
    m.Connectivity()
    
    return m 
 

def geo_poutre2D_ct_curve(rmin, rmax, ne_xi, ne_eta, degree_xi, degree_eta): 
    # Parametric space properties 
    p=2; q=1
    Xi = np.concatenate((np.repeat(0,p+1),np.repeat(1,p+1)))
    Eta = np.concatenate((np.repeat(0,q+1),np.repeat(1,q+1)))
    
    # Control points for a circular beam of radius rmin and rmax 
    xdmin = list(rmin*np.array([1,1,0])) 
    xdmax = list(rmax*np.array([1,1,0])) 
    ydmin = list(rmin*np.array([0,1,1])) 
    ydmax = list(rmax*np.array([0,1,1])) 
    x=np.array([xdmax,xdmin]).T  
    y=np.array([ydmax,ydmin]).T  
    w=np.array([[1,1/np.sqrt(2),1],[1,1/np.sqrt(2),1]]).T # Weights of the control points 
    
    ctrlPts = np.array([x,y,w])
  
    
    # Dictionary for the knot vector 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
     
    m = cim.Mesh(ctrlPts,np.array([p,q]),knot_vector)
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    m.KnotInsertion(ubar)  
 
    # Building connectivity 
    m.Connectivity()
    
    return m 

def disk(rmin,rmax,ne_xi,ne_eta,degree_xi,degree_eta):
    # Parametric space propreties
    if ne_xi %4 !=0:
        raise ValueError('Number of elements must be multiple of 4 for uniform subdivision')
     
    p=2; q=1
    Xi  = np.array([0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1])
    Eta = np.array([0,0,1,1])
    xdmax = list(rmax*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmax = list(rmax*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    xdmin =  list(rmin*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmin = list(rmin*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    c = 1/np.sqrt(2)
    w     =  list(np.array([1,c,1,c,1,c,1,c,1])) 
    
    x=np.array([xdmax,xdmin]).T  
    y=np.array([ydmax,ydmin]).T  
    w=np.array([w,w]).T # Weights of the control points 
    
    ctrlPts = np.array([x,y,w])
    
     # Dictionary for the knot vector 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
     
    m = cim.Mesh(ctrlPts,np.array([p,q]),knot_vector)
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = np.r_[ np.linspace(0,1/4,ne_xi//4+1)[1:-1], 
                    np.linspace(1/4,1/2,ne_xi//4+1)[1:-1] , 
                    np.linspace(1/2,3/4,ne_xi//4+1)[1:-1] , 
                    np.linspace(3/4,1,ne_xi//4+1)[1:-1]  ]
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    m.KnotInsertion(ubar)  
 
    # Building connectivity 
    m.Connectivity()
    
    return m  

def degenerated_disk(x0,y0,R,ne_xi,ne_eta,degree_xi,degree_eta):
    p = 2; q= 2 
    Xi  = np.array([0,0,0,1,1,1])
    Eta = np.array([0,0,0,1,1,1])
    x = np.array([[-1,-2, -1],
                  [ 0,  0,   0],
                  [ 1, 2,  1]])*R + x0
    y = np.array([[-1, 0, 1],
                  [-2, 0, 2],
                  [-1, 0, 1]])*R + y0 
    c = 1/np.sqrt(2)
    w = np.array([[1,c,1],
                  [c,1,c],
                  [1,c,1]])
    ctrlPts = np.array([x,y,w])
    
     # Dictionary for the knot vector 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
     
    m = cim.Mesh(ctrlPts,np.array([p,q]),knot_vector)        

    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta]))

    # Knot refinement 
    ubar=dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    m.KnotInsertion(ubar)  
 
    # Building connectivity 
    m.Connectivity()
    
    return m      

def cylinder2(x0,y0,z0,R,H,ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta):
    
    dxi   = 2 
    deta  = 2
    dzeta = 1
    Xi   = np.array([0,0,0,1,1,1])
    Eta  = np.array([0,0,0,1,1,1])
    Zeta = np.array([0,0,1,1])
    
    # Base of the cylinder 
    c = 1/np.sqrt(2)
    Xb = np.array([[-c,-2*c, -c],
                  [ 0,  0,   0],
                  [ c, 2*c,  c]])*R + x0
    Yb = np.array([[-c, 0, c],
                  [-2*c, 0, 2*c],
                  [-c, 0, c]])*R + y0 
    Wb = np.array([[1,c,1],
                  [c,1,c],
                  [1,c,1]])
    Z1 = np.ones(Xb.shape)
    
    X = np.array([Xb,Xb])
    Y = np.array([Yb,Yb])
    Z = np.array([Z1*z0,Z1*z0+H])
    W = np.array([Wb, Wb])
    
 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W, np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
    
    return m  


def quarterTorus(rmin,rmax,ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta):
    
    tr = (rmax -rmin)/2 # Torus tube radius 
    
    dxi   = 2 
    deta  = 1
    dzeta = 2
    Xi   = np.array([0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1])
    Eta  = np.array([0,0,1,1])
    Zeta = np.array([0,0,0,1,1,1])
    
    # Base of the cylinder 
    xdmax = list(tr*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmax = list(tr*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    xdmin = list(0*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmin = list(0*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    c = 1/np.sqrt(2)
    w     =  list(np.array([1,c,1,c,1,c,1,c,1])) 
    Xb = np.array([xdmax,xdmin]).T  
    Yb = np.array([ydmax,ydmin]).T 
    Wb = np.array([w,w]).T   # Weights of the control points 
    Zb = np.ones(Xb.shape)
    
    X1 = Xb  + (rmax+rmin)/2
    Y1 = Yb   
    Z1 = Zb  *0 
    
    Ry  = rotationMatrix(-np.pi/2, axis='y')
    C3  = np.c_[X1.ravel(), Y1.ravel(), Z1.ravel()].T
    C3r = Ry.dot(C3)
    X3 = C3r[0].reshape(X1.shape[0],-1)
    Y3 = C3r[1].reshape(X1.shape[0],-1)
    Z3 = C3r[2].reshape(X1.shape[0],-1)
 
    Ry  = rotationMatrix(-np.pi/4, axis='y')
    C2  = np.c_[X1.ravel(), Y1.ravel(), Z1.ravel()].T
    C2r = Ry.dot(C2)
    X2 = X1
    Y2 = C2r[1].reshape(X1.shape[0],-1)
    Z2 = Z3
    
    X = np.array([X1,X2,X3])
    Y = np.array([Y1,Y2,Y3])
    Z = np.array([Z1,Z2,Z3])
    W = np.array([Wb,Wb*c,Wb])
 
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W , np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = np.r_[ np.linspace(0,1/4,ne_xi//4+1)[1:-1], 
                    np.linspace(1/4,1/2,ne_xi//4+1)[1:-1] , 
                    np.linspace(1/2,3/4,ne_xi//4+1)[1:-1] , 
                    np.linspace(3/4,1,ne_xi//4+1)[1:-1]  ]
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
 
    
    return m 
        
    

def quarterTorus2(rmin,rmax,ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta):
    tr = (rmax -rmin)/2 # Torus tube radius 
    
    dxi   = 2 
    deta  = 2
    dzeta = 2
    Xi   = np.array([0,0,0,1,1,1])
    Eta  = np.array([0,0,0,1,1,1])
    Zeta = np.array([0,0,0,1,1,1])
    
   # Base of the cylinder 
    c = 1/np.sqrt(2)   
    Xb = np.array([[-c,-2*c, -c],
                  [ 0,  0,   0],
                  [ c, 2*c,  c]])*tr  
    Yb = np.array([[-c, 0, c],
                  [-2*c, 0, 2*c],
                  [-c, 0, c]])*tr  

    Wb = np.array([[1,c,1],
                  [c,1,c],
                  [1,c,1]])
    Zb = np.ones(Xb.shape)
    
    X1 = Xb  + (rmax+rmin)/2
    Y1 = Yb   
    Z1 = Zb  *0 
    
    Ry  = rotationMatrix(-np.pi/2, axis='y')
    C3  = np.c_[X1.ravel(), Y1.ravel(), Z1.ravel()].T
    C3r = Ry.dot(C3)
    X3 = C3r[0].reshape(X1.shape[0],-1)
    Y3 = C3r[1].reshape(X1.shape[0],-1)
    Z3 = C3r[2].reshape(X1.shape[0],-1)
 
    Ry  = rotationMatrix(-np.pi/4, axis='y')
    C2  = np.c_[X1.ravel(), Y1.ravel(), Z1.ravel()].T
    C2r = Ry.dot(C2)
    X2 = X1
    Y2 = C2r[1].reshape(X1.shape[0],-1)
    Z2 = Z3
    
    X = np.array([X1,X2,X3])
    Y = np.array([Y1,Y2,Y3])
    Z = np.array([Z1,Z2,Z3])
    W = np.array([Wb,Wb*c,Wb])

    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W, np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
    
    return m  
    

def geo_Cylinder_Solid(rmin,rmax,ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta): 
    
    dxi   = 2 
    deta  = 1
    dzeta = 1
    Xi   = np.array([0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1])
    Eta  = np.array([0,0,1,1])
    Zeta = np.array([0,0,1,1])
    
    # Base of the cylinder 
    xdmax = list(rmax*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmax = list(rmax*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    xdmin = list(rmin*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmin = list(rmin*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    c = 1/np.sqrt(2)
    w     =  list(np.array([1,c,1,c,1,c,1,c,1])) 
    Xb = np.array([xdmax,xdmin]).T  
    Yb = np.array([ydmax,ydmin]).T 
    Wb = np.array([w,w]).T   # Weights of the control points 
    Zb = np.ones(Xb.shape)
    
    X = np.array([Xb,Xb])
    Y = np.array([Yb,Yb])
    Z = np.array([Zb*0,Zb])
    W = np.array([Wb ,Wb ])
    
  
    
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W, np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = np.r_[ np.linspace(0,1/4,ne_xi//4+1)[1:-1], 
                    np.linspace(1/4,1/2,ne_xi//4+1)[1:-1] , 
                    np.linspace(1/2,3/4,ne_xi//4+1)[1:-1] , 
                    np.linspace(3/4,1,ne_xi//4+1)[1:-1]  ]
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
 
    
    return m 

def cuboid(xmin,xmax,ymin,ymax,zmin,zmax,ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta): 
    dxi   = 1 
    deta  = 1
    dzeta = 1
    Xi   = np.array([0,0,1,1])
    Eta  = np.array([0,0,1,1])
    Zeta = np.array([0,0,1,1])
    
    
    X = np.array([[[0,0],
                   [1,1.]],
                  [[0,0],
                   [1,1]]])*(xmax-xmin) + xmin  
    Y = np.array([[[0,1],
                   [0,1]],
                  [[0,1],
                   [0,1]]])*(ymax-ymin) + ymin  
    Z = np.array([[[0,0],
                   [0,0]],
                  [[1,1],
                   [1,1]]])*(zmax-zmin) + zmin  
    W = np.array([[[1,1],
                   [1,1]],
                  [[1,1],
                   [1,1]]]) 
    
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W, np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = 1/ne_xi*np.arange(1,ne_xi)
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
    
    return m 


def cilynder(x0,y0, rmin,rmax, H, ne_xi,ne_eta,ne_zeta,degree_xi,degree_eta,degree_zeta):
    
    dxi   = 2 
    deta  = 1
    dzeta = 1
    Xi   = np.array([0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1])
    Eta  = np.array([0,0,1,1])
    Zeta = np.array([0,0,1,1])
    
    # Base of the cylinder 
    xdmax = list(rmax*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmax = list(rmax*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    xdmin = list(rmin*np.array([0,1,1,1,0,-1,-1,-1,0]))
    ydmin = list(rmin*np.array([-1,-1,0,1,1,1,0,-1,-1]))
    c = 1/np.sqrt(2)
    w     =  list(np.array([1,c,1,c,1,c,1,c,1])) 
    Xb = np.array([xdmax,xdmin]).T  + x0 
    Yb = np.array([ydmax,ydmin]).T  + y0 
    Wb = np.array([w,w]).T   # Weights of the control points 
    Zb = np.ones(Xb.shape)*H
    
    X = np.array([Xb,Xb])
    Y = np.array([Yb,Yb])
    Z = np.array([Zb*0,Zb])
    W = np.array([Wb ,Wb ])
    
  
    
    knot_vector=dict()
    knot_vector[0]= Xi
    knot_vector[1]= Eta 
    knot_vector[2]= Zeta 
    
    m = cim.MeshVolume(X, Y, Z, W, np.array([dxi,deta,dzeta]), knot_vector)  
    
    # Degree elevation 
    m.DegElevation(np.array([degree_xi, degree_eta,degree_zeta]))
    
    # Knot refinement 
    ubar=dict()
    ubar[0] = np.r_[ np.linspace(0,1/4,ne_xi//4+1)[1:-1], 
                    np.linspace(1/4,1/2,ne_xi//4+1)[1:-1] , 
                    np.linspace(1/2,3/4,ne_xi//4+1)[1:-1] , 
                    np.linspace(3/4,1,ne_xi//4+1)[1:-1]  ]
    ubar[1] = 1/ne_eta*np.arange(1,ne_eta)
    ubar[2] = 1/ne_zeta*np.arange(1,ne_zeta)
    m.KnotInsertion(ubar)   
 
    
    return m    


 

def CylinderTetra(x0,y0,z0,R,h,lc):
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("cylinder")
    gmsh.model.geo.addPoint( x0, y0, z0, lc, 1)  
    gmsh.model.geo.addPoint( x0-R, y0, z0, lc, 2) 
    gmsh.model.geo.addPoint(x0+R, y0, z0, lc, 3) 
    gmsh.model.geo.addPoint( x0, y0, z0+R, lc, 4) 
    gmsh.model.geo.addPoint( x0, y0, z0-R, lc, 5) 
    gmsh.model.geo.addCircleArc(2, 1, 4, 1)   
    gmsh.model.geo.addCircleArc(4, 1, 3, 2) 
    gmsh.model.geo.addCircleArc(3, 1, 5, 3)  
    gmsh.model.geo.addCircleArc(5, 1, 2, 4)   
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)  
    gmsh.model.geo.addPlaneSurface([1],1)    
    gmsh.model.geo.extrude([(2,1)], 0, h, 0)       
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    
    nums,nodes,e=gmsh.model.mesh.getNodes()
    nodes=nodes.reshape((len(nums),3))
    
    nums,els = gmsh.model.mesh.getElementsByType(4)
    nnd=len(els)//len(nums)
    els=els.reshape((len(nums),nnd))-1
    
    m= cim.TetraMesh1(els,nodes)
    m.CleanMesh()
    return m 
 
def TetMeshBox(origin,box,lc=0.5):
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("box")
    w = box[0]
    h = box[1]
    L = box[2]
    gmsh.model.geo.addPoint( 0, 0 , 0, lc, 1)
    gmsh.model.geo.addPoint( w, 0 , 0, lc, 2)
    gmsh.model.geo.addPoint( w, h , 0, lc, 3)
    gmsh.model.geo.addPoint( 0, h , 0, lc, 4)    
    gmsh.model.geo.addLine(1,2,1)
    gmsh.model.geo.addLine(2,3,2)
    gmsh.model.geo.addLine(3,4,3)
    gmsh.model.geo.addLine(4,1,4)    
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addPlaneSurface([1],1)    
    gmsh.model.geo.extrude([(2, 1)], 0, 0, L)    
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    
    
    nums,nodes,e=gmsh.model.mesh.getNodes()
    nodes=nodes.reshape((len(nums),3))
    
    nums,els = gmsh.model.mesh.getElementsByType(4)
    nnd=len(els)//len(nums)
    els=els.reshape((len(nums),nnd))-1
    
    m = cim.TetraMesh1( els , nodes  )
    m.n[:,0] += origin[0]
    m.n[:,1] += origin[1]
    m.n[:,2] += origin[2]
    m.CleanMesh()
    return m 

    


    
def rotationMatrix(theta,axis):
    if axis=='x':
        return np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)], [0,np.sin(theta), np.cos(theta) ]]) 
    if axis=='y':
        return np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]])
    if axis=='z':
        return np.array([[np.cos(theta) , -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0], [0, 0, 1]])     





    




 
    
    

    
    
