#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:50:26 2021

@author: rouwane
"""
 

import scipy.sparse as sps
import numpy as np 
import meshio 
import pandas  
import triangles_to_VTK as tvtk 

def rotationMatrix(theta,axis):
    if axis=='x':
        return np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)], [0,np.sin(theta), np.cos(theta) ]]) 
    if axis=='y':
        return np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]])
    if axis=='z':
        return np.array([[np.cos(theta) , -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0], [0, 0, 1]])     

def TransfomationMatrix(tx,ty,tz,thetax,thetay,thetaz,sx,sy,sz):
    T = np.array([[1,0,0,tx],
                  [0,1,0,ty],
                  [0,0,1,tz],
                  [0,0,0,1]])
    S = np.array([[sx,0,0,0],
                  [0,sy,0,0],
                  [0,0,sz,0],
                  [0,0,0,1]])
    Rx = rotationMatrix(thetax, 'x')
    Ry = rotationMatrix(thetay, 'y')
    Rz = rotationMatrix(thetaz, 'z')
    Rnh = Rx@Ry@Rz
    R = np.vstack((Rnh,np.array([0,0,0])))
    R = np.hstack((R,np.array([[0,0,0,1]]).T)) 
    P = T@R@S 
    return P



def GetBinaryMatrix(index,size): 
    nnz = len(index)
    nnz_values = np.ones(nnz)
    D = sps.csr_matrix((nnz_values, (index,index)), shape=(size,size) )
    return D

def readRawImage(file,si,sj,sk,byteorder):
    #dtype = np.dtype(np.uint16)  
    dtype = np.dtype(np.uint8)     
    dtype = dtype.newbyteorder(byteorder) 
    pixf = np.fromfile(file, dtype=dtype) 
    pixf.shape = (si,sj,sk)  
    return pixf
 
def surfaceToOffFile(verts, faces, filename):
    """ 
    Surface (vertices and faces) to offset file 
    """ 
    with open(filename+".off", "w") as f: 
        f.write('OFF\n')
        f.write(str(verts.shape[0])+' '+str(faces.shape[0])+' '+str(0)+'\n\n')
        # loop over nodes 
        for i in range(verts.shape[0]):
            f.write(str(verts[i,0])+' '+str(verts[i,1])+' '+str(verts[i,2])+'\n')
        # loop over triangles 
        for i in range(faces.shape[0]):
            f.write('3  '+str(faces[i,0])+' '+str(faces[i,1])+' '+str(faces[i,2])+'\n')
        f.write(' ')
        
def exportMeshToVtk(filename, file_format, points, cells, output_fields_p=None, output_fields_c=None):
        """ 
        output_fields_p : dictionary of point data 
        output_fields_c : dictionary of cell data 
        """          
        mesh = meshio.Mesh( points, cells ) 
        if output_fields_p is not None: 
            mesh.point_data = output_fields_p 
        if output_fields_c is not None:
            mesh.cell_data = output_fields_c 
        mesh.write(filename+'.'+file_format) 
        print(filename+'.'+file_format + ' written')
 
def group_duplicate_index_v2(a):
    """
    Taken from 
    https://stackoverflow.com/questions/46629518/find-indices-of-duplicate-rows-in-pandas-dataframe
    """
    s = (np.max(a)+1)**np.arange(a.shape[1])
    sidx = a.dot(s).argsort()
    b = a[sidx]

    m = np.concatenate(([False], (b[1:] == b[:-1]).all(1), [False] ))
    idx = np.flatnonzero(m[1:] != m[:-1])
    index = np.arange(a.shape[0])
    I = index[sidx].tolist() 
    return [I[i:j] for i,j in zip(idx[::2],idx[1::2]+1)]        
 
 
def kronecker(A,B):
    """ test for C++ translation """ 
    C = np.zeros((A.shape[0]*B.shape[0],A.shape[1]*B.shape[1]))
    for iA in range(A.shape[0]):
        for iB in range(B.shape[0]):
            for jA in range(A.shape[1]):
                for jB in range(B.shape[1]):
                    C[iA*B.shape[0]+iB,jA*B.shape[1]+jB] = A[iA,jA]*B[iB,jB]
    return C
 
 
                    
        
