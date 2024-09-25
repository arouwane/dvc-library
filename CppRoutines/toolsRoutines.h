#ifndef TOOLS_ROUTINES_H
#define TOOLS_ROUTINES_H 
#include <iostream>
#include <cmath>
#include "memory.h"
#include "tools.h" 
#include "Eigen/Core"
#include <omp.h>



void EvaluateBspline3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v);  
					   
void EvaluateBsplineAndGradient3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v, 
					   double* dvdx, int s_dvdx, 
					   double* dvdy, int s_dvdy,
					   double* dvdz, int s_dvdz );

void EvaluateBsplineStructured3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v);  
					   
void EvaluateBsplineAndGradientStructured3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v, 
					   double* dvdx, int s_dvdx, 
					   double* dvdy, int s_dvdy,
					   double* dvdz, int s_dvdz );		
					   
 

void LocatePointsInTetraFE_Mesh1(double* x, int sx, 
                                double* y, int sy, 
                                double* z, int sz, 
                                int* e, int sie, int sje, 
                                double* n, int sin, int sjn, 
                                int* ie, int size_ie,
                                double eps );  


void LocatePointsInTetraFE_Mesh2(double* x, int sx, 
                                double* y, int sy, 
                                double* z, int sz, 
                                int* e, int sie, int sje, 
                                double* n, int sin, int sjn, 
                                int* ie, int size_ie,
                                double eps );  

void LocatePointsInTetraFE_Mesh3(double* x, int sx, 
                                 double* y, int sy, 
                                 double* z, int sz, 
                                 int* face_indices, int s_faces_indices, 
                                 int* e, int sie, int sje, 
                                 double* n, int sin, int sjn, 
                                 int* connFaces, int si_connFaces, int sj_connFaces,
                                 int* ie, int size_ie,
                                 double eps); 
                                                           
                                				    					   					   
#endif 