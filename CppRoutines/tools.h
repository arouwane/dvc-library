#ifndef TOOLS_H
#define TOOLS_H

#include "memory.h"
#include "Eigen/Core"
#include <stdlib.h>
#include <cmath>
#include <vector>


double GetGrayLevelOfPoint3D(double* image, int si, int sj, int sk,
                           double x, double y, double z, 
                           double xmin, double ymin, double zmin, double zmax,
                           double dx, double dy, double dz ); 

// Material properties
Eigen::MatrixXd hookeVolume(double E, double v);


// B-splines routines
int findspanUniformKnotVector(double* U, int sU, int deg, double l, double u);
int  findspan(int n, int p, double u, double *U);
void basisfuns(int i, double u, int p, double *U, double *N);
void dersbasisfuns(int d, double *k, double u, int s,int n, double **ders);
void bezierExtraction(double *U, int p, int m, int n_elems, double*** C);
double oneBasisFun(int i, double u, int p, double*U, int lenU);
Eigen::MatrixXi getNoelem(int nxi, int neta, int nzeta, int p, int q, int r);
 
 
// Gauss Integration routines
void gauleg(double x1, double x2, double x[], double w[], int n);
int GaussTriangle(std::vector<double> &xg, std::vector<double> &yg, std::vector<double> &wg, int n); 
Eigen::MatrixXd GaussTetrahedron(std::vector<double> &wg, int &npoints, int n);  


// B-spline basis function: analytical formula centered at x = 0
double CardinalBspline2(double x);  
double CardinalBspline3(double x); 
double CardinalBsplineDer2(double x); 
double CardinalBsplineDer3(double x);

double Evaluate3dBsplineOnOnePoint(double x, double y, double z,
                                   double* knotXi,  int sknotXi, 
					               double* knotEta, int sknotEta ,
					               double* knotZeta, int sknotZeta, 
					               int deg_xi, int deg_eta,int deg_zeta,
					               double* cp, int s_cp ); 



///////////////////////////  IMAGE TRILINEAR INTERPOLATION ROOTINES //////////////////////////////////////

double EvaluateTrilinearInterpolationOnOnePoint(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta,
        					        double* knotZeta, int sknotZeta,
					                double x, double y, double z ); 

 
void EvaluateTrilinearInterpolationAndGradientStructured(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v, 
            					    double* dvdx, 
            					    double* dvdy,
            					    double* dvdz); 
void EvaluateTrilinearInterpolationStructured(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v );             					    

int CheckCutCuboidTrilinearInterpolation(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
        					        double thrsh, 
        					        Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z); 	
        					        		        
///////////////////////////////////////////////////////////////////////////////////////////////


void interpolateLinearly(double xmin, double xmax, double lmin, double lmax, double &a, double &b);



///////////////////////////  IMAGE CBSPLINE2 INTERPOLATION ROOTINES //////////////////////////////////////

// Image evaluation with quadratic cardinal B-splines 

double EvaluateCBspline2OnOnePoint(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double xmax, double ymax, double zmax, 
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double x, double y, double z ); 

 
void EvaluateCBspline2AndGradientStructured(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double* x, int sx,
            			  		   double* y, int sy,
            			  		   double* z, int sz,   
            					   double* v, 
            					   double* dvdx, 
            					   double* dvdy,
            					   double* dvdz);             					            					        		        
///////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////  IMAGE CBSPLINE3 INTERPOLATION ROOTINES //////////////////////////////////////

// Image evaluation with quadratic cardinal B-splines 

double EvaluateCBspline3OnOnePoint(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double xmax, double ymax, double zmax, 
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double x, double y, double z ); 

 
void EvaluateCBspline3AndGradientStructured(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double* x, int sx,
            			  		   double* y, int sy,
            			  		   double* z, int sz,   
            					   double* v, 
            					   double* dvdx, 
            					   double* dvdy,
            					   double* dvdz);             					            					        		        
///////////////////////////////////////////////////////////////////////////////////////////////

 
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////


int PointIsInTetrahedron1(double* xn, double* yn, double*zn, double x, double y, double z, double eps ); 
int PointIsInTetrahedron2(double* xn, double* yn, double*zn, double x, double y, double z, double eps );

 
 
 

 
#endif
