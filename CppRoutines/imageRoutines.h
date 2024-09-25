#ifndef IMAGE_ROUTINES_H
#define IMAGE_ROUTINES_H 
#include <iostream>
#include <cmath>
#include "memory.h"
#include "tools.h" 
#include "toolsRoutines.h"
#include "Eigen/Core"
#include <omp.h>



void ComputeL2ProjLumpedCoefficients(double* image, int si, int sj, int sk, 
                           double xmin, double ymin, double zmin, double zmax,                     
                           double dx, double dy, double dz, 
    					   int deg_xi, int deg_eta,int deg_zeta,
					       double* knotXi,  int sknotXi, 
					       double* knotEta, int sknotEta ,
					       double* knotZeta, int sknotZeta,  
					       double* lsc, int s_lsc  );  
					       
void ComputeL2ProjLumpedCoefficientsSumFact(double* image, int si, int sj, int sk, 
                           double xmin, double ymin, double zmin, double zmax,                     
                           double dx, double dy, double dz, 
    					   int deg_xi, int deg_eta,int deg_zeta,
					       double* knotXi,  int sknotXi, 
					       double* knotEta, int sknotEta ,
					       double* knotZeta, int sknotZeta,  
					       double* lsc, int s_lsc  ); 					       
					       


void EvaluateTrilinearInterpolation(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
        					        double* v, int s_v );   

void EvaluateTrilinearInterpolationAndGradient(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v, int s_v, 
            					    double* dvdx, int s_dvdx, 
            					    double* dvdy, int s_dvdy,
            					    double* dvdz, int s_dvdz );	 
        					         
void EvaluateTrilinearInterpolationAndGradientStructured(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v, int s_v, 
            					    double* dvdx, int s_dvdx, 
            					    double* dvdy, int s_dvdy,
            					    double* dvdz, int s_dvdz );	 

void EvaluateCardBspline2(double* image, int si, int sj, int sk,
                          double xmin, double ymin, double zmin,
                          double dx, double dy, double dz,  
                          double* xc, int sxc, 
                          double* yc, int syc, 
                          double* zc, int szc, 
					      double* x, int sx,
			  		      double* y, int sy,
			  		      double* z, int sz,   
					      double* v, int s_v );                                   
void EvaluateCardBspline3(double* image, int si, int sj, int sk,
                          double xmin, double ymin, double zmin,
                          double dx, double dy, double dz,  
                          double* xc, int sxc, 
                          double* yc, int syc, 
                          double* zc, int szc, 
					      double* x, int sx,
			  		      double* y, int sy,
			  		      double* z, int sz,   
					      double* v, int s_v );     
					      
void EvaluateCardBsplineAndGradient2(double* image, int si, int sj, int sk,
                                               double xmin, double ymin, double zmin,
                                               double dx, double dy, double dz,  
                                               double* xc, int sxc, 
                                               double* yc, int syc, 
                                               double* zc, int szc, 
                    					       double* x, int sx,
                    			  		       double* y, int sy,
                    			  		       double* z, int sz,   
                    					       double* v, int s_v,
					                   		   double* dvdx, int s_dvdx, 
            					               double* dvdy, int s_dvdy,
            					               double* dvdz, int s_dvdz );	  
void EvaluateCardBsplineAndGradient3(double* image, int si, int sj, int sk,
                                               double xmin, double ymin, double zmin,
                                               double dx, double dy, double dz,  
                                               double* xc, int sxc, 
                                               double* yc, int syc, 
                                               double* zc, int szc, 
                    					       double* x, int sx,
                    			  		       double* y, int sy,
                    			  		       double* z, int sz,   
                    					       double* v, int s_v,
					                   		   double* dvdx, int s_dvdx, 
            					               double* dvdy, int s_dvdy,
            					               double* dvdz, int s_dvdz );					          					    

void EvaluateCardBsplineAndGradient2Structured(double* image, int si, int sj, int sk,
                                               double xmin, double ymin, double zmin,
                                               double dx, double dy, double dz,  
                                               double* xc, int sxc, 
                                               double* yc, int syc, 
                                               double* zc, int szc, 
                    					       double* x, int sx,
                    			  		       double* y, int sy,
                    			  		       double* z, int sz,   
                    					       double* v, int s_v,
					                   		   double* dvdx, int s_dvdx, 
            					               double* dvdy, int s_dvdy,
            					               double* dvdz, int s_dvdz );	  
void EvaluateCardBsplineAndGradient3Structured(double* image, int si, int sj, int sk,
                                               double xmin, double ymin, double zmin,
                                               double dx, double dy, double dz,  
                                               double* xc, int sxc, 
                                               double* yc, int syc, 
                                               double* zc, int szc, 
                    					       double* x, int sx,
                    			  		       double* y, int sy,
                    			  		       double* z, int sz,   
                    					       double* v, int s_v,
					                   		   double* dvdx, int s_dvdx, 
            					               double* dvdy, int s_dvdy,
            					               double* dvdz, int s_dvdz );
 
					               

void GetMeanImageAndStdOnMesh_TrilinearInterp(double* f, int sif, int sjf, int skf,
                                              double* knotXiImage, int sknotXiImage, 
                                              double* knotEtaImage, int sknotEtaImage, 
                                              double* knotZetaImage, int sknotZetaImage, 
                                              int deg_xi, int deg_eta, int deg_zeta, 
                                              double* knotXi, int sknotXi, 
                                              double* knotEta, int sknotEta, 
                                              double* knotZeta, int sknotZeta,
                                              int nipex, int nipey, int nipez,
                                              double* xg, int sxg, 
                                              double* yg, int syg, 
                                              double* zg, int szg,                                               
                                              double* fmean, int sfmean, 
                                              double* fstd, int sfstd, 
                                              double* dyne, int sdyne, 
                                              double* fip, int sfip, 
                                              double* dfipdx, int sdfipdx, 
                                              double* dfipdy, int sdfipdy, 
                                              double* dfipdz, int sdfipdz); 
 
                                         
void GetMeanImageAndStdOnMesh_CBspline2(double* f, int sif, int sjf, int skf,
                                        double xmin, double ymin, double zmin,
                                        double dx, double dy, double dz,  
                                        double* xc, int sxc, 
                                        double* yc, int syc, 
                                        double* zc, int szc, 
                                        int deg_xi, int deg_eta, int deg_zeta, 
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta,
                                        int nipex, int nipey, int nipez,
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg,                                               
                                        double* fmean, int sfmean, 
                                        double* fstd, int sfstd, 
                                        double* dyne, int sdyne, 
                                        double* fip, int sfip, 
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz);                                         

void GetMeanImageAndStdOnMesh_CBspline3(double* f, int sif, int sjf, int skf,
                                        double xmin, double ymin, double zmin,
                                        double dx, double dy, double dz,  
                                        double* xc, int sxc, 
                                        double* yc, int syc, 
                                        double* zc, int szc, 
                                        int deg_xi, int deg_eta, int deg_zeta, 
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta,
                                        int nipex, int nipey, int nipez,
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg,                                               
                                        double* fmean, int sfmean, 
                                        double* fstd, int sfstd, 
                                        double* dyne, int sdyne, 
                                        double* fip, int sfip, 
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz);     
                                        
                                              
                                              
void GetMeanImageAndStdOnMesh_thrsh_TrilinearInterp(double* f, int sif, int sjf, int skf,
                                              double thrsh, 
                                              double* knotXiImage, int sknotXiImage, 
                                              double* knotEtaImage, int sknotEtaImage, 
                                              double* knotZetaImage, int sknotZetaImage, 
                                              int deg_xi, int deg_eta, int deg_zeta, 
                                              double* knotXi, int sknotXi, 
                                              double* knotEta, int sknotEta, 
                                              double* knotZeta, int sknotZeta,
                                              int nipex, int nipey, int nipez,
                                              double* xg, int sxg, 
                                              double* yg, int syg, 
                                              double* zg, int szg,                                               
                                              double* fmean, int sfmean, 
                                              double* fstd, int sfstd, 
                                              double* fip, int sfip, 
                                              double* dfipdx, int sdfipdx, 
                                              double* dfipdy, int sdfipdy, 
                                              double* dfipdz, int sdfipdz); 
                                              
                                              
void GetC8MeshFromVoxelsTrilinearInterp(double* image, int si, int sj, int sk, 
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double thrsh, 
                                        double* knotXi, int sknotXi,
                                        double* knotEta, int sknotEta,  
                                        double* knotZeta, int sknotZeta,  
                                        int* elementMask, int se );                                                   
                                				    					   


                                                

void GetMeanImageAndStdOnFE_Mesh_TrilinearInterp(double* f, int sif, int sjf, int skf,
                                              double* knotXiImage, int sknotXiImage, 
                                              double* knotEtaImage, int sknotEtaImage, 
                                              double* knotZetaImage, int sknotZetaImage, 
                                              int* e, int sie, int sje,
                                              double* n, int sin, int sjn, 
                                              double* N, int siN, int sjN,                                               
                                              double* fmean, int sfmean, 
                                              double* fstd, int sfstd, 
                                              double* dyne, int sdyne, 
                                              double* fip, int sfip, 
                                              double* dfipdx, int sdfipdx, 
                                              double* dfipdy, int sdfipdy, 
                                              double* dfipdz, int sdfipdz); 

void GetMeanImageAndStdOnFE_Mesh_CBspline3(double* f, int sif, int sjf, int skf,
                                              double xmin, double ymin, double zmin, 
                                              double dx, double dy, double dz,
                                              double* xc, int sxc, 
                                              double* yc, int syc, 
                                              double* zc, int szc,   
                                              int* e, int sie, int sje,
                                              double* n, int sin, int sjn, 
                                              double* N, int siN, int sjN,                                               
                                              double* fmean, int sfmean, 
                                              double* fstd, int sfstd, 
                                              double* dyne, int sdyne, 
                                              double* fip, int sfip, 
                                              double* dfipdx, int sdfipdx, 
                                              double* dfipdy, int sdfipdy,       
                                              double* dfipdz, int sdfipdz);                            
                                                                                    
                                          
void GetMeanImageAndStdOnFE_Mesh_L2ProjLumped(double* cp, int scp,
                                                double* knotXiImage, int sknotXiImage, 
                                                double* knotEtaImage, int sknotEtaImage, 
                                                double* knotZetaImage, int sknotZetaImage, 
                                                int deg_xi, int deg_eta, int deg_zeta, 
                                                int* e, int sie, int sje,
                                                double* n, int sin, int sjn, 
                                                double* N, int siN, int sjN,                                               
                                                double* fmean, int sfmean, 
                                                double* fstd, int sfstd, 
                                                double* dyne, int sdyne, 
                                                double* fip, int sfip, 
                                                double* dfipdx, int sdfipdx, 
                                                double* dfipdy, int sdfipdy, 
                                                double* dfipdz, int sdfipdz);  
                                                
                                               
 					   
#endif 
