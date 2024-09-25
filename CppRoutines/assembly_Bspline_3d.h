#ifndef ASSEMBLY_BSPLINE_3D_H
#define ASSEMBLY_BSPLINE_3D_H 
#include <iostream>
#include <cmath>
#include "memory.h"
#include "tools.h" 
#include "meshUtils.h"
#include "octreeDecomposition.h"
#include "toolsRoutines.h"
#include "Eigen/Core"
#include <omp.h>

void ComputeIntegratedDiffTensors( double* knotXi,  int sknotXi,
                                   double* knotEta, int sknotEta,
                                   double* knotZeta, int sknotZeta,
                                   int deg_xi, int deg_eta, int deg_zeta ,
                                   int nipx, int nipy, int nipz,
                                   int nex, int ney, int nez, 
                                   double* dxi_dxi, int s_dxi_dxi, 
                                   double* dxi_xi,  int s_dxi_xi, 
                                   double* xi_dxi,  int s_xi_dxi, 
                                   double* xi_xi,   int s_xi_xi, 
                                   double* deta_deta, int s_deta_deta,   
                                   double* deta_eta, int s_deta_eta,   
                                   double* eta_deta, int s_eta_deta,    
                                   double* eta_eta, int s_eta_eta,   
                                   double* dzeta_dzeta, int s_dzeta_dzeta,  
                                   double* dzeta_zeta, int s_dzeta_zeta,        
                                   double* zeta_dzeta, int s_zeta_dzeta,      
                                   double* zeta_zeta, int s_zeta_zeta,  
                                   double* eta_eta_dxi_dxi, int s_eta_eta_dxi_dxi,  
                                   double* eta_deta_dxi_xi, int s_eta_deta_dxi_xi,  
                                   double* deta_eta_xi_dxi, int s_deta_eta_xi_dxi,   
                                   double* deta_deta_xi_xi, int s_deta_deta_xi_xi );  

void Laplacian(int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz); 

                                    


 

#endif