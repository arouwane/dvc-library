#ifndef ASSEMBLY_ROUTINES_H
#define ASSEMBLY_ROUTINES_H 
#include <iostream>
#include <cmath>
#include "memory.h"
#include "tools.h" 
#include "meshUtils.h"
#include "octreeDecomposition.h"
#include "toolsRoutines.h"
#include "Eigen/Core"
#include <omp.h>


void GetBsplineFunctionsMatrixStructured(double* x, int sx,
                                         double* y, int sy, 
                                         double* z, int sz,
                                         double* knotXi, int sknotXi,
                                         double* knotEta, int sknotEta,
                                         double* knotZeta, int sknotZeta,
                                         int deg_xi, int deg_eta, int deg_zeta, 
                                         double* valuesN, int s_valuesN,
                                         int* indexI, int s_indexI,
                                         int* indexJ, int s_indexJ);  
                                         
void GetBsplineFunctionsAndDerivativesMatrixStructured(double* x, int sx,
                                         double* y, int sy, 
                                         double* z, int sz,
                                         double* knotXi, int sknotXi,
                                         double* knotEta, int sknotEta,
                                         double* knotZeta, int sknotZeta,
                                         int deg_xi, int deg_eta, int deg_zeta, 
                                         double* valuesN, int s_valuesN, 
                                         double* valuesdNdx, int s_valuesdNdx,
                                         double* valuesdNdy, int s_valuesdNdy, 
                                         double* valuesdNdz, int s_valuesdNdz, 
                                         int* indexI, int s_indexI,
                                         int* indexJ, int s_indexJ);  

std::vector<double> FcmIntegrationTrilinearInterp(double* image, int si, int sj, int sk,
                                                  double* knotXiImage, int sknotXiImage, 
                                                  double* knotEtaImage, int sknotEtaImage, 
                                                  double* knotZetaImage, int sknotZetaImage, 
                                                  double thrsh, 
                                                  int deg_xi, int deg_eta, int deg_zeta, 
                                                  double* knotXi, int sknotXi, 
                                                  double* knotEta, int sknotEta, 
                                                  double* knotZeta, int sknotZeta, 
                                                  int lvlmax ); 

void FcmTrilinearInterpStiffnessAndBfIntegral(double* image, int si, int sj, int sk,
                                             double* knotXiImage, int sknotXiImage, 
                                             double* knotEtaImage, int sknotEtaImage, 
                                             double* knotZetaImage, int sknotZetaImage, 
                                             double thrsh, 
                                             int deg_xi, int deg_eta, int deg_zeta, 
                                             double* knotXi, int sknotXi,
                                             double* knotEta, int sknotEta,  
                                             double* knotZeta, int sknotZeta,
                                             int lvlmax, 
                                             double E, double nu, 
                                             int* indexI, int nnzI, 
                                             int* indexJ, int nnzJ,
                                             double* nnz_values, int nnz, 
                                             double* intBf, int sintBf ) ;                                                   


void FcmTrilinearInterpStiffnessAndBfIntegralParallel(double* image, int si, int sj, int sk,
                                             double* knotXiImage, int sknotXiImage, 
                                             double* knotEtaImage, int sknotEtaImage, 
                                             double* knotZetaImage, int sknotZetaImage, 
                                             double thrsh, 
                                             int deg_xi, int deg_eta, int deg_zeta, 
                                             double* knotXi, int sknotXi,
                                             double* knotEta, int sknotEta,  
                                             double* knotZeta, int sknotZeta,
                                             int lvlmax, 
                                             double E, double nu, 
                                             int* indexI, int nnzI, 
                                             int* indexJ, int nnzJ,
                                             double* nnz_values, int nnz, 
                                             double* intBf, int sintBf ) ; 
                                             

void Laplacian_Structured(int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz);
                          
void Laplacian_Structured_Parallel(int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz);                          

/*
void HomogeneousStiffness(double E, double nu, 
                          int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz);    */ 
                          

void L2Projection(int deg_xi, int deg_eta, int deg_zeta,
                  double* knotXi, int sknotXi,
                  double* knotEta, int sknotEta,
                  double* knotZeta, int sknotZeta,
                  int deg_xi2, int deg_eta2, int deg_zeta2,
                  double* knotXi2, int sknotXi2,
                  double* knotEta2, int sknotEta2,
                  double* knotZeta2, int sknotZeta2,
                  double* U2, int sU2, 
                  int* indexI, int nnzI, 
                  int* indexJ, int nnzJ,
                  double* nnz_values, int nnz,
                  double* rhs, int ndof);                                      


                                                  
std::vector<int> VoxelIntegrationThreshold(double* fip, int sfip,
                                                             double thrsh, 
                                                             int nipex, int nipey, int nipez, 
                                                             int deg_xi, int deg_eta, int deg_zeta, 
                                                             double* knotXi, int sknotXi, 
                                                             double* knotEta, int sknotEta, 
                                                             double* knotZeta, int sknotZeta ); 


std::vector<int> VoxelIntegrationMask(double* maskip, int smaskip,
                                                             int nipex, int nipey, int nipez, 
                                                             int deg_xi, int deg_eta, int deg_zeta, 
                                                             double* knotXi, int sknotXi, 
                                                             double* knotEta, int sknotEta, 
                                                             double* knotZeta, int sknotZeta ); 
                                                                                                                          
                                                             
                                                                    


void DVC_LHS_thrsh(double* dfipdx, int sdfipdx, 
                    double* dfipdy, int sdfipdy, 
                    double* dfipdz, int sdfipdz, 
                    int* ipIndices, int sipIndices, 
                    int deg_xi, int deg_eta, int deg_zeta, 
                    double* knotXi, int sknotXi, 
                    double* knotEta, int sknotEta, 
                    double* knotZeta, int sknotZeta, 
                    int nipex, int nipey, int nipez , 
                    double* Nxi, int siNxi, int sjNxi, 
                    double* Neta, int siNeta, int sjNeta, 
                    double* Nzeta, int siNzeta, int sjNzeta, 
                    int* indexI, int nnzI, 
                    int* indexJ, int nnzJ, 
                    double* nnz_values, int nnz);                                                                                                                


double DVC_RHS_thrsh_TrilinearInterp(double* g, int sgi, int sgj,int sgk, 
                                   double* knotXiImage, int sknotXiImage, 
                                   double* knotEtaImage, int sknotEtaImage, 
                                   double* knotZetaImage, int sknotZetaImage,                                    
                                   double* fip,    int sfip,                                         
                                   double* dfipdx, int sdfipdx, 
                                   double* dfipdy, int sdfipdy, 
                                   double* dfipdz, int sdfipdz, 
                                   int* ipIndices, int sipIndices, 
                                   int deg_xi, int deg_eta, int deg_zeta,
                                   double* knotXi, int sknotXi, 
                                   double* knotEta, int sknotEta, 
                                   double* knotZeta, int sknotZeta, 
                                   int* NOELEM, int siNOELEM, int sjNOELEM,  
                                   int nipex, int nipey, int nipez ,
                                   double* xg, int sxg, 
                                   double* yg, int syg, 
                                   double* zg, int szg,                                    
                                   double* Nxi, int siNxi, int sjNxi, 
                                   double* Neta, int siNeta, int sjNeta, 
                                   double* Nzeta, int siNzeta, int sjNzeta,
                                   double* U, int sU, 
                                   double* rhs, int ndof ); 



double DVC_RHS_ZN_thrsh_TrilinearInterp(double* g, int sgi, int sgj,int sgk, 
                                   double* knotXiImage, int sknotXiImage, 
                                   double* knotEtaImage, int sknotEtaImage, 
                                   double* knotZetaImage, int sknotZetaImage,                                    
                                   double* fip,    int sfip,                                         
                                   double* dfipdx, int sdfipdx, 
                                   double* dfipdy, int sdfipdy, 
                                   double* dfipdz, int sdfipdz, 
                                   double* fmeane, int sfmeane, 
                                   double* fstde,  int sfstde, 
                                   int* ipIndices, int sipIndices, 
                                   int deg_xi, int deg_eta, int deg_zeta,
                                   double* knotXi, int sknotXi, 
                                   double* knotEta, int sknotEta, 
                                   double* knotZeta, int sknotZeta, 
                                   int* NOELEM, int siNOELEM, int sjNOELEM,  
                                   int nipex, int nipey, int nipez ,
                                   double* xg, int sxg, 
                                   double* yg, int syg, 
                                   double* zg, int szg,                                    
                                   double* Nxi, int siNxi, int sjNxi, 
                                   double* Neta, int siNeta, int sjNeta, 
                                   double* Nzeta, int siNzeta, int sjNzeta,
                                   double* U, int sU, 
                                   double* rhs, int ndof ); 
                                   
                                   
                                   

void DVC_LHS_Structured(double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz, 
                                        int deg_xi, int deg_eta, int deg_zeta, 
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int nipex, int nipey, int nipez , 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        int* indexI, int nnzI, 
                                        int* indexJ, int nnzJ, 
                                        double* nnz_values, int nnz); 

void DVC_LHS_Structured_Parallel(double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz, 
                                        int deg_xi, int deg_eta, int deg_zeta, 
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int nipex, int nipey, int nipez , 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        int* indexI, int nnzI, 
                                        int* indexJ, int nnzJ, 
                                        double* nnz_values, int nnz); 


                                       
double DVC_RHS_Structured_TrilinearInterp(double* g, int sgi, int sgj,int sgk,
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double* fip,    int sfip,                                         
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz, 
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int* NOELEM, int siNOELEM, int sjNOELEM,  
                                        int nipex, int nipey, int nipez , 
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg, 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        double* U, int sU, 
                                        double* rhs, int ndof );  

double DVC_RHS_Structured_TrilinearInterp_Parallel(double* g, int sgi, int sgj,int sgk,
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double* fip,    int sfip,                                         
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz, 
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int* NOELEM, int siNOELEM, int sjNOELEM, 
                                        int nipex, int nipey, int nipez , 
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg, 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        double* U, int sU, 
                                        double* rhs, int ndof );  
                                                                    
 

void DVC_RHS_ZN_Structured_TrilinearInterp(double* g, int sgi, int sgj,int sgk,
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double* fip,    int sfip, 
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz,
                                        double* fmeane, int sfmeane, 
                                        double* fstde,  int sfstde, 
                                        double* dyne,   int sdyne,   
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta,
                                        int* NOELEM, int siNOELEM, int sjNOELEM,   
                                        int nipex, int nipey, int nipez ,
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg, 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        double* U, int sU, 
                                        double* rhs, int ndof,
                                        double* elementRes, int selementRes );     


void GLR_Structured_TrilinearInterp( double* fip, int sfip, 
                                     double* g, int sgi, int sgj,int sgk,
                                     double* knotXiImage, int sknotXiImage, 
                                     double* knotEtaImage, int sknotEtaImage, 
                                     double* knotZetaImage, int sknotZetaImage, 
                                     int deg_xi, int deg_eta, int deg_zeta,
                                     double* knotXi, int sknotXi, 
                                     double* knotEta, int sknotEta, 
                                     double* knotZeta, int sknotZeta, 
                                     int nipex, int nipey, int nipez ,
                                     double* xg, int sxg, 
                                     double* yg, int syg, 
                                     double* zg, int szg, 
                                     double* Nxi, int siNxi, int sjNxi, 
                                     double* Neta, int siNeta, int sjNeta, 
                                     double* Nzeta, int siNzeta, int sjNzeta, 
                                     double* U, int sU, 
                                     double* glr, int sglr ) ; 


void GLR_Structured_TrilinearInterp_Parallel( double* fip, int sfip, 
                                     double* g, int sgi, int sgj,int sgk,
                                     double* knotXiImage, int sknotXiImage, 
                                     double* knotEtaImage, int sknotEtaImage, 
                                     double* knotZetaImage, int sknotZetaImage, 
                                     int deg_xi, int deg_eta, int deg_zeta,
                                     double* knotXi, int sknotXi, 
                                     double* knotEta, int sknotEta, 
                                     double* knotZeta, int sknotZeta, 
                                     int nipex, int nipey, int nipez ,
                                     double* xg, int sxg, 
                                     double* yg, int syg, 
                                     double* zg, int szg, 
                                     double* Nxi, int siNxi, int sjNxi, 
                                     double* Neta, int siNeta, int sjNeta, 
                                     double* Nzeta, int siNzeta, int sjNzeta, 
                                     double* U, int sU, 
                                     double* glr, int sglr ) ; 
                                     
                                        
void GLR_ZN_Structured_TrilinearInterp( double* fip, int sfip, 
                                        double* fmeane, int sfmeane, 
                                        double* fstde,  int sfstde,  
                                        double* g, int sgi, int sgj,int sgk,
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int nipex, int nipey, int nipez ,
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg, 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        double* U, int sU, 
                                        double* glr, int sglr );  
 
                                        
                                        
                                 
 
                                       
/*                                           
void DVC_RHS_Structured_TrilinearInterp(double* g, int sgi, int sgj,int sgk,
                                        double xminIm, double xmaxIm, double yminIm, double ymaxIm, double zminIm, double zmaxIm, 
                                        double dx, double dy, double dz, 
                                        double* fip,    int sfip, 
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz, 
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta, 
                                        int nipex, int nipey, int nipez , 
                                        double* U, int sU, 
                                        double* rhs, int ndof ); */
 
                                       
double DVC_RHS_Structured_CBspline2(double* g, int sgi, int sgj,int sgk,
                                    double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                                    double dx, double dy, double dz, 
                                    double* xc, int sxc, 
                                    double* yc, int syc, 
                                    double* zc, int szc,            
                                    double* fip,    int sfip,                                         
                                    double* dfipdx, int sdfipdx, 
                                    double* dfipdy, int sdfipdy, 
                                    double* dfipdz, int sdfipdz, 
                                    int deg_xi, int deg_eta, int deg_zeta,
                                    double* knotXi, int sknotXi, 
                                    double* knotEta, int sknotEta, 
                                    double* knotZeta, int sknotZeta, 
                                    int* NOELEM, int siNOELEM, int sjNOELEM,  
                                    int nipex, int nipey, int nipez , 
                                    double* xg, int sxg, 
                                    double* yg, int syg, 
                                    double* zg, int szg, 
                                    double* Nxi, int siNxi, int sjNxi, 
                                    double* Neta, int siNeta, int sjNeta, 
                                    double* Nzeta, int siNzeta, int sjNzeta, 
                                    double* U, int sU, 
                                    double* rhs, int ndof );  
 

void GLR_Structured_CBspline2( double* fip, int sfip, 
                               double* g, int sgi, int sgj,int sgk,
                                double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                                double dx, double dy, double dz, 
                                double* xc, int sxc, 
                                double* yc, int syc, 
                                double* zc, int szc,            
                                int deg_xi, int deg_eta, int deg_zeta,
                                double* knotXi, int sknotXi, 
                                double* knotEta, int sknotEta, 
                                double* knotZeta, int sknotZeta, 
                                int nipex, int nipey, int nipez ,
                                double* xg, int sxg, 
                                double* yg, int syg, 
                                double* zg, int szg, 
                                double* Nxi, int siNxi, int sjNxi, 
                                double* Neta, int siNeta, int sjNeta, 
                                double* Nzeta, int siNzeta, int sjNzeta, 
                                double* U, int sU, 
                                double* glr, int sglr ); 
                                
double DVC_RHS_Structured_CBspline3(double* g, int sgi, int sgj,int sgk,
                                    double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                                    double dx, double dy, double dz, 
                                    double* xc, int sxc, 
                                    double* yc, int syc, 
                                    double* zc, int szc,            
                                    double* fip,    int sfip,                                         
                                    double* dfipdx, int sdfipdx, 
                                    double* dfipdy, int sdfipdy, 
                                    double* dfipdz, int sdfipdz, 
                                    int deg_xi, int deg_eta, int deg_zeta,
                                    double* knotXi, int sknotXi, 
                                    double* knotEta, int sknotEta, 
                                    double* knotZeta, int sknotZeta, 
                                    int* NOELEM, int siNOELEM, int sjNOELEM,  
                                    int nipex, int nipey, int nipez , 
                                    double* xg, int sxg, 
                                    double* yg, int syg, 
                                    double* zg, int szg, 
                                    double* Nxi, int siNxi, int sjNxi, 
                                    double* Neta, int siNeta, int sjNeta, 
                                    double* Nzeta, int siNzeta, int sjNzeta, 
                                    double* U, int sU, 
                                    double* rhs, int ndof );  
 

void GLR_Structured_CBspline3( double* fip, int sfip, 
                               double* g, int sgi, int sgj,int sgk,
                                double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                                double dx, double dy, double dz, 
                                double* xc, int sxc, 
                                double* yc, int syc, 
                                double* zc, int szc,            
                                int deg_xi, int deg_eta, int deg_zeta,
                                double* knotXi, int sknotXi, 
                                double* knotEta, int sknotEta, 
                                double* knotZeta, int sknotZeta, 
                                int nipex, int nipey, int nipez ,
                                double* xg, int sxg, 
                                double* yg, int syg, 
                                double* zg, int szg, 
                                double* Nxi, int siNxi, int sjNxi, 
                                double* Neta, int siNeta, int sjNeta, 
                                double* Nzeta, int siNzeta, int sjNzeta, 
                                double* U, int sU, 
                                double* glr, int sglr );                                 

void DVC_LHS_FE(int* e, int sie, int sje,
                   double* n, int sin, int sjn,
                   int* conn, int siconn, int sjconn,  
                   double* N, int siN, int sjN, 
                   double* dNdxi, int sidNdxi, int sjdNdxi, 
                   double* dNdeta, int sidNdeta, int sjdNdeta, 
                   double* dNdzeta, int sidNdzeta, int sjdNdzeta,
                   double* w, int sw, 
                   double* dfipdx, int sdfipdx, 
                   double* dfipdy, int sdfipdy, 
                   double* dfipdz, int sdfipdz, 
                   int* indexI, int nnzI, 
                   int* indexJ, int nnzJ, 
                   double* nnz_values, int nnz); 

double DVC_RHS_FE_TrilinearInterp(double* g, int sgi, int sgj, int sgk,
                double* knotXiImage, int sknotXiImage, 
                double* knotEtaImage, int sknotEtaImage, 
                double* knotZetaImage, int sknotZetaImage, 
                double* fip,    int sfip,                                         
                double* dfipdx, int sdfipdx, 
                double* dfipdy, int sdfipdy, 
                double* dfipdz, int sdfipdz, 
                int* e, int sie, int sje,
                double* n, int sin, int sjn, 
                int* conn, int siconn, int sjconn, 
                double* N, int siN, int sjN, 
                double* dNdxi, int sidNdxi, int sjdNdxi, 
                double* dNdeta, int sidNdeta, int sjdNdeta, 
                double* dNdzeta, int sidNdzeta, int sjdNdzeta,
                double* w, int sw, 
                double* U, int sU, 
                double* rhs, int ndof ) ;

void GLR_FE_TrilinearInterp(double* g, int sgi, int sgj,int sgk,
            double* knotXiImage, int sknotXiImage, 
            double* knotEtaImage, int sknotEtaImage, 
            double* knotZetaImage, int sknotZetaImage, 
            double* fip, int sfip, 
            int* e, int sie, int sje,
            double* n, int sin, int sjn, 
            int* conn, int siconn, int sjconn, 
            double* N, int siN, int sjN, 
            double* w, int sw, 
            double* U, int sU, 
            double* glr, int sglr ) ; 

void Gophi_FEMesh_TrilinearInterp(double*g, int sgi, int sgj, int sgk, 
                                  double* knotXiImage, int sknotXiImage, 
                                  double* knotEtaImage, int sknotEtaImage, 
                                  double* knotZetaImage, int sknotZetaImage, 
                                  double* xi, int sxi, 
                                  double* eta, int seta, 
                                  double* zeta, int szeta, 
                                  int* ie, int size_ie, 
                                  int* e, int sie, int sje, 
                                  double*n, int sin, int sjn, 
                                  int* conn, int siconn, int sjconn, 
                                  double* U, int sU, 
                                  double* glr, int sglr );   
                                  

 

double DVC_RHS_FE_CBspline3(double* g, int sgi, int sgj, int sgk,
                            double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                            double dx, double dy, double dz, 
                            double* xc, int sxc, 
                            double* yc, int syc, 
                            double* zc, int szc, 
                            double* fip,    int sfip,                                         
                            double* dfipdx, int sdfipdx, 
                            double* dfipdy, int sdfipdy, 
                            double* dfipdz, int sdfipdz, 
                            int* e, int sie, int sje,
                            double* n, int sin, int sjn, 
                            int* conn, int siconn, int sjconn, 
                            double* N, int siN, int sjN, 
                            double* dNdxi, int sidNdxi, int sjdNdxi, 
                            double* dNdeta, int sidNdeta, int sjdNdeta, 
                            double* dNdzeta, int sidNdzeta, int sjdNdzeta,
                            double* w, int sw, 
                            double* U, int sU, 
                            double* rhs, int ndof ) ;

void GLR_FE_CBspline3(double* g, int sgi, int sgj, int sgk,
                            double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                            double dx, double dy, double dz, 
                            double* xc, int sxc, 
                            double* yc, int syc, 
                            double* zc, int szc, 
                            double* fip, int sfip, 
                            int* e, int sie, int sje,
                            double* n, int sin, int sjn, 
                            int* conn, int siconn, int sjconn, 
                            double* N, int siN, int sjN, 
                            double* w, int sw, 
                            double* U, int sU, 
                            double* glr, int sglr ) ; 

void Gophi_FEMesh_CBspline3(double* g, int sgi, int sgj, int sgk,
                            double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
                            double dx, double dy, double dz, 
                            double* xc, int sxc, 
                            double* yc, int syc, 
                            double* zc, int szc, 
                            double* xi, int sxi, 
                            double* eta, int seta, 
                            double* zeta, int szeta, 
                            int* ie, int size_ie, 
                            int* e, int sie, int sje, 
                            double*n, int sin, int sjn, 
                            int* conn, int siconn, int sjconn, 
                            double* U, int sU, 
                            double* glr, int sglr );  



double DVC_RHS_FE_L2ProjLumped(double* lsc_g, int s_lsc_g,
                               double* knotXiImage_g, int sknotXiImage_g, 
                               double* knotEtaImage_g, int sknotEtaImage_g, 
                               double* knotZetaImage_g, int sknotZetaImage_g, 
                               int deg_xi_g, int deg_eta_g, int deg_zeta_g,
                               double* fip,    int sfip,                                         
                               double* dfipdx, int sdfipdx, 
                               double* dfipdy, int sdfipdy, 
                               double* dfipdz, int sdfipdz, 
                               int* e, int sie, int sje,
                               double* n, int sin, int sjn, 
                               int* conn, int siconn, int sjconn, 
                               double* N, int siN, int sjN, 
                               double* dNdxi, int sidNdxi, int sjdNdxi, 
                               double* dNdeta, int sidNdeta, int sjdNdeta, 
                               double* dNdzeta, int sidNdzeta, int sjdNdzeta,
                               double* w, int sw, 
                               double* U, int sU, 
                               double* rhs, int ndof ) ;
                               

void GLR_FE_L2ProjLumped(double* lsc_g, int s_lsc_g,
                         double* knotXiImage_g, int sknotXiImage_g, 
                         double* knotEtaImage_g, int sknotEtaImage_g, 
                         double* knotZetaImage_g, int sknotZetaImage_g, 
                         int deg_xi_g, int deg_eta_g, int deg_zeta_g,
                         double* fip, int sfip, 
                         int* e, int sie, int sje,
                         double* n, int sin, int sjn, 
                         int* conn, int siconn, int sjconn, 
                         double* N, int siN, int sjN, 
                         double* w, int sw, 
                         double* U, int sU, 
                         double* glr, int sglr ) ; 

void Gophi_FEMesh_L2ProjLumped(double* lsc_g, int s_lsc_g,
                               double* knotXiImage_g, int sknotXiImage_g, 
                               double* knotEtaImage_g, int sknotEtaImage_g, 
                               double* knotZetaImage_g, int sknotZetaImage_g, 
                               int deg_xi_g, int deg_eta_g, int deg_zeta_g,
                               double* xi, int sxi, 
                               double* eta, int seta, 
                               double* zeta, int szeta, 
                               int* ie, int size_ie, 
                               int* e, int sie, int sje, 
                               double*n, int sin, int sjn, 
                               int* conn, int siconn, int sjconn, 
                               double* U, int sU, 
                               double* glr, int sglr );  
                                                           
                          
                                               
               
                  
void Laplacian_FE(int* e, int sie, int sje,
                    double* n, int sin, int sjn, 
                    int* conn, int siconn, int sjconn, 
                    double* N, int siN, int sjN, 
                    double* dNdxi, int sidNdxi, int sjdNdxi, 
                    double* dNdeta, int sidNdeta, int sjdNdeta, 
                    double* dNdzeta, int sidNdzeta, int sjdNdzeta, 
                    double* w, int sw,                  
                    int* indexI, int nnzI, 
                    int* indexJ, int nnzJ, 
                    double* nnz_values, int nnz); 
                    
                    
void Stiffness_FE(  double E, double nu, 
                    int* e, int sie, int sje,
                    double* n, int sin, int sjn, 
                    int* conn, int siconn, int sjconn, 
                    double* N, int siN, int sjN, 
                    double* dNdxi, int sidNdxi, int sjdNdxi, 
                    double* dNdeta, int sidNdeta, int sjdNdeta, 
                    double* dNdzeta, int sidNdzeta, int sjdNdzeta, 
                    double* w, int sw,                  
                    int* indexI, int nnzI, 
                    int* indexJ, int nnzJ, 
                    double* nnz_values, int nnz); 


                                                       


                                                                 

#endif
