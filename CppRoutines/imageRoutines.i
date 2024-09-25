%module imageRoutines
%{
#define SWIG_FILE_WITH_INIT
#include "imageRoutines.h"
%}

%include "numpy.i"


//ComputeL2ProjLumpedCoefficients		       
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* lsc, int s_lsc) }

//ComputeL2ProjLumpedCoefficientsSumFact		       
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* lsc, int s_lsc) }
 

//EvaluateTrilinearInterpolation
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }


//EvaluateTrilinearInterpolationAndGradient
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }


//EvaluateTrilinearInterpolationAndGradientStructured 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }



//EvaluateCardBspline2
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }



//EvaluateCardBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }



//EvaluateCardBsplineAndGradient2 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }
 

            					               
//EvaluateCardBsplineAndGradient3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }            					                 
 
            					               



//EvaluateCardBsplineAndGradient2Structured
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }



//EvaluateCardBsplineAndGradient3Structured
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc,int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc,int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc,int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }



//GetMeanImageAndStdOnMesh_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) }  
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }        



// GetMeanImageAndStdOnMesh_CBspline2
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) } 
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)} 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) }                                         
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }   


// GetMeanImageAndStdOnMesh_CBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) } 
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)} 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) }                                         
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }  

                                              

                                                                   
 
//GetMeanImageAndStdOnMesh_thrsh_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }  



// GetC8MeshFromVoxelsTrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* elementMask, int se) }  


 
 

// GetMeanImageAndStdOnFE_Mesh_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }   



// GetMeanImageAndStdOnFE_Mesh_CBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* f, int sif, int sjf, int skf) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }  

 
// GetMeanImageAndStdOnFE_Mesh_L2ProjLumped
%apply (double* IN_ARRAY1, int DIM1) { (double* cp, int scp) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fmean, int sfmean) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fstd, int sfstd) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dyne, int sdyne) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }   

  
                                             



 
%init %{
    import_array();
%}

%include "imageRoutines.h"

 
