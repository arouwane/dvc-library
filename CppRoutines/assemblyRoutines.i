%module assemblyRoutines 
%{
#define SWIG_FILE_WITH_INIT
#include "assemblyRoutines.h"
%}

%include "numpy.i"
%include "std_vector.i"
 
//GetBsplineFunctionsMatrixStructured
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* valuesN, int s_valuesN) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int s_indexI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int s_indexJ) }



//GetBsplineFunctionsAndDerivativesMatrixStructured 
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* valuesN, int s_valuesN) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* valuesdNdx, int s_valuesdNdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* valuesdNdy, int s_valuesdNdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* valuesdNdz, int s_valuesdNdz) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int s_indexI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int s_indexJ) }



// FcmIntegrationTrilinearInterp 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }



// FcmTrilinearInterpStiffnessAndBfIntegralParallel
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* intBf, int sintBf) }



// FcmTrilinearInterpStiffnessAndBfIntegral
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* image, int si, int sj, int sk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* intBf, int sintBf) }

 
 
// Laplacian_Structured
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 

// HomogeneousStiffness
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 
 

// L2Projection 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi2,int sknotXi2) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta2,int sknotEta2) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta2,int sknotZeta2) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U2, int sU2) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 
 
 

// VoxelIntegrationThreshold
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }



// VoxelIntegrationMask
%apply (double* IN_ARRAY1, int DIM1) { (double* maskip, int smaskip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }

 


// DVC_LHS_thrsh
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ipIndices, int sipIndices) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 



// DVC_RHS_thrsh_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ipIndices, int sipIndices) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 



// DVC_RHS_ZN_thrsh_TrilinearInterp 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fmeane, int sfmeane) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fstde,  int sfstde) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ipIndices, int sipIndices) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) }



 
                                                                  
// DVC_LHS_Structured
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 


// DVC_LHS_Structured_Parallel
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 




 
 
// DVC_RHS_Structured_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 




// DVC_RHS_Structured_TrilinearInterp_Parallel
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 
 
 
 
// DVC_RHS_ZN_Structured_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fmeane, int sfmeane) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fstde,  int sfstde) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dyne,   int sdyne) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* elementRes, int selementRes) } 







// GLR_Structured_TrilinearInterp
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) } 



// GLR_Structured_TrilinearInterp_Parallel
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) } 



// GLR_ZN_Structured_TrilinearInterp
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* fmeane, int sfmeane) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* fstde,  int sfstde) } 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) } 


// DVC_RHS_Structured_TrilinearInterp
//%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
//%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
//%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
//%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
//%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// DVC_RHS_Structured_CBspline2
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// GLR_Structured_CBspline2
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) } 

 
//  DVC_RHS_Structured_CBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* NOELEM, int siNOELEM, int sjNOELEM ) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// GLR_Structured_CBspline3 
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi, int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta, int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta, int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) {(double* xg, int sxg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* yg, int syg)}
%apply (double* IN_ARRAY1, int DIM1) {(double* zg, int szg)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nxi, int siNxi, int sjNxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Neta, int siNeta, int sjNeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* Nzeta, int siNzeta, int sjNzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }    



// DVC_LHS_FE
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 


// DVC_RHS_FE_TrilinearInterp 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// GLR_FE_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }   



// Gophi_FEMesh_TrilinearInterp
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage, int sknotXiImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage, int sknotEtaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage, int sknotZetaImage) }
%apply (double* IN_ARRAY1, int DIM1) { (double* xi, int sxi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* eta, int seta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zeta, int szeta) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ie, int size_ie ) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }     



// DVC_RHS_FE_CBspline3 
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// GLR_FE_CBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }   



// Gophi_FEMesh_CBspline3
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) { (double* g, int sgi, int sgj, int sgk) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* xc, int sxc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* yc, int syc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zc, int szc) }
%apply (double* IN_ARRAY1, int DIM1) { (double* xi, int sxi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* eta, int seta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zeta, int szeta) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ie, int size_ie ) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }     

 

// DVC_RHS_FE_L2ProjLumped
%apply (double* IN_ARRAY1, int DIM1) { (double* lsc_g, int s_lsc_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage_g, int sknotXiImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage_g, int sknotEtaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage_g, int sknotZetaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdx, int sdfipdx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdy, int sdfipdy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* dfipdz, int sdfipdz) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* rhs, int ndof) } 


// GLR_FE_L2ProjLumped
%apply (double* IN_ARRAY1, int DIM1) { (double* lsc_g, int s_lsc_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage_g, int sknotXiImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage_g, int sknotEtaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage_g, int sknotZetaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* fip, int sfip) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) } 
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }   



// Gophi_FEMesh_L2ProjLumped
%apply (double* IN_ARRAY1, int DIM1) { (double* lsc_g, int s_lsc_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXiImage_g, int sknotXiImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEtaImage_g, int sknotEtaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZetaImage_g, int sknotZetaImage_g) }
%apply (double* IN_ARRAY1, int DIM1) { (double* xi, int sxi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* eta, int seta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* zeta, int szeta) }
%apply (int* IN_ARRAY1, int DIM1) { (int* ie, int size_ie ) } 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY1, int DIM1) { (double* U, int sU) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* glr, int sglr) }  

 
 
// Laplacian_FE 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 


// Stiffness_FE 
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* conn, int siconn, int sjconn) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* N, int siN, int sjN) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdxi, int sidNdxi, int sjdNdxi) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdeta, int sidNdeta, int sjdNdeta) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* dNdzeta, int sidNdzeta, int sjdNdzeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* w, int sw) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexI, int nnzI) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* indexJ, int nnzJ) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* nnz_values, int nnz) } 
 
 


 
 


 
 
// Return value std::vector<double>
%template(VectDouble) std::vector<double>;
// Return value std::vector<int>
%template(VectInt) std::vector<int>; 

%init %{
    import_array();
%}

%include "assemblyRoutines.h" 




