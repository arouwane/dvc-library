%module assembly_Bspline_3d
%{
#define SWIG_FILE_WITH_INIT
#include "assembly_Bspline_3d.h"
%}

%include "numpy.i"
 
 
// ComputeIntegratedDiffTensors
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dxi_dxi, int s_dxi_dxi) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dxi_xi,  int s_dxi_xi) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* xi_dxi,  int s_xi_dxi) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* xi_xi,   int s_xi_xi) }    
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* deta_deta, int s_deta_deta) }   
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* deta_eta, int s_deta_eta) }       
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* eta_deta, int s_eta_deta) }        
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* eta_eta, int s_eta_eta) }        
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dzeta_dzeta, int s_dzeta_dzeta) }   
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dzeta_zeta, int s_dzeta_zeta) }        
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* zeta_dzeta, int s_zeta_dzeta) }        
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* zeta_zeta, int s_zeta_zeta) }         
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* eta_eta_dxi_dxi, int s_eta_eta_dxi_dxi) } 
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* eta_deta_dxi_xi, int s_eta_deta_dxi_xi) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* deta_eta_xi_dxi, int s_deta_eta_xi_dxi) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* deta_deta_xi_xi, int s_deta_deta_xi_xi) }
 



 %init %{
    import_array();
%}
 
%include "assembly_Bspline_3d.h" 

