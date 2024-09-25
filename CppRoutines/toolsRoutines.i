%module toolsRoutines
%{
#define SWIG_FILE_WITH_INIT
#include "toolsRoutines.h"
%}

%include "numpy.i"



//EvaluateBspline3D
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* cp, int s_cp) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }



//EvaluateBsplineAndGradient3D
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* cp, int s_cp) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }



//EvaluateBsplineStructured3D 
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* cp, int s_cp) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }



//EvaluateBsplineAndGradientStructured3D
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotXi,int sknotXi) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotEta,int sknotEta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* knotZeta,int sknotZeta) }
%apply (double* IN_ARRAY1, int DIM1) { (double* cp, int s_cp) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* v, int s_v) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdx, int s_dvdx) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdy, int s_dvdy) }
%apply (double* ARGOUT_ARRAY1, int DIM1) { (double* dvdz, int s_dvdz) }


   
    
//LocatePointsInTetraFE_Mesh1
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* ie, int size_ie ) }     
 

//LocatePointsInTetraFE_Mesh2
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* ie, int size_ie ) }    


//LocatePointsInTetraFE_Mesh3
%apply (double* IN_ARRAY1, int DIM1) { (double* x, int sx) }
%apply (double* IN_ARRAY1, int DIM1) { (double* y, int sy) }
%apply (double* IN_ARRAY1, int DIM1) { (double* z, int sz) }
%apply (int* IN_ARRAY1, int DIM1) { (int* face_indices, int s_faces_indices) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* e, int sie, int sje) }
%apply (double* IN_ARRAY2, int DIM1, int DIM2) { (double* n, int sin, int sjn) }
%apply (int* IN_ARRAY2, int DIM1, int DIM2) { (int* connFaces, int si_connFaces, int sj_connFaces) }
%apply (int* ARGOUT_ARRAY1, int DIM1) { (int* ie, int size_ie ) }   

 
                                                                                                         

 
%init %{
    import_array();
%}

%include "toolsRoutines.h"

 
