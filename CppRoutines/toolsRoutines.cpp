#include "toolsRoutines.h"

void EvaluateBspline3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v)
{
	int spanx, spany, spanz ;
	int nbf_xi   = sknotXi   -1 - deg_xi   ;
	int nbf_eta  = sknotEta  -1 - deg_eta  ;
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	double* Nx = array_1d(deg_xi+1);
	double* Ny = array_1d(deg_eta+1);
	double* Nz = array_1d(deg_zeta+1);
	
    double s ;

    
	for ( int p = 0; p < sx; p++ )
	{		
		spanx = findspan(nxi,   deg_xi   ,  x[p], knotXi);
		spany = findspan(neta,  deg_eta  ,  y[p], knotEta);
		spanz = findspan(nzeta, deg_zeta ,  z[p], knotZeta);
		
		basisfuns(spanx, x[p], deg_xi,   knotXi  , Nx);
		basisfuns(spany, y[p], deg_eta,  knotEta , Ny);
		basisfuns(spanz, z[p], deg_zeta, knotZeta, Nz);
		
		s=0;
		
		for ( int k =0; k  < deg_zeta + 1 ; k++ )
		{
		  for ( int j=0; j < deg_eta + 1; j++ ) 
			{
				for ( int i =0; i < deg_xi + 1; i++)
				{
					s = s + Nz[k]*Ny[j]*Nx[i]*cp[ spanx-deg_xi+i +(spany-deg_eta+j)*nbf_xi + (spanz-deg_zeta+k)*nbf_xi*nbf_eta ];
				}
			}
		}		
		v[p] = s; 
	}    
	delete[] Nx;
	delete[] Ny;
	delete[] Nz;    
}


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
					   double* dvdz, int s_dvdz )   
{
	int spanx, spany, spanz ;
	int nbf_xi   = sknotXi   -1 - deg_xi   ;
	int nbf_eta  = sknotEta  -1 - deg_eta  ;
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	
	double** Nx = array_2d(2,deg_xi+1); 
	double** Ny = array_2d(2,deg_eta+1); 
	double** Nz = array_2d(2,deg_zeta+1); 
	
	
	double s  ;  // Value of the leve-set volume 
	double sxx ;  // Gradient of the volume in x direction 
	double syy ;  // Gradient of the volume in y direction 
	double szz ;  // Gradient of the volume in z direction 
	
 
    
    double cpc; 
    
	for ( int p = 0; p < sx; p++ )
	{
				
		spanx = findspan(nxi,   deg_xi   ,  x[p], knotXi)   ;  
		spany = findspan(neta,  deg_eta  ,  y[p], knotEta)  ;   
		spanz = findspan(nzeta, deg_zeta ,  z[p], knotZeta) ;  
		
		dersbasisfuns(deg_xi,   knotXi,    x[p], spanx,1, Nx);
		dersbasisfuns(deg_eta,  knotEta ,  y[p], spany,1, Ny); 
		dersbasisfuns(deg_zeta, knotZeta , z[p], spanz,1, Nz); 
		
		s=0;
		sxx=0; 
		syy=0;
		szz=0; 
		
		for ( int k =0; k  < deg_zeta + 1 ; k++ )
		{
		  for ( int j=0; j < deg_eta + 1; j++ ) 
			{
				for ( int i =0; i < deg_xi + 1; i++)
				{
    				cpc = cp[ spanx-deg_xi+i +(spany-deg_eta+j)*nbf_xi + (spanz-deg_zeta+k)*nbf_xi*nbf_eta ];
					s   = s   + Nz[0][k]*Ny[0][j]*Nx[0][i]*cpc; 
					sxx = sxx + Nz[0][k]*Ny[0][j]*Nx[1][i]*cpc;
					syy = syy + Nz[0][k]*Ny[1][j]*Nx[0][i]*cpc;  
					szz = szz + Nz[1][k]*Ny[0][j]*Nx[0][i]*cpc;  
				}
			}
		}		
		v[p]    = s; 
		dvdx[p] = sxx ; 
		dvdy[p] = syy ; 
		dvdz[p] = szz ; 
	}    
	delete_array_2d(Nx); 
	delete_array_2d(Ny);   
	delete_array_2d(Nz);  
} 

void EvaluateBsplineStructured3D(double* x, int sx,
			  		   double* y, int sy,
			  		   double* z, int sz, 
                       double* knotXi,  int sknotXi, 
					   double* knotEta, int sknotEta ,
					   double* knotZeta, int sknotZeta, 
					   int deg_xi, int deg_eta,int deg_zeta,
					   double* cp, int s_cp,
					   double* v, int s_v)
{
	int* Spanx = new int [sx];
	int* Spany = new int [sy];
	int* Spanz = new int [sz];  
	
	int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	double** Nx = array_2d(sx, deg_xi   + 1);
	double** Ny = array_2d(sy, deg_eta  + 1);
	double** Nz = array_2d(sz, deg_zeta + 1);	
	
	double s  ;  // Value of the leve-set volume 
	
	// First evaluating univariate basis functions
	for ( int k=0; k < sx; k++)
	{
		Spanx[k] = findspan(nxi, deg_xi, x[k], knotXi);
		basisfuns(Spanx[k], x[k], deg_xi, knotXi, Nx[k]);
	}
	for ( int k=0; k < sy; k++)
	{
		Spany[k] = findspan(neta, deg_eta, y[k], knotEta);
		basisfuns(Spany[k], y[k], deg_eta, knotEta, Ny[k]);
	}
	for ( int k=0; k < sz; k++)
	{
		Spanz[k] = findspan(nzeta, deg_zeta, z[k], knotZeta);
		basisfuns(Spanz[k], z[k], deg_zeta, knotZeta, Nz[k]);
	}
	int m = 0 ; 
	for ( int r= 0; r < sz; r++)
	{
		for ( int q =0; q < sy; q++)
		{
			for ( int p =0; p < sx; p++)
			{
				s=0;
				for ( int k =0; k  < deg_zeta +1 ; k++ )
				{
		    		for ( int j=0; j < deg_eta + 1; j++ )
					{
						for ( int i =0; i < deg_xi + 1; i++)
						{
							s = s + Nz[r][k]*Ny[q][j]*Nx[p][i]*cp[ Spanx[p] - deg_xi + i  + (Spany[q] - deg_eta +j)*nbf_xi + (Spanz[r]-deg_zeta+k)*nbf_xi*nbf_eta    ];
						}
					}
				}
				v[m] = s ;
				m++;
			}
		}
	}
	delete_array_2d(Nx);
	delete_array_2d(Ny);
	delete_array_2d(Nz);
	delete[] Spanx;
	delete[] Spany;
	delete[] Spanz;
}		


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
					   double* dvdz, int s_dvdz ) 	
{
	int* Spanx = new int [sx];
	int* Spany = new int [sy];
	int* Spanz = new int [sz];  
	
	int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	
	double** Nx1d = array_2d(2,deg_xi+1); 
	double** Ny1d = array_2d(2,deg_eta+1); 
	double** Nz1d = array_2d(2,deg_zeta+1); 
	double** Nx   = array_2d(sx, deg_xi   + 1);
	double** Ny   = array_2d(sy, deg_eta  + 1);
	double** Nz   = array_2d(sz, deg_zeta + 1);
	double** dNxdx = array_2d(sx, deg_xi   + 1);
	double** dNydy = array_2d(sy, deg_eta  + 1);
	double** dNzdz = array_2d(sz, deg_zeta + 1);	
	
	
	double s  ;  // Value of the leve-set volume 
	double sxx ;  // Gradient of the volume in x direction 
	double syy ;  // Gradient of the volume in y direction 
	double szz ;  // Gradient of the volume in z direction  
 
    
    double cpc; 
 
	// First evaluating univariate basis functions
	for ( int k=0; k < sx; k++)
	{
		Spanx[k] = findspan(nxi, deg_xi, x[k], knotXi);
		dersbasisfuns(deg_xi,   knotXi,   x[k], Spanx[k],1, Nx1d);
		for ( int l =0 ; l < deg_xi+1; l++ )
		{
			Nx[k][l]    = Nx1d[0][l] ;
			dNxdx[k][l] = Nx1d[1][l] ;
		}	
	}
	for ( int k=0; k < sy; k++)
	{
		Spany[k] = findspan(neta, deg_eta, y[k], knotEta);
		dersbasisfuns(deg_eta,   knotEta,   y[k], Spany[k],1, Ny1d);
		for ( int l =0 ; l < deg_eta+1; l++ )
		{
			Ny[k][l]    = Ny1d[0][l] ;
			dNydy[k][l] = Ny1d[1][l] ;
		}	
		
	}
	for ( int k=0; k < sz; k++)
	{
		Spanz[k] = findspan(nzeta, deg_zeta, z[k], knotZeta);
		dersbasisfuns(deg_zeta,   knotZeta,   z[k], Spanz[k],1, Nz1d);
		for ( int l =0 ; l < deg_zeta+1; l++ )
		{                                    
			Nz[k][l]    = Nz1d[0][l] ;      
			dNzdz[k][l] = Nz1d[1][l] ;    
		}	
	} 
	int m = 0 ; 
	for ( int r= 0; r < sz; r++)
	{
		for ( int q =0; q < sy; q++)
		{
			for ( int p =0; p < sx; p++)
			{
				s=0;
				sxx=0;
				syy=0;
				szz=0; 
				for ( int k =0; k  < deg_zeta +1 ; k++ )
				{
		    		for ( int j=0; j < deg_eta + 1; j++ )
					{
						for ( int i =0; i < deg_xi + 1; i++)
						{
    						cpc = cp[ Spanx[p] - deg_xi + i  + (Spany[q] - deg_eta +j)*nbf_xi + (Spanz[r]-deg_zeta+k)*nbf_xi*nbf_eta];
							s   = s   + Nz[r][k]*Ny[q][j]*Nx[p][i]*cpc ; 
							sxx = sxx + Nz[r][k]*Ny[q][j]*dNxdx[p][i]*cpc ; 
							syy = syy + Nz[r][k]*dNydy[q][j]*Nx[p][i]*cpc ; 
							szz = szz + dNzdz[r][k]*Ny[q][j]*Nx[p][i]*cpc ;  
						}
					}
				}       
				v[m]    = s   ;
				dvdx[m] = sxx ;  
				dvdy[m] = syy ; 
				dvdz[m] = szz ; 
				m++ ;  
			}
		}
	}
	delete_array_2d(Nx1d);
	delete_array_2d(Ny1d);
	delete_array_2d(Nz1d);
	delete_array_2d(Nx);
	delete_array_2d(Ny);
	delete_array_2d(Nz);
	delete_array_2d(dNxdx);
	delete_array_2d(dNydy);
	delete_array_2d(dNzdz);	
	delete[] Spanx;
	delete[] Spany;
	delete[] Spanz;
}

 

void LocatePointsInTetraFE_Mesh1(double* x, int sx, 
                           double* y, int sy, 
                           double* z, int sz, 
                           int* e, int sie, int sje, 
                           double* n, int sin, int sjn, 
                           int* ie, int size_ie,
                           double eps )  
{
    #pragma omp parallel 
    {
        #pragma omp for 
            for (int ip=0; ip< sx; ip++)
            {
                ie[ip] = -1 ; 
                for (int i=0; i< sie; i++)
                {
                    double xn[4] = {n[0+e[0+i*sje]*sjn],n[0+e[1+i*sje]*sjn],n[0+e[2+i*sje]*sjn],n[0+e[3+i*sje]*sjn]} ;
                    double yn[4] = {n[1+e[0+i*sje]*sjn],n[1+e[1+i*sje]*sjn],n[1+e[2+i*sje]*sjn],n[1+e[3+i*sje]*sjn]} ; 
                    double zn[4] = {n[2+e[0+i*sje]*sjn],n[2+e[1+i*sje]*sjn],n[2+e[2+i*sje]*sjn],n[2+e[3+i*sje]*sjn]} ; 
                    int in = PointIsInTetrahedron1(xn, yn, zn, x[ip], y[ip], z[ip], eps );
                    if (in==1)
                    {
                        ie[ip] = i ;
                        break ; 
                    } 
                }
                
            }
    }
}

void LocatePointsInTetraFE_Mesh2(double* x, int sx, 
                           double* y, int sy, 
                           double* z, int sz, 
                           int* e, int sie, int sje, 
                           double* n, int sin, int sjn, 
                           int* ie, int size_ie ,
                           double eps)  
{
    #pragma omp parallel 
    {
        #pragma omp for 
            for (int ip=0; ip< sx; ip++)
            {
                ie[ip] = -1 ; 
                for (int i=0; i< sie; i++)
                {
                    double xn[4] = {n[0+e[0+i*sje]*sjn],n[0+e[1+i*sje]*sjn],n[0+e[2+i*sje]*sjn],n[0+e[3+i*sje]*sjn]} ;
                    double yn[4] = {n[1+e[0+i*sje]*sjn],n[1+e[1+i*sje]*sjn],n[1+e[2+i*sje]*sjn],n[1+e[3+i*sje]*sjn]} ; 
                    double zn[4] = {n[2+e[0+i*sje]*sjn],n[2+e[1+i*sje]*sjn],n[2+e[2+i*sje]*sjn],n[2+e[3+i*sje]*sjn]} ; 
                    int in = PointIsInTetrahedron2(xn, yn, zn, x[ip], y[ip], z[ip], eps );
                    if (in==1)
                    {
                        ie[ip] = i ;
                        break ; 
                    } 
                }
                
            }
    }
}

void LocatePointsInTetraFE_Mesh3(double* x, int sx, 
                                 double* y, int sy, 
                                 double* z, int sz, 
                                 int* face_indices, int s_faces_indices, 
                                 int* e, int sie, int sje, 
                                 double* n, int sin, int sjn, 
                                 int* connFaces, int si_connFaces, int sj_connFaces,
                                 int* ie, int size_ie,
                                 double eps   )
{
    #pragma omp parallel 
    {
        #pragma omp for 
        for (int ip=0; ip<sx; ip++)
        {
            ie[ip] = -1 ;
            int fe ; // element index of face 
            // Loop on the two face candidates 
            for (int i=0; i<2; i++)
            {
                // Getting the element supporting the face
                fe = connFaces[i+face_indices[ip]*2] ; 
                if (fe != -1 )
                {
                    double xn[4] = {n[0+e[0+fe*sje]*sjn],n[0+e[1+fe*sje]*sjn],n[0+e[2+fe*sje]*sjn],n[0+e[3+fe*sje]*sjn]} ;
                    double yn[4] = {n[1+e[0+fe*sje]*sjn],n[1+e[1+fe*sje]*sjn],n[1+e[2+fe*sje]*sjn],n[1+e[3+fe*sje]*sjn]} ; 
                    double zn[4] = {n[2+e[0+fe*sje]*sjn],n[2+e[1+fe*sje]*sjn],n[2+e[2+fe*sje]*sjn],n[2+e[3+fe*sje]*sjn]} ; 
                    int in = PointIsInTetrahedron2(xn, yn, zn, x[ip], y[ip], z[ip], eps );
                    if (in==1)
                    {
                        ie[ip] = fe ;
                    } 
                }
                if (ie[ip]!=-1)
                {
                    break ; 
                }      
            }
        }        
    }
}



//     int eps = 1.e-15 ; 
//     #pragma omp parallel 
//     {
//         #pragma omp for 
//             for (int ip=0; ip< sx; ip++)
//             {
//                 double A[16] ;  
//                 double u[4] = { x[ip] +eps , y[ip] +eps, z[ip] +eps, 1  }; 
//                 double l[4]  ; 
//                 double invA[16] ; 
//                 double detA ;
//                 
//                 A[12] = 1 ;
//                 A[13] = 1 ; 
//                 A[14] = 1 ; 
//                 A[15] = 1 ; 
//                 
//                 ie[ip] = -1 ; 
//                 for (int i=0; i< sie; i++)
//                 {
//                     // We look if the points is in the current element 
//                     
//                     A[0] = n[0+e[0+i*sje]*sjn] ; 
//                     A[4] = n[1+e[0+i*sje]*sjn] ; 
//                     A[8] = n[2+e[0+i*sje]*sjn] ; 
//                     
//                     A[1] = n[0+e[1+i*sje]*sjn] ; 
//                     A[5] = n[1+e[1+i*sje]*sjn] ; 
//                     A[9] = n[2+e[1+i*sje]*sjn] ;
//                     
//                     A[2]  = n[0+e[2+i*sje]*sjn] ; 
//                     A[6]  = n[1+e[2+i*sje]*sjn] ; 
//                     A[10] = n[2+e[2+i*sje]*sjn] ;
//                     
//                     A[3]  = n[0+e[3+i*sje]*sjn] ;
//                     A[7]  = n[1+e[3+i*sje]*sjn] ;
//                     A[11] = n[2+e[3+i*sje]*sjn] ; 
//                     
//                     /******************************************************************/  
//                     // Inverting the (4x4) matrix 
//                     // Solving the system Al=u (Convex hull proprety all l_i must be greater than 0)
//                     invA[0] = A[5]  * A[10] * A[15] - 
//                              A[5]  * A[11] * A[14] - 
//                              A[9]  * A[6]  * A[15] + 
//                              A[9]  * A[7]  * A[14] +
//                              A[13] * A[6]  * A[11] - 
//                              A[13] * A[7]  * A[10];
//         
//                     invA[4] = -A[4]  * A[10] * A[15] + 
//                               A[4]  * A[11] * A[14] + 
//                               A[8]  * A[6]  * A[15] - 
//                               A[8]  * A[7]  * A[14] - 
//                               A[12] * A[6]  * A[11] + 
//                               A[12] * A[7]  * A[10];
//         
//                     invA[8] = A[4]  * A[9] * A[15] - 
//                              A[4]  * A[11] * A[13] - 
//                              A[8]  * A[5] * A[15] + 
//                              A[8]  * A[7] * A[13] + 
//                              A[12] * A[5] * A[11] - 
//                              A[12] * A[7] * A[9];
//         
//                     invA[12] = -A[4]  * A[9] * A[14] + 
//                                A[4]  * A[10] * A[13] +
//                                A[8]  * A[5] * A[14] - 
//                                A[8]  * A[6] * A[13] - 
//                                A[12] * A[5] * A[10] + 
//                                A[12] * A[6] * A[9];
//         
//                     invA[1] = -A[1]  * A[10] * A[15] + 
//                               A[1]  * A[11] * A[14] + 
//                               A[9]  * A[2] * A[15] - 
//                               A[9]  * A[3] * A[14] - 
//                               A[13] * A[2] * A[11] + 
//                               A[13] * A[3] * A[10];
//         
//                     invA[5] = A[0]  * A[10] * A[15] - 
//                              A[0]  * A[11] * A[14] - 
//                              A[8]  * A[2] * A[15] + 
//                              A[8]  * A[3] * A[14] + 
//                              A[12] * A[2] * A[11] - 
//                              A[12] * A[3] * A[10];
//         
//                     invA[9] = -A[0]  * A[9] * A[15] + 
//                               A[0]  * A[11] * A[13] + 
//                               A[8]  * A[1] * A[15] - 
//                               A[8]  * A[3] * A[13] - 
//                               A[12] * A[1] * A[11] + 
//                               A[12] * A[3] * A[9];
//         
//                     invA[13] = A[0]  * A[9] * A[14] - 
//                               A[0]  * A[10] * A[13] - 
//                               A[8]  * A[1] * A[14] + 
//                               A[8]  * A[2] * A[13] + 
//                               A[12] * A[1] * A[10] - 
//                               A[12] * A[2] * A[9];
//         
//                     invA[2] = A[1]  * A[6] * A[15] - 
//                              A[1]  * A[7] * A[14] - 
//                              A[5]  * A[2] * A[15] + 
//                              A[5]  * A[3] * A[14] + 
//                              A[13] * A[2] * A[7] - 
//                              A[13] * A[3] * A[6];
//         
//                     invA[6] = -A[0]  * A[6] * A[15] + 
//                               A[0]  * A[7] * A[14] + 
//                               A[4]  * A[2] * A[15] - 
//                               A[4]  * A[3] * A[14] - 
//                               A[12] * A[2] * A[7] + 
//                               A[12] * A[3] * A[6];
//         
//                     invA[10] = A[0]  * A[5] * A[15] - 
//                               A[0]  * A[7] * A[13] - 
//                               A[4]  * A[1] * A[15] + 
//                               A[4]  * A[3] * A[13] + 
//                               A[12] * A[1] * A[7] - 
//                               A[12] * A[3] * A[5];
//         
//                     invA[14] = -A[0]  * A[5] * A[14] + 
//                                A[0]  * A[6] * A[13] + 
//                                A[4]  * A[1] * A[14] - 
//                                A[4]  * A[2] * A[13] - 
//                                A[12] * A[1] * A[6] + 
//                                A[12] * A[2] * A[5];
//         
//                     invA[3] = -A[1] * A[6] * A[11] + 
//                               A[1] * A[7] * A[10] + 
//                               A[5] * A[2] * A[11] - 
//                               A[5] * A[3] * A[10] - 
//                               A[9] * A[2] * A[7] + 
//                               A[9] * A[3] * A[6];
//         
//                     invA[7] = A[0] * A[6] * A[11] - 
//                              A[0] * A[7] * A[10] - 
//                              A[4] * A[2] * A[11] + 
//                              A[4] * A[3] * A[10] + 
//                              A[8] * A[2] * A[7] - 
//                              A[8] * A[3] * A[6];
//         
//                     invA[11] = -A[0] * A[5] * A[11] + 
//                                A[0] * A[7] * A[9] + 
//                                A[4] * A[1] * A[11] - 
//                                A[4] * A[3] * A[9] - 
//                                A[8] * A[1] * A[7] + 
//                                A[8] * A[3] * A[5];
//         
//                     invA[15] = A[0] * A[5] * A[10] - 
//                               A[0] * A[6] * A[9] - 
//                               A[4] * A[1] * A[10] + 
//                               A[4] * A[2] * A[9] + 
//                               A[8] * A[1] * A[6] - 
//                               A[8] * A[2] * A[5];
//         
//                     detA = A[0] * invA[0] + A[1] * invA[4] + A[2] * invA[8] + A[3] * invA[12];    
//                     
//                     l[0] = (invA[0]*u[0]  + invA[1]*u[1]  + invA[2]*u[2]  + invA[3]*u[3])/detA ; 
//                     l[1] = (invA[4]*u[0]  + invA[5]*u[1]  + invA[6]*u[2]  + invA[7]*u[3])/detA ; 
//                     l[2] = (invA[8]*u[0]  + invA[9]*u[1]  + invA[10]*u[2] + invA[11]*u[3])/detA ; 
//                     l[3] = (invA[12]*u[0] + invA[13]*u[1] + invA[14]*u[2] + invA[15]*u[3])/detA ;         
//                     /******************************************************************/
//         
//                     // Check if the point is in the element 
//                     if ( l[0]>0 && l[1]>0 && l[2]>0 && l[3]>0 )
//                     {
//                         // In this case : the point (x[ip],y[ip],z[ip]) belongs to the element 
//                         ie[ip] = i ; 
//                         break ;     
//                     }
//                 }
//             }        
//     }

 

                            