#include "tools.h"

 
void EvaluateTrilinearInterpolationAndGradientStructured(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v, 
            					    double* dvdx, 
            					    double* dvdy,
            					    double* dvdz) 
{

	int deg_xi   = 1; 
	int deg_eta  = 1;
	int deg_zeta = 1; 
	
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
    int I,J,K; 


 
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
        					I = Spany[q]-deg_eta+j ;
        					J = sj -1 - (Spanz[r]-deg_zeta+k) ; 
        					K = Spanx[p]-deg_xi+i ;
        					cpc = image[K+J*sk+I*sk*sj];    
							s   = s   + Nz[r][k]*Ny[q][j]*Nx[p][i]*cpc    ;   
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

void EvaluateTrilinearInterpolationStructured(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
            					    double* v )
{

	int deg_xi   = 1; 
	int deg_eta  = 1;
	int deg_zeta = 1; 
	
	int* Spanx = new int [sx];
	int* Spany = new int [sy];
	int* Spanz = new int [sz];  
	
	int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
 
	double** Nx   = array_2d(sx, deg_xi   + 1);
	double** Ny   = array_2d(sy, deg_eta  + 1);
	double** Nz   = array_2d(sz, deg_zeta + 1);
 
	double s  ;  // Value of the leve-set volume 
 
 
    double cpc ;   
    int I,J,K  ; 
 
 
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
        					I = Spany[q]-deg_eta+j ;
        					J = sj -1 - (Spanz[r]-deg_zeta+k) ; 
        					K = Spanx[p]-deg_xi+i ;
        					cpc = image[K+J*sk+I*sk*sj];    
							s   = s   + Nz[r][k]*Ny[q][j]*Nx[p][i]*cpc ;   
						} 
					}
				}       
				v[m]    = s   ;
				m++ ;  
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



double EvaluateTrilinearInterpolationOnOnePoint(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double x, double y, double z )   
{

	double xc = x ; 
	double yc = y ; 
	double zc = z ; 
	
	if  ( x > knotXi[sknotXi-1] ) 
  	{
      	xc = knotXi[sknotXi-1] ;
  	}
  	else if ( x < knotXi[0] ) 
  	{
  		xc = knotXi[0]; 
  	}
  	if  ( y > knotEta[sknotEta-1] ) 
  	{
      	yc = knotEta[sknotEta-1] ;
  	}
  	else if ( y < knotEta[0] ) 
  	{
  		yc = knotEta[0]; 
  	}
  	if  ( z > knotZeta[sknotZeta-1] ) 
  	{
      	zc = knotZeta[sknotZeta-1] ;
  	}
  	else if ( z < knotZeta[0] ) 
  	{
  		zc = knotZeta[0]; 
  	}
  	
	
	
	int spanx, spany, spanz ;
	int deg_xi   = 1; 
	int deg_eta  = 1;
	int deg_zeta = 1; 
	
	int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	double* Nx = new double[deg_xi+1]   ; 
	double* Ny = new double[deg_eta+1]  ; 
	double* Nz = new double[deg_zeta+1] ; 
	
    double s ;
    
    int I,J,K; 
	
	spanx = findspan(nxi,   deg_xi   ,  xc, knotXi);
	spany = findspan(neta,  deg_eta  ,  yc, knotEta);
	spanz = findspan(nzeta, deg_zeta ,  zc, knotZeta);
		
	basisfuns(spanx, xc, deg_xi,   knotXi  , Nx);
	basisfuns(spany, yc, deg_eta,  knotEta , Ny);
	basisfuns(spanz, zc, deg_zeta, knotZeta, Nz);
		
	s=0;
		
	for ( int k =0; k  < deg_zeta + 1 ; k++ )
	{
    	J = sj -1 - (spanz-deg_zeta+k) ; 
    	for ( int j=0; j < deg_eta + 1; j++ ) 
		{
    		I = spany-deg_eta+j ;
    		for ( int i =0; i < deg_xi + 1; i++)
			{
        		K = spanx-deg_xi+i ;
				s = s + Nz[k]*Ny[j]*Nx[i]*image[K+J*sk+I*sk*sj];    
			}
		}
	}		
	delete[] Nx;
	delete[] Ny;
	delete[] Nz;
	return s; 
}


double Evaluate3dBsplineOnOnePoint(double x, double y, double z,
                                   double* knotXi,  int sknotXi, 
					               double* knotEta, int sknotEta ,
					               double* knotZeta, int sknotZeta, 
					               int deg_xi, int deg_eta,int deg_zeta,
					               double* cp, int s_cp) 
{
	double xc = x ; 
	double yc = y ; 
	double zc = z ; 
	
  	if  ( x >= knotXi[sknotXi-1] ) 
  	{
      	xc = knotXi[sknotXi-1] ;
  	}
  	else if ( x <= knotXi[0] ) 
  	{
  		xc = knotXi[0]; 
  	}
  	if  ( y >= knotEta[sknotEta-1] ) 
  	{
      	yc = knotEta[sknotEta-1] ;
  	}
  	else if ( y <= knotEta[0] ) 
  	{
  		yc = knotEta[0]; 
  	}
  	if  ( z >= knotZeta[sknotZeta-1] ) 
  	{
      	zc = knotZeta[sknotZeta-1] ;
  	}
  	else if ( z <= knotZeta[0] ) 
  	{
  		zc = knotZeta[0]; 
  	}
  	
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
	
	
	spanx = findspan(nxi,   deg_xi   ,  xc, knotXi);
	spany = findspan(neta,  deg_eta  ,  yc, knotEta);
	spanz = findspan(nzeta, deg_zeta ,  zc, knotZeta);
		
	basisfuns(spanx, xc, deg_xi,   knotXi  , Nx);
	basisfuns(spany, yc, deg_eta,  knotEta , Ny);
	basisfuns(spanz, zc, deg_zeta, knotZeta, Nz);
		
    double s=0;
		
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
	return s; 
}

  					                  


int CheckCutCuboidTrilinearInterpolation(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
        					        double thrsh, 
        					        Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z)					   
{
    int deg_xi   = 1; 
	int deg_eta  = 1;
	int deg_zeta = 1; 
	int* Spanx = new int [x.size()];
	int* Spany = new int [y.size()];
	int* Spanz = new int [z.size()];  
	
	int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;

	int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
	
	double** Nx = array_2d(x.size(), deg_xi   + 1);
	double** Ny = array_2d(y.size(), deg_eta  + 1);
	double** Nz = array_2d(z.size(), deg_zeta + 1);	
	
	double s  ;  // Value of the gray-level field 
	int I,J,K; 
	int countZero=0;
	
	// First evaluating univariate basis functions
	for ( int k=0; k < x.size(); k++)
	{
		Spanx[k] = findspan(nxi, deg_xi, x[k], knotXi);
		basisfuns(Spanx[k], x[k], deg_xi, knotXi, Nx[k]);
	}
	for ( int k=0; k < y.size(); k++)
	{
		Spany[k] = findspan(neta, deg_eta, y[k], knotEta);
		basisfuns(Spany[k], y[k], deg_eta, knotEta, Ny[k]);
	}
	for ( int k=0; k < z.size(); k++)
	{
		Spanz[k] = findspan(nzeta, deg_zeta, z[k], knotZeta);
		basisfuns(Spanz[k], z[k], deg_zeta, knotZeta, Nz[k]);
	}  
	// Face 1  x=xmin (y-z plane) 
	int kx = 0; 
	for (int kz=0; kz<z.size(); kz++)
	{
        for (int ky=0; ky<y.size(); ky++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				}
    			}
        	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}
        	   
        }
    }
	// Face 2  x=xmax
	kx = x.size()-1; 
	for (int kz=0; kz<z.size(); kz++)
	{
        for (int ky=0; ky<y.size(); ky++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				}
    			}
        	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}
        	   
        }
    }
	// Face 3  y=ymin 
	int ky = 0; 
	for (int kz=0; kz<z.size(); kz++)
	{
        for (int kx=0; kx<x.size(); kx++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				}
    			}
        	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}
        	   
        }
    }
	// Face 4  y=ymax 
	ky = y.size()-1; 
	for (int kz=0; kz<z.size(); kz++)
	{
        for (int kx=0; kx<x.size(); kx++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				}
    			}
        	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}
        	   
        }
    }
	// Face 5  z=zmin  
	int kz = 0; 
	for (int ky=0; ky<y.size(); ky++)
	{
        for (int kx=0; kx<x.size(); kx++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{ 
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				}
    			}  
         	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}
        }
    }
	// Face 6  z=zmax 
	kz = z.size()-1 ;   
	for (int ky=0; ky<y.size(); ky++)
	{
        for (int kx=0; kx<x.size(); kx++)
        {
            s = 0; 
            // We determine if the point is in the domain 
            for ( int k =0; k  < deg_zeta + 1 ; k++ )
        	{
            	J = sj -1 -(Spanz[kz]-deg_zeta+k) ; 
            	for ( int j=0; j < deg_eta + 1; j++ ) 
    			{
        			I = Spany[ky]-deg_eta+j ;
    				for ( int i =0; i < deg_xi + 1; i++)
    				{ 
        				K = Spanx[kx]-deg_xi+i ;	
    					s = s + Nz[kz][k]*Ny[ky][j]*Nx[kx][i]*image[K+J*sk+I*sk*sj]; 
    				} 
    			}    
         	}	
        	if ( s-thrsh < 0 )
        	{
            	countZero++; 
        	}      	   
        }
    }
 return countZero;    
}



void interpolateLinearly(double xmin, double xmax, double lmin, double lmax, double &a, double &b)
{
  a = (lmax-lmin)/(xmax-xmin) ;
  b = (lmin*xmax - lmax*xmin)/(xmax-xmin) ;
}	
	

 
					               
					               
double CardinalBspline2(double x)
{
    double shift = 1.5;
	double x1 = x + shift ;
	double x2 = x1*x1 ;

	if  (  (0 < x1)   &&  ( x1 < 1) )
	{
    	return 0.5*x2;
	}
    else if  ( (1 <= x1)  &&  (x1 < 2)) 
    {
		return 0.5+(x1-1)-((x1-1)*(x1-1));
	}
    else if  ((2 <= x1) && (x1 < 3))
    {
        return 0.5-(x1-2)+(0.5*(x1-2)*(x1-2));
	}
    else 
    {
        return 0. ; 
	}
}
double CardinalBsplineDer2(double x)
{
    double shift = 1.5;
	double x1 = x + shift ;
    if  ( (0 < x1)  &&  (x1 < 1)) 
    {
    	return x1;
    }
    else if  ((1 <= x1) &&  (x1 < 2)) 
    {
    	return 1-2*(x1-1);
    }
    else if  ((2 <= x1) && (x1 < 3)) 
    {
        return  -1 + x1-2;
    }
    else
    {
        return 0. ;
    }
}
double CardinalBspline3(double x)
{
    double shift = 2;
	double x1 = x + shift ;
	double x2 = x1*x1 ;
	double x3 = x1*x1*x1 ;
  	if ((0 < x1) && (x1 < 1)) {
  		return x3/6. ;
  	}
  	else if ((1 <= x1)  &&  (x1 < 2)) {
  		return -0.5*x3+2*x2-2*x1+2./3;
  	}
  	else if ((2 <= x1)  &&  (x1 < 3)) {
  		return 0.5*x3-4*x2+10*x1-22./3;
  	}
  	else if ((3 <= x1)  &&  (x1 < 4)) {
  		return -x3/6.+2*x2-8*x1+32./3 ;
  	}
  	else {
  		return 0. ;
  	}	  	
}
double CardinalBsplineDer3(double x)
{
    double shift = 2;
	double x1 = x + shift ;
	double x2 = x1*x1 ;
	if ((0 < x1) &&  (x1 < 1)) 
	{
		return 0.5*x2 ;
	}
	else if ((1 <= x1) &&  (x1 < 2)) 
	{
		return -3./2*x2+4*x1-2;
	}
	else if ((2 <= x1) &&  (x1 < 3)) 
	{
		return 3./2*x2-8*x1+10;
	}
	else if ((3 <= x1) &&  (x1 < 4)) 
	{
		return -0.5*x2+4*x1-8 ;
	}
	else 
	{
		return 0. ;
	}
}

Eigen::MatrixXi getNoelem(int nxi, int neta, int nzeta, int p, int q, int r)
{
    int nbf_elem    =  (p+1)*(q+1)*(r+1) ; 
    int n_elems     =  nxi*neta*nzeta ; 
    Eigen::MatrixXi noelem(n_elems,nbf_elem); 
    int ie=0; 
    int t; 
    for (int k=0; k<nzeta; k++ )
    {
        for (int j=0; j<neta; j++)
        {
            for (int i=0; i<nxi; i++)
            {
                t=0; 
                for (int kk=0; kk<r+1; kk++)
                {
                    for (int jj=0; jj<q+1; jj++)
                    {
                        for (int ii=0; ii<p+1; ii++)
                        {
                            // xi + nxi*eta + nxi*neta*zeta 
                            noelem(ie,t) = i+ii+(j+jj)*(nxi+p)+(k+kk)*(nxi+p)*(neta+q) ; 
                            t++; 
                        }        
                    }
                }
                ie++; 
            }
        }
    }
    return noelem;  
}



Eigen::MatrixXd hookeVolume(double E, double v)
{
    double f = E/((1+v)*(1-2*v)) ; 
    Eigen::MatrixXd hooke(6,6); 
    hooke.fill(0) ; 
    hooke(0,0) = (1-v)*f ; 
    hooke(0,1) = v*f ; 
    hooke(0,2) = v*f ;
    hooke(1,0) = v*f ; 
    hooke(1,1) = (1-v)*f ;
    hooke(1,2) = v*f  ;
    hooke(2,0) = v*f  ;
    hooke(2,1) = v*f  ;
    hooke(2,2) = (1-v)*f ;
    hooke(3,3) = 0.5*(1-2*v)*f ;
    hooke(4,4) = 0.5*(1-2*v)*f ;
    hooke(5,5) = 0.5*(1-2*v)*f ;
    return hooke; 
}


// Get the gray-level of a physical point (x,y,z)
double GetGrayLevelOfPoint3D(double* image, int si, int sj, int sk,
                           double x, double y, double z, 
                           double xmin, double ymin, double zmin, double zmax,
                           double dx, double dy, double dz )
{
        //double u = y ; 
        //double v = zmin + zmax - z ; 
        //double w = x ; 
        int i = int(std::floor((y-ymin)/dy)) ;
        int j = int(std::floor((zmax-z)/dz)) ;
        int k = int(std::floor((x-xmin)/dx)) ;   
        return image[k+j*sk+i*sk*sj] ;  
}


// B-splines routines

// Gets the knot span on a uniform knot vector
// Each knot span is of length l
// Works only for C^{-1} knot vectors at the boundaries
int findspanUniformKnotVector(double* U, int sU, int deg, double l, double u)
{
		if (u==U[sU-deg-1]) return(sU-deg-2);
		return    std::floor( (u-U[0])/l )+  deg  ;
}

// Find the knot span of the parametric point u.
//
// INPUT:
//
//   n - number of control points - 1
//   p - spline degree
//   u - parametric point
//   U - knot sequence
//
// RETURN:
//
//   s - knot span
// Source: Alex Wiltschko - Github
// Algorithm A2.1 from 'The NURBS BOOK' pg68.
int findspan(int n, int p, double u, double *U)
{
  int low, high, mid;

  // special case
  if (u == U[n+1]) return(n);

  // do binary search
  low = p;
  high = n + 1;
  mid = (low + high) / 2;
  while (u < U[mid] || u >= U[mid+1])
  {
    if (u < U[mid])
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2;
  }

  return(mid);
}

// Returns the evaluation of basis function i at point u
// NURBS Book
double oneBasisFun(int i, double u, int p, double*U, int lenU)
{
	int m  = lenU -1;
	double Nip;
	double saved;
	double temp, Uleft, Uright;
	// Verifying the local support proprety of the Bsplines
	if ( u < U[i] || u > U[i+p+1] ){
		Nip = 0;
		return Nip;
	}
	// Special case evaluation on the first and last knots
	if  ( ( i==0 && u==U[0]) || (i==m-p-1 && u==U[m]) ) {
		Nip = 1;
		return Nip;
	}
	double* N = array_1d(p+1);
    for ( int j =0; j<=p ; j++ )
	{
       	if (   ( u >= U[i+j] ) && ( u < U[i+j+1] ) )
		{
			N[j] = 1;
		}
		else
		{
			N[j] = 0 ;
		}
	}
	//  Initialize zeroth-degree functions
	for ( int k = 1; k <=p; k++ )
	{
		if ( N[0] == 0 )
		{
			saved = 0 ;
		}
		else
		{
			saved=((u-U[i])*N[0])/(U[i+k]-U[i]);
		}
		for ( int j=0; j <= p-k; j++)
		{
			Uleft=U[i+j+1];
            Uright=U[i+j+k+1];
            if  ( N[j+1]==0. )
			{
				N[j]=saved ;
                saved=0.;
			}
			else
			{
			  	temp=N[j+1]/(Uright-Uleft) ;
                N[j]=saved+(Uright-u)*temp ;
                saved=(u-Uleft)*temp ;
			}

		}
	}

	Nip = N[0];
	delete[] N;
    return Nip;

}



// Basis Function.
//
// INPUT:
//
//   i - knot span  ( from FindSpan() )
//   u - parametric point
//   p - spline degree
//   U - knot sequence
//
// OUTPUT:
//
//   N - Basis functions vector[p+1]
//  Source: Alex Wiltschko - Github
// Algorithm A2.2 from 'The NURBS BOOK' pg70.

void basisfuns(int i, double u, int p, double *U, double *N)
{
  int j,r;
  double saved, temp;

  // work space
  double *left  = array_1d(p+1);
  double *right = array_1d(p+1);

  N[0] = 1.0;
  for (j = 1; j <= p; j++)
  {
    left[j]  = u - U[i+1-j];
    right[j] = U[i+j] - u;
    saved = 0.0;

    for (r = 0; r < j; r++)
    {
      temp = N[r] / (right[r+1] + left[j-r]);
      N[r] = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    }

    N[j] = saved;
  }

  delete[] left;
  delete[] right;
}


// Compute Non-zero basis functions and their derivatives.
//
// INPUT:
//
//   d  - spline degree         integer
//   k  - knot sequence         double  vector(nk)
//   u  - parametric point      double
//   s  - knot span             integer
//   n  - number of derivatives integer
//
// OUTPUT:
//
//   dN -  Basis functions      double  matrix(n+1,d+1)
//         and derivatives upto the nth derivative (n < d)
//
// Algorithm A2.3 from 'The NURBS BOOK' pg72.

 

void dersbasisfuns(int d, double *k, double u, int s,int n, double **ders)
{
  int i,j,r,s1,s2,rk,pk,j1,j2;
  double temp, saved, der;
  double **ndu, **a, *left, *right;
 
  ndu = array_2d(d+1, d+1);
  a = array_2d(2,d+1);
  left  = array_1d(d+1);
  right = array_1d(d+1);


    ndu[0][0] = 1.0;
    
    for (j=1; j<=d; j++) {
        left[j] = u - k[s+1-j];
        right[j] = k[s+j]-u;
        saved = 0.0;
        for (r=0; r<j; r++) {
            ndu[j][r] = right[r+1] + left[j-r];
            temp = ndu[r][j-1]/ndu[j][r];
    
            ndu[r][j] = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        ndu[j][j] = saved;
    }
    
    for (j=0; j<=d; j++)
        ders[0][j] = ndu[j][d];
    
    for (r=0; r<=d; r++) {
        s1 = 0; s2 = 1;
        a[0][0] = 1.0;
    
        for (i = 1; i <= n; i++) {
            der = 0.0;
            rk = r-i;  pk = d-i;
    
            if (r >= i) {
                a[s2][0] = a[s1][0] / ndu[pk+1][rk];
                der = a[s2][0] * ndu[rk][pk];
            }  
            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r-1 <= pk)
                j2 = i-1;
            else
                j2 = d-r;
    
            for (j = j1; j <= j2; j++) {
                a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
                der += a[s2][j] * ndu[rk+j][pk];
            }
            if (r <= pk) {
                a[s2][i] = -a[s1][i-1] / ndu[pk+1][r];
                der += a[s2][i] * ndu[r][pk];
            }
            ders[i][r] = der;
            j = s1; s1 = s2; s2 = j;
        }        
    }
      
    r = d;
    for (i = 1; i <= n; i++) {
        for (j = 0; j <= d; j++)
            ders[i][j] *= r;
        r *= d-i;
    }

    delete_array_2d(ndu);
    delete_array_2d(a);
    delete[] left;
    delete[] right;

}


void bezierExtraction(double *U, int p, int m, int n_elems, double*** C)
{

//    print_3d_array(operators,n_elems,p+1,p+1);
//    identity(operators[0],p+1);
//    print_3d_array(operators,n_elems,p+1,p+1);
    int a    = p+1;
    int b    = a+1;
    int nb   = 0;
    int mult;
    int i,j,k,r,save,s,l;
    double numer;
    double alpha;
    double * alphas = array_1d(p);
    identity(C[nb],p+1);
    while (b < m)  {
        if (nb == n_elems -1)
            break;
        identity(C[nb+1],p+1); // Initialize the next extraction operator
        i = b;
        // Count multiplicity of the knot at location b.
        while (b < m && U[b] == U[b-1]) {
            b = b+1; }
        mult = b-i+1;
        if (mult < p ) {
            //Use (10) to compute the alphas.
            numer = U[b-1]-U[a-1];
            zero_1d_array(alphas,p);
            for (j=p; j> mult; j--){
                alphas[j-mult-1] = numer / (U[a+j-1]-U[a-1]);
            }
            r = p-mult;
            //Update the matrix coefficients for r new knots
            for (j=1;j<r+1;j++){
                save = r-j+1;
                s = mult + j;
                for ( k=p+1; k>s; k--){
                     alpha = alphas[k-s-1];
//                     The following line corresponds to (9).
                     for (l=0; l<p+1;l++){
                           C[nb][l][k-1] = alpha*C[nb][l][k-1] + (1.0-alpha)*C[nb][l][k-2];
                     }
                     }
                if (b<m) {
                    //  Update overlapping coefficients of the next operator.
                    for (l=0; l<=j; l++){
                        C[nb+1][l+save-1][save-1] = C[nb][l-j+p][p];
                    }
                }
             }
             nb = nb +1;
             if (b<m){
                // Update indices for the next operator.
                a = b;
                b = b +1;
             }

            }
        }
    delete[] alphas;
}



// Gauss Integration routine
// Taken from numerical recipes in C  book
void gauleg(double x1, double x2, double x[], double w[], int n)
{
int m,j,i;
double z1,z,xm,xl,pp,p3,p2,p1;
m=(n+1)/2;
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
double pi  = std::acos(-1);
	for (i=1;i<=m;i++)
	{
		z=std::cos(pi*(i-0.25)/(n+0.5));
		do{
			 p1=1.0;
			 p2=0.0;
			 for (j=1;j<=n;j++)
			  {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			  }
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (std::abs(z-z1) > 1.e-15);
		x[i-1]=xm-xl*z;
		x[n+1-i-1]=xm+xl*z;
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i-1]=w[i-1];
	}
}
// Gauss integration on the reference triangle in [0,1]
int GaussTriangle(std::vector<double> &xg, std::vector<double> &yg, std::vector<double> &wg, int n)
{
  int npoints;
  if (n==1)
  {
     // Deg = 1 : 1 point
     double x[] = { 1./3 }  ;
     double y[] = { 1./3 }  ;
     double w[] = { 1./2 }  ;
     npoints = 1;
     xg.insert(xg.end(), x, x + npoints );
     yg.insert(yg.end(), y, y + npoints );
     wg.insert(wg.end(), w, w + npoints );
     return npoints;
  }
  else if (n==2)
  {
    // Deg = 2 : 3 points
    double x[] = { 1./6, 4./6, 1./6 } ;
    double y[] = { 1./6, 1./6, 4./6 } ;
    double w[] = { 1./6, 1./6, 1./6 } ;
    npoints = 3;
    xg.insert(xg.end(), x, x + npoints );
    yg.insert(yg.end(), y, y + npoints );
    wg.insert(wg.end(), w, w + npoints );
    return npoints;
  }
  else if (n==3)
  {
    // Deg = 3 : 4 points
    double x[] = { 1./3, 1./5, 3./5, 1./5 } ;
    double y[] = { 1./3, 1./5, 1./5, 3./5 } ;
    double w[] = { -27./96, 25./96, 25./96, 25./96};
    npoints = 4;
    xg.insert(xg.end(), x, x + npoints );
    yg.insert(yg.end(), y, y + npoints );
    wg.insert(wg.end(), w, w + npoints );
    return npoints;
  }
  else if(n==4)
  {
    // Deg 4 : 6 points
    double a = 0.445948490915965 ;
    double b = 0.091576213509771 ;
    double c = 0.111690794839005 ;
    double d = 0.054975871827661 ;
    double x[] = { a, 1-2*a, a, b, 1-2*b, b } ;
    double y[] = { a, a, 1-2*a, b, b, 1-2*b } ;
    double w[] = { c, c, c, d, d, d };
    npoints = 6;
    xg.insert(xg.end(), x, x + npoints );
    yg.insert(yg.end(), y, y + npoints );
    wg.insert(wg.end(), w, w + npoints );
    return npoints;
  }
  else if (n==5)
  {
    // Deg 5 : 7 points
    double a =  (6.+std::sqrt(15))/21. ;
    double b = 4./7-a;
    double c = (155.+std::sqrt(15))/2400. ;
    double d = 31./240-c ;
    double x[] = { 1./3, a, 1-2*a, a, b, 1-2*b, b}  ;
    double y[] = { 1./3, a, a, 1-2*a, b,b, 1-2*b } ;
    double w[] = { 9./80, c,c,c, d,d,d } ;
    npoints = 7;
    xg.insert(xg.end(), x, x + npoints );
    yg.insert(yg.end(), y, y + npoints );
    wg.insert(wg.end(), w, w + npoints );
    return npoints;
  }
  else if (n==6)
  {
    // Deg = 6 : 12 points
    double x[] = {  0.5014265  ,
                  0.2492867  ,
                  0.2492867  ,
                  0.8738220  ,
                  0.06308901 ,
                  0.06308901 ,
                  0.6365025  ,
                  0.6365025  ,
                  0.05314505 ,
                  0.05314505 ,
                  0.3103525  ,
                  0.3103525 } ;
    double y[] =  { 0.2492867  ,
                  0.5014265  ,
                  0.2492867  ,
                  0.06308901  ,
                  0.8738220  ,
                  0.06308901 ,
                  0.05314505 ,
                  0.3103525 ,
                  0.6365025 ,
                  0.3103525 ,
                  0.6365025 ,
                  0.05314505 } ;
    double w[] = {  0.05839314 ,
                  0.05839314 ,
                  0.05839314 ,
                  0.02542245 ,
                  0.02542245 ,
                  0.02542245 ,
                  0.04142554 ,
                  0.04142554 ,
                  0.04142554 ,
                  0.04142554 ,
                  0.04142554 ,
                  0.04142554 } ;
      npoints = 12;
      xg.insert(xg.end(), x, x + npoints );
      yg.insert(yg.end(), y, y + npoints );
      wg.insert(wg.end(), w, w + npoints );
      return npoints;
  }
  else if (n==7)
  {
    // Deg 7 : 13 points
    double x[] =  { 0.3333333,
              0.4793081  ,
              0.2603460  ,
              0.2603460  ,
              0.8697398  ,
              0.06513010 ,
              0.06513010 ,
              0.6384442  ,
              0.6384442  ,
              0.04869032 ,
              0.04869032 ,
              0.3128655  ,
              0.3128655  };
    double y[] = { 0.3333333  ,
              0.2603460  ,
              0.4793081  ,
              0.2603460  ,
              0.06513010 ,
              0.8697398  ,
              0.06513010 ,
              0.04869032 ,
              0.3128655  ,
              0.6384442  ,
              0.3128655  ,
              0.6384442  ,
              0.04869032 } ;
    double w[] = { -0.07478502 ,
              0.08780763 ,
              0.08780763 ,
              0.08780763 ,
              0.02667362 ,
              0.02667362 ,
              0.02667362 ,
              0.03855688 ,
              0.03855688 ,
              0.03855688 ,
              0.03855688 ,
              0.03855688 ,
              0.03855688 } ;
      npoints = 13;
      xg.insert(xg.end(), x, x + npoints );
      yg.insert(yg.end(), y, y + npoints );
      wg.insert(wg.end(), w, w + npoints );
      return npoints;
  }
  else {
    exit (EXIT_FAILURE);
  }
}

Eigen::MatrixXd GaussTetrahedron(std::vector<double> &wg, int &npoints, int n)
{
    // Returns barycentric coordinates of the Gauss integration points 
    if (n==1)
    {
        npoints = 1; 
        double c[npoints][4] = {{1./4,1./4,1./4,1./4}};
        double w[npoints] = {1.} ;
        wg.insert(wg.end(), w, w + npoints );  
        Eigen::MatrixXd bc(npoints,4);
        for (int i=0; i<npoints; i++)
        {
            for (int j=0; j<4; j++)
            {
                bc(i,j)= c[i][j]; 
            }        
        }
        return bc; 
        
    }
    else if (n==2)
    {
        npoints = 4; 
        double w[npoints] = {1./4,1./4,1./4,1./4}; 
        wg.insert(wg.end(), w, w + npoints ); 
        Eigen::MatrixXd bc(npoints,4); 
        double a = (5.-std::sqrt(5))/20. ; 
        double c[npoints][4] ={{1-3*a,a,a,a},
                               {a,1-3*a,a,a},
                               {a,a,1-3*a,a},
                               {a,a,a,1-3*a}}; 
        for (int i=0; i<npoints; i++)
        {
            for (int j=0; j<4; j++)
            {
                bc(i,j)= c[i][j]; 
            }        
        }
        return bc;                                
    }
    else if (n==3)
    {
        npoints = 5; 
        double w[npoints] = {-4./5, 9./20, 9./20, 9./20,9./20};
        wg.insert(wg.end(), w, w + npoints ); 
        Eigen::MatrixXd bc(npoints,4); 
        double c[npoints][4] = {{1./4,1./4,1./4,1./4},
                                {1./2,1./6,1./6,1./6},
                                {1./6,1./2,1./6,1./6},
                                {1./6,1./6,1./2,1./6},
                                {1./6,1./6,1./6,1./2}} ;
        for (int i=0; i<npoints; i++)
        {
            for (int j=0; j<4; j++)
            {
                bc(i,j)= c[i][j]; 
            }        
        }
        return bc;                                 
    
    }
    else if (n==4 || n==5)
    {
        npoints = 15; 
        double a1 = (7.-std::sqrt(15))/34.  ; 
        double a2 = (7.+std::sqrt(15))/34.  ; 
        double a = (10.-2*std::sqrt(5))/40. ; 
        double w1 = (2665.+14*std::sqrt(15))/37800. ; 
        double w2 = (2665.-14*std::sqrt(15))/37800. ; 
        double w  = 10./189. ; 
        double wG[npoints] = {16./135., w1,w1,w1,w1, w2,w2,w2,w2, w,w,w,w,w,w};
        wg.insert(wg.end(), wG, wG + npoints ); 
        Eigen::MatrixXd bc(npoints,4);   
        double c[npoints][4] = {{1./4,1./4,1./4,1./4},
                                {1-2*a1,a1,a1,a1},
                                {a1,1-2*a1,a1,a1},
                                {a1,a1,1-2*a1,a1},
                                {a1,a1,a1,1-2*a1},
                                {1-2*a2,a2,a2,a2},
                                {a2,1-2*a2,a2,a2},
                                {a2,a2,1-2*a2,a2},
                                {a2,a2,a2,1-2*a2},
                                {1./2-a,1./2-a,a,a},
                                {1.2-a,a,1./2-a,a},
                                {1./2-a,a,a,1./2-a},
                                {a,1./2-a,1./2-a,a},
                                {a,a,1./2-a,1./2-a},
                                {a,1./2-a,1./2-a,a}};
        for (int i=0; i<npoints; i++)
        {
            for (int j=0; j<4; j++)
            {
                bc(i,j)= c[i][j]; 
            }        
        }
        return bc;                                  
                                                    
    }
    else 
    {
      exit (EXIT_FAILURE);
    }    

}

    





///////////////////////////  IMAGE CBSPLINE2 INTERPOLATION ROOTINES //////////////////////////////////////

// Image evaluation with quadratic cardinal B-splines 

double EvaluateCBspline2OnOnePoint(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double xmax, double ymax, double zmax, 
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double x, double y, double z ) 
{

    double xcopy = x ;  
    double ycopy = y ; 
    double zcopy = z ; 
    
    if ( x < xmin + 3*dx/2 )
    {
        xcopy = xmin + 3*dx/2 ;    
    }
    else if (x > xmax - 3*dx/2)
    {
        xcopy = xmax - 3*dx/2 ; 
    }
    
    if ( y < ymin + 3*dy/2 )
    {
        ycopy = ymin + 3*dy/2 ;    
    }
    else if (y > ymax - 3*dy/2)
    {
        ycopy = ymax - 3*dy/2 ; 
    }
    
    if ( z < zmin + 3*dz/2 )
    {
        zcopy = zmin + 3*dz/2 ;    
    }
    else if (z > zmax - 3*dz/2)
    {
        zcopy = zmax - 3*dz/2 ; 
    }        
   
    int tx,ty,tz; 
    int I,J,K; 
    double s; 

	tx = std::floor( (xcopy-(xmin+dx/2))/dx ) ;
	ty = std::floor( (ycopy-(ymin+dy/2))/dy ) ;
	tz = std::floor( (zcopy-(zmin+dz/2))/dz ) ;
	s = 0; 
	for (int k=tz-1; k<=tz+2; k++)
	{
  		J = sj-1-k; 
  		for (int j=ty-1; j<=ty+2; j++)
      	{
          	I = j; 
          	for (int i=tx-1; i<=tx+2; i++)
          	{
              	K = i; 
              	s = s + image[K+J*sk+I*sk*sj]*CardinalBspline2((xcopy-xc[i])/dx)*CardinalBspline2((ycopy-yc[j])/dy)*CardinalBspline2((zcopy-zc[k])/dz); 
          	}
      	}
	}
	return s ; 
}

 
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
            					    double* dvdz)
{
    double** Bx      = array_2d(sx,4);
    double** dBxdx   = array_2d(sx,4);
    int* spanx = new int[sx];
    double** By      = array_2d(sy,4);
    double** dBydy   = array_2d(sy,4);
    int* spany = new int[sy];
    double** Bz   = array_2d(sz,4);
    double** dBzdz   = array_2d(sz,4); 
    int* spanz = new int[sz];     
	double gl; 
	int I,J,K; 
	// First evaluating univariate basis functions
	int k; 
	for ( int i=0; i < sx; i++)
	{
		spanx[i] = std::floor((x[i]-(xmin+dx/2))/dx) ;
		k = 0 ;
		for ( int j= spanx[i] -1; j<= spanx[i] + 2 ; j++ )
		{
			Bx[i][k]    = CardinalBspline2((x[i]-xc[j])/dx); 
			dBxdx[i][k] = CardinalBsplineDer2((x[i]-xc[j])/dx); 
			k+=1;
		}
	}  
	for ( int i=0; i < sy; i++)
	{
		spany[i] = std::floor((y[i]-(ymin+dy/2))/dy) ;
		k = 0 ;
		for ( int j= spany[i] -1; j<= spany[i] + 2 ; j++ )
		{
			By[i][k]    = CardinalBspline2((y[i]-yc[j])/dy); 
			dBydy[i][k] = CardinalBsplineDer2((y[i]-yc[j])/dy); 
			k+=1;
		}
	}  
	for ( int i=0; i < sz; i++)
	{
		spanz[i] = std::floor((z[i]-(zmin+dz/2))/dz) ;
		k = 0 ;
		for ( int j= spanz[i] -1; j<= spanz[i] + 2 ; j++ )
		{
			Bz[i][k]    = CardinalBspline2((z[i]-zc[j])/dz); 
			dBzdz[i][k] = CardinalBsplineDer2((z[i]-zc[j])/dz); 
			k+=1;
		}
	}
	int m = 0 ; 
	int li,lj,lk; 
	double s  ;  // Value of the volume 
	double sxx ;  // Gradient of the volume in x direction 
	double syy ;  // Gradient of the volume in y direction 
	double szz ;  // Gradient of the volume in z direction  
	
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
				
				lk = 0; 
				for ( int k =spanz[r]-1; k  <= spanz[r]+2 ; k++ )
				{
    				J = sj-1-k; 
    				lj =0; 
		    		for ( int j=spany[q]-1; j <= spany[q]+2; j++ )
					{
    					I = j;
    					li =0;  
						for ( int i =spanx[p]-1; i <= spanx[p]+2; i++)
						{
    						K = i;  
        					gl = image[K+J*sk+I*sk*sj];   
        					s   = s   + Bx[p][li]*By[q][lj]*Bz[r][lk]*gl; 
        					sxx = sxx + dBxdx[p][li]*By[q][lj]*Bz[r][lk]*gl;
        					syy = syy + Bx[p][li]*dBydy[q][lj]*Bz[r][lk]*gl;  
							szz = szz + Bx[p][li]*By[q][lj]*dBzdz[r][lk]*gl; 
                            li++;
						} 
    					lj++; 
					}
					lk++; 
				}       
				v[m]    = s   ;
				dvdx[m] = sxx/dx ;  
				dvdy[m] = syy/dy ; 
				dvdz[m] = szz/dz ; 
				m++ ;  
			}
		}
	}
	delete[] spanx;
	delete[] spany;
	delete[] spanz;
	delete_array_2d(Bx);
	delete_array_2d(dBxdx);
	delete_array_2d(By);
	delete_array_2d(dBydy);
	delete_array_2d(Bz);
	delete_array_2d(dBzdz);
} 
        					            					        		        
///////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////  IMAGE CBSPLINE3 INTERPOLATION ROOTINES //////////////////////////////////////

// Image evaluation with cubic cardinal B-splines 

double EvaluateCBspline3OnOnePoint(double* image, int si, int sj, int sk, 
                                   double xmin, double ymin, double zmin,
                                   double xmax, double ymax, double zmax, 
                                   double dx, double dy, double dz,  
        					       double* xc,  int sxc, 
        					       double* yc,  int syc,
        					       double* zc,  int szc,
					               double x, double y, double z ) 
{

    double xcopy = x ;  
    double ycopy = y ; 
    double zcopy = z ; 
    
    if ( x < xmin + 3*dx/2 )
    {
        xcopy = xmin + 3*dx/2 ;    
    }
    else if (x > xmax - 3*dx/2)
    {
        xcopy = xmax - 3*dx/2 ; 
    }
    
    if ( y < ymin + 3*dy/2 )
    {
        ycopy = ymin + 3*dy/2 ;    
    }
    else if (y > ymax - 3*dy/2)
    {
        ycopy = ymax - 3*dy/2 ; 
    }
    
    if ( z < zmin + 3*dz/2 )
    {
        zcopy = zmin + 3*dz/2 ;    
    }
    else if (z > zmax - 3*dz/2)
    {
        zcopy = zmax - 3*dz/2 ; 
    }        
   
    int tx,ty,tz; 
    int I,J,K; 
    double s; 

	tx = std::floor( (xcopy-(xmin+dx/2))/dx ) ;
	ty = std::floor( (ycopy-(ymin+dy/2))/dy ) ;
	tz = std::floor( (zcopy-(zmin+dz/2))/dz ) ;
	s = 0; 
	for (int k=tz-1; k<=tz+2; k++)
	{
  		J = sj-1-k; 
  		for (int j=ty-1; j<=ty+2; j++)
      	{
          	I = j; 
          	for (int i=tx-1; i<=tx+2; i++)
          	{
              	K = i; 
              	s = s + image[K+J*sk+I*sk*sj]*CardinalBspline3((xcopy-xc[i])/dx)*CardinalBspline3((ycopy-yc[j])/dy)*CardinalBspline3((zcopy-zc[k])/dz); 
          	}
      	}
	}
	return s ; 
}

 
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
            					    double* dvdz)
{
    double** Bx      = array_2d(sx,4);
    double** dBxdx   = array_2d(sx,4);
    int* spanx = new int[sx];
    double** By      = array_2d(sy,4);
    double** dBydy   = array_2d(sy,4);
    int* spany = new int[sy];
    double** Bz   = array_2d(sz,4);
    double** dBzdz   = array_2d(sz,4); 
    int* spanz = new int[sz];     
	double gl; 
	int I,J,K; 
	// First evaluating univariate basis functions
	int k; 
	for ( int i=0; i < sx; i++)
	{
		spanx[i] = std::floor((x[i]-(xmin+dx/2))/dx) ;
		k = 0 ;
		for ( int j= spanx[i] -1; j<= spanx[i] + 2 ; j++ )
		{
			Bx[i][k]    = CardinalBspline3((x[i]-xc[j])/dx); 
			dBxdx[i][k] = CardinalBsplineDer3((x[i]-xc[j])/dx); 
			k+=1;
		}
	}  
	for ( int i=0; i < sy; i++)
	{
		spany[i] = std::floor((y[i]-(ymin+dy/2))/dy) ;
		k = 0 ;
		for ( int j= spany[i] -1; j<= spany[i] + 2 ; j++ )
		{
			By[i][k]    = CardinalBspline3((y[i]-yc[j])/dy); 
			dBydy[i][k] = CardinalBsplineDer3((y[i]-yc[j])/dy); 
			k+=1;
		}
	}  
	for ( int i=0; i < sz; i++)
	{
		spanz[i] = std::floor((z[i]-(zmin+dz/2))/dz) ;
		k = 0 ;
		for ( int j= spanz[i] -1; j<= spanz[i] + 2 ; j++ )
		{
			Bz[i][k]    = CardinalBspline3((z[i]-zc[j])/dz); 
			dBzdz[i][k] = CardinalBsplineDer3((z[i]-zc[j])/dz); 
			k+=1;
		}
	}
	int m = 0 ; 
	int li,lj,lk; 
	double s  ;  // Value of the volume 
	double sxx ;  // Gradient of the volume in x direction 
	double syy ;  // Gradient of the volume in y direction 
	double szz ;  // Gradient of the volume in z direction  
	
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
				
				lk = 0; 
				for ( int k =spanz[r]-1; k  <= spanz[r]+2 ; k++ )
				{
    				J = sj-1-k; 
    				lj =0; 
		    		for ( int j=spany[q]-1; j <= spany[q]+2; j++ )
					{
    					I = j;
    					li =0;  
						for ( int i =spanx[p]-1; i <= spanx[p]+2; i++)
						{
    						K = i;  
        					gl = image[K+J*sk+I*sk*sj];   
        					s   = s   + Bx[p][li]*By[q][lj]*Bz[r][lk]*gl; 
        					sxx = sxx + dBxdx[p][li]*By[q][lj]*Bz[r][lk]*gl;
        					syy = syy + Bx[p][li]*dBydy[q][lj]*Bz[r][lk]*gl;  
							szz = szz + Bx[p][li]*By[q][lj]*dBzdz[r][lk]*gl; 
                            li++;
						} 
    					lj++; 
					}
					lk++; 
				}       
				v[m]    = s   ;
				dvdx[m] = sxx/dx ;  
				dvdy[m] = syy/dy ; 
				dvdz[m] = szz/dz ; 
				m++ ;  
			}
		}
	}
	delete[] spanx;
	delete[] spany;
	delete[] spanz;
	delete_array_2d(Bx);
	delete_array_2d(dBxdx);
	delete_array_2d(By);
	delete_array_2d(dBydy);
	delete_array_2d(Bz);
	delete_array_2d(dBzdz);
} 
        					            					        		        
///////////////////////////////////////////////////////////////////////////////////////////////
int PointIsInTetrahedron1(double* xn, double* yn, double*zn, double x, double y, double z, double eps )
{
    // Returns 1 if the point is inside the tetrahedron 
    // Returns 0 if the point is outside the tetrahedron  
    // PAY CAUTION WITH EPS
    double A[16] ;  
    double u[4] = { x +eps , y+eps, z+eps, 1  }; 
    double l[4]  ; 
    double invA[16] ; 
    double detA ;
    
    A[12] = 1 ;
    A[13] = 1 ; 
    A[14] = 1 ; 
    A[15] = 1 ;     

    A[0] =  xn[0] ; 
    A[4] =  yn[0] ; 
    A[8] =  zn[0] ;  
    
    A[1] =  xn[1] ; 
    A[5] =  yn[1] ; 
    A[9] =  zn[1] ;
    
    A[2]  =  xn[2] ; 
    A[6]  =  yn[2] ; 
    A[10] =  zn[2] ;
    
    A[3]  =  xn[3] ;
    A[7]  =  yn[3] ;
    A[11] =  zn[3] ; 
    
    /******************************************************************/  
    // Inverting the (4x4) matrix 
    // Solving the system Al=u (Convex hull proprety all l_i must be greater than 0)
    invA[0] = A[5]  * A[10] * A[15] - 
             A[5]  * A[11] * A[14] - 
             A[9]  * A[6]  * A[15] + 
             A[9]  * A[7]  * A[14] +
             A[13] * A[6]  * A[11] - 
             A[13] * A[7]  * A[10];

    invA[4] = -A[4]  * A[10] * A[15] + 
              A[4]  * A[11] * A[14] + 
              A[8]  * A[6]  * A[15] - 
              A[8]  * A[7]  * A[14] - 
              A[12] * A[6]  * A[11] + 
              A[12] * A[7]  * A[10];

    invA[8] = A[4]  * A[9] * A[15] - 
             A[4]  * A[11] * A[13] - 
             A[8]  * A[5] * A[15] + 
             A[8]  * A[7] * A[13] + 
             A[12] * A[5] * A[11] - 
             A[12] * A[7] * A[9];

    invA[12] = -A[4]  * A[9] * A[14] + 
               A[4]  * A[10] * A[13] +
               A[8]  * A[5] * A[14] - 
               A[8]  * A[6] * A[13] - 
               A[12] * A[5] * A[10] + 
               A[12] * A[6] * A[9];

    invA[1] = -A[1]  * A[10] * A[15] + 
              A[1]  * A[11] * A[14] + 
              A[9]  * A[2] * A[15] - 
              A[9]  * A[3] * A[14] - 
              A[13] * A[2] * A[11] + 
              A[13] * A[3] * A[10];

    invA[5] = A[0]  * A[10] * A[15] - 
             A[0]  * A[11] * A[14] - 
             A[8]  * A[2] * A[15] + 
             A[8]  * A[3] * A[14] + 
             A[12] * A[2] * A[11] - 
             A[12] * A[3] * A[10];

    invA[9] = -A[0]  * A[9] * A[15] + 
              A[0]  * A[11] * A[13] + 
              A[8]  * A[1] * A[15] - 
              A[8]  * A[3] * A[13] - 
              A[12] * A[1] * A[11] + 
              A[12] * A[3] * A[9];

    invA[13] = A[0]  * A[9] * A[14] - 
              A[0]  * A[10] * A[13] - 
              A[8]  * A[1] * A[14] + 
              A[8]  * A[2] * A[13] + 
              A[12] * A[1] * A[10] - 
              A[12] * A[2] * A[9];

    invA[2] = A[1]  * A[6] * A[15] - 
             A[1]  * A[7] * A[14] - 
             A[5]  * A[2] * A[15] + 
             A[5]  * A[3] * A[14] + 
             A[13] * A[2] * A[7] - 
             A[13] * A[3] * A[6];

    invA[6] = -A[0]  * A[6] * A[15] + 
              A[0]  * A[7] * A[14] + 
              A[4]  * A[2] * A[15] - 
              A[4]  * A[3] * A[14] - 
              A[12] * A[2] * A[7] + 
              A[12] * A[3] * A[6];

    invA[10] = A[0]  * A[5] * A[15] - 
              A[0]  * A[7] * A[13] - 
              A[4]  * A[1] * A[15] + 
              A[4]  * A[3] * A[13] + 
              A[12] * A[1] * A[7] - 
              A[12] * A[3] * A[5];

    invA[14] = -A[0]  * A[5] * A[14] + 
               A[0]  * A[6] * A[13] + 
               A[4]  * A[1] * A[14] - 
               A[4]  * A[2] * A[13] - 
               A[12] * A[1] * A[6] + 
               A[12] * A[2] * A[5];

    invA[3] = -A[1] * A[6] * A[11] + 
              A[1] * A[7] * A[10] + 
              A[5] * A[2] * A[11] - 
              A[5] * A[3] * A[10] - 
              A[9] * A[2] * A[7] + 
              A[9] * A[3] * A[6];

    invA[7] = A[0] * A[6] * A[11] - 
             A[0] * A[7] * A[10] - 
             A[4] * A[2] * A[11] + 
             A[4] * A[3] * A[10] + 
             A[8] * A[2] * A[7] - 
             A[8] * A[3] * A[6];

    invA[11] = -A[0] * A[5] * A[11] + 
               A[0] * A[7] * A[9] + 
               A[4] * A[1] * A[11] - 
               A[4] * A[3] * A[9] - 
               A[8] * A[1] * A[7] + 
               A[8] * A[3] * A[5];

    invA[15] = A[0] * A[5] * A[10] - 
              A[0] * A[6] * A[9] - 
              A[4] * A[1] * A[10] + 
              A[4] * A[2] * A[9] + 
              A[8] * A[1] * A[6] - 
              A[8] * A[2] * A[5];

    detA = A[0] * invA[0] + A[1] * invA[4] + A[2] * invA[8] + A[3] * invA[12];    
    
    l[0] = (invA[0]*u[0]  + invA[1]*u[1]  + invA[2]*u[2]  + invA[3]*u[3])/detA ; 
    l[1] = (invA[4]*u[0]  + invA[5]*u[1]  + invA[6]*u[2]  + invA[7]*u[3])/detA ; 
    l[2] = (invA[8]*u[0]  + invA[9]*u[1]  + invA[10]*u[2] + invA[11]*u[3])/detA ; 
    l[3] = (invA[12]*u[0] + invA[13]*u[1] + invA[14]*u[2] + invA[15]*u[3])/detA ;         
    /******************************************************************/
    if ( l[0]>0 && l[1]>0 && l[2]>0 && l[3]>0 )
    {
        return 1 ; 
    }
    else
    {
        return 0 ; 
    }
}

int PointIsInTetrahedron2(double* xn, double* yn, double*zn, double x, double y, double z, double eps )
{
    // Evaluating the zero level-set equation at the point 
    double nx,ny,nz; 
    double x1,y1,z1,x2,y2,z2,x3,y3,z3 ; 
    double phi ;  // Level-set 
    // Each face (triangle) is defined by 3 nodes
    // ATTENTION: ALL ELEMENTS MUST HAVE THE SAME ORIENTATION 
    // IN ADDITION: THE ORIENTATION MUST BE CHECKED BEFORE USING THIS ROUTINE 
    int i1[4] = {1,0,0,0} ; 
    int i2[4] = {2,3,1,2} ; 
    int i3[4] = {3,2,3,1} ; 
    phi = - 1 ; 
    int i = 0 ; 
    while (phi<eps && i<4)
    {
        x1 = xn[i1[i]] ; y1 = yn[i1[i]] ; z1 = zn[i1[i]] ; 
        x2 = xn[i2[i]] ; y2 = yn[i2[i]] ; z2 = zn[i2[i]] ; 
        x3 = xn[i3[i]] ; y3 = yn[i3[i]] ; z3 = zn[i3[i]] ; 
        
        nx = (y2-y1)*(z3-z2) - (y3-y2)*(z2-z1) ; 
        ny = (x3-x2)*(z2-z1) - (x2-x1)*(z3-z2) ; 
        nz = (x2-x1)*(y3-y2) - (x3-x2)*(y2-y1) ; 

        phi = nx*(x-x2)+ny*(y-y2)+nz*(z-z2) ;
        i++ ;      
    }
    if (phi<eps)
    {
        return 1 ; 
    }
    else 
    {
        return 0 ; 
    }   
}









