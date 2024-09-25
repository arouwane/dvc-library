#include "imageRoutines.h"

void ComputeL2ProjLumpedCoefficients(double* image, int si, int sj, int sk, 
                           double xmin, double ymin, double zmin, double zmax,                     
                           double dx, double dy, double dz, 
    					   int deg_xi, int deg_eta,int deg_zeta,
					       double* knotXi,  int sknotXi, 
					       double* knotEta, int sknotEta ,
					       double* knotZeta, int sknotZeta,  
					       double* lsc, int s_lsc  )  
{
  		int nbg_xi   = deg_xi   + 1;  // Number of Gauss points in xi direction 
  		int nbg_eta  = deg_eta  + 1;  // Number of Gauss points in eta direction 
  		int nbg_zeta = deg_zeta + 1;  // Number of Gauss points in zeta direction 
    	double* xig   = array_1d(nbg_xi); // Arrays of gauss points and weights  
		double* wgxi  = array_1d(nbg_xi); 
		double* etag  = array_1d(nbg_eta); 
		double* wgeta = array_1d(nbg_eta);   
		double* zetag = array_1d(nbg_zeta);  
		double* wgzeta = array_1d(nbg_zeta); 
		gauleg(-1,1,xig,wgxi,nbg_xi);
		gauleg(-1,1,etag,wgeta,nbg_eta);
		gauleg(-1,1,zetag,wgzeta,nbg_zeta);  	
		// For this method the number of elements is equal to the number of voxels 	
		int nbf_xi   =  sk   +  deg_xi  ; // Number of basis functions in xi direction   
        int nbf_eta  =  si   +  deg_eta ; // Number of basis functions in eta direction 
		int nbf_zeta =  sj   +  deg_zeta ; // Number of basis functions in zeta direction
		// Note that the size of lsc is s_lsc = nbf_xi * nbf_eta *nbf_zeta
		// Evaluation in 1D dimension due to the tensor product nature 		
		int nbge_xi   = nbg_xi*nbg_xi ;
		int nbge_eta  = nbg_eta*nbg_eta ;
		int nbge_zeta = nbg_zeta*nbg_zeta ; 		
		Eigen::ArrayXXd NxiArray    = Eigen::ArrayXXd::Zero(nbf_xi,nbge_xi) ; 
		Eigen::ArrayXXd NetaArray   = Eigen::ArrayXXd::Zero(nbf_eta,nbge_eta); 
		Eigen::ArrayXXd NzetaArray  = Eigen::ArrayXXd::Zero(nbf_zeta,nbge_zeta);
        // The vector xigArray will contain the gauss point coordinates
        // Very important here the physical space starts from 0 
        // If knotXi[0] <0, the arrays xigArray should be initialised with a constant outside the domain boundary 
        // And the test must be performed when computing num and denom 
        // The whole vector is initialized with a value outside the domain ( for example knotXi[0]-1)
        //Eigen::ArrayXXd xigArray    = Eigen::ArrayXXd::Constant(nbf_xi, nbge_xi,  knotXi[0]   - 1); 
        //Eigen::ArrayXXd etagArray   = Eigen::ArrayXXd::Constant(nbf_eta,nbge_eta, knotEta[0]  - 1); 
        //Eigen::ArrayXXd zetagArray  = Eigen::ArrayXXd::Constant(nbf_zeta,nbge_zeta, knotZeta[0] - 1); 
        Eigen::ArrayXXd xigArray     = Eigen::ArrayXXd::Zero(nbf_xi, nbge_xi); 
        Eigen::ArrayXXd etagArray    = Eigen::ArrayXXd::Zero(nbf_eta,nbge_eta); 
        Eigen::ArrayXXd zetagArray   = Eigen::ArrayXXd::Zero(nbf_zeta,nbge_zeta); 
		Eigen::ArrayXXd wxigArray    = Eigen::ArrayXXd::Zero(nbf_xi,nbge_xi)  ; 
		Eigen::ArrayXXd wetagArray   = Eigen::ArrayXXd::Zero(nbf_eta,nbge_eta);
		Eigen::ArrayXXd wzetagArray  = Eigen::ArrayXXd::Zero(nbf_zeta,nbge_zeta);
		Eigen::ArrayXi  counter_xi   = Eigen::ArrayXi::Zero(nbf_xi); 
		Eigen::ArrayXi  counter_eta  = Eigen::ArrayXi::Zero(nbf_eta);
		Eigen::ArrayXi  counter_zeta = Eigen::ArrayXi::Zero(nbf_zeta); 
		int span ; 
		double xg,yg,zg; 
		double *Nxi   = new double[nbg_xi]; 
		double *Neta  = new double[nbg_eta];
		double *Nzeta = new double[nbg_zeta]; 
		int conn ; 
 	
		// xi   direction
		for ( int exi = 0 ; exi < sk ; exi++ ) 
		{
			span = exi + deg_xi ; 
			for ( int ii=0; ii < nbg_xi; ii++ ) 
			{
				xg = knotXi[span] + 0.5*dx*(xig[ii]+1); 
				basisfuns(span, xg, deg_xi, knotXi, Nxi); 
				for ( int l =0 ; l < nbg_xi; l++ ) 
				{
					conn = span - deg_xi + l; 
					NxiArray  ( conn , counter_xi(conn) ) = Nxi[l]; 
					xigArray  ( conn , counter_xi(conn) ) = xg; 
					wxigArray ( conn , counter_xi(conn) ) = wgxi[ii]; 
					counter_xi(conn) = counter_xi(conn) + 1 ;  
				} 
			}
		}
		// eta  direction 
		for ( int eeta = 0; eeta < si ; eeta ++) 
		{
			span = eeta + deg_eta; 
			for ( int jj = 0; jj < nbg_eta; jj ++ )
			{
				yg = knotEta[span] + 0.5*dy*(etag[jj]+1); 
				basisfuns(span, yg, deg_eta, knotEta, Neta); 
				for ( int l=0; l<nbg_eta; l++ )
				{
					conn = span - deg_eta + l; 
					NetaArray( conn, counter_eta(conn) ) = Neta[l]; 
					etagArray( conn, counter_eta(conn) ) = yg; 
					wetagArray( conn, counter_eta(conn) ) = wgeta[jj]; 
					counter_eta(conn) = counter_eta(conn) + 1 ;
				} 
			}
		}
		// zeta direction 
		for ( int ezeta = 0; ezeta < sj ; ezeta ++) 
		{
			span = ezeta + deg_zeta; 
			for ( int kk = 0; kk < nbg_zeta; kk ++ )
			{
				zg = knotZeta[span] + 0.5*dz*(zetag[kk]+1); 
				basisfuns(span, zg, deg_zeta, knotZeta, Nzeta); 
				for ( int l=0; l<nbg_zeta; l++ )
				{
					conn = span - deg_zeta + l; 
					NzetaArray( conn, counter_zeta(conn) ) = Nzeta[l]; 
					zetagArray( conn, counter_zeta(conn) ) = zg; 
					wzetagArray( conn, counter_zeta(conn) ) = wgzeta[kk]; 
					counter_zeta(conn) = counter_zeta(conn) + 1 ;
				} 
			}
		}		 		
 
 
        double mes = dx*dy*dz;  // Constant measure of an element = voxel volume 
        
        # pragma omp parallel 
        {
            #pragma omp for 
        		for ( int kbf=0; kbf<nbf_zeta; kbf++  )
        		{
        			for ( int jbf=0; jbf<nbf_eta; jbf++ ) 
        			{
        				for (int ibf =0; ibf < nbf_xi; ibf++ )
        				{	
        					// Current basis function
        					double num = 0; 
        					double denom = 0; 
        					double N ; 
        					double wg ; 
        					double g ;
        					double x,y,z ;  
        		            //loop over elements supported by this function 
        					for ( int kk = 0; kk < nbge_zeta; kk ++ ) 
        					{
        						for ( int jj = 0 ; jj < nbge_eta; jj ++ )
        						{
        							for (int ii = 0; ii < nbge_xi; ii ++)
        							{
        								x = xigArray(ibf,ii); 
        								y = etagArray(jbf,jj);
        								z = zetagArray(kbf,kk);  
        					            //if (  (x !=knotXi[0]-1) && (y !=knotEta[0] - 1) && (z !=knotZeta[0] - 1) )
                                        //{
                                          N  = NxiArray(ibf,ii)*NetaArray(jbf,jj)*NzetaArray(kbf,kk); 
                                          wg = wxigArray(ibf,ii)*wetagArray(jbf,jj)*wzetagArray(kbf,kk); 
                                          g  = GetGrayLevelOfPoint3D(image,si,sj,sk,x,y,z,xmin,ymin,zmin,zmax,dx,dy,dz ); 
                                        //}
                                        //else 
                                        //{
                                        //  N =  0;
                                        //  wg = 0;
                                        //  g  = 0;
                                        //}                                
                                        num = num + wg*g*N*mes/8;
                                        denom = denom + wg*N*mes/8;								
        							} 
        						} 
        					}
        					lsc[ibf + jbf*nbf_xi + kbf*nbf_xi*nbf_eta ] = num/denom ;   
        				} 
        			}
        		}
        } 
		delete[] xig;
		delete[] etag;
		delete[] zetag; 
		delete[] wgxi; 
		delete[] wgeta; 
		delete[] wgzeta; 
		delete[] Nxi; 
		delete[] Neta; 
		delete[] Nzeta; 			 
}


void ComputeL2ProjLumpedCoefficientsSumFact(double* image, int si, int sj, int sk, 
                           double xmin, double ymin, double zmin, double zmax,                     
                           double dx, double dy, double dz, 
    					   int deg_xi, int deg_eta,int deg_zeta,
					       double* knotXi,  int sknotXi, 
					       double* knotEta, int sknotEta ,
					       double* knotZeta, int sknotZeta,  
					       double* lsc, int s_lsc  )  
{
  		int nbg_xi   = std::ceil(  (deg_xi+1)/2   );  // Number of Gauss points in xi direction 
  		int nbg_eta  = std::ceil(  (deg_eta+1)/2  ); // Number of Gauss points in eta direction 
  		int nbg_zeta = std::ceil(  (deg_zeta+1)/2 );  // Number of Gauss points in zeta direction 
    	double* xig   = array_1d(nbg_xi); // Arrays of gauss points and weights  
		double* wgxi  = array_1d(nbg_xi); 
		double* etag  = array_1d(nbg_eta); 
		double* wgeta = array_1d(nbg_eta);   
		double* zetag = array_1d(nbg_zeta);  
		double* wgzeta = array_1d(nbg_zeta); 
		gauleg(-1,1,xig,wgxi,nbg_xi);
		gauleg(-1,1,etag,wgeta,nbg_eta);
		gauleg(-1,1,zetag,wgzeta,nbg_zeta);  	
		// For this method the number of elements is equal to the number of voxels 	
		int ne_xi   = sk ; 
		int ne_eta  = si ;  
		int ne_zeta = sj ;  
		int nbf_xi   =  ne_xi    +  deg_xi  ; // Number of basis functions in xi direction   
        int nbf_eta  =  ne_eta   +  deg_eta ; // Number of basis functions in eta direction 
		int nbf_zeta =  ne_zeta  +  deg_zeta ; // Number of basis functions in zeta direction
		
		int nbf_elem_xi   = deg_xi   + 1 ; 
		int nbf_elem_eta  = deg_eta  + 1 ; 
		int nbf_elem_zeta = deg_zeta + 1 ;  
		// Note that the size of lsc is s_lsc = nbf_xi * nbf_eta *nbf_zeta
		// Evaluation in 1D dimension due to the tensor product nature 		
 
		
		Eigen::ArrayXXd intNxiArray   = Eigen::ArrayXXd::Zero(nbf_xi  , nbf_elem_xi   ) ;   
		Eigen::ArrayXXd intNetaArray  = Eigen::ArrayXXd::Zero(nbf_eta , nbf_elem_eta ) ;  
		Eigen::ArrayXXd intNzetaArray = Eigen::ArrayXXd::Zero(nbf_zeta, nbf_elem_zeta ) ;   
		
		Eigen::ArrayXd  elem_coord_xi   = Eigen::ArrayXd::Zero(ne_xi)   ; 
		Eigen::ArrayXd  elem_coord_eta  = Eigen::ArrayXd::Zero(ne_eta)  ;
		Eigen::ArrayXd  elem_coord_zeta = Eigen::ArrayXd::Zero(ne_zeta) ; 
		
		Eigen::ArrayXi  elemsPerBf_xi   = Eigen::ArrayXi::Zero(nbf_xi) ; 
		Eigen::ArrayXi  elemsPerBf_eta  = Eigen::ArrayXi::Zero(nbf_eta) ; 
		Eigen::ArrayXi  elemsPerBf_zeta = Eigen::ArrayXi::Zero(nbf_zeta) ; 
		
		Eigen::ArrayXXi  elemsOfBf_xi   = Eigen::ArrayXXi::Zero(nbf_xi  , nbf_elem_xi   ) ;   
		Eigen::ArrayXXi  elemsOfBf_eta  = Eigen::ArrayXXi::Zero(nbf_eta  , nbf_elem_eta ) ;  
		Eigen::ArrayXXi  elemsOfBf_zeta = Eigen::ArrayXXi::Zero(nbf_zeta  , nbf_elem_zeta ) ;   
		
		
		int span ; 
		double xg,yg,zg; 
		double *Nxi   = new double[nbf_elem_xi]; 
		double *Neta  = new double[nbf_elem_eta];
		double *Nzeta = new double[nbf_elem_zeta]; 
		int conn ; 
				
		// xi   direction
		// Loop over univariate elements 
		for ( int exi = 0 ; exi < ne_xi ; exi++ ) 
		{
			span = exi + deg_xi ; 
			elem_coord_xi(exi) = knotXi[span] + 0.5*dx*(xig[0]+1);
			// Loop over univariate gauss integration points  
			for ( int ii=0; ii < nbg_xi; ii++ ) 
			{
        		xg = knotXi[span] + 0.5*dx*(xig[ii]+1); 	
        		basisfuns(span, xg, deg_xi, knotXi, Nxi); 
				// Loop over univarite basis functions 
				for ( int l =0 ; l < nbf_elem_xi; l++ ) 
				{
					conn = span - deg_xi + l; // Indice of the basis function  
					intNxiArray(conn, elemsPerBf_xi(conn) )  += wgxi[ii]* Nxi[l]*dx/2 ; 
				}
			}			
			for ( int l =0 ; l < nbf_elem_xi; l++ ) 
			{
    			conn = span - deg_xi + l ; 
    			elemsOfBf_xi(conn,elemsPerBf_xi(conn)) = exi ;
    			elemsPerBf_xi(conn) += 1 ; 
			}
		}
		// eta  direction 
		// Loop over univariate elements 
		for ( int eeta = 0 ; eeta < ne_eta ; eeta++ ) 
		{
			span = eeta + deg_eta ; 
			elem_coord_eta(eeta) = knotEta[span] + 0.5*dy*(etag[0]+1);
			// Loop over univariate gauss integration points  
			for ( int ii=0; ii < nbg_eta; ii++ ) 
			{
        		yg = knotEta[span] + 0.5*dy*(etag[ii]+1); 	
        		basisfuns(span, yg, deg_eta, knotEta, Neta); 
				// Loop over univarite basis functions 
				for ( int l =0 ; l < nbf_elem_eta; l++ ) 
				{
					conn = span - deg_eta + l; // Indice of the basis function  
					intNetaArray(conn, elemsPerBf_eta(conn) )  += wgeta[ii]* Neta[l]*dy/2 ; 
				}
			}			
			for ( int l =0 ; l < nbf_elem_eta; l++ ) 
			{
    			conn = span - deg_eta + l ; 
    			elemsOfBf_eta(conn,elemsPerBf_eta(conn)) = eeta ;
    			elemsPerBf_eta(conn) += 1 ; 
			}
		}
		// zeta direction 
		// Loop over univariate elements 
		for ( int ezeta = 0 ; ezeta < ne_zeta ; ezeta++ ) 
		{
			span = ezeta + deg_zeta ; 
			elem_coord_zeta(ezeta) = knotZeta[span] + 0.5*dz*(zetag[0]+1);
			// Loop over univariate gauss integration points  
			for ( int ii=0; ii < nbg_zeta; ii++ ) 
			{
        		zg = knotZeta[span] + 0.5*dz*(zetag[ii]+1); 	
        		basisfuns(span, zg, deg_zeta, knotZeta, Nzeta); 
				// Loop over univarite basis functions 
				for ( int l =0 ; l < nbf_elem_zeta; l++ ) 
				{
					conn = span - deg_zeta + l; // Indice of the basis function  
					intNzetaArray(conn, elemsPerBf_zeta(conn) )  += wgzeta[ii]* Nzeta[l]*dz/2 ; 
				}
			}			
			for ( int l =0 ; l < nbf_elem_zeta; l++ ) 
			{
    			conn = span - deg_zeta + l ; 
    			elemsOfBf_zeta(conn,elemsPerBf_zeta(conn)) = ezeta ;
    			elemsPerBf_zeta(conn) += 1 ; 
			}
		}		 		
      
        # pragma omp parallel 
        {
            #pragma omp for 
        		for ( int kbf=0; kbf<nbf_zeta; kbf++  )
        		{
        			for ( int jbf=0; jbf<nbf_eta; jbf++ ) 
        			{
        				for (int ibf =0; ibf < nbf_xi; ibf++ )
        				{	
        					// Current basis function
        					double num = 0   ; 
        					double denom = 0 ; 
        					double g ; 
        					double intProd_zy ;
        					double intProd_xyz ;  
        					double x,y,z ; 
        					
        					// Loop over univariate elements supported by this function 
        					for (int kk = 0; kk < elemsPerBf_zeta(kbf) ; kk++)
        					{
            					z = elem_coord_zeta( elemsOfBf_zeta(kbf,kk) ); 
            					for (int jj=0; jj< elemsPerBf_eta(jbf) ; jj++)
            					{
                					y = elem_coord_eta( elemsOfBf_eta(jbf,jj) ); 
                					intProd_zy = intNzetaArray(kbf,kk) * intNetaArray(jbf,jj) ; 
                					for (int ii=0; ii< elemsPerBf_xi(ibf) ; ii++)
                					{
                            			x = elem_coord_xi( elemsOfBf_xi(ibf,ii) ); 
        								g  = GetGrayLevelOfPoint3D(image,si,sj,sk,x,y,z,xmin,ymin,zmin,zmax,dx,dy,dz ); 
        								intProd_xyz = intProd_zy * intNxiArray(ibf,ii) ; 
        								num   += g*intProd_xyz  ;
        								denom += intProd_xyz ; 
                					}
            					}
        					}
        					lsc[ibf + jbf*nbf_xi + kbf*nbf_xi*nbf_eta ] = num/denom ;   
        				} 
        			}
        		}
        } 
		delete[] xig;
		delete[] etag;
		delete[] zetag; 
		delete[] wgxi; 
		delete[] wgeta; 
		delete[] wgzeta; 
		delete[] Nxi; 
		delete[] Neta; 
		delete[] Nzeta; 			 
}




void EvaluateTrilinearInterpolation(double* image, int si, int sj, int sk, 
        					        double* knotXi,  int sknotXi, 
        					        double* knotEta, int sknotEta ,
        					        double* knotZeta, int sknotZeta,
					                double* x, int sx,
            			  		    double* y, int sy,
            			  		    double* z, int sz,   
        					        double* v, int s_v )
{
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
	
	double* Nx = array_1d(deg_xi+1);
	double* Ny = array_1d(deg_eta+1);
	double* Nz = array_1d(deg_zeta+1);
	
    double s ;
    
    int I,J,K; 
    
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
					I = spany-deg_eta+j ;
					J = sj -1 - (spanz-deg_zeta+k) ; 
					K = spanx-deg_xi+i ;
					s = s + Nz[k]*Ny[j]*Nx[i]*image[K+J*sk+I*sk*sj];    
				}
			}
		}		
		v[p] = s; 
	}    
	delete[] Nx;
	delete[] Ny;
	delete[] Nz;    
}

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
            					    double* dvdz, int s_dvdz )        					        
{
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
	
	double** Nx = array_2d(2,deg_xi+1); 
	double** Ny = array_2d(2,deg_eta+1); 
	double** Nz = array_2d(2,deg_zeta+1); 
	
    double s,sxx,syy,szz;
    double gl; 
    
    int I,J,K; 
    
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
					I = spany-deg_eta+j ;
					J = sj -1 - (spanz-deg_zeta+k) ; 
					K = spanx-deg_xi+i ;
					gl = image[K+J*sk+I*sk*sj]; 
					s   = s   + Nz[0][k]*Ny[0][j]*Nx[0][i]*gl;  
					sxx = sxx + Nz[0][k]*Ny[0][j]*Nx[1][i]*gl;  
					syy = syy + Nz[0][k]*Ny[1][j]*Nx[0][i]*gl;  
					szz = szz + Nz[1][k]*Ny[0][j]*Nx[0][i]*gl;  
					  
				}
			}
		}		
		v[p] = s;
		dvdx[p] = sxx; 
		dvdy[p] = syy; 
		dvdz[p] = szz;  
		
	}    
	delete_array_2d(Nx);
	delete_array_2d(Ny);
	delete_array_2d(Nz);  
}



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
            					    double* dvdz, int s_dvdz )
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

void EvaluateCardBspline2(double* image, int si, int sj, int sk,
                          double xmin, double ymin, double zmin,
                          double dx, double dy, double dz,  
                          double* xc, int sxc, 
                          double* yc, int syc, 
                          double* zc, int szc, 
					      double* x, int sx,
			  		      double* y, int sy,
			  		      double* z, int sz,   
					      double* v, int s_v )
{
    int tx,ty,tz; 
    int I,J,K; 
    double s; 
	for ( int p = 0; p < sx; p++ )
	{
		tx = std::floor( (x[p]-(xmin+dx/2))/dx ) ;
		ty = std::floor( (y[p]-(ymin+dy/2))/dy ) ;
		tz = std::floor( (z[p]-(zmin+dz/2))/dz ) ;
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
                	s = s + image[K+J*sk+I*sk*sj]*CardinalBspline2((x[p]-xc[i])/dx)*CardinalBspline2((y[p]-yc[j])/dy)*CardinalBspline2((z[p]-zc[k])/dz); 
            	}
        	}
		}
		v[p] = s; 
	}
}
void EvaluateCardBspline3(double* image, int si, int sj, int sk,
                          double xmin, double ymin, double zmin,
                          double dx, double dy, double dz,  
                          double* xc, int sxc, 
                          double* yc, int syc, 
                          double* zc, int szc, 
					      double* x, int sx,
			  		      double* y, int sy,
			  		      double* z, int sz,   
					      double* v, int s_v )
{
    int tx,ty,tz; 
    int I,J,K; 
    double s; 
	for ( int p = 0; p < sx; p++ )
	{
		tx = std::floor( (x[p]-(xmin+dx/2))/dx ) ;
		ty = std::floor( (y[p]-(ymin+dy/2))/dy ) ;
		tz = std::floor( (z[p]-(zmin+dz/2))/dz ) ;
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
                	s = s + image[K+J*sk+I*sk*sj]*CardinalBspline3((x[p]-xc[i])/dx)*CardinalBspline3((y[p]-yc[j])/dy)*CardinalBspline3((z[p]-zc[k])/dz); 
            	}
        	}
		}
		v[p] = s; 
	}
}
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
					                 double* dvdz, int s_dvdz )
{
    int tx,ty,tz; 
    int I,J,K; 
    double s,sxx,syy,szz; 
    double gl; 
	for ( int p = 0; p < sx; p++ )
	{
		tx = std::floor( (x[p]-(xmin+dx/2))/dx ) ;
		ty = std::floor( (y[p]-(ymin+dy/2))/dy ) ;
		tz = std::floor( (z[p]-(zmin+dz/2))/dz ) ;
		s = 0; 
		sxx = 0;
		syy = 0; 
		szz = 0; 
		for (int k=tz-1; k<=tz+2; k++)
		{
    		J = sj-1-k; 
    		for (int j=ty-1; j<=ty+2; j++)
        	{
            	I = j; 
            	for (int i=tx-1; i<=tx+2; i++)
            	{
                	K = i; 
                	gl = image[K+J*sk+I*sk*sj]; 
                	s   =  s  + gl*CardinalBspline2((x[p]-xc[i])/dx)*CardinalBspline2((y[p]-yc[j])/dy)*CardinalBspline2((z[p]-zc[k])/dz);
                	sxx = sxx + gl*CardinalBsplineDer2((x[p]-xc[i])/dx)*CardinalBspline2((y[p]-yc[j])/dy)*CardinalBspline2((z[p]-zc[k])/dz); 
                	syy = syy + gl*CardinalBspline2((x[p]-xc[i])/dx)*CardinalBsplineDer2((y[p]-yc[j])/dy)*CardinalBspline2((z[p]-zc[k])/dz);
                	szz = szz + gl*CardinalBspline2((x[p]-xc[i])/dx)*CardinalBspline2((y[p]-yc[j])/dy)*CardinalBsplineDer2((z[p]-zc[k])/dz);
            	}
        	}
		}
		v[p] = s; 
		dvdx[p] = sxx/dx; 
		dvdy[p] = syy/dy; 
		dvdz[p] = szz/dz; 
	}
}
            					                 
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
					               double* dvdz, int s_dvdz )	
{
    int tx,ty,tz; 
    int I,J,K; 
    double s,sxx,syy,szz; 
    double gl; 
	for ( int p = 0; p < sx; p++ )
	{
		tx = std::floor( (x[p]-(xmin+dx/2))/dx ) ;
		ty = std::floor( (y[p]-(ymin+dy/2))/dy ) ;
		tz = std::floor( (z[p]-(zmin+dz/2))/dz ) ;
		s = 0; 
		sxx = 0;
		syy = 0; 
		szz = 0; 
		for (int k=tz-1; k<=tz+2; k++)
		{
    		J = sj-1-k; 
    		for (int j=ty-1; j<=ty+2; j++)
        	{
            	I = j; 
            	for (int i=tx-1; i<=tx+2; i++)
            	{
                	K = i; 
                	gl = image[K+J*sk+I*sk*sj]; 
                	s   =  s  + gl*CardinalBspline3((x[p]-xc[i])/dx)*CardinalBspline3((y[p]-yc[j])/dy)*CardinalBspline3((z[p]-zc[k])/dz);
                	sxx = sxx + gl*CardinalBsplineDer3((x[p]-xc[i])/dx)*CardinalBspline3((y[p]-yc[j])/dy)*CardinalBspline3((z[p]-zc[k])/dz); 
                	syy = syy + gl*CardinalBspline3((x[p]-xc[i])/dx)*CardinalBsplineDer3((y[p]-yc[j])/dy)*CardinalBspline3((z[p]-zc[k])/dz);
                	szz = szz + gl*CardinalBspline3((x[p]-xc[i])/dx)*CardinalBspline3((y[p]-yc[j])/dy)*CardinalBsplineDer3((z[p]-zc[k])/dz);
            	}
        	}
		}
		v[p] = s; 
		dvdx[p] = sxx/dx; 
		dvdy[p] = syy/dy; 
		dvdz[p] = szz/dz; 
	}
}            					               
            					               

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
            					               double* dvdz, int s_dvdz )
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
            					               double* dvdz, int s_dvdz )
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
                                              double* dfipdz, int sdfipdz) 
{
    // Returns the mean and standard deviation of the gray-levels on each element 
    // Returns also the image evaluation at the integration points and the gradient vector 
    // Order (integration points of each element)

    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipe = nipex*nipey*nipez ; 
 
    
    // Coordinates of 1d integration points on one element
    double* xge = new double[nipex]; 
    double* yge = new double[nipey];
    double* zge = new double[nipez];
    
    double* fe    = new double[nipe] ; 
    double* dfdxe = new double[nipe] ; 
    double* dfdye = new double[nipe] ; 
    double* dfdze = new double[nipe] ;  
    
    double fm ; // Mean value of gray-level 
    double fms; // Mean value of the square of gray-level  

     //Getting the coordinates of the integration points 
     /*
    for (int i=0; i<nip_xi ; i++)
    {
        xg[i] = knotXi[0]+mes_xi/2+i*mes_xi; 
    }
    for (int i=0; i<nip_eta ; i++)
    {
        yg[i] = knotEta[0]+mes_eta/2+i*mes_eta; 
    }    
    for (int i=0; i<nip_zeta; i++)
    {
        zg[i] = knotZeta[0]+mes_zeta/2+i*mes_zeta;        
    }*/
    
    int ecounter=0 ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int kk=k*nipez, ig1dz=0; kk<(k+1)*nipez && ig1dz< nipez; kk++,ig1dz++)
        {
            zge[ig1dz] = zg[kk] ; 
        }
        for (int j=0; j<ne_eta; j++)
        {
            for (int jj=j*nipey, ig1dy=0; jj<(j+1)*nipey && ig1dy< nipey; jj++,ig1dy++)
            {
                yge[ig1dy] = yg[jj] ;
            }               
            for (int i=0; i<ne_xi; i++)
            {
                for (int ii=i*nipex, ig1dx=0; ii<(i+1)*nipex && ig1dx< nipex; ii++,ig1dx++)
                {
                    xge[ig1dx] = xg[ii] ;
                }
                
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                
                EvaluateTrilinearInterpolationAndGradientStructured(f,sif,sjf,skf,
                knotXiImage, sknotXiImage, knotEtaImage, sknotEtaImage, knotZetaImage, sknotZetaImage, 
                xge, nipex, yge, nipey, zge, nipez, fe, dfdxe, dfdye, dfdze); 
                
                fm  = 0. ; 
                fms = 0. ; 
                for (int ip=0; ip<nipe; ip++)
                {
                    fm  = fm + fe[ip] ; 
                    fms = fms + std::pow(fe[ip],2) ;  
                    
                    fip[ip + ecounter*nipe]    = fe[ip]    ;
                    dfipdx[ip + ecounter*nipe] = dfdxe[ip] ;  
                    dfipdy[ip + ecounter*nipe] = dfdye[ip] ; 
                    dfipdz[ip + ecounter*nipe] = dfdze[ip] ; 
                } 
                dyne[ecounter]  = *std::max_element(fe,fe+nipe) - *std::min_element(fe,fe+nipe) ; 
                fmean[ecounter] = fm/nipe;  
                fstd[ecounter]  = std::sqrt( fms/nipe - std::pow(fmean[ecounter],2) ); 
                ecounter++; 			    
            }
        }
        
    }
   
    delete[] xge ; 
    delete[] yge ; 
    delete[] zge ; 
    
    delete[] fe    ; 
    delete[] dfdxe ;  
    delete[] dfdye ;  
    delete[] dfdze ;   
        
        
} 
 

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
                                        double* dfipdz, int sdfipdz)
{
    // Returns the mean and standard deviation of the gray-levels on each element 
    // Returns also the image evaluation at the integration points and the gradient vector 
    // Order (integration points of each element)

    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipe = nipex*nipey*nipez ; 
 
    
    // Coordinates of 1d integration points on one element
    double* xge = new double[nipex]; 
    double* yge = new double[nipey];
    double* zge = new double[nipez];
    
    double* fe    = new double[nipe] ; 
    double* dfdxe = new double[nipe] ; 
    double* dfdye = new double[nipe] ; 
    double* dfdze = new double[nipe] ;  
    
    double fm ; // Mean value of gray-level 
    double fms; // Mean value of the square of gray-level  

    int ecounter=0 ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int kk=k*nipez, ig1dz=0; kk<(k+1)*nipez && ig1dz< nipez; kk++,ig1dz++)
        {
            zge[ig1dz] = zg[kk] ; 
        }
        for (int j=0; j<ne_eta; j++)
        {
            for (int jj=j*nipey, ig1dy=0; jj<(j+1)*nipey && ig1dy< nipey; jj++,ig1dy++)
            {
                yge[ig1dy] = yg[jj] ;
            }               
            for (int i=0; i<ne_xi; i++)
            {
                for (int ii=i*nipex, ig1dx=0; ii<(i+1)*nipex && ig1dx< nipex; ii++,ig1dx++)
                {
                    xge[ig1dx] = xg[ii] ;
                }
                
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                
                EvaluateCBspline2AndGradientStructured(f,sif,sjf,skf,
                                                       xmin, ymin, zmin, dx, dy, dz, 
                                                       xc, sxc, yc, syc, zc, szc, 
                                                       xge, nipex, yge, nipey, zge, nipez, fe, dfdxe, dfdye, dfdze); 
                fm  = 0. ; 
                fms = 0. ; 
                for (int ip=0; ip<nipe; ip++)
                {
                    fm  = fm + fe[ip] ; 
                    fms = fms + std::pow(fe[ip],2) ;  
                    
                    fip[ip + ecounter*nipe]    = fe[ip]    ;
                    dfipdx[ip + ecounter*nipe] = dfdxe[ip] ;  
                    dfipdy[ip + ecounter*nipe] = dfdye[ip] ; 
                    dfipdz[ip + ecounter*nipe] = dfdze[ip] ; 
                } 
                dyne[ecounter]  = *std::max_element(fe,fe+nipe) - *std::min_element(fe,fe+nipe) ; 
                fmean[ecounter] = fm/nipe;  
                fstd[ecounter]  = std::sqrt( fms/nipe - std::pow(fmean[ecounter],2) ); 
                ecounter++; 			    
            }
        }
        
    }
   
    delete[] xge ; 
    delete[] yge ; 
    delete[] zge ; 
    
    delete[] fe    ; 
    delete[] dfdxe ;  
    delete[] dfdye ;  
    delete[] dfdze ;  

}


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
                                        double* dfipdz, int sdfipdz)
{
    // Returns the mean and standard deviation of the gray-levels on each element 
    // Returns also the image evaluation at the integration points and the gradient vector 
    // Order (integration points of each element)

    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipe = nipex*nipey*nipez ; 
 
    
    // Coordinates of 1d integration points on one element
    double* xge = new double[nipex]; 
    double* yge = new double[nipey];
    double* zge = new double[nipez];
    
    double* fe    = new double[nipe] ; 
    double* dfdxe = new double[nipe] ; 
    double* dfdye = new double[nipe] ; 
    double* dfdze = new double[nipe] ;  
    
    double fm ; // Mean value of gray-level 
    double fms; // Mean value of the square of gray-level  

    int ecounter=0 ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int kk=k*nipez, ig1dz=0; kk<(k+1)*nipez && ig1dz< nipez; kk++,ig1dz++)
        {
            zge[ig1dz] = zg[kk] ; 
        }
        for (int j=0; j<ne_eta; j++)
        {
            for (int jj=j*nipey, ig1dy=0; jj<(j+1)*nipey && ig1dy< nipey; jj++,ig1dy++)
            {
                yge[ig1dy] = yg[jj] ;
            }               
            for (int i=0; i<ne_xi; i++)
            {
                for (int ii=i*nipex, ig1dx=0; ii<(i+1)*nipex && ig1dx< nipex; ii++,ig1dx++)
                {
                    xge[ig1dx] = xg[ii] ;
                }
                
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                
                EvaluateCBspline3AndGradientStructured(f,sif,sjf,skf,
                                                       xmin, ymin, zmin, dx, dy, dz, 
                                                       xc, sxc, yc, syc, zc, szc, 
                                                       xge, nipex, yge, nipey, zge, nipez, fe, dfdxe, dfdye, dfdze); 
                fm  = 0. ; 
                fms = 0. ; 
                for (int ip=0; ip<nipe; ip++)
                {
                    fm  = fm + fe[ip] ; 
                    fms = fms + std::pow(fe[ip],2) ;  
                    
                    fip[ip + ecounter*nipe]    = fe[ip]    ;
                    dfipdx[ip + ecounter*nipe] = dfdxe[ip] ;  
                    dfipdy[ip + ecounter*nipe] = dfdye[ip] ; 
                    dfipdz[ip + ecounter*nipe] = dfdze[ip] ; 
                } 
                dyne[ecounter]  = *std::max_element(fe,fe+nipe) - *std::min_element(fe,fe+nipe) ; 
                fmean[ecounter] = fm/nipe;  
                fstd[ecounter]  = std::sqrt( fms/nipe - std::pow(fmean[ecounter],2) ); 
                ecounter++; 			    
            }
        }
        
    }
   
    delete[] xge ; 
    delete[] yge ; 
    delete[] zge ; 
    
    delete[] fe    ; 
    delete[] dfdxe ;  
    delete[] dfdye ;  
    delete[] dfdze ;  

}

                                        

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
                                              double* dfipdz, int sdfipdz) 
{
    // Returns the mean and standard deviation of the gray-levels on each element 
    // Returns also the image evaluation at the integration points and the gradient vector 
    // Order (integration points of each element)

    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
    
 
    int nipe = nipex*nipey*nipez ; 
 
    
    // Coordinates of 1d integration points on one element
    double* xge = new double[nipex]; 
    double* yge = new double[nipey];
    double* zge = new double[nipez];
    
    double* fe    = new double[nipe] ; 
    double* dfdxe = new double[nipe] ; 
    double* dfdye = new double[nipe] ; 
    double* dfdze = new double[nipe] ;  
    
    double fm ; // Mean value of gray-level 
    double fms; // Mean value of the square of gray-level  
    int nv ; // Number of keeped voxels per element  
 
    
    int ecounter=0 ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int kk=k*nipez, ig1dz=0; kk<(k+1)*nipez && ig1dz< nipez; kk++,ig1dz++)
        {
            zge[ig1dz] = zg[kk] ; 
        }
        for (int j=0; j<ne_eta; j++)
        {
            for (int jj=j*nipey, ig1dy=0; jj<(j+1)*nipey && ig1dy< nipey; jj++,ig1dy++)
            {
                yge[ig1dy] = yg[jj] ;
            }               
            for (int i=0; i<ne_xi; i++)
            {
                for (int ii=i*nipex, ig1dx=0; ii<(i+1)*nipex && ig1dx< nipex; ii++,ig1dx++)
                {
                    xge[ig1dx] = xg[ii] ;
                }
                
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                
                EvaluateTrilinearInterpolationAndGradientStructured(f,sif,sjf,skf,
                knotXiImage, sknotXiImage, knotEtaImage, sknotEtaImage, knotZetaImage, sknotZetaImage, 
                xge, nipex, yge, nipey, zge, nipez, fe, dfdxe, dfdye, dfdze); 
                
                fm  = 0. ; 
                fms = 0. ; 
                nv = 0 ; 
                for (int ip=0; ip<nipe; ip++)
                {
                    if (fe[ip]>thrsh)
                    {
                        fm  = fm + fe[ip] ; 
                        fms = fms + std::pow(fe[ip],2) ; 
                        nv++ ;                         
                    } 
                    
                    fip[ip + ecounter*nipe]    = fe[ip]    ;
                    dfipdx[ip + ecounter*nipe] = dfdxe[ip] ;  
                    dfipdy[ip + ecounter*nipe] = dfdye[ip] ; 
                    dfipdz[ip + ecounter*nipe] = dfdze[ip] ; 
                } 
                if (nv !=0)
                {
                    fmean[ecounter] = fm/nv;  
                    fstd[ecounter]  = std::sqrt( fms/nv - std::pow(fmean[ecounter],2) );                     
                }
                else 
                {
                    fmean[ecounter] =  0 ;  
                    fstd[ecounter]  =  0 ;     
                }

                ecounter++; 			    
            }
        }
        
    }
   
    delete[] xge ; 
    delete[] yge ; 
    delete[] zge ; 
    
    delete[] fe    ; 
    delete[] dfdxe ;  
    delete[] dfdye ;  
    delete[] dfdze ;   
        
        
}      

void GetC8MeshFromVoxelsTrilinearInterp(double* image, int si, int sj, int sk, 
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double thrsh, 
                                        double* knotXi, int sknotXi,
                                        double* knotEta, int sknotEta,  
                                        double* knotZeta, int sknotZeta,  
                                        int* elementMask, int se )
{
    int deg_xi   = 1 ; 
    int deg_eta  = 1 ;
    int deg_zeta = 1 ;

    int nbf_xi_Image   = sknotXiImage   -1 - deg_xi   ; 
	int nbf_eta_Image  = sknotEtaImage  -1 - deg_eta  ;  
	int nbf_zeta_Image = sknotZetaImage -1 - deg_zeta ;
	
    int nbf_xi    = sknotXi    -1 - deg_xi   ; 
	int nbf_eta   = sknotEta   -1 - deg_eta  ;  
	int nbf_zeta  = sknotZeta  -1 - deg_zeta ;    
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    
    int nxi   = nbf_xi_Image   - 1 ;  
	int neta  = nbf_eta_Image  - 1 ;
	int nzeta = nbf_zeta_Image - 1 ; 
 
 
	
	double mes_xi   = knotXi[deg_xi+1] - knotXi[deg_xi]         ; 
    double mes_eta  = knotEta[deg_eta+1] - knotEta[deg_eta]     ; 
    double mes_zeta = knotZeta[deg_zeta+1] - knotZeta[deg_zeta] ; 
    
    double  xg   ; 
    double  yg   ; 
    double  zg   ; 
    
    int* SpanxImage = new int [ne_xi]   ;
	int* SpanyImage = new int [ne_eta]  ;
	int* SpanzImage = new int [ne_zeta] ;  
 
	
	double** Nx  = array_2d(ne_xi, deg_xi+1)   ; 
	double** Ny  = array_2d(ne_eta,deg_eta+1)  ; 
	double** Nz  = array_2d(ne_zeta,deg_zeta+1) ; 
    

    //Getting the coordinates of the univariate voxels centers  
    for (int i=0; i<ne_xi ; i++)
    {
        xg  = knotXi[0]+mes_xi/2+i*mes_xi; 
        SpanxImage[i] = findspan(nxi, deg_xi, xg , knotXiImage);
        basisfuns(SpanxImage[i], xg , deg_xi, knotXiImage, Nx[i]);
    }
    for (int i=0; i<ne_eta ; i++)
    {
        yg  = knotEta[0]+mes_eta/2+i*mes_eta;
        SpanyImage[i] = findspan(neta, deg_eta, yg , knotEtaImage);
        basisfuns(SpanyImage[i], yg , deg_eta, knotEtaImage, Ny[i]);
    }    
    for (int i=0; i<ne_zeta; i++)
    {
        zg  = knotZeta[0]+mes_zeta/2+i*mes_zeta;   
        SpanzImage[i] = findspan(nzeta, deg_zeta, zg , knotZetaImage); 
        basisfuns(SpanzImage[i], zg , deg_zeta, knotZetaImage, Nz[i]);    
    } 
    int I,J,K; 
    double s; 
    int ecounter=0 ; 
 
    for ( int r= 0; r < ne_zeta; r++)
	{
		for ( int q =0; q < ne_eta; q++)
		{
			for ( int p =0; p < ne_xi; p++)
			{
    			s=0; 
			    for (int k=0; k<deg_zeta+1 ;k++)
                {  
                    for (int j=0; j<deg_eta+1; j++)
                    {   
                        for (int i=0; i<deg_xi+1; i++)
                        {
                            // Evaluating the gray-level of the voxel 
                            I = SpanyImage[q]-deg_eta+j ;
                    		J = sj -1 - (SpanzImage[r]-deg_zeta+k) ; 
                    		K = SpanxImage[p]-deg_xi+i ;
                    		s   = s   + Nz[r][k]*Ny[q][j]*Nx[p][i]*image[K+J*sk+I*sk*sj]; 
                        }
                    }
                }
                if (s>thrsh)
                {
                    elementMask[ecounter] = 1 ;    
                }
                else
                {
                    elementMask[ecounter] = 0 ; 
                }
                ecounter++;    	
			}
        }
    }
 
  
    delete[] SpanxImage ;
    delete[] SpanyImage ; 
    delete[] SpanzImage ; 
    
    delete_array_2d(Nx);
	delete_array_2d(Ny);
	delete_array_2d(Nz);  
 
}   


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
                                              double* dfipdz, int sdfipdz)
{
    int nip = siN ;  // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                 double* fe    = new double[nip] ; 
                 double* dfdxe = new double[nip] ; 
                 double* dfdye = new double[nip] ; 
                 double* dfdze = new double[nip] ; 
                 
                 // Mapping the integration points from the reference element 
                 double* x = new double[nip] ;
                 double* y = new double[nip] ; 
                 double* z = new double[nip] ; 
                 
                 for (int ip=0; ip<nip ; ip++)
                 {
                     x[ip] = 0 ; 
                     y[ip] = 0 ; 
                     z[ip] = 0 ; 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         x[ip] += N[ibf + ip*sjN] * n[0+e[ibf+i*sje]*sjn] ; 
                         y[ip] += N[ibf + ip*sjN] * n[1+e[ibf+i*sje]*sjn] ; 
                         z[ip] += N[ibf + ip*sjN] * n[2+e[ibf+i*sje]*sjn] ; 
                     } 
                 }
                 
                // Evaluate image and gradient on each voxel of the element 
                 EvaluateTrilinearInterpolationAndGradient(f, sif, sjf, skf, 
                 knotXiImage, sknotXiImage, knotEtaImage, sknotEtaImage, knotZetaImage, sknotZetaImage, 
                 x, nip, y, nip, z, nip, fe, nip, dfdxe, nip, dfdye, nip,  dfdze, nip ); 
                     
                 // Loop on the integration points 
                 double fm  = 0. ; 
                 double fms = 0. ;                  
                 for (int ip=0; ip<nip ; ip++)
                 {
                     fm  = fm + fe[ip] ; 
                     fms = fms + std::pow(fe[ip],2) ;  
                     fip[ip + i*nip]    = fe[ip]    ;
                     dfipdx[ip + i*nip] = dfdxe[ip] ;  
                     dfipdy[ip + i*nip] = dfdye[ip] ; 
                     dfipdz[ip + i*nip] = dfdze[ip] ;    
                 } 
                 dyne[i]  = *std::max_element(fe,fe+nip) - *std::min_element(fe,fe+nip) ; 
                 fmean[i] = fm/nip;  
                 fstd[i]  = std::sqrt( fms/nip - std::pow(fmean[i],2) );    
                 
                 delete[] fe    ; 
                 delete[] dfdxe ;  
                 delete[] dfdye ;  
                 delete[] dfdze ; 
                 
                 delete[] x ; 
                 delete[] y ; 
                 delete[] z ;                      				       
            }
    }  
}

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
                                              double* dfipdz, int sdfipdz)
{
    int nip = siN ;  // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                 double* fe    = new double[nip] ; 
                 double* dfdxe = new double[nip] ; 
                 double* dfdye = new double[nip] ; 
                 double* dfdze = new double[nip] ; 
                 
                 // Mapping the integration points from the reference element 
                 double* x = new double[nip] ;
                 double* y = new double[nip] ; 
                 double* z = new double[nip] ; 
                 
                 for (int ip=0; ip<nip ; ip++)
                 {
                     x[ip] = 0 ; 
                     y[ip] = 0 ; 
                     z[ip] = 0 ; 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         x[ip] += N[ibf + ip*sjN] * n[0+e[ibf+i*sje]*sjn] ; 
                         y[ip] += N[ibf + ip*sjN] * n[1+e[ibf+i*sje]*sjn] ; 
                         z[ip] += N[ibf + ip*sjN] * n[2+e[ibf+i*sje]*sjn] ; 
                     } 
                 }
                 
                // Evaluate image and gradient on each voxel of the element 
                EvaluateCardBsplineAndGradient3(f,sif,sjf,skf,
                                       xmin, ymin, zmin, dx, dy, dz, 
                                       xc, sxc, yc, syc, zc, szc, 
                                      x, nip, y, nip, z, nip, fe, nip, dfdxe, nip, dfdye, nip,  dfdze, nip ); 
 
                 // Loop on the integration points 
                 double fm  = 0. ; 
                 double fms = 0. ;                  
                 for (int ip=0; ip<nip ; ip++)
                 {
                     fm  = fm + fe[ip] ; 
                     fms = fms + std::pow(fe[ip],2) ;  
                     fip[ip + i*nip]    = fe[ip]    ;
                     dfipdx[ip + i*nip] = dfdxe[ip] ;  
                     dfipdy[ip + i*nip] = dfdye[ip] ; 
                     dfipdz[ip + i*nip] = dfdze[ip] ;    
                 } 
                 dyne[i]  = *std::max_element(fe,fe+nip) - *std::min_element(fe,fe+nip) ; 
                 fmean[i] = fm/nip;  
                 fstd[i]  = std::sqrt( fms/nip - std::pow(fmean[i],2) );    
                 
                 delete[] fe    ; 
                 delete[] dfdxe ;  
                 delete[] dfdye ;  
                 delete[] dfdze ; 
                 
                 delete[] x ; 
                 delete[] y ; 
                 delete[] z ;                      				       
            }
    }  
}
 
 
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
                                                double* dfipdz, int sdfipdz)  
{
    int nip = siN ;  // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                  double* fe    = new double[nip] ; 
                  double* dfdxe = new double[nip] ; 
                  double* dfdye = new double[nip] ; 
                  double* dfdze = new double[nip] ; 
                  
                  // Mapping the integration points from the reference element 
                  double* x = new double[nip] ;
                  double* y = new double[nip] ; 
                  double* z = new double[nip] ; 
                  
                  for (int ip=0; ip<nip ; ip++)
                  {
                      x[ip] = 0 ; 
                      y[ip] = 0 ; 
                      z[ip] = 0 ; 
                      for (int ibf=0; ibf < nbf_elem; ibf++)
                      {
                          x[ip] += N[ibf + ip*sjN] * n[0+e[ibf+i*sje]*sjn] ; 
                          y[ip] += N[ibf + ip*sjN] * n[1+e[ibf+i*sje]*sjn] ; 
                          z[ip] += N[ibf + ip*sjN] * n[2+e[ibf+i*sje]*sjn] ; 
                      } 
                  }
                  
                  // Evaluate image and gradient on each voxel of the element                  
                  EvaluateBsplineAndGradient3D(x, nip, y, nip, z, nip,
                                              knotXiImage, sknotXiImage, 
                                              knotEtaImage, sknotEtaImage, 
                                              knotZetaImage, sknotZetaImage, 
                                              deg_xi, deg_eta, deg_zeta, 
                                              cp, scp, 
                                              fe, nip, dfdxe, nip, dfdye, nip,  dfdze, nip ); 
  
                  // Loop on the integration points 
                  double fm  = 0. ; 
                  double fms = 0. ;                  
                  for (int ip=0; ip<nip ; ip++)
                  {
                      fm  = fm + fe[ip] ; 
                      fms = fms + std::pow(fe[ip],2) ;  
                      fip[ip + i*nip]    = fe[ip]    ;
                      dfipdx[ip + i*nip] = dfdxe[ip] ;  
                      dfipdy[ip + i*nip] = dfdye[ip] ; 
                      dfipdz[ip + i*nip] = dfdze[ip] ;    
                  } 
                  dyne[i]  = *std::max_element(fe,fe+nip) - *std::min_element(fe,fe+nip) ; 
                  fmean[i] = fm/nip;  
                  fstd[i]  = std::sqrt( fms/nip - std::pow(fmean[i],2) );    
                  
                  delete[] fe    ; 
                  delete[] dfdxe ;  
                  delete[] dfdye ;  
                  delete[] dfdze ; 
                  
                  delete[] x ; 
                  delete[] y ; 
                  delete[] z ;                      				       
            }
    }  
} 
                                               

                                        
                                               					               
            					                 




