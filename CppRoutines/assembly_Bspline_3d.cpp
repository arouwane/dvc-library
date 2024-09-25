#include "assembly_Bspline_3d.h"

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
                                   double* deta_deta_xi_xi, int s_deta_deta_xi_xi ) 
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_xi_elem   = deg_xi+1 ; 
	int nbf_eta_elem  = deg_eta+1 ;
	int nbf_zeta_elem = deg_zeta+1 ;  
 
	int nelems_xi   = nbf_xi   - deg_xi   ;
	int nelems_eta  = nbf_eta  - deg_eta  ; 
	int nelems_zeta = nbf_zeta - deg_zeta ;  

	double* e_xi   = new double[nelems_xi+1]   ;
	double* e_eta  = new double[nelems_eta+1]  ;
	double* e_zeta = new double[nelems_zeta+1] ; 

	for (int i=0;i< nelems_xi+1; i++)
	{
		e_xi[i] = knotXi[i+deg_xi];
	}
    for(int i=0; i< nelems_eta+1; i++)
    {
		e_eta[i] = knotEta[i+deg_eta];
	}
    for(int i=0; i< nelems_zeta+1; i++)
    {
		e_zeta[i] = knotZeta[i+deg_zeta];
	} 
	
	double* xig     = new double[nipx];
    double* etag    = new double[nipy];
    double* zetag   = new double[nipz]; 
    double* wg_xi   = new double[nipx]; 
    double* wg_eta  = new double[nipy]; 
    double* wg_zeta = new double[nipz]; 
    
    gauleg(-1,1,xig,wg_xi,nipx);    
	gauleg(-1,1,etag,wg_eta,nipy); 
    gauleg(-1,1,zetag,wg_zeta,nipz); 
    
    double xi_min   ; 
    double eta_min  ; 
    double zeta_min ; 
 
    double mes_xi   =  (e_xi[1]   - e_xi[0] )   / nex ; 
    double mes_eta  =  (e_eta[1]  - e_eta[0])   / ney ; 
    double mes_zeta =  (e_zeta[1] - e_zeta[0])  / nez ; 
    
    double xi_p ; 
    double eta_p ; 
    double zeta_p ; 
    
    double** Nxi         = array_2d(2,nbf_xi_elem);
	double** Neta        = array_2d(2,nbf_eta_elem);
	double** Nzeta       = array_2d(2,nbf_zeta_elem);
    
    // Initializing integrated operators with zero values 
    for (int i=0; i < s_dxi_dxi; i++)
    {
        dxi_dxi[i] = 0 ; 
        dxi_xi[i] = 0 ; 
        xi_dxi[i] = 0 ; 
        xi_xi[i] = 0 ; 
    }
    for (int i=0; i<s_deta_deta; i++)
    {
        deta_deta[i] = 0 ; 
        deta_eta[i] = 0 ; 
        eta_deta[i] = 0 ; 
        eta_eta[i] = 0 ; 
    }
    for (int i=0; i<s_dzeta_dzeta; i++)
    {
        dzeta_dzeta[i] = 0 ; 
        dzeta_zeta[i] = 0 ;        
        zeta_dzeta[i] = 0 ;     
        zeta_zeta[i] = 0 ;  
    }
 
    // Integration in xi direction 
    for (int i=0; i< nelems_xi; i++)
    {
        // Loop over integration elements of the current 1d element 
        for (int j=0; j<nex; j++)
        {
            xi_min = e_xi[i] + j*mes_xi ; 
            for (int ip=0; ip < nipx; ip++)
            {
                 xi_p = xi_min  + 0.5*(xig[ip]+1)*mes_xi ;
                 dersbasisfuns(deg_xi,  knotXi  ,  xi_p, i+deg_xi, 1, Nxi);
                 // Perform the tensor product                  
                 for (int ibf=0; ibf<nbf_xi_elem ; ibf++)
                 {
                     for (int jbf=0; jbf<nbf_xi_elem; jbf++)
                     {
                         dxi_dxi[jbf + ibf*nbf_xi_elem + (j+i*nex)*nbf_xi_elem*nbf_xi_elem]  += wg_xi[ip]*0.5*mes_xi* Nxi[1][ibf]*Nxi[1][jbf] ; 
                         dxi_xi[jbf + ibf*nbf_xi_elem + (j+i*nex)*nbf_xi_elem*nbf_xi_elem]   += wg_xi[ip]*0.5*mes_xi*Nxi[1][ibf]*Nxi[0][jbf] ; 
                         xi_dxi[jbf + ibf*nbf_xi_elem + (j+i*nex)*nbf_xi_elem*nbf_xi_elem]   += wg_xi[ip]*0.5*mes_xi* Nxi[0][ibf]*Nxi[1][jbf] ;
                         xi_xi[ jbf + ibf*nbf_xi_elem + (j+i*nex)*nbf_xi_elem*nbf_xi_elem]   += wg_xi[ip]*0.5*mes_xi* Nxi[0][ibf]*Nxi[0][jbf] ;                         
                     }
                 } 
            }
        }
    }
    // Integration in eta direction 
    for (int i=0; i< nelems_eta; i++)
    {
        // Loop over integration elements of the current 1d element 
        for (int j=0; j<ney; j++)
        {
            eta_min = e_eta[i] + j*mes_eta ; 
            for (int ip=0; ip < nipy; ip++)
            {
                 eta_p = eta_min  + 0.5*(etag[ip]+1)*mes_eta ;
                 dersbasisfuns(deg_eta,  knotEta  ,  eta_p, i+deg_eta, 1, Neta);
                 // Perform the tensor product 
                 for (int ibf=0; ibf<nbf_eta_elem ; ibf++)
                 {
                     for (int jbf=0; jbf<nbf_eta_elem; jbf++)
                     {
                         deta_deta[jbf + ibf*nbf_eta_elem + (j+i*ney)*nbf_eta_elem*nbf_eta_elem] += wg_eta[ip]*0.5*mes_eta* Neta[1][ibf]*Neta[1][jbf] ; 
                         deta_eta[jbf + ibf*nbf_eta_elem + (j+i*ney)*nbf_eta_elem*nbf_eta_elem]  += wg_eta[ip]*0.5*mes_eta* Neta[1][ibf]*Neta[0][jbf] ; 
                         eta_deta[jbf + ibf*nbf_eta_elem + (j+i*ney)*nbf_eta_elem*nbf_eta_elem]  += wg_eta[ip]*0.5*mes_eta* Neta[0][ibf]*Neta[1][jbf] ;
                         eta_eta[jbf + ibf*nbf_eta_elem + (j+i*ney)*nbf_eta_elem*nbf_eta_elem]   += wg_eta[ip]*0.5*mes_eta* Neta[0][ibf]*Neta[0][jbf] ;
                     }
                 }     
            }
        }
    }    

    // Integration in zeta direction 
    for (int i=0; i< nelems_zeta; i++)
    {
        // Loop over integration elements of the current 1d element 
        for (int j=0; j<nez; j++)
        {
            zeta_min = e_zeta[i] + j*mes_zeta ; 
            for (int ip=0; ip < nipz; ip++)
            {
                 zeta_p = zeta_min  + 0.5*(zetag[ip]+1)*mes_zeta ;
                 dersbasisfuns(deg_zeta,  knotZeta  ,  zeta_p, i+deg_zeta, 1, Nzeta);
                 // Perform the tensor product 
                 for (int ibf=0; ibf<nbf_zeta_elem ; ibf++)
                 {
                     for (int jbf=0; jbf<nbf_zeta_elem; jbf++)
                     {
                         dzeta_dzeta[jbf + ibf*nbf_zeta_elem + (j+i*nez)*nbf_zeta_elem*nbf_zeta_elem] += wg_zeta[ip]*0.5*mes_zeta* Nzeta[1][ibf]*Nzeta[1][jbf] ; 
                         dzeta_zeta[jbf + ibf*nbf_zeta_elem + (j+i*nez)*nbf_zeta_elem*nbf_zeta_elem]  += wg_zeta[ip]*0.5*mes_zeta* Nzeta[1][ibf]*Nzeta[0][jbf] ; 
                         zeta_dzeta[jbf + ibf*nbf_zeta_elem + (j+i*nez)*nbf_zeta_elem*nbf_zeta_elem]  += wg_zeta[ip]*0.5*mes_zeta* Nzeta[0][ibf]*Nzeta[1][jbf] ;
                         zeta_zeta[jbf + ibf*nbf_zeta_elem + (j+i*nez)*nbf_zeta_elem*nbf_zeta_elem]   += wg_zeta[ip]*0.5*mes_zeta* Nzeta[0][ibf]*Nzeta[0][jbf] ;
                     }
                 }     
            }
        }
    }
 
    int ei;
    int ej; 
    int iie ; 
    int nex_tot = nelems_xi*nex ; 
    int nbf_elem2d = nbf_xi_elem*nbf_eta_elem ; 
    
    
    // Integration in (xi,eta) plane by kronecker product 
    // Loop over 2d basis functions in the plane 
    for (int j=0; j< nelems_eta; j++ )
    {
        for (int i=0; i<nelems_xi; i++)
        {
            int I; 
            int Ixi; 
            int Ieta; 
            // Loop over integration elements 
            for (int ji=0; ji < ney ; ji++)
            {
                ej = ji + j*ney ;
                for (int ii=0; ii<nex; ii++)
                {
                    ei = ii + i*nex ; 
                    iie = ei + ej*nex_tot ;
                    // Kronecker product of the univariate integrated operators
                    // Kronkecker rootine for AxB 
                    for (int ieta =0; ieta< nbf_eta_elem  ; ieta++ )
                    {
                        for (int ixi=0 ; ixi<nbf_xi_elem ; ixi++)
                        {
                            for (int jeta=0; jeta < nbf_eta_elem; jeta++ )
                            {
                                for (int jxi=0; jxi < nbf_xi_elem ; jxi++ )
                                {      
                                    I    = jeta*nbf_xi_elem+jxi + (ieta*nbf_xi_elem+ixi)*nbf_elem2d + iie*nbf_elem2d*nbf_elem2d ;
                                    Ixi  = jxi  + ixi*nbf_xi_elem + ei*nbf_xi_elem*nbf_xi_elem ; 
                                    Ieta = jeta + ieta*nbf_eta_elem + ej*nbf_eta_elem*nbf_eta_elem ; 
     
                                    eta_eta_dxi_dxi[I] =  eta_eta[Ieta]*dxi_dxi[Ixi]  ; 
                                    eta_deta_dxi_xi[I] =  eta_deta[Ieta]*dxi_xi[Ixi]  ;                  
                                    deta_eta_xi_dxi[I] =  deta_eta[Ieta]*xi_dxi[Ixi] ;              
                                    deta_deta_xi_xi[I] =  deta_deta[Ieta]*xi_xi[Ixi] ;                                                                                                                                                                                              
                                }
                            }
                        }
                    }
                }
            }       
        }
    } 
    // Clearing memory 
    delete[] e_xi   ; 
    delete[] e_eta  ; 
    delete[] e_zeta ; 
    delete[] xig ; 
    delete[] etag ; 
    delete[] zetag ; 
    delete[] wg_xi; 
    delete[] wg_eta; 
    delete[] wg_zeta;   
    delete_array_2d(Nxi);
	delete_array_2d(Neta);
	delete_array_2d(Nzeta);     
}

void Laplacian(int deg_xi, int deg_eta, int deg_zeta,
               double* knotXi, int sknotXi, 
               double* knotEta, int sknotEta, 
               double* knotZeta, int sknotZeta, 
               int* NOELEM, int siNOELEM, int sjNOELEM, 
               
               double* zeta_zeta 
               double* dzeta_dzeta
               double* eta_eta_dxi_dxi 
               double* deta_deta_xi_xi
               double* eta_eta_xi_xi
 
               
               
               int* indexI, int nnzI, 
               int* indexJ, int nnzJ,
               double* nnz_values, int nnz)
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int nbf_xi_elem   = deg_xi+1 ; 
	int nbf_eta_elem  = deg_eta+1 ;
	int nbf_zeta_elem = deg_zeta+1 ;  
 
	
	int nbf_elem = nbf_xi_elem * nbf_eta_elem * nbf_zeta_elem ; 
	int nbf_elem2d =  nbf_xi_elem * nbf_eta_elem ; 
 
	int nelems_xi   = nbf_xi   - deg_xi   ;
	int nelems_eta  = nbf_eta  - deg_eta  ; 
	int nelems_zeta = nbf_zeta - deg_zeta ;  
	
	for (int k=0; k< n_elems_zeta ; k++)
	{
    	for (int j=0; j< n_elems_eta; j++)
    	{
        	for (int i=0; i< n_elems_xi ; i++)
        	{
                int ecounter = i +j*n_elems_xi + k*n_elems_xi*n_elems_eta ;
                double** Leb = array_2d(nbf_elem, nbf_elem) ; // Block of the elementary Laplacian Operator 
                //double** zeta_zeta_eta_eta_dxi_dxi = array_2d(nbf_elem,nbf_elem);
                //double** zeta_zeta_eta_eta_dxi_xi  = array_2d(nbf_elem,nbf_elem);
                //double** zeta_dzeta_eta_eta_dxi_xi = array_2d(nbf_elem,nbf_elem);
                //double** zeta_zeta_deta_eta_xi_dxi = array_2d(nbf_elem,nbf_elem);
                //double** zeta_zeta_deta_deta_xi_xi = array_2d(nbf_elem,nbf_elem);
                //double** zeta_dzeta_deta_eta_xi_xi = array_2d(nbf_elem,nbf_elem);
                //double** dzeta_zeta_eta_eta_xi_dxi = array_2d(nbf_elem,nbf_elem);
                //double** dzeta_zeta_eta_deta_xi_xi = array_2d(nbf_elem,nbf_elem);
                //double** dzeta_dzeta_eta_eta_xi_xi = array_2d(nbf_elem,nbf_elem);

                double** dNdx_dNdx = array_2d(nbf_elem,nbf_elem);
                //double** dNdx_dNdy  = array_2d(nbf_elem,nbf_elem);
                //double** dNdx_dNdz = array_2d(nbf_elem,nbf_elem);
                //double** dNdy_dNdx = array_2d(nbf_elem,nbf_elem);
                double** dNdy_dNdy = array_2d(nbf_elem,nbf_elem);
                //double** dNdy_dNdz = array_2d(nbf_elem,nbf_elem);
                //double** dNdz_dNdx = array_2d(nbf_elem,nbf_elem);
                //double** dNdz_dNdy = array_2d(nbf_elem,nbf_elem);
                double** dNdz_dNdz = array_2d(nbf_elem,nbf_elem);
                                
                
                int Is; 
                int Izeta; 
                
                int ei ;
                int ej ; 
                
                for (int ibf=0; ibf < nbf_elem ; ibf ++)
                {
                    for (int jbf=0; jbf < nbf_elem; jbf++)
                    {
                        Leb[ibf][jbf] = 0 ; 
                    }
                }
                // Loop over integration elements in z direction 
                for (int ki=0; ki<nei_zeta; ki++)
                {
                    // Loop over 2d integration elements (of the z slice)
                    for (int ji=0; ji<nei_eta; ji++)
                    {
                        ej = ji + j*nei_eta ; 
                        for (int ii=0; ii<nei_xi ; ii++)
                        {
                            ei = ii + i*nei_xi 
                            iie = ei + ej*nei_xi_tot
                            
                            // Perform kronecker product to get the trivariate elementary (integration operators)
                            for(int izeta=0; izeta<nbf_zeta_elem; izeta++)
                            {
                                for (int is=0; is< nbf_elem2d; is++ )
                                {
                                    for (int jzeta=0; jzeta<nbf_zeta_elem; jzeta++)
                                    {
                                        for (int js=0; js<nbf_elem2d; js++)
                                        {
                                            Is    =  js + is*nbf_elem2d + iie*nbf_elem ; 
                                            Izeta =  jzeta + izeta*nbf_zeta_elem ;  
                                            dNdx_dNdx[izeta*nbf_elem2d+is][jzeta*nbf_elem2d+js] =  zeta_zeta[Izeta]*eta_eta_dxi_dxi[Is] ; 
                                            dNdy_dNdy[izeta*nbf_elem2d+is][jzeta*nbf_elem2d+js] =  zeta_zeta[Izeta]*deta_deta_xi_xi[Is] ;   
                                            dNdz_dNdz[izeta*nbf_elem2d+is][jzeta*nbf_elem2d+js] =  dzeta_dzeta[Izeta]*eta_eta_xi_xi[Is] ;                                    
                                        }
                                    }
                                }
                            }  
                            // Filling the elementary differential matrix 
                            for (int ibf=0; ibf<nbf_elem;ibf++)
                            {
                                for (int jbf=0; jbf<bf_elem;jbf++)
                                {
                                    Leb[ibf][jbf] +== dNdx_dNdx[ibf][jbf] + dNdy_dNdy[ibf][jbf] + dNdz_dNdz[ibf][jbf] ;   
                                }
                            }
                            
                        }
                    }
                    
                }
                // Elementary contribution to the global matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM[ibf + ecounter* sjNOELEM] ; 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = NOELEM[jbf + ecounter*sjNOELEM] ;
    					
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = Leb[ibf][jbf];
    		            sp_count ++ ;  	
    		            
		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = Leb[ibf][jbf];
    		            sp_count ++ ;
    		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Leb[ibf][jbf];
    		            sp_count ++ ; 
                    }
                }                
            	delete_array_2d(Leb); 
            	delete_array_2d(dNdx_dNdx);
            	delete_array_2d(dNdy_dNdy);
            	delete_array_2d(dNdz_dNdz);
        	}
    	}
	}
	

 
                       
                       




} 
 
    
    


