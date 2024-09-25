#include "assemblyRoutines.h"

 
 
void GetBsplineFunctionsMatrixStructured(double* x, int sx,
                                         double* y, int sy, 
                                         double* z, int sz,
                                         double* knotXi, int sknotXi,
                                         double* knotEta, int sknotEta,
                                         double* knotZeta, int sknotZeta,
                                         int deg_xi, int deg_eta, int deg_zeta, 
                                         double* valuesN, int s_valuesN,
                                         int* indexI, int s_indexI,
                                         int* indexJ, int s_indexJ)
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
	
	int l=0; 
	int indexPoint = 0; 
	for ( int r= 0; r < sz; r++)
	{
		for ( int q =0; q < sy; q++)
		{
			for ( int p =0; p < sx; p++)
			{
            	for ( int k =0; k  < deg_zeta +1 ; k++ )
				{
		    		for ( int j=0; j < deg_eta + 1; j++ )
					{
						for ( int i =0; i < deg_xi + 1; i++)
						{
    						valuesN[l] = Nz[r][k]*Ny[q][j]*Nx[p][i] ; 
    						indexI[l]  = indexPoint ; 
    						indexJ[l]   = Spanx[p]-deg_xi+i +(Spany[q]-deg_eta+j)*nbf_xi + (Spanz[r]-deg_zeta+k)*nbf_xi*nbf_eta ; 
    						l++; 
						}
    			    }
			    }
			    indexPoint++; 
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
                                         int* indexJ, int s_indexJ)
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
	
	int l=0; 
	int indexPoint = 0; 
	for ( int r= 0; r < sz; r++)
	{
		for ( int q =0; q < sy; q++)
		{
			for ( int p =0; p < sx; p++)
			{
            	for ( int k =0; k  < deg_zeta +1 ; k++ )
				{
		    		for ( int j=0; j < deg_eta + 1; j++ )
					{
						for ( int i =0; i < deg_xi + 1; i++)
						{
    						valuesN[l] = Nz[r][k]*Ny[q][j]*Nx[p][i] ; 
    						valuesdNdx[l] = Nz[r][k]*Ny[q][j]*dNxdx[p][i]; 
    						valuesdNdy[l] = Nz[r][k]*dNydy[q][j]*Nx[p][i]; 
    						valuesdNdz[l] = dNzdz[r][k]*Ny[q][j]*Nx[p][i]; 							
    						indexI[l]     = indexPoint ;  
    						indexJ[l]     = Spanx[p]-deg_xi+i +(Spany[q]-deg_eta+j)*nbf_xi + (Spanz[r]-deg_zeta+k)*nbf_xi*nbf_eta ; 
    						l++; 
						}
    			    }
			    }
			    indexPoint++; 
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


/*
std::vector<double> FcmIntegration_K_L_TrilinearInterp(double* image, int si, int sj, int sk,
                                                  double xminIm, double xmaxIm, double yminIm, double ymaxIm, double zminIm, double zmaxIm, 
                                                  double dx, double dy, double dz,  
                                                  double thrsh, 
                                                  int deg_xi, int deg_eta, int deg_zeta, 
                                                  double* knotXi, int sknotXi, 
                                                  double* knotEta, int sknotEta, 
                                                  double* knotZeta, int sknotZeta, 
                                                  int lvlmax )
{



}     */                                              
                                                  

std::vector<double> FcmIntegrationTrilinearInterp(double* image, int si, int sj, int sk,
                                                  double* knotXiImage, int sknotXiImage, 
                                                  double* knotEtaImage, int sknotEtaImage, 
                                                  double* knotZetaImage, int sknotZetaImage,
                                                  double thrsh, 
                                                  int deg_xi, int deg_eta, int deg_zeta, 
                                                  double* knotXi, int sknotXi, 
                                                  double* knotEta, int sknotEta, 
                                                  double* knotZeta, int sknotZeta, 
                                                  int lvlmax ) 
{
    // Mesh parameters 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	//int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	int nbf_xi_elem   = deg_xi+1 ; 
	int nbf_eta_elem  = deg_eta+1 ;
	int nbf_zeta_elem = deg_zeta+1 ;  
	int nbf_elem = nbf_xi_elem*nbf_eta_elem*nbf_zeta_elem ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
	
	double* e_xi   = new double[ne_xi+1]   ;
	double* e_eta  = new double[ne_eta+1]  ;
	double* e_zeta = new double[ne_zeta+1] ; 
	for (int i=0;i< ne_xi+1; i++)
	{
		e_xi[i] = knotXi[i+deg_xi];
	}
    for(int i=0; i< ne_eta+1; i++)
    {
		e_eta[i] = knotEta[i+deg_eta];
	}
    for(int i=0; i< ne_zeta+1; i++)
    {
		e_zeta[i] = knotZeta[i+deg_zeta];
	} 
 
 
	// Gauss Integration parameters 
    int nbg_xi   = deg_xi   +1; 
    int nbg_eta  = deg_eta  +1; 
    int nbg_zeta = deg_zeta +1; 
    double* xig  = new double[nbg_xi];
    double* etag = new double[nbg_eta];
    double* zetag = new double[nbg_zeta]; 
    double* wg_xi   = new double[nbg_xi]; 
    double* wg_eta  = new double[nbg_eta]; 
    double* wg_zeta = new double[nbg_zeta]; 
    double* wg = new double[nbg_xi*nbg_eta*nbg_zeta]; 
    gauleg(-1,1,xig,wg_xi,nbg_xi);    
	gauleg(-1,1,etag,wg_eta,nbg_eta); 
    gauleg(-1,1,zetag,wg_zeta,nbg_zeta); 
    int t=0 ; 
    for (int kk=0; kk<nbg_zeta; kk++)
    {
        for (int jj=0; jj<nbg_eta; jj++)
        {
            for (int ii=0; ii<nbg_xi; ii++)
            {
                wg[t] = wg_xi[ii]*wg_eta[jj]*wg_zeta[kk] ; 
                t++; 
            }
        }
    }
    
    double* Cxig   = new double[nbg_xi]; 
    double* Cetag  = new double[nbg_eta];
    double* Czetag = new double[nbg_zeta];
    
    // Integration on a tetrahedron 
    std::vector<double> GaussTetraW ;  
    int nbgTetrahedron ; 
    Eigen::MatrixXd GaussTetraBc = GaussTetrahedron(GaussTetraW, nbgTetrahedron, std::max(deg_xi,std::max(deg_eta,deg_zeta)) ) ; 
    
    // Elementary operators
    Eigen::MatrixXd Ke(3*nbf_elem,3*nbf_elem);   // Elementary stiffness matrix 
    Eigen::MatrixXd Bc(6,3*nbf_elem); // Elementary differential matrix  
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
    Eigen::MatrixXd hooke = hookeVolume(0,0); 
 
    

    // For tri-linear image evaluation 
//    double* NxImage = new double[2];      
//    double* NyImage = new double[2];
//    double* NzImage = new double[2]; 
    
    
    
    // Memory allocation for basis functions routines output 
    double** bfOutputXi   = array_2d(2,nbg_xi);
	double** bfOutputEta  = array_2d(2,nbg_eta);
	double** bfOutputZeta = array_2d(2,nbg_zeta);
	
	double** Nxi         = array_2d(nbg_xi,nbf_xi_elem);
	double** dNxidxi     = array_2d(nbg_xi,nbf_xi_elem);
	double** Neta        = array_2d(nbg_eta,nbf_eta_elem);
	double** dNetadeta   = array_2d(nbg_eta,nbf_eta_elem);
	double** Nzeta       = array_2d(nbg_zeta,nbf_zeta_elem);
	double** dNzetadzeta = array_2d(nbg_zeta,nbf_zeta_elem);
	double* N       = new double[nbf_elem];  
	double* dNdxi   = new double[nbf_elem];  
	double* dNdeta  = new double[nbf_elem]; 
	double* dNdzeta = new double[nbf_elem]; 
	
	
    std::vector<double> integrationPoints; 
 
    MeshParameters meshParam(image, si, sj, sk,
                            knotXiImage, sknotXiImage, 
                            knotEtaImage, sknotEtaImage,
                            knotZetaImage, sknotZetaImage,
                            thrsh, 
                            deg_xi, deg_eta, deg_zeta, 
                            knotXi, sknotXi, knotEta, sknotEta, knotZeta, sknotZeta, 
                            nbg_xi, nbg_eta, nbg_zeta, 
                            xig, etag, zetag, wg, 
                            Cxig, Cetag, Czetag, 
                            bfOutputXi, bfOutputEta, bfOutputZeta,
                            Nxi, dNxidxi, Neta, dNetadeta, Nzeta, dNzetadzeta, 
                            N, dNdxi, dNdeta, dNdzeta,  
                            nbgTetrahedron, GaussTetraBc, GaussTetraW, 
                            &Ke, &Bc, 
                            NOELEM, hooke);  
         
      
    double xmin, xmax, ymin, ymax, zmin, zmax ; 
    // Loop over elements 
    for (int k=0; k<ne_zeta ;k++)
    {
        zmin = e_zeta[k]; 
        zmax = e_zeta[k+1];
        for (int j=0; j<ne_eta; j++)
        {
            ymin = e_eta[j]  ; 
            ymax = e_eta[j+1]; 
            for (int i=0; i<ne_xi; i++)
            {
                xmin = e_xi[i]   ;
                xmax = e_xi[i+1] ; 
                Cell c(xmin, xmax, ymin, ymax, zmin, zmax); 

                // Run the octree algorithm on the element 
                c.DecomposeTesselationTrilinearInterp(lvlmax, &meshParam, &integrationPoints);   
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
    delete[] wg ; 
    delete[] Cxig; 
    delete[] Cetag; 
    delete[] Czetag; 
    //delete[] NxImage ; 
    //delete[] NyImage ; 
    //delete[] NzImage ;
    
    delete_array_2d(bfOutputXi);
	delete_array_2d(bfOutputEta); 
    delete_array_2d(bfOutputZeta); 
    
    delete_array_2d(Nxi);
	delete_array_2d(dNxidxi);
	delete_array_2d(Neta);
	delete_array_2d(dNetadeta);
	delete_array_2d(Nzeta);
	delete_array_2d(dNzetadzeta); 
	delete[] N ; 
	delete[] dNdxi   ; 
	delete[] dNdeta  ; 
	delete[] dNdzeta ; 
	
    return integrationPoints ; 
} 	

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
                                             double* intBf, int sintBf ) 
{

    
    // Mesh parameters 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	int nbf_xi_elem   = deg_xi+1 ; 
	int nbf_eta_elem  = deg_eta+1 ;
	int nbf_zeta_elem = deg_zeta+1 ;  
	int nbf_elem = nbf_xi_elem*nbf_eta_elem*nbf_zeta_elem ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
	
	
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 

	
	double* e_xi   = new double[ne_xi+1]   ;
	double* e_eta  = new double[ne_eta+1]  ;
	double* e_zeta = new double[ne_zeta+1] ; 
	for (int i=0;i< ne_xi+1; i++)
	{
		e_xi[i] = knotXi[i+deg_xi];
	}
    for(int i=0; i< ne_eta+1; i++)
    {
		e_eta[i] = knotEta[i+deg_eta];
	}
    for(int i=0; i< ne_zeta+1; i++)
    {
		e_zeta[i] = knotZeta[i+deg_zeta];
	} 
 
 
	// Gauss Integration parameters 
    int nbg_xi   = deg_xi   +1; 
    int nbg_eta  = deg_eta  +1; 
    int nbg_zeta = deg_zeta +1; 
    double* xig  = new double[nbg_xi];
    double* etag = new double[nbg_eta];
    double* zetag = new double[nbg_zeta]; 
    double* wg_xi   = new double[nbg_xi]; 
    double* wg_eta  = new double[nbg_eta]; 
    double* wg_zeta = new double[nbg_zeta]; 
    double* wg = new double[nbg_xi*nbg_eta*nbg_zeta]; 
    gauleg(-1,1,xig,wg_xi,nbg_xi);    
	gauleg(-1,1,etag,wg_eta,nbg_eta); 
    gauleg(-1,1,zetag,wg_zeta,nbg_zeta); 
    int t=0 ; 
    for (int kk=0; kk<nbg_zeta; kk++)
    {
        for (int jj=0; jj<nbg_eta; jj++)
        {
            for (int ii=0; ii<nbg_xi; ii++)
            {
                wg[t] = wg_xi[ii]*wg_eta[jj]*wg_zeta[kk] ; 
                t++; 
            }
        }
    }
   
    double* Cxig   = new double[nbg_xi]; 
    double* Cetag  = new double[nbg_eta];
    double* Czetag = new double[nbg_zeta];
    
    
    // Integration on a tetrahedron 
    std::vector<double> GaussTetraW ;  
    int nbgTetrahedron ; 
    Eigen::MatrixXd GaussTetraBc = GaussTetrahedron(GaussTetraW, nbgTetrahedron, std::max(deg_xi,std::max(deg_eta,deg_zeta)) ) ; 
 
     
    // Elementary operators
    Eigen::MatrixXd Ke(3*nbf_elem,3*nbf_elem);   // Elementary stiffness matrix 
    Eigen::MatrixXd Bc(6,3*nbf_elem); // Elementary differential matrix  
    Bc.fill(0); 
    Eigen::MatrixXd hooke = hookeVolume(E, nu) ; 
 
    
    // Initializing the integral of the basis functions 
    for (int i=0; i<sintBf; i++)
    {
        intBf[i] = 0 ; 
    }
    

    
    // Memory allocation for basis functions routines output 
    double** bfOutputXi   = array_2d(2,nbg_xi);
	double** bfOutputEta  = array_2d(2,nbg_eta);
	double** bfOutputZeta = array_2d(2,nbg_zeta);
	
	double** Nxi         = array_2d(nbg_xi,nbf_xi_elem);
	double** dNxidxi     = array_2d(nbg_xi,nbf_xi_elem);
	double** Neta        = array_2d(nbg_eta,nbf_eta_elem);
	double** dNetadeta   = array_2d(nbg_eta,nbf_eta_elem);
	double** Nzeta       = array_2d(nbg_zeta,nbf_zeta_elem);
	double** dNzetadzeta = array_2d(nbg_zeta,nbf_zeta_elem);
	double* N       = new double[nbf_elem];  
	double* dNdxi   = new double[nbf_elem];  
	double* dNdeta  = new double[nbf_elem]; 
	double* dNdzeta = new double[nbf_elem]; 
 
 
    MeshParameters meshParam(image, si, sj, sk,
                            knotXiImage, sknotXiImage, 
                            knotEtaImage, sknotEtaImage,
                            knotZetaImage, sknotZetaImage,
                            thrsh, 
                            deg_xi, deg_eta, deg_zeta, 
                            knotXi, sknotXi, knotEta, sknotEta, knotZeta, sknotZeta, 
                            nbg_xi, nbg_eta, nbg_zeta, 
                            xig, etag, zetag, wg, 
                            Cxig, Cetag, Czetag, 
                            bfOutputXi, bfOutputEta, bfOutputZeta,
                            Nxi, dNxidxi, Neta, dNetadeta, Nzeta, dNzetadzeta, 
                            N, dNdxi, dNdeta, dNdzeta,  
                            nbgTetrahedron, GaussTetraBc, GaussTetraW, 
                            &Ke, &Bc, 
                            NOELEM, hooke);   
 
    
    double xmin, xmax, ymin, ymax, zmin, zmax ; 
    // Loop over elements 
    int I,J; 
    int sp_count = 0 ;  
    int ecounter = 0 ; // index of element 
    for (int k=0; k<ne_zeta ;k++)
    {
        zmin = e_zeta[k]; 
        zmax = e_zeta[k+1];
        for (int j=0; j<ne_eta; j++)
        {
            ymin = e_eta[j]  ; 
            ymax = e_eta[j+1]; 
            for (int i=0; i<ne_xi; i++)
            {
                //std::cout << "Element " << ecounter << std::endl ; 
                xmin = e_xi[i]   ; 
                xmax = e_xi[i+1] ; 
                
                // Initializing elementary stiffness matrix
                Ke.fill(0); 
                
                Cell c(xmin, xmax, ymin, ymax, zmin, zmax);    

                // Run the octree algorithm on the element  
                c.DecomposeTesselationTrilinearInterpAssembly(ecounter, i+deg_xi, j+deg_eta, k+deg_zeta, lvlmax, &meshParam, intBf); 
                
                // Elementary contribution to the global Stiffness matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM(ecounter,ibf); 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = NOELEM(ecounter,jbf); 
    					
    		            indexI[sp_count] = I ; 
    		            indexJ[sp_count] = J ; 
    		            nnz_values[sp_count] = Ke(ibf,jbf); 
    		            sp_count ++ ; 
    		            
    		            indexI[sp_count] = I+nbf ; 
    		            indexJ[sp_count] = J ; 
    		            nnz_values[sp_count] = Ke(ibf+nbf_elem,jbf); 
    		            sp_count ++ ; 
    		            
		                indexI[sp_count] = I ;
    		            indexJ[sp_count] = J +nbf;
    		            nnz_values[sp_count] = Ke(ibf,jbf+nbf_elem);
    		            sp_count ++ ;
    		            
		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = Ke(ibf+nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ;
    		                   
    		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = Ke(ibf+2*nbf_elem,jbf);
    		            sp_count ++ ;    		            
    		            
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Ke(ibf,jbf+2*nbf_elem);
    		            sp_count ++ ;
    		            
    		            indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Ke(ibf+nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ;   		            
 
     		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] =  Ke(ibf+2*nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ; 
    		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Ke(ibf+2*nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ; 
                    }
                }
                ecounter++;      		 
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
    delete[] wg ; 
    delete[] Cxig; 
    delete[] Cetag; 
    delete[] Czetag; 
 
    
    delete_array_2d(bfOutputXi);
	delete_array_2d(bfOutputEta); 
    delete_array_2d(bfOutputZeta); 
    
    delete_array_2d(Nxi);
	delete_array_2d(dNxidxi);
	delete_array_2d(Neta);
	delete_array_2d(dNetadeta);
	delete_array_2d(Nzeta);
	delete_array_2d(dNzetadzeta); 
	delete[] N ; 
	delete[] dNdxi   ; 
	delete[] dNdeta  ; 
	delete[] dNdzeta ; 
 
} 	


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
                                             double* intBf, int sintBf )
{

    // Mesh parameters 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	int nbf_xi_elem   = deg_xi+1   ; 
	int nbf_eta_elem  = deg_eta+1  ;
	int nbf_zeta_elem = deg_zeta+1 ;  
	int nbf_elem = nbf_xi_elem*nbf_eta_elem*nbf_zeta_elem ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
	
	
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 

	
	double* e_xi   = new double[ne_xi+1]   ;
	double* e_eta  = new double[ne_eta+1]  ;
	double* e_zeta = new double[ne_zeta+1] ; 
	for (int i=0;i< ne_xi+1; i++)
	{
		e_xi[i] = knotXi[i+deg_xi];
	}
    for(int i=0; i< ne_eta+1; i++)
    {
		e_eta[i] = knotEta[i+deg_eta];
	}
    for(int i=0; i< ne_zeta+1; i++)
    {
		e_zeta[i] = knotZeta[i+deg_zeta];
	} 
 
 
	// Gauss Integration parameters 
    int nbg_xi    = deg_xi   +1 ; 
    int nbg_eta   = deg_eta  +1 ; 
    int nbg_zeta  = deg_zeta +1 ; 
    double* xig   = new double[nbg_xi]    ;
    double* etag  = new double[nbg_eta]   ;
    double* zetag = new double[nbg_zeta] ; 
    double* wg_xi   = new double[nbg_xi]; 
    double* wg_eta  = new double[nbg_eta]; 
    double* wg_zeta = new double[nbg_zeta]; 
    double* wg = new double[nbg_xi*nbg_eta*nbg_zeta]; 
    gauleg(-1,1,xig,wg_xi,nbg_xi);    
	gauleg(-1,1,etag,wg_eta,nbg_eta); 
    gauleg(-1,1,zetag,wg_zeta,nbg_zeta); 
    int t=0 ; 
    for (int kk=0; kk<nbg_zeta; kk++)
    {
        for (int jj=0; jj<nbg_eta; jj++)
        {
            for (int ii=0; ii<nbg_xi; ii++)
            {
                wg[t] = wg_xi[ii]*wg_eta[jj]*wg_zeta[kk] ; 
                t++; 
            }
        }
    }
    // Integration on a tetrahedron 
    std::vector<double> GaussTetraW ;  
    int nbgTetrahedron ; 
    Eigen::MatrixXd GaussTetraBc = GaussTetrahedron(GaussTetraW, nbgTetrahedron, std::max(deg_xi,std::max(deg_eta,deg_zeta)) ) ; 
    // Material properties 
    Eigen::MatrixXd hooke = hookeVolume(E, nu) ; 

    // Initializing the integral of the basis functions 
    for (int i=0; i<sintBf; i++)
    {
        intBf[i] = 0 ; 
    }
    // Defining the shared parameters object 
    // These variables will be ony read by the threads 
    SharedParameters globalParam(image, si, sj, sk,
                                knotXiImage, sknotXiImage, 
                                knotEtaImage, sknotEtaImage,
                                knotZetaImage, sknotZetaImage,
                                thrsh, 
                                deg_xi, deg_eta, deg_zeta, 
                                knotXi, sknotXi, knotEta, sknotEta, knotZeta, sknotZeta, 
                                nbg_xi, nbg_eta, nbg_zeta, 
                                xig, etag, zetag, wg,  
                                NOELEM, 
                                nbgTetrahedron, GaussTetraBc, GaussTetraW, 
                                hooke ); 


    // Loop over elements 
    
    //int sp_count = 0 ;  
    //int ecounter = 0 ; // index of element 
    #pragma omp parallel 
    {
        #pragma omp for  
        for (int k=0; k<ne_zeta ;k++)
        {
            double zmin = e_zeta[k]; 
            double zmax = e_zeta[k+1];
            for (int j=0; j<ne_eta; j++)
            {
                double ymin = e_eta[j]  ; 
                double ymax = e_eta[j+1]; 
                for (int i=0; i<ne_xi; i++)
                {
                    //std::cout << "Element " << ecounter << std::endl ; 
                    
                    int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ; 
                    
                    double xmin = e_xi[i]   ; 
                    double xmax = e_xi[i+1] ; 
                    
                    
                    // Elementary operators 
                    Eigen::MatrixXd Ke(3*nbf_elem,3*nbf_elem);   // Elementary stiffness matrix 
                    Eigen::MatrixXd Bc(6,3*nbf_elem); // Elementary differential matrix  
                    // Initializing elementary stiffness matrix
                    Ke.fill(0); 
                    Bc.fill(0); 
    
                    // Memory allocation for basis functions routines output 
                    double* Cxig   = new double[nbg_xi]; 
                    double* Cetag  = new double[nbg_eta];
                    double* Czetag = new double[nbg_zeta];
        
                    double** bfOutputXi   = array_2d(2,nbg_xi);
                	double** bfOutputEta  = array_2d(2,nbg_eta);
                	double** bfOutputZeta = array_2d(2,nbg_zeta);
                	
                	double** Nxi         = array_2d(nbg_xi,nbf_xi_elem);
                	double** dNxidxi     = array_2d(nbg_xi,nbf_xi_elem);
                	double** Neta        = array_2d(nbg_eta,nbf_eta_elem);
                	double** dNetadeta   = array_2d(nbg_eta,nbf_eta_elem);
                	double** Nzeta       = array_2d(nbg_zeta,nbf_zeta_elem);
                	double** dNzetadzeta = array_2d(nbg_zeta,nbf_zeta_elem);
                	double* N       = new double[nbf_elem];  
                	double* dNdxi   = new double[nbf_elem];  
                	double* dNdeta  = new double[nbf_elem]; 
                	double* dNdzeta = new double[nbf_elem];  
                	
                	
            	    PrivateParameters elemParam( Cxig, Cetag, Czetag,  
                    bfOutputXi, bfOutputEta, bfOutputZeta, 
                    Nxi, dNxidxi, Neta, dNetadeta, Nzeta, dNzetadzeta,  
                    N, dNdxi, dNdeta, dNdzeta,
                    &Ke, &Bc );  
                           
                    // Defining the integration cell: starting from the element 
                    Cell c(xmin, xmax, ymin, ymax, zmin, zmax);    
    
                    // Run the octree algorithm on the element  
                    // This routines fills Ke 
                    c.DecomposeTesselationTrilinearInterpAssemblyParallel(ecounter, i+deg_xi, j+deg_eta, k+deg_zeta, lvlmax, 
                                                                  &globalParam, &elemParam, intBf );  
                                                                                  
                    // Assembly 
                    // Elementary contribution to the global Stiffness matrix 
                    int I,J; 
                    int sg = 9*nbf_elem*nbf_elem*ecounter  ;
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        I = NOELEM(ecounter,ibf); 
                        for (int jbf=0; jbf<nbf_elem; jbf++)
                        {
        					J = NOELEM(ecounter,jbf); 
        					
        		            indexI[sg] = I ; 
        		            indexJ[sg] = J ; 
        		            nnz_values[sg] = Ke(ibf,jbf); 
        		            sg ++ ; 
        		            
        		            indexI[sg] = I+nbf ; 
        		            indexJ[sg] = J ; 
        		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf); 
        		            sg ++ ; 
        		            
    		                indexI[sg] = I ;
        		            indexJ[sg] = J +nbf;
        		            nnz_values[sg] = Ke(ibf,jbf+nbf_elem);
        		            sg ++ ;
        		            
    		            	indexI[sg] = I+nbf ;
        		            indexJ[sg] = J+nbf ;
        		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf+nbf_elem);
        		            sg ++ ;
        		                   
        		            indexI[sg] = I+2*nbf ;
        		            indexJ[sg] = J ;
        		            nnz_values[sg] = Ke(ibf+2*nbf_elem,jbf);
        		            sg ++ ;    		            
        		            
        		            indexI[sg] = I ;
        		            indexJ[sg] = J+2*nbf ;
        		            nnz_values[sg] =  Ke(ibf,jbf+2*nbf_elem);
        		            sg ++ ;
        		            
        		            indexI[sg] = I+nbf ;
        		            indexJ[sg] = J+2*nbf ;
        		            nnz_values[sg] =  Ke(ibf+nbf_elem,jbf+2*nbf_elem);
        		            sg ++ ;   		            
     
         		            indexI[sg] = I+2*nbf ;
        		            indexJ[sg] = J+nbf ;
        		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+nbf_elem);
        		            sg ++ ; 
        		               		            
    		                indexI[sg] = I+2*nbf ;
        		            indexJ[sg] = J+2*nbf ;
        		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+2*nbf_elem);
        		            sg ++ ; 
                        }
                    }
                    ecounter++;   
                    // Cleaning memory allocated by the element 
                    delete[] Cxig   ; 
                    delete[] Cetag  ; 
                    delete[] Czetag ;
                    delete_array_2d(bfOutputXi);
                	delete_array_2d(bfOutputEta); 
                    delete_array_2d(bfOutputZeta); 
        
                    delete_array_2d(Nxi);
                	delete_array_2d(dNxidxi);
                	delete_array_2d(Neta);
                	delete_array_2d(dNetadeta);
                	delete_array_2d(Nzeta);
                	delete_array_2d(dNzetadzeta); 
                	delete[] N ; 
                	delete[] dNdxi   ; 
                	delete[] dNdeta  ; 
                	delete[] dNdzeta ;                     		 
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
    delete[] wg ;     
}

void Laplacian_Structured(int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz)
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
    int nbg_xi     = deg_xi   + 1 ; 
    int nbg_eta    = deg_eta  + 1 ; 
    int nbg_zeta   = deg_zeta + 1 ; 
    double* xig    = new double[nbg_xi]   ; 
    double* wgxi   = new double[nbg_xi]   ; 
    double* etag   = new double[nbg_eta]  ; 
    double* wgeta  = new double[nbg_eta]  ;
    double* zetag  = new double[nbg_zeta] ; 
    double* wgzeta = new double[nbg_zeta] ;
    gauleg(-1,1,xig,wgxi,nbg_xi); 
    gauleg(-1,1,etag,wgeta,nbg_eta);  
    gauleg(-1,1,zetag,wgzeta,nbg_zeta);
    
    int nbg_xi_tot   =  nbg_xi*ne_xi     ;
	int nbg_eta_tot  =  nbg_eta*ne_eta   ;
	int nbg_zeta_tot =  nbg_zeta*ne_zeta ; 
	
    double* xg = new double[nbg_xi_tot]   ; 
    double* yg = new double[nbg_eta_tot]  ; 
    double* zg = new double[nbg_zeta_tot] ;
    double* wx = new double[nbg_xi_tot]   ;
	double* wy = new double[nbg_eta_tot]  ; 
	double* wz = new double[nbg_zeta_tot] ; 
	double* mesxg = new double[nbg_xi_tot] ; 
	double* mesyg = new double[nbg_eta_tot] ; 
	double* meszg = new double[nbg_zeta_tot] ;   
	
	double **Nxi   = array_2d(2,nbg_xi)   ;
	double **Neta  = array_2d(2,nbg_eta)  ;
	double **Nzeta = array_2d(2,nbg_zeta) ; 
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	double** Nxig     = array_2d(nbg_xi_tot, nbg_xi);
	double** dNdxig   = array_2d(nbg_xi_tot, nbg_xi);
	double** Netag    = array_2d(nbg_eta_tot, nbg_eta);
	double** dNdetag  = array_2d(nbg_eta_tot, nbg_eta);
	double** Nzetag   = array_2d(nbg_zeta_tot, nbg_zeta);
	double** dNdzetag = array_2d(nbg_zeta_tot, nbg_zeta);    
	
	
	double* dNdxi   = new double[nbf_elem]; 
	double* dNdeta  = new double[nbf_elem];
	double* dNdzeta = new double[nbf_elem];
	
	
	double mes;
	int k=0;
	// Loop on 1D xi elements
	for ( int sxi = deg_xi ; sxi < sknotXi - deg_xi -1 ; sxi++ )
	{
		mes = knotXi[sxi+1] - knotXi[sxi];
		for ( int ii=0; ii < nbg_xi; ii++ )
		{
			xg[k]    = knotXi[sxi] + 0.5*(xig[ii]+1)*mes ;
			wx[k]    = wgxi[ii];
			mesxg[k] = mes;
			dersbasisfuns(deg_xi,  knotXi  ,  xg[k], sxi, 1, Nxi);
			for ( int l =0 ; l < nbg_xi; l++ )
			{
				Nxig[k][l]   = Nxi[0][l] ;
				dNdxig[k][l] = Nxi[1][l] ;
			}
			k++;
		}
	}
	k = 0;
	// Loop on 1D eta elements
	for  ( int seta = deg_eta ; seta < sknotEta - deg_eta -1 ; seta++ )
	{
		mes = knotEta[seta+1] - knotEta[seta];
		for ( int jj = 0; jj < nbg_eta; jj ++ )
		{
			yg[k] = knotEta[seta] + 0.5*(etag[jj]+1)*mes;
			wy[k] = wgeta[jj];
			mesyg[k] = mes;
			dersbasisfuns(deg_eta,  knotEta  ,  yg[k], seta, 1, Neta);
			for ( int l =0 ; l < nbg_eta; l++ )
			{
				Netag[k][l]   = Neta[0][l] ;
				dNdetag[k][l] = Neta[1][l] ;
			}
			k++;
		}
	}
	k = 0 ; 
	// Loop on 1D zeta elements 
	for  ( int szeta = deg_zeta ; szeta < sknotZeta - deg_zeta -1 ; szeta++ )
	{
		mes = knotZeta[szeta+1] - knotZeta[szeta];
		for ( int kk = 0; kk < nbg_zeta; kk ++ )
		{
			zg[k] = knotZeta[szeta] + 0.5*(zetag[kk]+1)*mes;
			wz[k] = wgzeta[kk];
			meszg[k] = mes;
			dersbasisfuns(deg_zeta,  knotZeta  ,  zg[k], szeta, 1, Nzeta);
			for ( int l =0 ; l < nbg_zeta; l++ )
			{
				Nzetag[k][l]   = Nzeta[0][l] ;
				dNdzetag[k][l] = Nzeta[1][l] ;
			}
			k++;
		}
	}	
	
	// Elementary operators
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
    Eigen::MatrixXd Leb(nbf_elem,nbf_elem);  // Block of the elementary Laplacian operator 
    Eigen::MatrixXd dNeb(3,nbf_elem);       // Block elementary differential operator 
     
    
    
    // Assembly of the matrix 
    int I,J;
    int sp_count = 0 ;  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {        
            for (int i=0; i<ne_xi; i++)
            {
                // Current element (i,j,k)
                Leb.fill(0) ;  
                // Loop on the integration points 
                for (int kk=k*nbg_zeta; kk<(k+1)*nbg_zeta; kk++)
                {
                    for (int jj=j*nbg_eta; jj<(j+1)*nbg_eta; jj++)
                    {
                        for (int ii=i*nbg_xi; ii<(i+1)*nbg_xi; ii++)
                        {
                            // Getting the trivariate basis functions on the Gauss integration point 
                            bfc = 0 ; 
                            for (int kbf=0; kbf < nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf < nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf < nbf_elem_xi; ibf++)
                                    {
                                        dNdxi[bfc]   = Nzetag[kk][kbf]*Netag[jj][jbf]*dNdxig[ii][ibf] ;   
                                        dNdeta[bfc]  = Nzetag[kk][kbf]*dNdetag[jj][jbf]*Nxig[ii][ibf] ;   
                                        dNdzeta[bfc] = dNdzetag[kk][kbf]*Netag[jj][jbf]*Nxig[ii][ibf] ; 
                                        bfc++; 
                                    }
                                }
                            }  
                        
                            // Setting the elementary differential operator 
                            for ( int ibf =0; ibf < nbf_elem; ibf++)
            				{
                    			dNeb(0,ibf)  = dNdxi[ibf] ; 
                				dNeb(1,ibf)  = dNdeta[ibf] ;  
                				dNeb(2,ibf)  = dNdzeta[ibf] ; 
            				}
                            Leb = Leb + dNeb.transpose()*dNeb*(wx[ii]*wy[jj]*wz[kk]*mesxg[ii]*mesyg[jj]*meszg[kk]/8) ; 
                        }
                    }
                }
                // Elementary contribution to the global matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM(ecounter,ibf); 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = NOELEM(ecounter,jbf);
    					
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = Leb(ibf,jbf);
    		            sp_count ++ ;  	
    		            
		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = Leb(ibf,jbf);
    		            sp_count ++ ;
    		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Leb(ibf,jbf);
    		            sp_count ++ ; 
                    }
                }
                ecounter++;                            
            }
        }
    }
 
 
    delete[] xig   ;  
    delete[] wgxi  ; 
    delete[] etag  ; 
    delete[] wgeta ;
    delete[] zetag ; 
    delete[] wgzeta ;
    
    
    delete[] xg  ; 
    delete[] yg  ; 
    delete[] zg  ;
    delete[] wx  ;
	delete[] wy  ; 
	delete[] wz  ; 
	delete[] mesxg ; 
	delete[] mesyg ; 
	delete[] meszg ;   
	
	delete_array_2d(Nxi); 
	delete_array_2d(Neta);  
	delete_array_2d(Nzeta);  
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	delete_array_2d(Nxig)     ;
	delete_array_2d(dNdxig)   ;
	delete_array_2d(Netag)    ;
	delete_array_2d(dNdetag)  ;
	delete_array_2d(Nzetag)   ;
	delete_array_2d(dNdzetag) ;    
		
	delete[] dNdxi    ; 
	delete[] dNdeta   ;
	delete[] dNdzeta  ;
   
}

void Laplacian_Structured_Parallel(int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz)
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
    int nbg_xi     = deg_xi   + 1 ; 
    int nbg_eta    = deg_eta  + 1 ; 
    int nbg_zeta   = deg_zeta + 1 ; 
    double* xig    = new double[nbg_xi]   ; 
    double* wgxi   = new double[nbg_xi]   ; 
    double* etag   = new double[nbg_eta]  ; 
    double* wgeta  = new double[nbg_eta]  ;
    double* zetag  = new double[nbg_zeta] ; 
    double* wgzeta = new double[nbg_zeta] ;
    gauleg(-1,1,xig,wgxi,nbg_xi); 
    gauleg(-1,1,etag,wgeta,nbg_eta);  
    gauleg(-1,1,zetag,wgzeta,nbg_zeta);
    
    int nbg_xi_tot   =  nbg_xi*ne_xi     ;
	int nbg_eta_tot  =  nbg_eta*ne_eta   ;
	int nbg_zeta_tot =  nbg_zeta*ne_zeta ; 
	
    double* xg = new double[nbg_xi_tot]   ; 
    double* yg = new double[nbg_eta_tot]  ; 
    double* zg = new double[nbg_zeta_tot] ;
    double* wx = new double[nbg_xi_tot]   ;
	double* wy = new double[nbg_eta_tot]  ; 
	double* wz = new double[nbg_zeta_tot] ; 
	double* mesxg = new double[nbg_xi_tot] ; 
	double* mesyg = new double[nbg_eta_tot] ; 
	double* meszg = new double[nbg_zeta_tot] ;   
	
	double **Nxi   = array_2d(2,nbg_xi)   ;
	double **Neta  = array_2d(2,nbg_eta)  ;
	double **Nzeta = array_2d(2,nbg_zeta) ; 
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	double** Nxig     = array_2d(nbg_xi_tot, nbg_xi);
	double** dNdxig   = array_2d(nbg_xi_tot, nbg_xi);
	double** Netag    = array_2d(nbg_eta_tot, nbg_eta);
	double** dNdetag  = array_2d(nbg_eta_tot, nbg_eta);
	double** Nzetag   = array_2d(nbg_zeta_tot, nbg_zeta);
	double** dNdzetag = array_2d(nbg_zeta_tot, nbg_zeta);    
	
 	
	double mes;
	int k=0;
	// Loop on 1D xi elements
	for ( int sxi = deg_xi ; sxi < sknotXi - deg_xi -1 ; sxi++ )
	{
		mes = knotXi[sxi+1] - knotXi[sxi];
		for ( int ii=0; ii < nbg_xi; ii++ )
		{
			xg[k]    = knotXi[sxi] + 0.5*(xig[ii]+1)*mes ;
			wx[k]    = wgxi[ii];
			mesxg[k] = mes;
			dersbasisfuns(deg_xi,  knotXi  ,  xg[k], sxi, 1, Nxi);
			for ( int l =0 ; l < nbg_xi; l++ )
			{
				Nxig[k][l]   = Nxi[0][l] ;
				dNdxig[k][l] = Nxi[1][l] ;
			}
			k++;
		}
	}
	k = 0;
	// Loop on 1D eta elements
	for  ( int seta = deg_eta ; seta < sknotEta - deg_eta -1 ; seta++ )
	{
		mes = knotEta[seta+1] - knotEta[seta];
		for ( int jj = 0; jj < nbg_eta; jj ++ )
		{
			yg[k] = knotEta[seta] + 0.5*(etag[jj]+1)*mes;
			wy[k] = wgeta[jj];
			mesyg[k] = mes;
			dersbasisfuns(deg_eta,  knotEta  ,  yg[k], seta, 1, Neta);
			for ( int l =0 ; l < nbg_eta; l++ )
			{
				Netag[k][l]   = Neta[0][l] ;
				dNdetag[k][l] = Neta[1][l] ;
			}
			k++;
		}
	}
	k = 0 ; 
	// Loop on 1D zeta elements 
	for  ( int szeta = deg_zeta ; szeta < sknotZeta - deg_zeta -1 ; szeta++ )
	{
		mes = knotZeta[szeta+1] - knotZeta[szeta];
		for ( int kk = 0; kk < nbg_zeta; kk ++ )
		{
			zg[k] = knotZeta[szeta] + 0.5*(zetag[kk]+1)*mes;
			wz[k] = wgzeta[kk];
			meszg[k] = mes;
			dersbasisfuns(deg_zeta,  knotZeta  ,  zg[k], szeta, 1, Nzeta);
			for ( int l =0 ; l < nbg_zeta; l++ )
			{
				Nzetag[k][l]   = Nzeta[0][l] ;
				dNdzetag[k][l] = Nzeta[1][l] ;
			}
			k++;
		}
	}	
 
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
 
    
    // Assembly of the matrix 
    #pragma omp parallel 
    {
        #pragma omp for 
            for (int k=0; k<ne_zeta ;k++)
            {   
                for (int j=0; j<ne_eta; j++)
                {        
                    for (int i=0; i<ne_xi; i++)
                    {
                        int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ;    
                        Eigen::MatrixXd Leb(nbf_elem,nbf_elem);  // Block of the elementary Laplacian operator 
                        Eigen::MatrixXd dNeb(3,nbf_elem);       // Block elementary differential operator 
                        double* dNdxi   = new double[nbf_elem]; 
                    	double* dNdeta  = new double[nbf_elem];
                    	double* dNdzeta = new double[nbf_elem];
                        // Current element (i,j,k)
                        Leb.fill(0) ;  
                        // Loop on the integration points 
                        int bfc ; 
                        for (int kk=k*nbg_zeta; kk<(k+1)*nbg_zeta; kk++)
                        {
                            for (int jj=j*nbg_eta; jj<(j+1)*nbg_eta; jj++)
                            {
                                for (int ii=i*nbg_xi; ii<(i+1)*nbg_xi; ii++)
                                {
                                    // Getting the trivariate basis functions on the Gauss integration point 
                                    bfc = 0 ; 
                                    for (int kbf=0; kbf < nbf_elem_zeta; kbf++)
                                    {
                                        for (int jbf=0; jbf < nbf_elem_eta; jbf++)
                                        {
                                            for (int ibf=0; ibf < nbf_elem_xi; ibf++)
                                            {
                                                dNdxi[bfc]   = Nzetag[kk][kbf]*Netag[jj][jbf]*dNdxig[ii][ibf] ;   
                                                dNdeta[bfc]  = Nzetag[kk][kbf]*dNdetag[jj][jbf]*Nxig[ii][ibf] ;   
                                                dNdzeta[bfc] = dNdzetag[kk][kbf]*Netag[jj][jbf]*Nxig[ii][ibf] ; 
                                                bfc++; 
                                            }
                                        }
                                    }  
                                
                                // Setting the elementary differential operator 
                                for ( int ibf =0; ibf < nbf_elem; ibf++)
                				{
                        			dNeb(0,ibf)  = dNdxi[ibf] ; 
                    				dNeb(1,ibf)  = dNdeta[ibf] ;  
                    				dNeb(2,ibf)  = dNdzeta[ibf] ; 
                				}
                                Leb = Leb + dNeb.transpose()*dNeb*(wx[ii]*wy[jj]*wz[kk]*mesxg[ii]*mesyg[jj]*meszg[kk]/8) ; 
                                }
                            }
                        }
                        int I,J ; 
                        int sg = 3*nbf_elem*nbf_elem*ecounter ; 
                        // Elementary contribution to the global matrix 
                        for (int ibf=0; ibf<nbf_elem; ibf++)
                        {
                            I = NOELEM(ecounter,ibf); 
                            for (int jbf=0; jbf<nbf_elem; jbf++)
                            {
            					J = NOELEM(ecounter,jbf);
            					
            		            indexI[sg] = I ;
            		            indexJ[sg] = J ;
            		            nnz_values[sg] = Leb(ibf,jbf);
            		            sg ++ ;  	
            		            
        		            	indexI[sg] = I+nbf ;
            		            indexJ[sg] = J+nbf ;
            		            nnz_values[sg] = Leb(ibf,jbf);
            		            sg ++ ;
            		               		            
        		                indexI[sg] = I+2*nbf ;
            		            indexJ[sg] = J+2*nbf ;
            		            nnz_values[sg] =  Leb(ibf,jbf);
            		            sg ++ ; 
                            }
                        }
                        delete[] dNdxi    ; 
                    	delete[] dNdeta   ;
                    	delete[] dNdzeta  ;                         
                    }
                }
            }                
    }
 
    delete[] xig   ;  
    delete[] wgxi  ; 
    delete[] etag  ; 
    delete[] wgeta ;
    delete[] zetag ; 
    delete[] wgzeta ;
    
    
    delete[] xg  ; 
    delete[] yg  ; 
    delete[] zg  ;
    delete[] wx  ;
	delete[] wy  ; 
	delete[] wz  ; 
	delete[] mesxg ; 
	delete[] mesyg ; 
	delete[] meszg ;   
	
	delete_array_2d(Nxi); 
	delete_array_2d(Neta);  
	delete_array_2d(Nzeta);  
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	delete_array_2d(Nxig)     ;
	delete_array_2d(dNdxig)   ;
	delete_array_2d(Netag)    ;
	delete_array_2d(dNdetag)  ;
	delete_array_2d(Nzetag)   ;
	delete_array_2d(dNdzetag) ;    

}


void HomogeneousStiffness(double E, double nu, 
                          int deg_xi, int deg_eta, int deg_zeta,
                          double* knotXi, int sknotXi, 
                          double* knotEta, int sknotEta, 
                          double* knotZeta, int sknotZeta, 
                          int* indexI, int nnzI, 
                          int* indexJ, int nnzJ,
                          double* nnz_values, int nnz)  
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
    int nbg_xi     = deg_xi   + 1 ; 
    int nbg_eta    = deg_eta  + 1 ; 
    int nbg_zeta   = deg_zeta + 1 ; 
    double* xig    = new double[nbg_xi]   ; 
    double* wgxi   = new double[nbg_xi]   ; 
    double* etag   = new double[nbg_eta]  ; 
    double* wgeta  = new double[nbg_eta]  ;
    double* zetag  = new double[nbg_zeta] ; 
    double* wgzeta = new double[nbg_zeta] ;
    gauleg(-1,1,xig,wgxi,nbg_xi); 
    gauleg(-1,1,etag,wgeta,nbg_eta);  
    gauleg(-1,1,zetag,wgzeta,nbg_zeta);
    
    int nbg_xi_tot   =  nbg_xi*ne_xi     ;
	int nbg_eta_tot  =  nbg_eta*ne_eta   ;
	int nbg_zeta_tot =  nbg_zeta*ne_zeta ; 
	
    double* xg = new double[nbg_xi_tot]   ; 
    double* yg = new double[nbg_eta_tot]  ; 
    double* zg = new double[nbg_zeta_tot] ;
    double* wx = new double[nbg_xi_tot]   ;
	double* wy = new double[nbg_eta_tot]  ; 
	double* wz = new double[nbg_zeta_tot] ; 
	double* mesxg = new double[nbg_xi_tot] ; 
	double* mesyg = new double[nbg_eta_tot] ; 
	double* meszg = new double[nbg_zeta_tot] ;   
	
	double **Nxi   = array_2d(2,nbg_xi)   ;
	double **Neta  = array_2d(2,nbg_eta)  ;
	double **Nzeta = array_2d(2,nbg_zeta) ; 
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	double** Nxig     = array_2d(nbg_xi_tot, nbg_xi);
	double** dNdxig   = array_2d(nbg_xi_tot, nbg_xi);
	double** Netag    = array_2d(nbg_eta_tot, nbg_eta);
	double** dNdetag  = array_2d(nbg_eta_tot, nbg_eta);
	double** Nzetag   = array_2d(nbg_zeta_tot, nbg_zeta);
	double** dNdzetag = array_2d(nbg_zeta_tot, nbg_zeta);    
	
 	
	double mes;
	int k=0;
	// Loop on 1D xi elements
	for ( int sxi = deg_xi ; sxi < sknotXi - deg_xi -1 ; sxi++ )
	{
		mes = knotXi[sxi+1] - knotXi[sxi];
		for ( int ii=0; ii < nbg_xi; ii++ )
		{
			xg[k]    = knotXi[sxi] + 0.5*(xig[ii]+1)*mes ;
			wx[k]    = wgxi[ii];
			mesxg[k] = mes;
			dersbasisfuns(deg_xi,  knotXi  ,  xg[k], sxi, 1, Nxi);
			for ( int l =0 ; l < nbg_xi; l++ )
			{
				Nxig[k][l]   = Nxi[0][l] ;
				dNdxig[k][l] = Nxi[1][l] ;
			}
			k++;
		}
	}
	k = 0;
	// Loop on 1D eta elements
	for  ( int seta = deg_eta ; seta < sknotEta - deg_eta -1 ; seta++ )
	{
		mes = knotEta[seta+1] - knotEta[seta];
		for ( int jj = 0; jj < nbg_eta; jj ++ )
		{
			yg[k] = knotEta[seta] + 0.5*(etag[jj]+1)*mes;
			wy[k] = wgeta[jj];
			mesyg[k] = mes;
			dersbasisfuns(deg_eta,  knotEta  ,  yg[k], seta, 1, Neta);
			for ( int l =0 ; l < nbg_eta; l++ )
			{
				Netag[k][l]   = Neta[0][l] ;
				dNdetag[k][l] = Neta[1][l] ;
			}
			k++;
		}
	}
	k = 0 ; 
	// Loop on 1D zeta elements 
	for  ( int szeta = deg_zeta ; szeta < sknotZeta - deg_zeta -1 ; szeta++ )
	{
		mes = knotZeta[szeta+1] - knotZeta[szeta];
		for ( int kk = 0; kk < nbg_zeta; kk ++ )
		{
			zg[k] = knotZeta[szeta] + 0.5*(zetag[kk]+1)*mes;
			wz[k] = wgzeta[kk];
			meszg[k] = mes;
			dersbasisfuns(deg_zeta,  knotZeta  ,  zg[k], szeta, 1, Nzeta);
			for ( int l =0 ; l < nbg_zeta; l++ )
			{
				Nzetag[k][l]   = Nzeta[0][l] ;
				dNdzetag[k][l] = Nzeta[1][l] ;
			}
			k++;
		}
	}	
 
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
	Eigen::MatrixXd hooke = hookeVolume(E, nu) ; 
 
    
    // Assembly of the matrix 
    #pragma omp parallel 
    {
        #pragma omp for 
            for (int k=0; k<ne_zeta ;k++)
            {   
                for (int j=0; j<ne_eta; j++)
                {        
                    for (int i=0; i<ne_xi; i++)
                    {
                        int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ;    
                        
                        Eigen::MatrixXd Ke(3*nbf_elem,3*nbf_elem);   // Elementary stiffness matrix 
                        Eigen::MatrixXd Be(6,3*nbf_elem);   // Elementary differential matrix 
                
                        double* dNdxi   = new double[nbf_elem]; 
                    	double* dNdeta  = new double[nbf_elem];
                    	double* dNdzeta = new double[nbf_elem];
                        // Current element (i,j,k)
                        Ke.fill(0) ;
                        Be.fill(0) ;   
                        // Loop on the integration points 
                        int bfc ; 
                        for (int kk=k*nbg_zeta; kk<(k+1)*nbg_zeta; kk++)
                        {
                            for (int jj=j*nbg_eta; jj<(j+1)*nbg_eta; jj++)
                            {
                                for (int ii=i*nbg_xi; ii<(i+1)*nbg_xi; ii++)
                                {
                                    // Getting the trivariate basis functions on the Gauss integration point 
                                    bfc = 0 ; 
                                    for (int kbf=0; kbf < nbf_elem_zeta; kbf++)
                                    {
                                        for (int jbf=0; jbf < nbf_elem_eta; jbf++)
                                        {
                                            for (int ibf=0; ibf < nbf_elem_xi; ibf++)
                                            {
                                                dNdxi[bfc]   = Nzetag[kk][kbf]*Netag[jj][jbf]*dNdxig[ii][ibf] ;   
                                                dNdeta[bfc]  = Nzetag[kk][kbf]*dNdetag[jj][jbf]*Nxig[ii][ibf] ;   
                                                dNdzeta[bfc] = dNdzetag[kk][kbf]*Netag[jj][jbf]*Nxig[ii][ibf] ; 
                                                bfc++; 
                                            }
                                        }
                                    }  
                                // Setting the elemental basis function matrix 
                                for ( int ibf =0; ibf < nbf_elem; ibf++)
             					{
                     				Be( 0, ibf )           = dNdxi[ibf];
                    				Be( 1, ibf+nbf_elem)   = dNdeta[ibf] ;
                    				Be( 2, ibf+2*nbf_elem) = dNdzeta[ibf] ;
                    				
                    				Be( 3, ibf)            = dNdeta[ibf] ;
                    				Be( 3, ibf+nbf_elem)   = dNdxi[ibf];
                    				
                    				Be( 4, ibf)            = dNdzeta[ibf] ;
                    				Be( 4, ibf+2*nbf_elem) = dNdxi[ibf];  
                    				
                    				Be( 5, ibf+nbf_elem)   = dNdzeta[ibf] ;
                    				Be( 5, ibf+2*nbf_elem) = dNdeta[ibf] ; 
             					}                                
                                // Setting the elementary differential operator
                                Ke  = Ke + Be.transpose() * hooke * Be *(wx[ii]*wy[jj]*wz[kk]*mesxg[ii]*mesyg[jj]*meszg[kk]/8) ;
                                }
                            }
                        }
                        int I,J ; 
                        int sg = 9*nbf_elem*nbf_elem*ecounter ; 
                        // Elementary contribution to the global matrix 
                        for (int ibf=0; ibf<nbf_elem; ibf++)
                        {
                            I = NOELEM(ecounter,ibf); 
                            for (int jbf=0; jbf<nbf_elem; jbf++)
                            {
            					J = NOELEM(ecounter,jbf);
            					
             		            indexI[sg] = I ;
             		            indexJ[sg] = J ;
             		            nnz_values[sg] = Ke(ibf,jbf);
             		            sg ++ ;
             		            
             		            indexI[sg] = I+nbf ;
             		            indexJ[sg] = J ;
             		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf);
             		            sg ++ ;
             		            
         		                indexI[sg] = I ;
             		            indexJ[sg] = J +nbf;
             		            nnz_values[sg] = Ke(ibf,jbf+nbf_elem);
             		            sg ++ ;
             		            
         		            	indexI[sg] = I+nbf ;
             		            indexJ[sg] = J+nbf ;
             		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf+nbf_elem);
             		            sg ++ ;
             		                   
             		            indexI[sg] = I+2*nbf ;
             		            indexJ[sg] = J ;
             		            nnz_values[sg] = Ke(ibf+2*nbf_elem,jbf);
             		            sg ++ ;    		            
             		            
             		            indexI[sg] = I ;
             		            indexJ[sg] = J+2*nbf ;
             		            nnz_values[sg] =  Ke(ibf,jbf+2*nbf_elem);
             		            sg ++ ;
             		            
             		            indexI[sg] = I+nbf ;
             		            indexJ[sg] = J+2*nbf ;
             		            nnz_values[sg] =  Ke(ibf+nbf_elem,jbf+2*nbf_elem);
             		            sg ++ ;   		            
          
              		            indexI[sg] = I+2*nbf ;
             		            indexJ[sg] = J+nbf ;
             		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+nbf_elem);
             		            sg ++ ; 
             		               		            
         		                indexI[sg] = I+2*nbf ;
             		            indexJ[sg] = J+2*nbf ;
             		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+2*nbf_elem);
             		            sg ++ ; 
                            }
                        }
                        delete[] dNdxi    ; 
                    	delete[] dNdeta   ;
                    	delete[] dNdzeta  ;                         
                    }
                }
            }                
    }
 
    delete[] xig   ;  
    delete[] wgxi  ; 
    delete[] etag  ; 
    delete[] wgeta ;
    delete[] zetag ; 
    delete[] wgzeta ;
    
    
    delete[] xg  ; 
    delete[] yg  ; 
    delete[] zg  ;
    delete[] wx  ;
	delete[] wy  ; 
	delete[] wz  ; 
	delete[] mesxg ; 
	delete[] mesyg ; 
	delete[] meszg ;   
	
	delete_array_2d(Nxi); 
	delete_array_2d(Neta);  
	delete_array_2d(Nzeta);  
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	delete_array_2d(Nxig)     ;
	delete_array_2d(dNdxig)   ;
	delete_array_2d(Netag)    ;
	delete_array_2d(dNdetag)  ;
	delete_array_2d(Nzetag)   ;
	delete_array_2d(dNdzetag) ;    

}

              


              


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
                  double* rhs, int ndof)            
{
    // Projection of field 2 (defined as a B-spline displacement) ont the the nodes of the B-spline mesh 1 
    // This function returns the L2-projection system to be solved 
    // LHS: Sparse Mass matrix of the final mesh 
    // RHS: Integral of the projected displacement weighted by the B-spline basis functions  
    
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
	
    //int nbg_xi     = deg_xi   + 1 ; 
    //int nbg_eta    = deg_eta  + 1 ; 
    //int nbg_zeta   = deg_zeta + 1 ; 
    
    // As an L2 projection is performed 
    // We consider an integration rule that takes into account the regularity of both fields 
    
    int nbg_xi   = std::max(deg_xi,   deg_xi2)   + 1 ; 
    int nbg_eta  = std::max(deg_eta,  deg_eta2)  + 1 ; 
    int nbg_zeta = std::max(deg_zeta, deg_zeta2) + 1 ; 
    
    
    double* xig    = new double[nbg_xi]   ; 
    double* wgxi   = new double[nbg_xi]   ; 
    double* etag   = new double[nbg_eta]  ; 
    double* wgeta  = new double[nbg_eta]  ;
    double* zetag  = new double[nbg_zeta] ; 
    double* wgzeta = new double[nbg_zeta] ;
    
    gauleg(-1,1,xig,wgxi,nbg_xi); 
    gauleg(-1,1,etag,wgeta,nbg_eta);  
    gauleg(-1,1,zetag,wgzeta,nbg_zeta);
    
    int nbg_xi_tot   =  nbg_xi*ne_xi     ;
	int nbg_eta_tot  =  nbg_eta*ne_eta   ;
	int nbg_zeta_tot =  nbg_zeta*ne_zeta ; 
	
    double* xg = new double[nbg_xi_tot]   ; 
    double* yg = new double[nbg_eta_tot]  ; 
    double* zg = new double[nbg_zeta_tot] ;
    double* wx = new double[nbg_xi_tot]   ;
	double* wy = new double[nbg_eta_tot]  ; 
	double* wz = new double[nbg_zeta_tot] ; 
	double* mesxg = new double[nbg_xi_tot] ; 
	double* mesyg = new double[nbg_eta_tot] ; 
	double* meszg = new double[nbg_zeta_tot] ;   

	
	
    // Univariate basis functions at the 1D integration points 
	double** Nxig   = array_2d(nbg_xi_tot , nbf_elem_xi);
	double** Netag  = array_2d(nbg_eta_tot, nbf_elem_eta);
    double** Nzetag = array_2d(nbg_zeta_tot, nbf_elem_zeta); 
    //int spanx; 
    //int spany; 
    //int spanz; 
    //int nxi   = nbf_xi   - 1 ;  
	//int neta  = nbf_eta  - 1 ;
	//int nzeta = nbf_zeta - 1 ;
	
	// Info of the second mesh 
	int* Spanx2 = new int[nbg_xi_tot]   ; 
	int* Spany2 = new int[nbg_eta_tot]  ; 
	int* Spanz2 = new int[nbg_zeta_tot] ; 
	
	int nbf_xi2   = sknotXi2   -1 - deg_xi2   ; 
	int nbf_eta2  = sknotEta2  -1 - deg_eta2  ;  
	int nbf_zeta2 = sknotZeta2 -1 - deg_zeta2 ;
	int nbf2 = nbf_xi2*nbf_eta2*nbf_zeta2 ;  
	
	int nxi2   = nbf_xi2   - 1 ;  
	int neta2  = nbf_eta2  - 1 ;
	int nzeta2 = nbf_zeta2 - 1 ;
	
	int nbf_elem_xi2   = deg_xi2   + 1 ;
	int nbf_elem_eta2  = deg_eta2  + 1 ;
	int nbf_elem_zeta2 = deg_zeta2 + 1 ; 
	int nbf_elem2      = nbf_elem_xi2*nbf_elem_eta2*nbf_elem_zeta2;
	
	int ne_xi2   = nbf_xi2   - deg_xi2   ;
	int ne_eta2  = nbf_eta2  - deg_eta2  ; 
	int ne_zeta2 = nbf_zeta2 - deg_zeta2 ; 
	
	double** Nxig2   = array_2d(nbg_xi_tot  , nbf_elem_xi2);
	double** Netag2  = array_2d(nbg_eta_tot , nbf_elem_eta2);
    double** Nzetag2  = array_2d(nbg_zeta_tot, nbf_elem_zeta2); 
    
			
	
	double mes;
	int k=0;
	// Loop on 1D xi elements
	for ( int sxi = deg_xi ; sxi < sknotXi - deg_xi -1 ; sxi++ )
	{
		mes = knotXi[sxi+1] - knotXi[sxi];
		for ( int ii=0; ii < nbg_xi; ii++ )
		{
			xg[k]    = knotXi[sxi] + 0.5*(xig[ii]+1)*mes ;
			wx[k]    = wgxi[ii];
			mesxg[k] = mes;
			basisfuns(sxi, xg[k], deg_xi, knotXi, Nxig[k]);
			
			// Second mesh 
			Spanx2[k] = findspan(nxi2, deg_xi2, xg[k], knotXi2) ; 
			basisfuns(Spanx2[k], xg[k], deg_xi2, knotXi2, Nxig2[k]  ) ; 
			
			
			k++;
		}
	}
	k = 0;
	// Loop on 1D eta elements
	for ( int seta = deg_eta ; seta < sknotEta - deg_eta -1 ; seta++ )
	{
		mes = knotEta[seta+1] - knotEta[seta];
		for ( int ii=0; ii < nbg_eta; ii++ )
		{
			yg[k]    = knotEta[seta] + 0.5*(etag[ii]+1)*mes ;
			wy[k]    = wgeta[ii] ;
			mesyg[k] = mes ;
			basisfuns(seta, yg[k], deg_eta, knotEta, Netag[k]);
			
			// Second mesh 
			Spany2[k] = findspan(neta2, deg_eta2, yg[k], knotEta2) ; 
			basisfuns(Spany2[k], yg[k], deg_eta2, knotEta2, Netag2[k]  ) ;
			
						
			k++;
		}
	}
	k = 0;	
	// Loop on 1D eta elements
	for ( int szeta = deg_zeta ; szeta < sknotZeta - deg_zeta -1 ; szeta++ )
	{
		mes = knotZeta[szeta+1] - knotZeta[szeta];
		for ( int ii=0; ii < nbg_zeta; ii++ )
		{
			zg[k]    = knotZeta[szeta] + 0.5*(zetag[ii]+1)*mes ;
			wz[k]    = wgzeta[ii] ;
			meszg[k] = mes ;
			basisfuns(szeta, zg[k], deg_zeta, knotZeta, Nzetag[k]);
			
			// Second mesh 
			Spanz2[k] = findspan(nzeta2, deg_zeta2, zg[k], knotZeta2) ; 
			basisfuns(Spanz2[k], zg[k], deg_zeta2, knotZeta2, Nzetag2[k]  ) ;
			
			k++;
		}
	}

     
    double* N  = new double[nbf_elem]  ; 
    double* be = new double[3*nbf_elem]; // Elementary right hand side vector 
	// Elementary operators
	Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
    Eigen::MatrixXd Meb(nbf_elem,nbf_elem);  // Block of the elementary Mass operator  
    Eigen::MatrixXd Neb(nbf_elem,nbf_elem);       // Block elementary differential operator 
 
    double* N2 = new double[nbf_elem2] ; 
	Eigen::MatrixXi NOELEM2 = getNoelem(ne_xi2, ne_eta2, ne_zeta2, deg_xi2, deg_eta2, deg_zeta2 ); 
 
	
	for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
    double uxg ; // Displacement  ux of a Gauss point 
    double uyg ; // Displacement  uy of a Gauss point 
    double uzg ; // Displacement  uz of a Gauss point 
    double wgMes ; 

    
    // Assembly of the matrix 
    int I,J;
    int sp_count = 0 ;  
    int ecounter=0 ; // index of element in mesh 
    int ie ; // index of element in mesh 2 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {        
            for (int i=0; i<ne_xi; i++)
            {
                // Current element (i,j,k)
                Meb.fill(0) ;  
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
                
                
                // Loop on the integration points 
                for (int kk=k*nbg_zeta; kk<(k+1)*nbg_zeta; kk++)
                {
                    for (int jj=j*nbg_eta; jj<(j+1)*nbg_eta; jj++)
                    {
                        for (int ii=i*nbg_xi; ii<(i+1)*nbg_xi; ii++)
                        {
                            // Getting the trivariate basis functions on the Gauss integration point 
                            bfc = 0 ; 
                            for (int kbf=0; kbf < nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf < nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf < nbf_elem_xi; ibf++)
                                    {
                                        N[bfc]   = Nzetag[kk][kbf]*Netag[jj][jbf]*Nxig[ii][ibf] ;    
                                        bfc++; 
                                    }
                                }
                            }  
                        
                            // Setting the block elementary differential operator 
                            // N tensor N is the a block of the diagonal elemental mass matrix 
                            for ( int ibf =0; ibf < nbf_elem; ibf++)
            				{
                				for (int jbf=0; jbf < nbf_elem; jbf++)
                				{
        				             Neb(ibf,jbf)  = N[ibf]*N[jbf] ; 
                				}
            				}
            				wgMes = (wx[ii]*wy[jj]*wz[kk]*mesxg[ii]*mesyg[jj]*meszg[kk]/8) ; 
            				Meb   = Meb + Neb*wgMes ; 
            				
            				
            				
            				// Getting the trivariate basis functions of the projected mesh 
				            bfc = 0 ; 
                            for (int kbf=0; kbf < nbf_elem_zeta2; kbf++)
                            {
                                for (int jbf=0; jbf < nbf_elem_eta2; jbf++)
                                {
                                    for (int ibf=0; ibf < nbf_elem_xi2; ibf++)
                                    {
                                        N2[bfc] = Nzetag2[kk][kbf]*Netag2[jj][jbf]*Nxig2[ii][ibf] ;    
                                        bfc++ ; 
                                    }
                                }
                            }
                            
                            // Evaluating the displacement to be projected at the Gauss point    
                            ie = Spanx2[ii] -deg_xi2 + (Spany2[jj]-deg_eta2)*ne_xi2 + (Spanz2[kk]-deg_zeta2)*ne_xi2*ne_eta2 ;  
            				uxg = 0 ; 
            				uyg = 0 ;
            				uzg = 0 ; 
            			    
                            for (int ibf=0; ibf< nbf_elem2; ibf++)
                            {
                                uxg  +=  N2[ibf]*U2[NOELEM2(ie,ibf)]; 
                                uyg  +=  N2[ibf]*U2[NOELEM2(ie,ibf)+nbf2];
                                uzg  +=  N2[ibf]*U2[NOELEM2(ie,ibf)+2*nbf2];     
                            }
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            +=  wgMes*uxg*N[ibf] ; 
                                be[ibf+nbf_elem]   +=  wgMes*uyg*N[ibf] ; 
                                be[ibf+2*nbf_elem] +=  wgMes*uzg*N[ibf] ; 
                            }                                        				
                        }
                    }
                }
                
                // Elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM(ecounter,ibf)]       += be[ibf]; 
                    rhs[NOELEM(ecounter,ibf)+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM(ecounter,ibf)+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                // Elementary contribution to the global matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM(ecounter,ibf); 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = NOELEM(ecounter,jbf);
    					
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = Meb(ibf,jbf);
    		            sp_count ++ ;  	
    		            
		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = Meb(ibf,jbf);
    		            sp_count ++ ;
    		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  Meb(ibf,jbf);
    		            sp_count ++ ; 
                    }
                }
                ecounter++;                            
            }
        }
    }    
    
    delete[] xig   ;  
    delete[] wgxi  ; 
    delete[] etag  ; 
    delete[] wgeta ;
    delete[] zetag ; 
    delete[] wgzeta ;
    
    
    delete[] xg  ; 
    delete[] yg  ; 
    delete[] zg  ;
    delete[] wx  ;
	delete[] wy  ; 
	delete[] wz  ; 
	delete[] mesxg ; 
	delete[] mesyg ; 
	delete[] meszg ;   
	
    delete[] Spanx2 ; 
    delete[] Spany2 ; 
    delete[] Spanz2 ;
    
    delete[] be ; 
    delete[] N  ; 
    delete[] N2 ;
	 
	// Univariate basis functions and their derivatives at  1D gauss points
	delete_array_2d(Nxig)     ;
	delete_array_2d(Netag)    ;
	delete_array_2d(Nzetag)   ;
	delete_array_2d(Nxig2)     ;
	delete_array_2d(Netag2)    ;
	delete_array_2d(Nzetag2)   ;  
    
}
                          
                          
                          

std::vector<int> VoxelIntegrationThreshold(double* fip, int sfip,
                                                             double thrsh, 
                                                             int nipex, int nipey, int nipez, 
                                                             int deg_xi, int deg_eta, int deg_zeta, 
                                                             double* knotXi, int sknotXi, 
                                                             double* knotEta, int sknotEta, 
                                                             double* knotZeta, int sknotZeta )
{
    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
	std::vector<int> ipIndices;  
    std::vector<int> ipIndicesElem;  

    int vc=0; // voxel counter   
    int nv; // Number of taken voxels for an element 
    

    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {   
                // Current element 
                ipIndicesElem.clear();     
                
                // Segmentation 
                // Loop over the voxel integration points and select only those 
                // which are greater than the threshold value 
                nv = 0 ;               
                for (int kkg=k*nipez, kkl=0 ; kkg<(k+1)*nipez && kkl<nipez ; kkg++, kkl++ )
                {
                    for (int jjg=j*nipey, jjl=0 ; jjg<(j+1)*nipey && jjl<nipey  ; jjg++,jjl++ )
                    {
                        for (int iig=i*nipex, iil=0 ; iig<(i+1)*nipex && iil<nipex  ; iig++, iil++ ) 
                        {
                            if (fip[vc] > thrsh) 
                            {  
                                ipIndicesElem.push_back(iig); 
                                ipIndicesElem.push_back(jjg); 
                                ipIndicesElem.push_back(kkg);
                                ipIndicesElem.push_back(iil); 
                                ipIndicesElem.push_back(jjl); 
                                ipIndicesElem.push_back(kkl);                                  
                                nv++; 
                            } 
                            vc++; 
                        }
                    }
                }  
                // Output data structure (number of voxels of Element 1, i1,j1,k1, i2,j2,k2, .... )
                ipIndices.push_back(nv); 
                ipIndices.insert(ipIndices.end(), ipIndicesElem.begin(), ipIndicesElem.end() ); 	    
            }
        }
        
    }  
    return ipIndices;              
}


std::vector<int> VoxelIntegrationMask(double* maskip, int smaskip,
                                                             int nipex, int nipey, int nipez, 
                                                             int deg_xi, int deg_eta, int deg_zeta, 
                                                             double* knotXi, int sknotXi, 
                                                             double* knotEta, int sknotEta, 
                                                             double* knotZeta, int sknotZeta )
{
    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	
	std::vector<int> ipIndices;  
    std::vector<int> ipIndicesElem;  

    int vc=0; // voxel counter   
    int nv; // Number of taken voxels for an element 
    

    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {   
                // Current element 
                ipIndicesElem.clear();     
                
                // Segmentation 
                // Loop over the voxel integration points and select only those 
                // which belong to the mask 
                nv = 0 ;               
                for (int kkg=k*nipez, kkl=0 ; kkg<(k+1)*nipez && kkl<nipez ; kkg++, kkl++ )
                {
                    for (int jjg=j*nipey, jjl=0 ; jjg<(j+1)*nipey && jjl<nipey  ; jjg++,jjl++ )
                    {
                        for (int iig=i*nipex, iil=0 ; iig<(i+1)*nipex && iil<nipex  ; iig++, iil++ ) 
                        {
                            if (maskip[vc] ==1) 
                            {  
                                ipIndicesElem.push_back(iig); 
                                ipIndicesElem.push_back(jjg); 
                                ipIndicesElem.push_back(kkg);
                                ipIndicesElem.push_back(iil); 
                                ipIndicesElem.push_back(jjl); 
                                ipIndicesElem.push_back(kkl);                                  
                                nv++; 
                            } 
                            vc++; 
                        }
                    }
                }  
                // Output data structure (number of voxels of Element 1, i1,j1,k1, i2,j2,k2, .... )
                ipIndices.push_back(nv); 
                ipIndices.insert(ipIndices.end(), ipIndicesElem.begin(), ipIndicesElem.end() ); 	    
            }
        }
        
    }  
    return ipIndices;              
}



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
                                        double* nnz_values, int nnz)
{

    // Integration rule (rectangle uniform subdivision of each element)
 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;
	int ne = ne_xi*ne_eta*ne_zeta     ;    
 
 
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
 
    int nipeG = nipex*nipey*nipez ;  
 
 
    // Elementary operators
    Eigen::MatrixXd He(3*nbf_elem,3*nbf_elem);  
    Eigen::MatrixXd Ne(3,3*nbf_elem); 
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
    Ne.fill(0.);  
    Eigen::MatrixXd gradFgradFt(3,3) ;   // Tensor product of the gradient of gray-levels with itself 
    

    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    
    // Assembly of the Hessian matrix 
    int I,J;
    int Ig; 
    int sp_count = 0 ; 
    int ecounter=0 ; // index of element 
    int bfc ; 
    int nipe ;  // local number of integration points per element 
    int jmin ; 
    int jmax ; 
    int iig ; 
    int jjg ; 
    int kkg ;
    int iil ; 
    int jjl ; 
    int kkl ;  
    nipe = ipIndices[0] ; 
    if (nipe!=0)
    {
        jmin = 1 ; 
        jmax = jmin + 6*nipe; 
    } 
    else 
    {
        jmax = 1 ; 
    }   
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {        
            for (int i=0; i<ne_xi; i++)
            {
                // Current element (i,j,k)
                He.fill(0.); // Elementary hessian matrix initialization 
                
                // Loop over the integration points of the element 
                for (int ip=0; ip<nipe; ip++)
                {
                    iig = ipIndices[jmin+ip]   ;
                    jjg = ipIndices[jmin+ip+1] ;
                    kkg = ipIndices[jmin+ip+2] ;
                    iil = ipIndices[jmin+ip+3] ;
                    jjl = ipIndices[jmin+ip+4] ;
                    kkl = ipIndices[jmin+ip+5] ;
                    jmin+=5 ; 
                    
                    bfc = 0 ; 
                    // Getting the trivariate basis functions on the integration point 
                    for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                    {
                        for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                        {
                            for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                            {                                        
                                N[bfc] = Nzeta[kbf+kkg*sjNzeta]*Neta[jbf+jjg*sjNeta]*Nxi[ibf+iig*sjNxi] ; 
                                bfc++;  
                            }
                        }
                    } 
                    // Setting the elemental basis function matrix 
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
   					{
       					Ne(0,ibf)  = N[ibf] ; 
       					Ne(1,ibf+nbf_elem)   = N[ibf] ;  
       					Ne(2,ibf+2*nbf_elem) = N[ibf] ; 
   					}
					//Setting the gradient tensor product 
					Ig = iil +jjl*nipex + kkl*nipex*nipey + ecounter*nipeG ;  
                    gradFgradFt(0,0) = dfipdx[Ig]*dfipdx[Ig] ;     
                    gradFgradFt(0,1) = dfipdx[Ig]*dfipdy[Ig] ;
                    gradFgradFt(0,2) = dfipdx[Ig]*dfipdz[Ig] ;
                    gradFgradFt(1,0) = dfipdy[Ig]*dfipdx[Ig] ;
                    gradFgradFt(1,1) = dfipdy[Ig]*dfipdy[Ig] ;
                    gradFgradFt(1,2) = dfipdy[Ig]*dfipdz[Ig] ;
                    gradFgradFt(2,0) = dfipdz[Ig]*dfipdx[Ig] ;
                    gradFgradFt(2,1) = dfipdz[Ig]*dfipdy[Ig] ;
                    gradFgradFt(2,2) = dfipdz[Ig]*dfipdz[Ig] ;
                    // Adding to the voxel summation 
        			He = He + mes*(Ne.transpose()*gradFgradFt*Ne) ;                    
                }
                

                // Elementary contribution to the global Hessian matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM(ecounter,ibf); 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = NOELEM(ecounter,jbf);
    					
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf,jbf);
    		            sp_count ++ ;
    		            
    		            indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf+nbf_elem,jbf);
    		            sp_count ++ ;
    		            
		                indexI[sp_count] = I ;
    		            indexJ[sp_count] = J +nbf;
    		            nnz_values[sp_count] = He(ibf,jbf+nbf_elem);
    		            sp_count ++ ;
    		            
		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = He(ibf+nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ;
    		                   
    		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf+2*nbf_elem,jbf);
    		            sp_count ++ ;    		            
    		            
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf,jbf+2*nbf_elem);
    		            sp_count ++ ;
    		            
    		            indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf+nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ;   		            
 
     		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] =  He(ibf+2*nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ; 
    		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf+2*nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ; 
                    }
                }                
                if ( ecounter == ne-1 ) 
                {
                    break ; 
                } 
                nipe = ipIndices[jmax] ; 
                if (nipe==0)
                {
                    jmax += 1; 
                }
                else 
                {
                    jmin = jmax + 1 ; 
                    jmax = jmin + 6*nipe ; 
                }
                ecounter++;      					    
            }
        }
    }

    delete[] N; 
 
}   



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
                                   double* rhs, int ndof )
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	int ne = ne_xi*ne_eta*ne_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
    
    int nipeG = nipex*nipey*nipez ; 
    
    
 
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhi ; // g(x+u(x)) : advected gray-level 
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    double ssd = 0 ; // Gray-level sum of squared differences 
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
 

 
    // Assembly of the Hessian matrix 
    int Ig; 
    int ecounter=0 ; // index of element 
    int bfc ; 
    int nipe ;  // local number of integration points per element 
    int jmin; 
    int jmax; 
    int iig ; 
    int jjg ; 
    int kkg ;
    int iil ; 
    int jjl ; 
    int kkl ; 
    nipe = ipIndices[0] ; 
    if (nipe!=0)
    {
        jmin = 1 ; 
        jmax = jmin + 6*nipe; 
    } 
    else 
    {
        jmax = 1 ; 
    }   
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {
                            
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
               
                // Loop over the integration points (voxels) of the element
                for (int ip=0; ip<nipe; ip++)
                {
                    iig = ipIndices[jmin+ip]   ;
                    jjg = ipIndices[jmin+ip+1] ;
                    kkg = ipIndices[jmin+ip+2] ;
                    iil = ipIndices[jmin+ip+3] ;
                    jjl = ipIndices[jmin+ip+4] ;
                    kkl = ipIndices[jmin+ip+5] ;
                    jmin+=5 ;
                    
                    
                    //Getting trivariate basis functions at the integration point
                    bfc = 0 ;   
                    for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                    {
                        for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                        {
                            for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                            {
                                N[bfc] = Nzeta[kbf+kkg*sjNzeta]*Neta[jbf+jjg*sjNeta]*Nxi[ibf+iig*sjNxi] ;  
                                bfc++ ;  
                            }
                        }
                    }
                    
                    // Computing the displacement at the integration points 
                    // by the B-spline linear combination 
                    uxg = 0 ; 
                    uyg = 0 ; 
                    uzg = 0 ; 
                    for (int ibf=0; ibf< nbf_elem; ibf++)
                    {
                        
                        uxg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                        uyg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                        uzg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                    }
                    // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                    goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                            knotXiImage, sknotXiImage,   
                            knotEtaImage, sknotEtaImage,   
                            knotZetaImage, sknotZetaImage,   
                            xg[iig]+uxg, yg[jjg]+uyg ,zg[kkg]+uzg ); 
                    
                    Ig = iil +jjl*nipex + kkl*nipex*nipey + ecounter*nipeG ;
                    res = fip[Ig] - goPhi ;          
                    wgResdfdx =  mes*res*dfipdx[Ig] ;
                    wgResdfdy =  mes*res*dfipdy[Ig] ; 
                    wgResdfdz =  mes*res*dfipdz[Ig] ; 
                    
                    ssd += std::pow(res, 2) ; 
                    
                    // Setting the elemental right hand side vector 
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        be[ibf]            +=  wgResdfdx*N[ibf] ; 
                        be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                        be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ; 
                    }                     
                } 
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                if ( ecounter == ne-1 ) 
                {
                    break ; 
                } 
                nipe = ipIndices[jmax] ; 
                if (nipe==0)
                {
                    jmax += 1; 
                }
                else 
                {
                    jmin = jmax + 1 ; 
                    jmax = jmin + 6*nipe ; 
                }
                ecounter++;     
            }
        }
    }
 
    delete[] N  ; 
    delete[] be ; 
    
    return ssd ; 
}   


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
                                   double* rhs, int ndof )
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ; 
	int ne = ne_xi*ne_eta*ne_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
    
    int nipeG = nipex*nipey*nipez ; 
 
    
    // Elementary operations 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhiMean ; // Mean of g(x+u(x)) over the element 
    double goPhiStd  ; // Standard deviation of g(x+u(x)) over the element  
    double goPhiSquare ; // Square of g(x+u(x))
    double* goPhi = new double[nipeG] ; // g(x+u(x)) : advected gray-level is saved at the voxels of the element so that the mean and std will be computed  
    double** Nelem = array_2d(nipeG,nbf_elem); // Saving for the current element all the basis functions evaluated at the voxels  
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    double znssd = 0 ; // gray-level znssd 
    
    
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
 

 
    // Assembly of the Hessian matrix 
    int Ig; 
    int ecounter=0 ; // index of element 
    int bfc ; 
    int nipe ;  // local number of integration points per element 
    int jmin; 
    int jminOld; 
    int jmax; 
    int iig ; 
    int jjg ; 
    int kkg ;
    int iil ; 
    int jjl ; 
    int kkl ; 
    nipe = ipIndices[0] ; 
    if (nipe!=0)
    {
        jmin = 1 ; 
        jmax = jmin + 6*nipe; 
    } 
    else 
    {
        jmax = 1 ; 
    }   
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {
                            
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                

                // First evaluating goPhi and computing elemental mean and std of goPhi
                // Also saving the basis functions of the element 
                goPhiMean = 0 ; 
                goPhiSquare  = 0 ; 
                jminOld = jmin ; 
                for (int ip=0; ip<nipe; ip++)
                {
                
                    iig = ipIndices[jmin+ip]   ;
                    jjg = ipIndices[jmin+ip+1] ;
                    kkg = ipIndices[jmin+ip+2] ;
                    iil = ipIndices[jmin+ip+3] ;
                    jjl = ipIndices[jmin+ip+4] ;
                    kkl = ipIndices[jmin+ip+5] ;
                    jmin+=5 ;
                    
                    
                    //Getting trivariate basis functions at the integration point
                    bfc = 0 ;   
                    for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                    {
                        for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                        {
                            for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                            {
                                Nelem[ip][bfc] = Nzeta[kbf+kkg*sjNzeta]*Neta[jbf+jjg*sjNeta]*Nxi[ibf+iig*sjNxi] ;  
                                bfc++ ;  
                            }
                        }
                    } 
                    // Computing the displacement at the integration points 
                    // by the B-spline linear combination 
                    uxg = 0 ; 
                    uyg = 0 ; 
                    uzg = 0 ; 
                    for (int ibf=0; ibf< nbf_elem; ibf++)
                    {
                        uxg  +=  Nelem[ip][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                        uyg  +=  Nelem[ip][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                        uzg  +=  Nelem[ip][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                    }
                    // Saving goPhi at the current voxel = integration point 
                    goPhi[ip] = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                 knotXiImage, sknotXiImage, 
                                 knotEtaImage, sknotEtaImage,
                                 knotZetaImage, sknotZetaImage, 
                                 xg[iig]+uxg, yg[jjg]+uyg ,zg[kkg]+uzg );  
                    
                    goPhiMean = goPhiMean + goPhi[ip] ;  
                    goPhiSquare  = goPhiSquare  + std::pow(goPhi[ip],2) ;    
                }
                // Getting the mean and standard deviation of g(x+u(x)) over the element 
                goPhiMean = goPhiMean/nipe ;  
                goPhiStd  = std::sqrt( goPhiSquare/nipe - std::pow(goPhiMean,2) ) ; 
                // Loop again over the integration points for the assembly 
                jmin = jminOld ; 
                for (int ip=0; ip<nipe; ip++)
                {
                    iig = ipIndices[jmin+ip]   ;
                    jjg = ipIndices[jmin+ip+1] ;
                    kkg = ipIndices[jmin+ip+2] ;
                    iil = ipIndices[jmin+ip+3] ;
                    jjl = ipIndices[jmin+ip+4] ;
                    kkl = ipIndices[jmin+ip+5] ;
                    jmin+=5 ;
                    Ig = iil +jjl*nipex + kkl*nipex*nipey + ecounter*nipeG ; 
                    res = fip[Ig] - fmeane[ecounter] - (fstde[ecounter]/goPhiStd)*(goPhi[ip]-goPhiMean) ;  
                    wgResdfdx =  mes*res*dfipdx[Ig] ;
                    wgResdfdy =  mes*res*dfipdy[Ig] ; 
                    wgResdfdz =  mes*res*dfipdz[Ig] ; 
                    
                    znssd += std::pow(  (fip[Ig] - fmeane[ecounter])/fstde[ecounter]  -  (goPhi[ip]-goPhiMean)/goPhiStd   ,2) ; 
                    
                    // Setting the elemental right hand side vector 
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        be[ibf]            += wgResdfdx*Nelem[ip][ibf] ; 
                        be[ibf+nbf_elem]   += wgResdfdy*Nelem[ip][ibf] ; 
                        be[ibf+2*nbf_elem] += wgResdfdz*Nelem[ip][ibf] ; 
                    } 
                }      
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem ;  ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf] ; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem] ; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem] ; 
                }
                if ( ecounter == ne-1 ) 
                {
                    break ; 
                } 
                nipe = ipIndices[jmax] ; 
                if (nipe==0)
                {
                    jmax += 1; 
                }
                else 
                {
                    jmin = jmax + 1 ; 
                    jmax = jmin + 6*nipe ; 
                }
                ecounter++;     
            }
        }
    }
 
    delete[] goPhi ; 

    delete[] be ;

	delete_array_2d(Nelem); 
	
	return znssd ; 
    
}    


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
                                        double* nnz_values, int nnz) 
{    
    // Integration rule (rectangle uniform subdivision of each element)
  
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
 	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
 	int nbf_zeta = sknotZeta -1 - deg_zeta ;
 	
 	int nbf_elem_xi   = deg_xi  +1 ;
 	int nbf_elem_eta  = deg_eta +1 ;
 	int nbf_elem_zeta = deg_zeta+1 ; 
 	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
 	
 	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
 	
 	int ne_xi   = nbf_xi   - deg_xi   ;
 	int ne_eta  = nbf_eta  - deg_eta  ; 
 	int ne_zeta = nbf_zeta - deg_zeta ;  
  
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
  
  
    // Elementary operators
    Eigen::MatrixXd He(3*nbf_elem,3*nbf_elem);  
    Eigen::MatrixXd Ne(3,3*nbf_elem); 
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 
    Ne.fill(0.);  
    Eigen::MatrixXd gradFgradFt(3,3);// Tensor product of the gradient of gray-levels with itself 
    

    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    
    // Assembly of the Hessian matrix 
    int I,J;
    int sp_count = 0 ; 
    int ig=0 ; // index of the integration points  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {        
            for (int i=0; i<ne_xi; i++)
            {
                // Current element (i,j,k)
                He.fill(0.); // Elementary hessian matrix initialization 
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                
                //EvaluateTrilinearInterpolationAndGradientStructured(f,sif,sjf,skf,
                //knotXiImage, sknotXiImage, knotEtaImage, sknotEtaImage, knotZetaImage, sknotZetaImage, 
                //xge, nipex, yge, nipey, zge, nipez, fe, dfdxe, dfdye, dfdze);  
                
                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                            bfc = 0 ; 
                            //Getting trivariate basis functions on the gauss integration point 
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        // ibf + jbf*nbf_elem_xi + kbf*nbf_elem_xi*nbf_elem_eta
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ; 
                                        bfc++;  
                                    }
                                }
                            }
                            // Setting the elemental basis function matrix 
                            for ( int ibf =0; ibf < nbf_elem; ibf++)
         					{
             					Ne(0,ibf)  = N[ibf] ; 
             					Ne(1,ibf+nbf_elem)   = N[ibf] ;  
             					Ne(2,ibf+2*nbf_elem) = N[ibf] ; 
         					}
 					        //Setting the gradient tensor product 
                            gradFgradFt(0,0) = dfipdx[ig]*dfipdx[ig] ; 
                            gradFgradFt(0,1) = dfipdx[ig]*dfipdy[ig] ;
                            gradFgradFt(0,2) = dfipdx[ig]*dfipdz[ig] ;
                            gradFgradFt(1,0) = dfipdy[ig]*dfipdx[ig] ;
                            gradFgradFt(1,1) = dfipdy[ig]*dfipdy[ig] ;
                            gradFgradFt(1,2) = dfipdy[ig]*dfipdz[ig] ;
                            gradFgradFt(2,0) = dfipdz[ig]*dfipdx[ig] ;
                            gradFgradFt(2,1) = dfipdz[ig]*dfipdy[ig] ;
                            gradFgradFt(2,2) = dfipdz[ig]*dfipdz[ig] ;
                            // Adding to the voxel summation 
         					He = He + mes*(Ne.transpose()*gradFgradFt*Ne) ; 
                            ig++; 
                        }  
                    }        
                }
                // Elementary contribution to the global Hessian matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = NOELEM(ecounter,ibf); 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
     					J = NOELEM(ecounter,jbf);
     					
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf,jbf);
    		            sp_count ++ ;
    		            
    		            indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf+nbf_elem,jbf);
    		            sp_count ++ ;
    		            
		                indexI[sp_count] = I ;
    		            indexJ[sp_count] = J +nbf;
    		            nnz_values[sp_count] = He(ibf,jbf+nbf_elem);
    		            sp_count ++ ;
    		            
 		            	indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] = He(ibf+nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ;
     		                   
    		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J ;
    		            nnz_values[sp_count] = He(ibf+2*nbf_elem,jbf);
    		            sp_count ++ ;    		            
    		            
    		            indexI[sp_count] = I ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf,jbf+2*nbf_elem);
    		            sp_count ++ ;
    		            
    		            indexI[sp_count] = I+nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf+nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ;   		            
  
      		            indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+nbf ;
    		            nnz_values[sp_count] =  He(ibf+2*nbf_elem,jbf+nbf_elem);
    		            sp_count ++ ; 
     		               		            
		                indexI[sp_count] = I+2*nbf ;
    		            indexJ[sp_count] = J+2*nbf ;
    		            nnz_values[sp_count] =  He(ibf+2*nbf_elem,jbf+2*nbf_elem);
    		            sp_count ++ ; 
                    }
                }
                ecounter++;      					    
            }
        }
    }

    delete[] N; 
  
}                                                                             
                                        
                                   
 
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
                                        double* nnz_values, int nnz) 
{    
    // Integration rule (rectangle uniform subdivision of each element)
 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
 
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
 
    int nipeG = nipex*nipey*nipez ;  
 
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta ); 

    
 
   
     
//    #pragma omp parallel   
//     default(shared)  
//     private(i, ecounter, He, Ne, gradFgradFt, Ne, N, ig, Ig, kk, jj, ii, bfc, kbf, jbf, ibf, I, J ) 
//    #pragma omp parallel default(shared) private()

    // Assembly of the Hessian matrix 
    
    #pragma omp parallel 
    {
        #pragma omp for   
            for (int k=0; k<ne_zeta ;k++)
            {   
                //#pragma omp parallel shared(k,ne_zeta)
                {
                    //#pragma omp for   
                        for (int j=0; j<ne_eta; j++)
                        {  
                            // #pragma omp parallel 
                            {
                                // #pragma omp for
                                    for (int i=0; i<ne_xi; i++)
                                    {
                                        int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ;    
                                        Eigen::MatrixXd He(3*nbf_elem,3*nbf_elem);
                                        Eigen::MatrixXd Ne(3,3*nbf_elem); 
                                        Eigen::MatrixXd gradFgradFt(3,3);
                                        double* N = new double[nbf_elem] ;
                                        
                                        // Current element (i,j,k)
                                        Ne.fill(0.);  
                                        He.fill(0.); // Elementary hessian matrix initialization 
                                      
                                        // Loop over the integration points (voxels) of the element
                                        int ig = 0 ;
                                        int Ig ; 
                                        int bfc ; 
                                        for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                                        {
                                            for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                                            {
                                                for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                                                {
                                                    bfc = 0 ; 
                                                    //Getting trivariate basis functions on the gauss integration point 
                                                    for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                                                    {
                                                        for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                                        {
                                                            for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                                            {
                                                                // ibf + jbf*nbf_elem_xi + kbf*nbf_elem_xi*nbf_elem_eta
                                                                N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ; 
                                                                bfc++;  
                                                            }
                                                        }
                                                    }
                                                    // Setting the elemental basis function matrix 
                                                    for ( int ibf =0; ibf < nbf_elem; ibf++)
                                					{
                                    					Ne(0,ibf)  = N[ibf] ; 
                                    					Ne(1,ibf+nbf_elem)   = N[ibf] ;  
                                    					Ne(2,ibf+2*nbf_elem) = N[ibf] ; 
                                					}
                        					        //Setting the gradient tensor product 
                        					        Ig = ig + ecounter*nipeG ; 
                                                    gradFgradFt(0,0) = dfipdx[Ig]*dfipdx[Ig] ; 
                                                    gradFgradFt(0,1) = dfipdx[Ig]*dfipdy[Ig] ;
                                                    gradFgradFt(0,2) = dfipdx[Ig]*dfipdz[Ig] ;
                                                    gradFgradFt(1,0) = dfipdy[Ig]*dfipdx[Ig] ;
                                                    gradFgradFt(1,1) = dfipdy[Ig]*dfipdy[Ig] ; 
                                                    gradFgradFt(1,2) = dfipdy[Ig]*dfipdz[Ig] ;
                                                    gradFgradFt(2,0) = dfipdz[Ig]*dfipdx[Ig] ;
                                                    gradFgradFt(2,1) = dfipdz[Ig]*dfipdy[Ig] ;
                                                    gradFgradFt(2,2) = dfipdz[Ig]*dfipdz[Ig] ;
                                                    // Adding to the voxel summation 
                                					He = He + mes*(Ne.transpose()*gradFgradFt*Ne) ; 
                                					ig++ ; 
                                                }  
                                            }        
                                        }
                                        int I,J ; 
                                        int sg = 9*nbf_elem*nbf_elem*ecounter  ; // global indexing in the total sparse value vector  
                                        // Elementary contribution to the global Hessian matrix 
                                        for (int ibf=0; ibf<nbf_elem; ibf++)
                                        {
                                            I = NOELEM(ecounter,ibf); 
                                            for (int jbf=0; jbf<nbf_elem; jbf++)
                                            {
                            					J = NOELEM(ecounter,jbf);
                            					
                            		            indexI[sg] = I ;
                            		            indexJ[sg] = J ;
                            		            nnz_values[sg] = He(ibf,jbf);
                            		            sg ++ ;
                            		            
                            		            indexI[sg] = I+nbf ;
                            		            indexJ[sg] = J ;
                            		            nnz_values[sg] = He(ibf+nbf_elem,jbf);
                            		            sg ++ ;
                            		            
                        		                indexI[sg] = I ;
                            		            indexJ[sg] = J +nbf;
                            		            nnz_values[sg] = He(ibf,jbf+nbf_elem);
                            		            sg ++ ;
                            		            
                        		            	indexI[sg] = I+nbf ;
                            		            indexJ[sg] = J+nbf ;
                            		            nnz_values[sg] = He(ibf+nbf_elem,jbf+nbf_elem);
                            		            sg ++ ;
                            		                   
                            		            indexI[sg] = I+2*nbf ;
                            		            indexJ[sg] = J ;
                            		            nnz_values[sg] = He(ibf+2*nbf_elem,jbf);
                            		            sg ++ ;    		            
                            		            
                            		            indexI[sg] = I ;
                            		            indexJ[sg] = J+2*nbf ;
                            		            nnz_values[sg] =  He(ibf,jbf+2*nbf_elem);
                            		            sg ++ ;
                            		            
                            		            indexI[sg] = I+nbf ;
                            		            indexJ[sg] = J+2*nbf ;
                            		            nnz_values[sg] =  He(ibf+nbf_elem,jbf+2*nbf_elem);
                            		            sg ++ ;   		            
                         
                             		            indexI[sg] = I+2*nbf ;
                            		            indexJ[sg] = J+nbf ;
                            		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+nbf_elem);
                            		            sg ++ ; 
                            		               		            
                        		                indexI[sg] = I+2*nbf ;
                            		            indexJ[sg] = J+2*nbf ;
                            		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+2*nbf_elem);
                            		            sg ++ ; 
                                            }
                                        }   
                                        delete[] N;   					    
                                    }
                            }
                        }
                 }
            }  
    }                
}    
                                        
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
                                        double* rhs, int ndof ) 
{ 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
 
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhi ; // g(x+u(x)) : advected gray-level 
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    double ssd=0; 
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
 
    int ig=0 ; // index of the integration points  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
               
                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                                uyg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                                uzg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                    knotXiImage, sknotXiImage, 
                                    knotEtaImage, sknotEtaImage,
                                    knotZetaImage, sknotZetaImage, 
                                    xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                                    
                            res = fip[ig] - goPhi ; 
                            wgResdfdx =  mes*res*dfipdx[ig] ;
                            wgResdfdy =  mes*res*dfipdy[ig] ; 
                            wgResdfdz =  mes*res*dfipdz[ig] ; 
                            
                            ssd += std::pow(res, 2) ; 
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            +=  wgResdfdx*N[ibf] ; 
                                be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                                be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ; 
                            } 
                            ig++; 
                        }  
                    }        
                }
                
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                ecounter++;        
            }
        }
    }
 
    delete[] N  ; 
    delete[] be ; 
    
    return ssd*mes ;  
}   

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
                                        double* rhs, int ndof ) 
{ 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
 
    // Elementary operations 
    
    int nipeG = nipex*nipey*nipez ;  
    
    
    double ssd=0; 
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }

    #pragma omp parallel 
    {
        #pragma omp for reduction(+:ssd)
            for (int k=0; k<ne_zeta ;k++)
            {   
                for (int j=0; j<ne_eta; j++)
                {             
                    for (int i=0; i<ne_xi; i++)
                    {             
                        // Current element (i,j,k)
                        double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
                        double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
                        double uxg ; // Displacement  ux of the voxel 
                        double uyg ; // Displacement  uy of the voxel
                        double uzg ; // Displacement  uz of the voxel
                        double res ; // f-goPhi at the integration point 
                        double goPhi ; // g(x+u(x)) : advected gray-level 
                        double wgResdfdx;  
                        double wgResdfdy;
                        double wgResdfdz;
                        int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ;    
                        double Nzeta_x_Neta ; 
                        
                        // Setting elementary right hand side vector to 0
                        for (int idofe=0; idofe<3*nbf_elem; idofe++)
                        {
                            be[idofe] = 0 ; 
                        } 
                        
                       
                        // Loop over the integration points (voxels) of the element
                        int ig = 0 ; 
                        int Ig ; 
                        int bfc ; 
                        for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                        {
                            for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                            {
                                for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                                {
                                
                                    //Getting trivariate basis functions at the integration point
                                    bfc = 0 ;   
                                    for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                                    {
                                        for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                        {
                                            Nzeta_x_Neta = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta] ; 
                                            for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                            {
                                                N[bfc] = Nzeta_x_Neta*Nxi[ibf+ii*sjNxi] ;  
                                                bfc++ ;  
                                            }
                                        }
                                    }
                                    // Computing the displacement at the integration points 
                                    // by the B-spline linear combination 
                                    uxg = 0 ; 
                                    uyg = 0 ; 
                                    uzg = 0 ; 
                                    for (int ibf=0; ibf< nbf_elem; ibf++)
                                    {
                                        uxg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                                        uyg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                                        uzg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                                    }
                                    // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                                    goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                            knotXiImage, sknotXiImage, 
                                            knotEtaImage, sknotEtaImage,
                                            knotZetaImage, sknotZetaImage, 
                                            xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                                    
                                    // Getting the gray-level index 
                                    Ig = ig + ecounter*nipeG ; 
                                            
                                    res = fip[Ig] - goPhi ; 
                                    wgResdfdx =  mes*res*dfipdx[Ig] ;
                                    wgResdfdy =  mes*res*dfipdy[Ig] ; 
                                    wgResdfdz =  mes*res*dfipdz[Ig] ; 
                                    
                                    ssd += std::pow(res, 2) ; 
                                    
                                    // Setting the elemental right hand side vector 
                                    for (int ibf=0; ibf<nbf_elem; ibf++)
                                    {
                                        be[ibf]            +=  wgResdfdx*N[ibf] ; 
                                        be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                                        be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ; 
                                    } 
                                    ig++; 
                                }  
                            }        
                        }
                        #pragma omp critical 
                        {
                            // Adding the elementary contribution to the right hand side 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                                rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                                rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                            }
                        }
                        delete[] N  ; 
                        delete[] be ;      
                    }
                }
            }        
    }            
    return ssd*mes ;  
}   

 
void DVC_RHS_ZN_Structured_TrilinearInterp(double* g, int sgi, int sgj,int sgk,
                                        double* knotXiImage, int sknotXiImage, 
                                        double* knotEtaImage, int sknotEtaImage, 
                                        double* knotZetaImage, int sknotZetaImage, 
                                        double* fip,    int sfip, 
                                        double* dfipdx, int sdfipdx, 
                                        double* dfipdy, int sdfipdy, 
                                        double* dfipdz, int sdfipdz,
                                        double* fmeane, int sfmeane, 
                                        double* fstde,  int sfstde , 
                                        double* dyne,   int sdyne,     
                                        int deg_xi, int deg_eta, int deg_zeta,
                                        double* knotXi, int sknotXi, 
                                        double* knotEta, int sknotEta, 
                                        double* knotZeta, int sknotZeta,
                                        int* NOELEM, int siNOELEM, int sjNOELEM,  
                                        int nipex, int nipey, int nipez,
                                        double* xg, int sxg, 
                                        double* yg, int syg, 
                                        double* zg, int szg, 
                                        double* Nxi, int siNxi, int sjNxi, 
                                        double* Neta, int siNeta, int sjNeta, 
                                        double* Nzeta, int siNzeta, int sjNzeta, 
                                        double* U, int sU, 
                                        double* rhs, int ndof, 
                                        double* elementRes, int selementRes )
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;   
	
    int nipe = nipex*nipey*nipez ; 
    
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
    
    // Elementary operations 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhiMean ; // Mean of g(x+u(x)) over the element 
    double goPhiStd  ; // Standard deviation of g(x+u(x)) over the element  
    double goPhiSquare ; // Square of g(x+u(x))
    double* goPhi = new double[nipe] ; // g(x+u(x)) : advected gray-level is saved at the voxels of the element so that the mean and std will be computed  
    double** Nelem = array_2d(nipe,nbf_elem); // Saving for the current element all the basis functions evaluated at the voxels  
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    
    double znssde = 0 ; // Elemental znssd 

    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
 
    int ig=0 ; // index of the integration points  
    int ige; // Index of the element integration point  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                // Firs evaluating goPhi and computing elemental mean and std of goPhi
                // Also saving the basis functions of the element 
                ige=0; // integration point (or voxel) counter 
                goPhiMean = 0 ; 
                goPhiSquare  = 0 ; 
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ; // Basis function counter 
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        // ibf + jbf*nbf_elem_xi + kbf*nbf_elem_xi*nbf_elem_eta
                                        Nelem[ige][bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  Nelem[ige][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                                uyg  +=  Nelem[ige][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                                uzg  +=  Nelem[ige][ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                            }
                            // Saving goPhi at the current voxel = integration point 
                            goPhi[ige] = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                         knotXiImage, sknotXiImage, 
                                         knotEtaImage, sknotEtaImage,
                                         knotZetaImage, sknotZetaImage, 
                                         xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                            
                            goPhiMean = goPhiMean + goPhi[ige] ;  
                            goPhiSquare  = goPhiSquare  + std::pow(goPhi[ige],2) ;  
                            ige++ ;                             
                            
                        }
                    }
                }
                // Getting the mean and standard deviation of g(x+u(x)) over the element 
                goPhiMean = goPhiMean/nipe ;  
                goPhiStd  = std::sqrt( goPhiSquare/nipe - std::pow(goPhiMean,2) ); 
             
                // Loop again  over the integration points (voxels) for the assembly 
                ige = 0 ; 
                znssde = 0; 
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {   
                            res = fip[ig] - fmeane[ecounter] - (fstde[ecounter]/goPhiStd)*(goPhi[ige]-goPhiMean) ;    
                            wgResdfdx =  mes*res*dfipdx[ig] ;
                            wgResdfdy =  mes*res*dfipdy[ig] ; 
                            wgResdfdz =  mes*res*dfipdz[ig] ; 
                            
                            znssde += std::pow(  (fip[ig] - fmeane[ecounter])/fstde[ecounter]  -  (goPhi[ige]-goPhiMean)/goPhiStd   ,2) ; 
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            += wgResdfdx*Nelem[ige][ibf] ; 
                                be[ibf+nbf_elem]   += wgResdfdy*Nelem[ige][ibf] ; 
                                be[ibf+2*nbf_elem] += wgResdfdz*Nelem[ige][ibf] ; 
                            } 
                            ig++; 
                            ige++; 
                        }  
                    }        
                }
                
                // Setting the elementary residual 
                elementRes[ecounter] = znssde ;  
              
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                 ecounter++;    
            }
        }
    }

    delete[] goPhi ; 

    delete[] be ;

	delete_array_2d(Nelem); 
 
	 
}
 

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
                                     double* glr, int sglr ) 
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipx  = nipex*ne_xi  ; 
    int nipy  = nipey*ne_eta ; 
    int nipxy =  nipx * nipy ; 
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
    
    
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double goPhi ; // g(x+u(x)) : advected gray-level 
    
    int ig=0; 
    int ecounter=0; 
    int bfc; 
 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)

                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                    knotXiImage, sknotXiImage, 
                                    knotEtaImage, sknotEtaImage,
                                    knotZetaImage, sknotZetaImage, 
                                    xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                                    
                            glr[ii+jj*nipx+kk*nipxy] = fip[ig] - goPhi ; 
 
                            ig++; 
                        }  
                    }        
                }
                ecounter++;      
            }
        }
    }
    
    delete[] N;    
} 

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
                                     double* glr, int sglr ) 
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipx  = nipex*ne_xi  ; 
    int nipy  = nipey*ne_eta ; 
    int nipxy =  nipx * nipy ; 
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
 
    int nipeG = nipex*nipey*nipez ;  
   
 
    #pragma omp parallel 
    {
        for (int k=0; k<ne_zeta ;k++)
        {   
            for (int j=0; j<ne_eta; j++)
            {             
                for (int i=0; i<ne_xi; i++)
                {
                    // Current element (i,j,k)
                    
                    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
                    double uxg ; // Displacement  ux of the voxel 
                    double uyg ; // Displacement  uy of the voxel
                    double uzg ; // Displacement  uz of the voxel
                    double goPhi ; // g(x+u(x)) : advected gray-level 
            
                    int ecounter = i + j*ne_xi + k*ne_xi*ne_eta ;    
                    int ig = 0 ; 
                    int Ig ; 
                    int bfc ; 
                    // Loop over the integration points (voxels) of the element
                    
                    for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                    {
                        for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                        {
                            for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                            {
                                //Getting trivariate basis functions at the integration point
                                bfc = 0 ;   
                                for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                                {
                                    for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                    {
                                        for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                        {
                                            N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                            bfc++ ;  
                                        }
                                    }
                                }
                                // Computing the displacement at the integration points 
                                // by the B-spline linear combination 
                                uxg = 0 ; 
                                uyg = 0 ; 
                                uzg = 0 ; 
                                for (int ibf=0; ibf< nbf_elem; ibf++)
                                {
                                    uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                    uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                    uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];     
                                }
                                // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                                goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                        knotXiImage, sknotXiImage, 
                                        knotEtaImage, sknotEtaImage,
                                        knotZetaImage, sknotZetaImage, 
                                        xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                                        
                                Ig = ig + ecounter*nipeG ;         
                                glr[ii+jj*nipx+kk*nipxy] = fip[Ig] - goPhi ; 
     
                                ig++; 
                            }  
                        }        
                    }
                    delete[] N;      
                }
            }
        }   
    }
 
} 


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
                                        double* glr, int sglr ) 
{

    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipe = nipex*nipey*nipez ; 


    int nipx  = nipex*ne_xi  ; 
    int nipy  = nipey*ne_eta ; 
    int nipxy =  nipx * nipy ; 
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
    
    
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel

    double goPhiMean ; // Mean of g(x+u(x)) over the element 
    double goPhiStd  ; // Standard deviation of g(x+u(x)) over the element  
    double goPhiSquare ; // Square of g(x+u(x))
    double* goPhi = new double[nipe] ; // g(x+u(x)) : advected gray-level is saved at the voxels of the element so that the mean and std will be computed  
    
    int ige; 
    int ig=0; 
    int ecounter=0; 
    int bfc; 
    
    
 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)

                // Firs evaluating goPhi and computing elemental mean and std of goPhi
                // Also saving the basis functions of the element 
                ige=0; // integration point (or voxel) counter 
                goPhiMean = 0 ; 
                goPhiSquare  = 0 ; 
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                         //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];     
                            }
                            // Saving goPhi at the current voxel = integration point 
                            goPhi[ige] = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                         knotXiImage, sknotXiImage, 
                                         knotEtaImage, sknotEtaImage,
                                         knotZetaImage, sknotZetaImage, 
                                         xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                            
                            goPhiMean    = goPhiMean + goPhi[ige] ;  
                            goPhiSquare  = goPhiSquare  + std::pow(goPhi[ige],2) ;  
                            ige++ ;                             
                            
                        }
                    }
                }
                // Getting the mean and standard deviation of g(x+u(x)) over the element 
                goPhiMean = goPhiMean/nipe ;  
                goPhiStd  = std::sqrt( goPhiSquare/nipe - std::pow(goPhiMean,2) ); 
             
                // Loop again  over the integration points (voxels) for the assembly 
                ige = 0 ; 
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {   
                            glr[ii+jj*nipx+kk*nipxy] = (fip[ig] - fmeane[ecounter])/fstde[ecounter]  -  (goPhi[ige]-goPhiMean)/goPhiStd ; 
                            ig++; 
                            ige++; 
                        }  
                    }        
                }
                 ecounter++;    
            }
        }
    }

    delete[] goPhi ; 

    delete[] N ;
	 
}
                                        
                                          
 
					                
					                     
 

/*
std::vector<double> VoxelIntegrationThresholdTrilinearInterp(double* image, int si, int sj, int sk, 
                                                    double* knotXiImage, int sknotXiImage, 
                                                    double* knotEtaImage, int sknotEtaImage, 
                                                    double* knotZetaImage, int sknotZetaImage, 
                                                    double xminIm, double xmaxIm, double yminIm, double ymaxIm, double zminIm, double zmaxIm,
                                                    double dx, double dy, double dz, 
                                                    double thrsh, 
                                                    int nipex, int nipey, int nipez, 
                                                    double* xg, int sxg, 
                                                    double* yg, int syg, 
                                                    double* zg, int szg, 
                                                    int deg_xi, int deg_eta, int deg_zeta, 
                                                    double* knotXi, int sknotXi, 
                                                    double* knotEta, int sknotEta, 
                                                    double* knotZeta, int sknotZeta                                                    
                                                    ) 
{    
    // Mesh parameters 

    //Total 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
    
    int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
    
    int nip_xi   =  ne_xi*nipex   ; 
    int nip_eta  =  ne_eta*nipey  ; 
    int nip_zeta =  ne_zeta*nipez ;  
    int nipe = nipex*nipey*nipez ; 
 
	
	double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex         ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey     ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
	
 
    double* xg = new double[nip_xi]   ; 
    double* yg = new double[nip_eta]  ; 
    double* zg = new double[nip_zeta] ; 
    
    // Coordinates of 1d integration points on one element
    double* xge = new double[nipex]; 
    double* yge = new double[nipey];
    double* zge = new double[nipez];
    
    double* fe    = new double[nipe] ; 

 
    std::vector<double> ipIndices;  
    std::vector<double> ipIndicesElem;  

    int vc; // voxel counter   
    int nv; // Number of taken voxels 
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
                
                // Current element 
                
                
                // Evaluate image gradient on each voxel of element using tensor product 1d-Bsplines 
                EvaluateTrilinearInterpolationStructured(image, si,sj,sk,
                knotXiImage, sknotXiImage, knotEtaImage, sknotEtaImage, knotZetaImage, sknotZetaImage, 
                xge, nipex, yge, nipey, zge, nipez, fe ); 
                
                // Segmentation 
                // Loop over the voxel integration points and select only those 
                // which are greater than the threshold value 
                vc = 0 ;   
                nv = 0 ;               
                for (int kk=k*nipez, ig1dz=0; kk<(k+1)*nipez && ig1dz< nipez; kk++,ig1dz++)
                {
                    for (int jj=j*nipey, ig1dy=0; jj<(j+1)*nipey && ig1dy< nipey; jj++,ig1dy++)
                    {
                        for (int ii=i*nipex, ig1dx=0; ii<(i+1)*nipex && ig1dx< nipex; ii++,ig1dx++)
                        {
                            if (fe[vc] > thrsh)
                            {
                                ipIndicesElem.push_back(ii); 
                                ipIndicesElem.push_back(jj); 
                                ipIndicesElem.push_back(kk);  
                                nv++; 
                            } 
                            vc++; 
                        }
                    }
                }
                // Output data structure (number of voxels of Element 1, i1,j1,k1, i2,j2,k2, .... )
                ipIndices.push_back(nv); 
                ipIndices.insert(ipIndices.end(), ipIndicesElem.begin(), ipIndicesElem.end() ); 	    
            }
        }
        
    }
    
    delete[] knotXiImage   ; 
    delete[] knotEtaImage  ; 
    delete[] knotZetaImage ; 
    
    delete[] xg  ; 
    delete[] yg  ; 
    delete[] zg  ; 
    
    delete[] xge ; 
    delete[] yge ; 
    delete[] zge ; 
    
    delete[] fe  ; 
      
    return ipIndices;              
} 
*/ 


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
                                        double* rhs, int ndof ) 
{
    // Image parameters (x=xi=k, y=eta=i, z=zeta=j)
    Eigen::VectorXd pixCx = Eigen::VectorXd::LinSpaced(sgk,xminIm+dx/2,xmaxIm-dx/2); 
    Eigen::VectorXd pixCy = Eigen::VectorXd::LinSpaced(sgi,yminIm+dy/2,ymaxIm-dy/2); 
    Eigen::VectorXd pixCz = Eigen::VectorXd::LinSpaced(sgj,zminIm+dz/2,zmaxIm-dz/2); 
    int sknotXiImage   = sgk+2;
    int sknotEtaImage  = sgi+2;
    int sknotZetaImage = sgj+2;

    double* knotXiImage   = new double[sknotXiImage];  
    double* knotEtaImage  = new double[sknotEtaImage];
    double* knotZetaImage = new double[sknotZetaImage]; 
    
    knotXiImage[0]    = xminIm+dx/2 ; 
    knotXiImage[sknotXiImage-1] = xmaxIm-dx/2 ; 
    for (int i=1;i<sgk+1;i++)
    {
        knotXiImage[i] = pixCx(i-1); 
    } 
    knotEtaImage[0]    = yminIm+dy/2; 
    knotEtaImage[sknotEtaImage-1] = ymaxIm-dy/2;   
    for (int i=1;i<sgi+1;i++)
    {
        knotEtaImage[i] = pixCy(i-1); 
    }
    knotZetaImage[0]    = zminIm+dz/2;
    knotZetaImage[sknotZetaImage-1] = zmaxIm-dz/2;         
    for (int i=1; i<sgj+1; i++)
    {
        knotZetaImage[i] = pixCz(i-1); 
    }
    
 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi   + 1 ;
	int nbf_elem_eta  = deg_eta  + 1 ;
	int nbf_elem_zeta = deg_zeta + 1 ;  
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta ;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
	
 
	int nip_xi   =  ne_xi*nipex   ; 
    int nip_eta  =  ne_eta*nipey  ; 
    int nip_zeta =  ne_zeta*nipez ;  
    //int nipe = nipex*nipey*nipez ; 
    
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
    
    
    double* xg = new double[nip_xi]   ; 
    double* yg = new double[nip_eta]  ; 
    double* zg = new double[nip_zeta] ; 
    
 
 
     
    // Univariate basis functions at the 1D integration points 
	double** Nxig   = array_2d(nip_xi , nbf_elem_xi);
	double** Netag  = array_2d(nip_eta, nbf_elem_eta);
    double** Nzetag = array_2d(nip_zeta, nbf_elem_zeta); 
    int spanx; 
    int spany; 
    int spanz; 
    int nxi   = nbf_xi   - 1 ;  
	int neta  = nbf_eta  - 1 ;
	int nzeta = nbf_zeta - 1 ; 
    for (int i=0; i<nip_xi ; i++)
    {
        xg[i] = knotXi[0]+mes_xi/2+i*mes_xi; 
        spanx = findspan(nxi, deg_xi, xg[i], knotXi);
        basisfuns(spanx, xg[i], deg_xi, knotXi, Nxig[i]);
    }
    for (int i=0; i<nip_eta ; i++)
    {
        yg[i] = knotEta[0]+mes_eta/2+i*mes_eta; 
        spany = findspan(neta, deg_eta, yg[i], knotEta);
        basisfuns(spany, yg[i], deg_eta, knotEta, Netag[i]);
    }    
    for (int i=0; i<nip_zeta; i++)
    {
        zg[i] = knotZeta[0]+mes_zeta/2+i*mes_zeta; 
        spanz = findspan(nzeta, deg_zeta, zg[i], knotZeta);
        basisfuns(spanz, zg[i], deg_zeta, knotZeta, Nzetag[i]);        
    }  
    // Elementary operations 
    double* N  = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhi ; // g(x+u(x)) : advected gray-level 
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
    int ig=0 ; // index of the integration points  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
               
                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ; 
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        // ibf + jbf*nbf_elem_xi + kbf*nbf_elem_xi*nbf_elem_eta
                                        N[bfc] = Nzetag[kk][kbf]*Netag[jj][jbf]*Nxig[ii][ibf] ; 
                                        bfc++;     
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0. ; 
                            uyg = 0. ; 
                            uzg = 0. ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];    
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                                    knotXiImage, sknotXiImage, 
                                    knotEtaImage, sknotEtaImage,
                                    knotZetaImage, sknotZetaImage, 
                                    xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg ); 
                            
                                    
                            res = fip[ig] - goPhi ; 
                            wgResdfdx =  mes*res*dfipdx[ig] ;
                            wgResdfdy =  mes*res*dfipdy[ig] ; 
                            wgResdfdz =  mes*res*dfipdz[ig] ; 
                            
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            +=  wgResdfdx*N[ibf] ; 
                                be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                                be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ;  
                                 
                            } 
                            ig++; 
                        }  
                    }        
                }
                
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM(ecounter,ibf)]       += be[ibf]; 
                    rhs[NOELEM(ecounter,ibf)+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM(ecounter,ibf)+2*nbf] += be[ibf+2*nbf_elem];
                }
                ecounter++;        
            }
        }
    }
    
    delete[] knotXiImage   ; 
    delete[] knotEtaImage  ; 
    delete[] knotZetaImage ;  
    delete[] xg ; 
    delete[] yg ; 
    delete[] zg ; 
 
    delete[] N  ; 
    delete[] be ;
    
	delete_array_2d(Nxig);
	delete_array_2d(Netag);
	delete_array_2d(Nzetag);
 
} */

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
                                    double* rhs, int ndof )
{ 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
 
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhi ; // g(x+u(x)) : advected gray-level 
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    double ssd=0; 
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
 
    int ig=0 ; // index of the integration points  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
               
                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                                uyg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                                uzg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 					               
                            goPhi  = EvaluateCBspline2OnOnePoint(g,sgi,sgj,sgk,
                                     xmin, ymin, zmin, xmax, ymax, zmax, 
                                     dx, dy, dz, xc, sxc, yc, syc, zc, szc, 
                                     xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg );  
                                    
                            res = fip[ig] - goPhi ; 
                            wgResdfdx =  mes*res*dfipdx[ig] ;
                            wgResdfdy =  mes*res*dfipdy[ig] ; 
                            wgResdfdz =  mes*res*dfipdz[ig] ; 
                            
                            ssd += std::pow(res, 2) ; 
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            +=  wgResdfdx*N[ibf] ; 
                                be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                                be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ; 
                            } 
                            ig++; 
                        }  
                    }        
                }
                
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                ecounter++;        
            }
        }
    }
 
    delete[] N  ; 
    delete[] be ; 
    
    return ssd*mes ;  
}   
   
                                    


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
                                double* glr, int sglr )
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipx  = nipex*ne_xi  ; 
    int nipy  = nipey*ne_eta ; 
    int nipxy =  nipx * nipy ; 
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
    
    
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double goPhi ; // g(x+u(x)) : advected gray-level 
    
    int ig=0; 
    int ecounter=0; 
    int bfc; 
 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)

                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            goPhi  = EvaluateCBspline2OnOnePoint(g,sgi,sgj,sgk,
                                     xmin, ymin, zmin, xmax, ymax, zmax, 
                                     dx, dy, dz, xc, sxc, yc, syc, zc, szc, 
                                     xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg );                                      
                                    
                                    
                            glr[ii+jj*nipx+kk*nipxy] = fip[ig] - goPhi ; 
 
                            ig++; 
                        }  
                    }        
                }
                ecounter++;      
            }
        }
    }
    
    delete[] N;    
} 

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
                                    double* rhs, int ndof )
{ 
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
  
    double mes_xi   = (knotXi[deg_xi+1] - knotXi[deg_xi])/nipex ; 
    double mes_eta  = (knotEta[deg_eta+1] - knotEta[deg_eta])/nipey ; 
    double mes_zeta = (knotZeta[deg_zeta+1] - knotZeta[deg_zeta])/nipez ; 
    double mes = mes_xi*mes_eta*mes_zeta ; // Constant measure of an integration voxel approximative to 1
    //double volElem = mes*nipe;  // Volume of an element 
 
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double res ; // f-goPhi at the integration point 
    double goPhi ; // g(x+u(x)) : advected gray-level 
    double wgResdfdx;  
    double wgResdfdy;
    double wgResdfdz;
    
    double ssd=0; 
    
    // Assembly 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }
    
 
    int ig=0 ; // index of the integration points  
    int ecounter=0 ; // index of element 
    int bfc ; 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)
                
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                } 
                
               
                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]]; 
                                uyg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+nbf];
                                uzg  +=  N[ibf]*U[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 					               
                            goPhi  = EvaluateCBspline3OnOnePoint(g,sgi,sgj,sgk,
                                     xmin, ymin, zmin, xmax, ymax, zmax, 
                                     dx, dy, dz, xc, sxc, yc, syc, zc, szc, 
                                     xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg );  
                                    
                            res = fip[ig] - goPhi ; 
                            wgResdfdx =  mes*res*dfipdx[ig] ;
                            wgResdfdy =  mes*res*dfipdy[ig] ; 
                            wgResdfdz =  mes*res*dfipdz[ig] ; 
                            
                            ssd += std::pow(res, 2) ; 
                            
                            // Setting the elemental right hand side vector 
                            for (int ibf=0; ibf<nbf_elem; ibf++)
                            {
                                be[ibf]            +=  wgResdfdx*N[ibf] ; 
                                be[ibf+nbf_elem]   +=  wgResdfdy*N[ibf] ; 
                                be[ibf+2*nbf_elem] +=  wgResdfdz*N[ibf] ; 
                            } 
                            ig++; 
                        }  
                    }        
                }
                
                // Adding the elementary contribution to the right hand side 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]]       += be[ibf]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+nbf]   += be[ibf+nbf_elem]; 
                    rhs[NOELEM[ibf+ecounter*sjNOELEM]+2*nbf] += be[ibf+2*nbf_elem]; 
                }
                ecounter++;        
            }
        }
    }
 
    delete[] N  ; 
    delete[] be ; 
    
    return ssd*mes ;  
}   
   
                                    


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
                                double* glr, int sglr ) 
{
    int nbf_xi   = sknotXi   -1 - deg_xi   ; 
	int nbf_eta  = sknotEta  -1 - deg_eta  ;  
	int nbf_zeta = sknotZeta -1 - deg_zeta ;
	
	int nbf_elem_xi   = deg_xi  +1 ;
	int nbf_elem_eta  = deg_eta +1 ;
	int nbf_elem_zeta = deg_zeta+1 ; 
	int nbf_elem = nbf_elem_xi*nbf_elem_eta*nbf_elem_zeta;
	
	int nbf = nbf_xi*nbf_eta*nbf_zeta ;  
	
	int ne_xi   = nbf_xi   - deg_xi   ;
	int ne_eta  = nbf_eta  - deg_eta  ; 
	int ne_zeta = nbf_zeta - deg_zeta ;  
 
    int nipx  = nipex*ne_xi  ; 
    int nipy  = nipey*ne_eta ; 
    int nipxy =  nipx * nipy ; 
    
    Eigen::MatrixXi NOELEM = getNoelem(ne_xi, ne_eta, ne_zeta, deg_xi, deg_eta, deg_zeta );
    
    
    // Elementary operations 
    double* N = new double[nbf_elem] ; // Trivariate basis functions defined on a point 
    double uxg ; // Displacement  ux of the voxel 
    double uyg ; // Displacement  uy of the voxel
    double uzg ; // Displacement  uz of the voxel
    double goPhi ; // g(x+u(x)) : advected gray-level 
    
    int ig=0; 
    int ecounter=0; 
    int bfc; 
 
    for (int k=0; k<ne_zeta ;k++)
    {   
        for (int j=0; j<ne_eta; j++)
        {             
            for (int i=0; i<ne_xi; i++)
            {             
                // Current element (i,j,k)

                // Loop over the integration points (voxels) of the element
                
                for (int kk=k*nipez; kk<(k+1)*nipez; kk++)
                {
                    for (int jj=j*nipey; jj<(j+1)*nipey; jj++)
                    {
                        for (int ii=i*nipex; ii<(i+1)*nipex; ii++)
                        {
                        
                            //Getting trivariate basis functions at the integration point
                            bfc = 0 ;   
                            for (int kbf=0; kbf<nbf_elem_zeta; kbf++)
                            {
                                for (int jbf=0; jbf<nbf_elem_eta; jbf++)
                                {
                                    for (int ibf=0; ibf<nbf_elem_xi; ibf++)
                                    {
                                        N[bfc] = Nzeta[kbf+kk*sjNzeta]*Neta[jbf+jj*sjNeta]*Nxi[ibf+ii*sjNxi] ;  
                                        bfc++ ;  
                                    }
                                }
                            }
                            // Computing the displacement at the integration points 
                            // by the B-spline linear combination 
                            uxg = 0 ; 
                            uyg = 0 ; 
                            uzg = 0 ; 
                            for (int ibf=0; ibf< nbf_elem; ibf++)
                            {
                                uxg  +=  N[ibf]*U[NOELEM(ecounter,ibf)]; 
                                uyg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+nbf];
                                uzg  +=  N[ibf]*U[NOELEM(ecounter,ibf)+2*nbf];     
                            }
                            // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            goPhi  = EvaluateCBspline3OnOnePoint(g,sgi,sgj,sgk,
                                     xmin, ymin, zmin, xmax, ymax, zmax, 
                                     dx, dy, dz, xc, sxc, yc, syc, zc, szc, 
                                     xg[ii]+uxg, yg[jj]+uyg ,zg[kk]+uzg );                                      
                                    
                                    
                            glr[ii+jj*nipx+kk*nipxy] = fip[ig] - goPhi ; 
 
                            ig++; 
                        }  
                    }        
                }
                ecounter++;      
            }
        }
    }
    
    delete[] N;    
}                              


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
                   double* nnz_values, int nnz)                           
{
    int nbf_elem = sjN ; // Number of nodes per element = number of lagrange basis functions 
    int nip = siN ; // Number of integration points per element 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                Eigen::MatrixXd He(3*nbf_elem,3*nbf_elem);
                Eigen::MatrixXd Ne(3,3*nbf_elem); 
                Eigen::MatrixXd gradFgradFt(3,3);
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                int I,J ; 
                // Current element 
                Ne.fill(0.);  
                He.fill(0.); // Elementary hessian matrix initialization      
                
                // Loop over integration points   
                int Ip ; 
                for (int ip=0; ip< nip; ip++)
                {
                    // Setting the elemental basis function matrix 
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                     
                    // Setting the elementary basis functions matrix and 
                    // Computing the components of the jacobian of the isoparametric transformation 
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
   					{
                        I =  ibf+i*sje    ; 
                        J =  ibf + ip*sjN ; 
       					
       					Ne(0,ibf)            = N[J] ; 
       					Ne(1,ibf+nbf_elem)   = N[J] ;  
       					Ne(2,ibf+2*nbf_elem) = N[J] ; 
       					
       					dxdxi   +=  dNdxi[J]   * n[0+e[I]*sjn] ; 
       					dxdeta  +=  dNdeta[J]  * n[0+e[I]*sjn] ;   
       					dxdzeta +=  dNdzeta[J] * n[0+e[I]*sjn] ; 
       					
       					dydxi   +=  dNdxi[J]   * n[1+e[I]*sjn] ; 
       					dydeta  +=  dNdeta[J]  * n[1+e[I]*sjn] ; 
       					dydzeta +=  dNdzeta[J] * n[1+e[I]*sjn] ;
       					
       					dzdxi   +=  dNdxi[J]   * n[2+e[I]*sjn] ;
       					dzdeta  +=  dNdeta[J]  * n[2+e[I]*sjn] ;  
       					dzdzeta +=  dNdzeta[J] * n[2+e[I]*sjn] ; 	
   					}  
   					
   					// Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
   					
   					//Setting the gradient tensor product 
   					Ip = ip + i*nip ;    
                    gradFgradFt(0,0) = dfipdx[Ip]*dfipdx[Ip] ; 
                    gradFgradFt(0,1) = dfipdx[Ip]*dfipdy[Ip] ;
                    gradFgradFt(0,2) = dfipdx[Ip]*dfipdz[Ip] ;
                    gradFgradFt(1,0) = dfipdy[Ip]*dfipdx[Ip] ;
                    gradFgradFt(1,1) = dfipdy[Ip]*dfipdy[Ip] ; 
                    gradFgradFt(1,2) = dfipdy[Ip]*dfipdz[Ip] ;
                    gradFgradFt(2,0) = dfipdz[Ip]*dfipdx[Ip] ;
                    gradFgradFt(2,1) = dfipdz[Ip]*dfipdy[Ip] ;
                    gradFgradFt(2,2) = dfipdz[Ip]*dfipdz[Ip] ;
                    
                    // Adding to the voxel summation 
   					He = He + w[ip]*std::abs(detJ)*(Ne.transpose()*gradFgradFt*Ne) ;    										                     
                }
                int sg = 9*nbf_elem*nbf_elem* i ;
                // Elementary contribution to the global Hessian matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = e[ibf+i*sje] ; 
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
       					J = e[jbf+i*sje] ;
     
       		            indexI[sg] = conn[0+I*sjconn]   ;
       		            indexJ[sg] = conn[0+J*sjconn]   ;
       		            nnz_values[sg] = He(ibf,jbf);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn]  ;
       		            indexJ[sg] = conn[0+J*sjconn]  ;
       		            nnz_values[sg] = He(ibf+nbf_elem,jbf);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[0+I*sjconn]  ;
       		            indexJ[sg] = conn[1+J*sjconn]  ;
       		            nnz_values[sg] = He(ibf,jbf+nbf_elem);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn] ;
       		            indexJ[sg] = conn[1+J*sjconn] ;
       		            nnz_values[sg] = He(ibf+nbf_elem,jbf+nbf_elem);
       		            sg ++ ;
       		                   
       		            indexI[sg] = conn[2+I*sjconn]  ;
       		            indexJ[sg] = conn[0+J*sjconn]  ;
       		            nnz_values[sg] = He(ibf+2*nbf_elem,jbf);
       		            sg ++ ;    		            
       		            
       		            indexI[sg] = conn[0+I*sjconn]  ;
       		            indexJ[sg] = conn[2+J*sjconn]  ;
       		            nnz_values[sg] =  He(ibf,jbf+2*nbf_elem);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn]  ;
       		            indexJ[sg] = conn[2+J*sjconn]  ;
       		            nnz_values[sg] =  He(ibf+nbf_elem,jbf+2*nbf_elem);
       		            sg ++ ;   		            
    
    		            indexI[sg] = conn[2+I*sjconn]   ;
       		            indexJ[sg] = conn[1+J*sjconn]   ;
       		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+nbf_elem);
       		            sg ++ ; 
       		               		            
   		                indexI[sg] = conn[2+I*sjconn]   ;
       		            indexJ[sg] = conn[2+J*sjconn]   ;
       		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+2*nbf_elem);
       		            sg ++ ; 
                    }
               }                
            }
    }        
}

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
                double* rhs, int ndof )                 
{
    int nbf_elem = sjN ; 
    int nip = siN ; 
    
    double ssd=0; 
    
 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }  
     
    # pragma omp parallel 
    {
        #pragma omp for reduction(+:ssd)
            for (int i=0; i<sie; i++)
            {
                double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
                double res ; // f-goPhi at the integration point 
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double detJwgResdfdx;  
                double detJwgResdfdy;
                double detJwgResdfdz;   
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                int I,J ; 
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                }
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
                     
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                    
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                         
                         dxdxi   +=  dNdxi[J]   * n[0+e[I]*sjn] ; 
       					 dxdeta  +=  dNdeta[J]  * n[0+e[I]*sjn] ;   
       					 dxdzeta +=  dNdzeta[J] * n[0+e[I]*sjn] ; 
       				 	
       					 dydxi   +=  dNdxi[J]   * n[1+e[I]*sjn] ; 
       					 dydeta  +=  dNdeta[J]  * n[1+e[I]*sjn] ; 
       					 dydzeta +=  dNdzeta[J] * n[1+e[I]*sjn] ;
       				 	
       					 dzdxi   +=  dNdxi[J]   * n[2+e[I]*sjn] ;
       					 dzdeta  +=  dNdeta[J]  * n[2+e[I]*sjn] ;  
       					 dzdzeta +=  dNdzeta[J] * n[2+e[I]*sjn] ; 	
                     } 
                    // Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
                             
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                     goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                            knotXiImage, sknotXiImage, 
                            knotEtaImage, sknotEtaImage,
                            knotZetaImage, sknotZetaImage, 
                            xa, ya , za ); 
   					 Ip  = ip + i*nip ;    
                     res = fip[Ip] - goPhi ; 
                     detJwgResdfdx =  w[ip]*std::abs(detJ)*res*dfipdx[Ip] ;
                     detJwgResdfdy =  w[ip]*std::abs(detJ)*res*dfipdy[Ip] ; 
                     detJwgResdfdz =  w[ip]*std::abs(detJ)*res*dfipdz[Ip] ;  
                     
                     ssd += std::abs(detJ) * w[ip] * std::pow(res,2) ;   // Adding the contribution to the global SSD 
                     
                     // Setting the elemental right hand side vector 
                     for (int ibf=0; ibf<nbf_elem; ibf++)
                     {
                        J = ibf + ip*sjN ; 
                        be[ibf]            +=  detJwgResdfdx*N[J] ; 
                        be[ibf+nbf_elem]   +=  detJwgResdfdy*N[J] ; 
                        be[ibf+2*nbf_elem] +=  detJwgResdfdz*N[J] ; 
                     }                                            
                }
                # pragma omp critical 
                {
                    
                    // Adding the elementary contribution to the right hand side 
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        I = e[ibf+i*sje] ; 
                        rhs[conn[0+I*sjconn]]  += be[ibf] ;  
                        rhs[conn[1+I*sjconn]]  += be[ibf+nbf_elem] ; 
                        rhs[conn[2+I*sjconn]]  += be[ibf+2*nbf_elem];   
                    }                
                }
                delete[] be ;        
            } 
    }
    return ssd ; 
}

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
            double* glr, int sglr )  
{
    int nip      = siN ; // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
            for (int i=0; i<sie; i++)
            {
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                int I,J ; 
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
 
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                     } 
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                     goPhi = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                            knotXiImage, sknotXiImage, 
                            knotEtaImage, sknotEtaImage,
                            knotZetaImage, sknotZetaImage, 
                            xa, ya , za ); 
   					 Ip  = ip + i*nip ;    
                     glr[Ip] = fip[Ip] - goPhi ;                                      
                } 
            }     
    }
} 

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
                                  double* glr, int sglr ) 
{
    #pragma omp parallel 
    {
          #pragma omp for 
            for (int ip=0; ip< sxi; ip++)
            {
                double N[4] = { xi[ip], eta[ip], zeta[ip], 1-xi[ip]-eta[ip]-zeta[ip] } ; 
                double xa = 0 ; 
                double ya = 0 ; 
                double za = 0 ; 
                int I ; 
                // Computing the advected voxels 
                for (int ibf=0; ibf < 4; ibf++)
                {
                    I =  ibf+  ie[ip]*sje   ; 
                    xa += N[ibf] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                    ya += N[ibf] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                    za += N[ibf] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                } 
                // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                glr[ip] = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
                            knotXiImage, sknotXiImage, 
                            knotEtaImage, sknotEtaImage,
                            knotZetaImage, sknotZetaImage, 
                            xa, ya , za );                                             
            }   
    }
} 


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
                            double* rhs, int ndof )               
{
    int nbf_elem = sjN ; 
    int nip = siN ; 
    
    double ssd=0; 
    
 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }  
     
    # pragma omp parallel 
    {
        #pragma omp for reduction(+:ssd)
            for (int i=0; i<sie; i++)
            {
                double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
                double res ; // f-goPhi at the integration point 
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double detJwgResdfdx;  
                double detJwgResdfdy;
                double detJwgResdfdz;   
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                int I,J ; 
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                }
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
                     
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                    
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                         
                         dxdxi   +=  dNdxi[J]   * n[0+e[I]*sjn] ; 
       					 dxdeta  +=  dNdeta[J]  * n[0+e[I]*sjn] ;   
       					 dxdzeta +=  dNdzeta[J] * n[0+e[I]*sjn] ; 
       				 	
       					 dydxi   +=  dNdxi[J]   * n[1+e[I]*sjn] ; 
       					 dydeta  +=  dNdeta[J]  * n[1+e[I]*sjn] ; 
       					 dydzeta +=  dNdzeta[J] * n[1+e[I]*sjn] ;
       				 	
       					 dzdxi   +=  dNdxi[J]   * n[2+e[I]*sjn] ;
       					 dzdeta  +=  dNdeta[J]  * n[2+e[I]*sjn] ;  
       					 dzdzeta +=  dNdzeta[J] * n[2+e[I]*sjn] ; 	
                     } 
                    // Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
                             
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                            
                    goPhi  = EvaluateCBspline3OnOnePoint(g,sgi,sgj,sgk,
                                                    xmin, ymin, zmin,
                                                    xmax, ymax, zmax, 
                                                    dx, dy, dz, 
                                                    xc, sxc, yc, syc, zc, szc, 
                                                    xa, ya , za ); 

   					 Ip  = ip + i*nip ;    
                     res = fip[Ip] - goPhi ; 
                     detJwgResdfdx =  w[ip]*std::abs(detJ)*res*dfipdx[Ip] ;
                     detJwgResdfdy =  w[ip]*std::abs(detJ)*res*dfipdy[Ip] ; 
                     detJwgResdfdz =  w[ip]*std::abs(detJ)*res*dfipdz[Ip] ;  
                     
                     ssd += std::abs(detJ) * w[ip] * std::pow(res,2) ;   // Adding the contribution to the global SSD 
                     
                     // Setting the elemental right hand side vector 
                     for (int ibf=0; ibf<nbf_elem; ibf++)
                     {
                        J = ibf + ip*sjN ; 
                        be[ibf]            +=  detJwgResdfdx*N[J] ; 
                        be[ibf+nbf_elem]   +=  detJwgResdfdy*N[J] ; 
                        be[ibf+2*nbf_elem] +=  detJwgResdfdz*N[J] ; 
                     }                                            
                }
                # pragma omp critical 
                {
                    
                    // Adding the elementary contribution to the right hand side 
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        I = e[ibf+i*sje] ; 
                        rhs[conn[0+I*sjconn]]  += be[ibf] ;  
                        rhs[conn[1+I*sjconn]]  += be[ibf+nbf_elem] ; 
                        rhs[conn[2+I*sjconn]]  += be[ibf+2*nbf_elem];   
                    }                
                }
                delete[] be ;        
            } 
    }
    return ssd ; 
}

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
                            double* glr, int sglr ) 
{
    int nip      = siN ; // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
            for (int i=0; i<sie; i++)
            {
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                int I,J ; 
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
 
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                     } 
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                     goPhi = EvaluateCBspline3OnOnePoint(g,sgi,sgj,sgk,
                                                    xmin, ymin, zmin,
                                                    xmax, ymax, zmax, 
                                                    dx, dy, dz, 
                                                    xc, sxc, yc, syc, zc, szc, 
                                                    xa, ya , za ); 
                            
   					 Ip  = ip + i*nip ;    
                     glr[Ip] = fip[Ip] - goPhi ;                                      
                } 
            }     
    }
} 

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
                            double* glr, int sglr )
{
    #pragma omp parallel 
    {
          #pragma omp for 
            for (int ip=0; ip< sxi; ip++)
            {
                double N[4] = { xi[ip], eta[ip], zeta[ip], 1-xi[ip]-eta[ip]-zeta[ip] } ; 
                double xa = 0 ; 
                double ya = 0 ; 
                double za = 0 ; 
                int I ; 
                // Computing the advected voxels 
                for (int ibf=0; ibf < 4; ibf++)
                {
                    I =  ibf+  ie[ip]*sje   ; 
                    xa += N[ibf] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                    ya += N[ibf] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                    za += N[ibf] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                } 
                // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                glr[ip] = EvaluateCBspline3OnOnePoint(g,sgi,sgj,sgk,
                                                    xmin, ymin, zmin,
                                                    xmax, ymax, zmax, 
                                                    dx, dy, dz, 
                                                    xc, sxc, yc, syc, zc, szc, 
                                                    xa, ya , za );                                         
            }   
    }
} 


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
                               double* rhs, int ndof ) 
{
    int nbf_elem = sjN ; 
    int nip = siN ; 
    
    double ssd=0; 
    
 
    // Setting the right hand side to zero as we are adding contributions 
    for (int idof=0; idof < ndof ; idof++ )
    {
        rhs[idof] = 0 ; 
    }  
     
    # pragma omp parallel 
    {
        #pragma omp for reduction(+:ssd)
            for (int i=0; i<sie; i++)
            {
                double* be = new double[3*nbf_elem]; // Elemntary right hand side vector 
                double res ; // f-goPhi at the integration point 
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double detJwgResdfdx;  
                double detJwgResdfdy;
                double detJwgResdfdz;   
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                int I,J ; 
                // Setting elementary right hand side vector to 0
                for (int idofe=0; idofe<3*nbf_elem; idofe++)
                {
                    be[idofe] = 0 ; 
                }
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
                     
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                    
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                         
                         dxdxi   +=  dNdxi[J]   * n[0+e[I]*sjn] ; 
       					 dxdeta  +=  dNdeta[J]  * n[0+e[I]*sjn] ;   
       					 dxdzeta +=  dNdzeta[J] * n[0+e[I]*sjn] ; 
       				 	
       					 dydxi   +=  dNdxi[J]   * n[1+e[I]*sjn] ; 
       					 dydeta  +=  dNdeta[J]  * n[1+e[I]*sjn] ; 
       					 dydzeta +=  dNdzeta[J] * n[1+e[I]*sjn] ;
       				 	
       					 dzdxi   +=  dNdxi[J]   * n[2+e[I]*sjn] ;
       					 dzdeta  +=  dNdeta[J]  * n[2+e[I]*sjn] ;  
       					 dzdzeta +=  dNdzeta[J] * n[2+e[I]*sjn] ; 	
                     } 
                    // Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
                             
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z))                                          
                     goPhi = Evaluate3dBsplineOnOnePoint(xa,ya,za, 
                                                          knotXiImage_g, sknotXiImage_g, 
                                                          knotEtaImage_g, sknotEtaImage_g, 
                                                          knotZetaImage_g, sknotZetaImage_g, 
                                                          deg_xi_g, deg_eta_g, deg_zeta_g, 
                                                          lsc_g, s_lsc_g ) ; 
 
   					 Ip  = ip + i*nip ;    
                     res = fip[Ip] - goPhi ; 
                     detJwgResdfdx =  w[ip]*std::abs(detJ)*res*dfipdx[Ip] ;
                     detJwgResdfdy =  w[ip]*std::abs(detJ)*res*dfipdy[Ip] ; 
                     detJwgResdfdz =  w[ip]*std::abs(detJ)*res*dfipdz[Ip] ;  
                     
                     ssd += std::abs(detJ) * w[ip] * std::pow(res,2) ;   // Adding the contribution to the global SSD 
                     
                     // Setting the elemental right hand side vector 
                     for (int ibf=0; ibf<nbf_elem; ibf++)
                     {
                        J = ibf + ip*sjN ; 
                        be[ibf]            +=  detJwgResdfdx*N[J] ; 
                        be[ibf+nbf_elem]   +=  detJwgResdfdy*N[J] ; 
                        be[ibf+2*nbf_elem] +=  detJwgResdfdz*N[J] ; 
                     }                                            
                }
                # pragma omp critical 
                {
                    
                    // Adding the elementary contribution to the right hand side 
                    for (int ibf=0; ibf<nbf_elem; ibf++)
                    {
                        I = e[ibf+i*sje] ; 
                        rhs[conn[0+I*sjconn]]  += be[ibf] ;  
                        rhs[conn[1+I*sjconn]]  += be[ibf+nbf_elem] ; 
                        rhs[conn[2+I*sjconn]]  += be[ibf+2*nbf_elem];   
                    }                
                }
                delete[] be ;        
            } 
    }
    return ssd ; 
}     
                          

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
                         double* glr, int sglr ) 
{
    int nip      = siN ; // Number of integration points per element 
    int nbf_elem = sjN ; // Number of basis functions per element 
    # pragma omp parallel 
    {
            for (int i=0; i<sie; i++)
            {
                double goPhi ; // g(x+u(x)) : advected gray-level 
                double xa   ; // Advected voxels 
                double ya   ; 
                double za   ;             
                int I,J ; 
                int Ip ; 
                
                for (int ip=0; ip<nip ; ip++)
                {
                 
                    xa = 0 ; 
                    ya = 0 ; 
                    za = 0 ; 
 
                     // Computing the advected voxels 
                     // Computing the components of the jacobian of the isoparametric transformation 
                     for (int ibf=0; ibf < nbf_elem; ibf++)
                     {
                         I =  ibf+  i*sje   ; 
                         J =  ibf + ip*sjN  ; 
                         
                         xa += N[J] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                         ya += N[J] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                         za += N[J] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                     } 
                     // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                     goPhi = Evaluate3dBsplineOnOnePoint(xa,ya,za, 
                                                          knotXiImage_g, sknotXiImage_g, 
                                                          knotEtaImage_g, sknotEtaImage_g, 
                                                          knotZetaImage_g, sknotZetaImage_g, 
                                                          deg_xi_g, deg_eta_g, deg_zeta_g, 
                                                          lsc_g, s_lsc_g ) ; 
                            
   					 Ip  = ip + i*nip ;    
                     glr[Ip] = fip[Ip] - goPhi ;                                      
                } 
            }     
    }
}                          

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
                               double* glr, int sglr )  
{
    #pragma omp parallel 
    {
          #pragma omp for 
            for (int ip=0; ip< sxi; ip++)
            {
                double N[4] = { xi[ip], eta[ip], zeta[ip], 1-xi[ip]-eta[ip]-zeta[ip] } ; 
                double xa = 0 ; 
                double ya = 0 ; 
                double za = 0 ; 
                int I ; 
                // Computing the advected voxels 
                for (int ibf=0; ibf < 4; ibf++)
                {
                    I =  ibf+  ie[ip]*sje   ; 
                    xa += N[ibf] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
                    ya += N[ibf] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
                    za += N[ibf] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
                } 
                // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
                glr[ip] = Evaluate3dBsplineOnOnePoint(xa,ya,za, 
                                                      knotXiImage_g, sknotXiImage_g, 
                                                      knotEtaImage_g, sknotEtaImage_g, 
                                                      knotZetaImage_g, sknotZetaImage_g, 
                                                      deg_xi_g, deg_eta_g, deg_zeta_g, 
                                                      lsc_g, s_lsc_g ) ;                                       
            }   
    }
}                                
                               
                               


            
            
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
                    double* nnz_values, int nnz)
{
    int nbf_elem = sjN ; 
    int nip = siN ; 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                Eigen::MatrixXd Leb(nbf_elem,nbf_elem);  // Block of the elementary Laplacian operator 
                Eigen::MatrixXd dNeb(3,nbf_elem);       // Block elementary differential operator 
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                Leb.fill(0) ;  
                // Loop over integration points 
                for (int ip=0; ip<nip; ip++)
                {
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                    // Setting the elementary differential operator
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
        			{
       					dxdxi   +=  dNdxi[ibf + ip*sjN]   * n[0+e[ibf+i*sje]*sjn] ; 
       					dxdeta  +=  dNdeta[ibf + ip*sjN]  * n[0+e[ibf+i*sje]*sjn] ;   
       					dxdzeta +=  dNdzeta[ibf + ip*sjN] * n[0+e[ibf+i*sje]*sjn] ; 
       					
       					dydxi   +=  dNdxi[ibf + ip*sjN]   * n[1+e[ibf+i*sje]*sjn] ; 
       					dydeta  +=  dNdeta[ibf + ip*sjN]  * n[1+e[ibf+i*sje]*sjn] ; 
       					dydzeta +=  dNdzeta[ibf + ip*sjN] * n[1+e[ibf+i*sje]*sjn] ;
       					
       					dzdxi   +=  dNdxi[ibf + ip*sjN]   * n[2+e[ibf+i*sje]*sjn] ;
       					dzdeta  +=  dNdeta[ibf + ip*sjN]  * n[2+e[ibf+i*sje]*sjn] ;  
       					dzdzeta +=  dNdzeta[ibf + ip*sjN] * n[2+e[ibf+i*sje]*sjn] ; 	    					            			 
        			}
   					// Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
                    // Setting the block differential operator at the integration point 
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
        			{
                        dNeb(0,ibf)  = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dNdxi[ibf + ip*sjN] +  
                                       (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dNdeta[ibf + ip*sjN] + 
                                       ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dNdzeta[ibf + ip*sjN] ; 
                                       
                		dNeb(1,ibf)  = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dNdxi[ibf + ip*sjN] +
                                       ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dNdeta[ibf + ip*sjN] +
                                       (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dNdzeta[ibf + ip*sjN] ; 
                		
                		dNeb(2,ibf)  = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dNdxi[ibf + ip*sjN] +
                                       (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dNdeta[ibf + ip*sjN] +
                                       ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dNdzeta[ibf + ip*sjN] ;     					            			 
        			}      
			        Leb = Leb + dNeb.transpose()*dNeb*(w[ip]*std::abs(detJ)) ;                				 
                }
                int I,J ; 
                int sg = 3*nbf_elem*nbf_elem * i ; 
                // Elementary contribution to the global matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = e[ibf+i*sje] ;
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = e[jbf+i*sje];
    					
    		            indexI[sg] = conn[0+I*sjconn] ;
    		            indexJ[sg] = conn[0+J*sjconn] ;
    		            nnz_values[sg] = Leb(ibf,jbf);
    		            sg ++ ;  	
    		            
		            	indexI[sg] = conn[1+I*sjconn] ;
    		            indexJ[sg] = conn[1+J*sjconn] ;
    		            nnz_values[sg] = Leb(ibf,jbf);
    		            sg ++ ;
    		               		            
		                indexI[sg] = conn[2+I*sjconn] ;
    		            indexJ[sg] = conn[2+J*sjconn] ;
    		            nnz_values[sg] =  Leb(ibf,jbf);
    		            sg ++ ; 
                    }
                }                
            }
    }

}

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
                    double* nnz_values, int nnz)
{
    int nbf_elem = sjN ; 
    int nip = siN ; 
    Eigen::MatrixXd hooke = hookeVolume(E, nu) ; 
    # pragma omp parallel 
    {
        #pragma omp for 
            for (int i=0; i<sie; i++)
            {
                Eigen::MatrixXd Ke(3*nbf_elem,3*nbf_elem);   // Elementary stiffness matrix 
                Eigen::MatrixXd Be(6,3*nbf_elem);   // Elementary differential matrix 
 
                double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
                double detJ ; 
                double dNdx, dNdy, dNdz ; 
                Ke.fill(0) ;
                Be.fill(0) ;   
                // Loop over integration points 
                for (int ip=0; ip<nip; ip++)
                {
                    dxdxi   = 0 ; 
                    dxdeta  = 0 ; 
                    dxdzeta = 0 ; 
                    
                    dydxi   = 0 ; 
                    dydeta  = 0 ; 
                    dydzeta = 0 ; 
                    
                    dzdxi   =  0 ; 
                    dzdeta  =  0 ; 
                    dzdzeta =  0 ;  
                    // Setting the elementary differential operator
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
        			{
       					dxdxi   +=  dNdxi[ibf + ip*sjN]   * n[0+e[ibf+i*sje]*sjn] ; 
       					dxdeta  +=  dNdeta[ibf + ip*sjN]  * n[0+e[ibf+i*sje]*sjn] ;   
       					dxdzeta +=  dNdzeta[ibf + ip*sjN] * n[0+e[ibf+i*sje]*sjn] ; 
       					
       					dydxi   +=  dNdxi[ibf + ip*sjN]   * n[1+e[ibf+i*sje]*sjn] ; 
       					dydeta  +=  dNdeta[ibf + ip*sjN]  * n[1+e[ibf+i*sje]*sjn] ; 
       					dydzeta +=  dNdzeta[ibf + ip*sjN] * n[1+e[ibf+i*sje]*sjn] ;
       					
       					dzdxi   +=  dNdxi[ibf + ip*sjN]   * n[2+e[ibf+i*sje]*sjn] ;
       					dzdeta  +=  dNdeta[ibf + ip*sjN]  * n[2+e[ibf+i*sje]*sjn] ;  
       					dzdzeta +=  dNdzeta[ibf + ip*sjN] * n[2+e[ibf+i*sje]*sjn] ; 	    					            			 
        			}
   					// Determinant of the jacobian matrix 
   					detJ   = dxdxi*dydeta*dzdzeta + 
                             dydxi*dzdeta*dxdzeta + 
                             dzdxi*dxdeta*dydzeta - 
                             dzdxi*dydeta*dxdzeta - 
                             dydxi*dxdeta*dzdzeta - 
                             dxdxi*dzdeta*dydzeta ; 
                    // Setting the differential operator Be at the integration point 
                    for ( int ibf =0; ibf < nbf_elem; ibf++)
        			{
            			
                        dNdx  = ((dydeta*dzdzeta-dzdeta*dydzeta)/detJ)*dNdxi[ibf + ip*sjN] +  
                                       (-(dydxi*dzdzeta-dzdxi*dydzeta)/detJ)*dNdeta[ibf + ip*sjN] + 
                                       ((dydxi*dzdeta-dzdxi*dydeta)/detJ)*dNdzeta[ibf + ip*sjN] ; 
                                       
                		dNdy  = ((-(dxdeta*dzdzeta-dzdeta*dxdzeta))/detJ)*dNdxi[ibf + ip*sjN] +
                                       ((dxdxi*dzdzeta-dzdxi*dxdzeta)/detJ)*dNdeta[ibf + ip*sjN] +
                                       (-(dxdxi*dzdeta-dzdxi*dxdeta)/detJ)*dNdzeta[ibf + ip*sjN] ; 
                		
                		dNdz  = ((dxdeta*dydzeta-dydeta*dxdzeta)/detJ)*dNdxi[ibf + ip*sjN] +
                                       (-(dxdxi*dydzeta-dydxi*dxdzeta)/detJ)*dNdeta[ibf + ip*sjN] +
                                       ((dxdxi*dydeta-dydxi*dxdeta)/detJ)*dNdzeta[ibf + ip*sjN] ;   
                                       
                                    
                        Be( 0, ibf )           = dNdx;
        				Be( 1, ibf+nbf_elem)   = dNdy ;
        				Be( 2, ibf+2*nbf_elem) = dNdz ;
        				
        				Be( 3, ibf)           = dNdy ;
        				Be( 3, ibf+nbf_elem)  = dNdx;
        				
        				Be( 4, ibf)            = dNdz ;
        				Be( 4, ibf+2*nbf_elem) = dNdx;  
        				
        				Be( 5, ibf+nbf_elem)   = dNdz ;
        				Be( 5, ibf+2*nbf_elem) = dNdy ;   
            				  					            			 
        			}        
			        Ke = Ke + Be.transpose()*hooke*Be*(w[ip]*std::abs(detJ)) ;     
 
                }
                int I,J ; 
                int sg = 9*nbf_elem*nbf_elem * i ; 
                // Elementary contribution to the global matrix 
                for (int ibf=0; ibf<nbf_elem; ibf++)
                {
                    I = e[ibf+i*sje] ;
                    for (int jbf=0; jbf<nbf_elem; jbf++)
                    {
    					J = e[jbf+i*sje];
    					
       		            indexI[sg] = conn[0+I*sjconn]   ;
       		            indexJ[sg] = conn[0+J*sjconn]   ;
       		            nnz_values[sg] = Ke(ibf,jbf);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn]  ;
       		            indexJ[sg] = conn[0+J*sjconn]  ;
       		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[0+I*sjconn]  ;
       		            indexJ[sg] = conn[1+J*sjconn]  ;
       		            nnz_values[sg] = Ke(ibf,jbf+nbf_elem);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn] ;
       		            indexJ[sg] = conn[1+J*sjconn] ;
       		            nnz_values[sg] = Ke(ibf+nbf_elem,jbf+nbf_elem);
       		            sg ++ ;
       		                   
       		            indexI[sg] = conn[2+I*sjconn]  ;
       		            indexJ[sg] = conn[0+J*sjconn]  ;
       		            nnz_values[sg] = Ke(ibf+2*nbf_elem,jbf);
       		            sg ++ ;    		            
       		            
       		            indexI[sg] = conn[0+I*sjconn]  ;
       		            indexJ[sg] = conn[2+J*sjconn]  ;
       		            nnz_values[sg] =  Ke(ibf,jbf+2*nbf_elem);
       		            sg ++ ;
       		            
       		            indexI[sg] = conn[1+I*sjconn]  ;
       		            indexJ[sg] = conn[2+J*sjconn]  ;
       		            nnz_values[sg] =  Ke(ibf+nbf_elem,jbf+2*nbf_elem);
       		            sg ++ ;   		            
    
    		            indexI[sg] = conn[2+I*sjconn]   ;
       		            indexJ[sg] = conn[1+J*sjconn]   ;
       		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+nbf_elem);
       		            sg ++ ; 
       		               		            
   		                indexI[sg] = conn[2+I*sjconn]   ;
       		            indexJ[sg] = conn[2+J*sjconn]   ;
       		            nnz_values[sg] =  Ke(ibf+2*nbf_elem,jbf+2*nbf_elem);
       		            sg ++ ; 
       		            
                    }
                }                
            
            }
    }

}
 




// void Gophi_TrilinearInterpFEMesh(double* g, int sgi, int sgj, int sgk,
//                                  double* knotXiImage, int sknotXiImage, 
//                                  double* knotEtaImage, int sknotEtaImage, 
//                                  double* knotZetaImage, int sknotZetaImage, 
//                                  double* xi,  int sxi, 
//                                  double* eta, int seta, 
//                                  double* zeta, int szeta, 
//                                  int* ie, int size_ie,  
//                                  int* e, int sie, int sje, 
//                                  double* n, int sin, int sjn, 
//                                  int* conn, int siconn, int sjconn, 
//                                  double* U, int sU, 
//                                  double* glr, int sglr ) 
// {
//     #pragma omp parallel 
//     {
//          #pragma omp for 
//             for (int ip=0; ip< sxi; ip++)
//            {
//                  double N[4] = { xi[ip], eta[ip], zeta[ip], 1-xi[ip]-eta[ip]-zeta[ip] } ; 
//                  double xa = 0 ; 
//                  double ya = 0 ; 
//                  double za = 0 ; 
//                  int I ; 
//                  // Computing the advected voxels 
//                  // Computing the components of the jacobian of the isoparametric transformation 
//                  for (int ibf=0; ibf < 4; ibf++)
//                  {
//                      I =  ibf+  ie[ip]*sje   ; 
//                      xa += N[ibf] * ( n[0+e[I]*sjn] + U[conn[0+e[I]*sjconn]] ) ; 
//                      ya += N[ibf] * ( n[1+e[I]*sjn] + U[conn[1+e[I]*sjconn]] ) ; 
//                      za += N[ibf] * ( n[2+e[I]*sjn] + U[conn[2+e[I]*sjconn]] ) ; 
//                  } 
//                  // Gray-level residual f(x,y,z) - g(x+ux(x,y,z), y+uy(x,y,z), z+uz(x,y,z)) 
//                  glr[ip] = EvaluateTrilinearInterpolationOnOnePoint(g,sgi,sgj,sgk,
//                             knotXiImage, sknotXiImage, 
//                             knotEtaImage, sknotEtaImage,
//                             knotZetaImage, sknotZetaImage, 
//                             xa, ya , za ); 
//             }   
//     }
// 
// }
 
 

// void DVC_LHS_Tetra1(double* dfipdx, int sdfipdx, 
//                     double* dfipdy, int sdfipdy, 
//                     double* dfipdz, int sdfipdz,
//                     int* e, int sie, int sje, 
//                     int* conn, int siconn, int sjconn, 
//                     int* n, int sin, int sjn, 
//                     double* xi, int sxi, 
//                     double* eta, int seta, 
//                     double* zeta, int szeta, 
//                     double* w, int sw, 
//                     int* indexI, int nnzI, 
//                     int* indexJ, int nnzJ, 
//                     double* nnz_values, int nnz)
// {
//     int nbf_elem = 4 ;  // 4 fe nodes per element s
//     // Assembles Hessian correlation matrix 
//     // With the integration rule (x,y,z,w) defined by the user 
//     double* N       = array_2d(sx, nbf_elem) ; // 
//     double* dNdxi   = new double[nbf_elem] ;
//     double* dNdeta  = new double[nbf_elem] ; 
//     double* dNdzeta = new double[nbf_elem] ; 
//     
//     // Computing the basis functions at integration points 
//     for (int ig=0; ig< sx; ig++)
//     {   
//         N[ig][0] = xi[ig]   ; 
//         N[ig][1] = eta[ig]  ; 
//         N[ig][2] = zeta[ig] ; 
//         N[ig][3] = 1-xi[ig]-eta[ig]-zeta[ig] ; 
//     }
//     dNdxi[0] =  1  ; 
//     dNdxi[1] =  0  ; 
//     dNdxi[2] =  0  ; 
//     dNdxi[3] =  -1 ; 
//     
//     dNdeta[0] = 0  ; 
//     dNdeta[1] = 1  ; 
//     dNdeta[2] = 0  ; 
//     dNdeta[3] = -1  ; 
//     
//     dNdzeta[0] = 0  ; 
//     dNdzeta[1] = 0  ; 
//     dNdzeta[2] = 1  ; 
//     dNdzeta[3] = -1  ;   
//                 
//     //# pragma omp parallel 
//     {
//         //#pragma omp for 
//             for (int i=0; i<sie; i++)
//             {
//                 Eigen::MatrixXd He(3*nbf_elem,3*nbf_elem);
//                 Eigen::MatrixXd Ne(3,3*nbf_elem); 
//                 Eigen::MatrixXd gradFgradFt(3,3);
//                 double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4 ; 
//                 double dxdxi, dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta ;  
//                 
//                 // Current element 
//                 Ne.fill(0.);  
//                 He.fill(0.); // Elementary hessian matrix initialization      
//                 
//                 // Loop over integration points   
//                 int Ig ; 
//                 for (int ig=0; ig< sx; ig++)
//                 {
//                     // Seting basis functions and derivatives on the integration point 
//                     N[0] = xi[ig] ;   
//                     N[1] = eta[ig] ; 
//                     N[2] = zeta[ig] ; 
//                     N[3] = 1-xi[ig]-eta[ig]-zeta[ig] ; 
//                    
//                     
//                     // Setting the elemental basis function matrix 
//                     for ( int ibf =0; ibf < nbf_elem; ibf++)
//    					{
//        					Ne(0,ibf)  = N[ibf] ; 
//        					Ne(1,ibf+nbf_elem)   = N[ibf] ;  
//        					Ne(2,ibf+2*nbf_elem) = N[ibf] ; 
//    					}  
//    					//Setting the gradient tensor product 
//    					Ig = ig + i*sx ;    
//                     gradFgradFt(0,0) = dfipdx[Ig]*dfipdx[Ig] ; 
//                     gradFgradFt(0,1) = dfipdx[Ig]*dfipdy[Ig] ;
//                     gradFgradFt(0,2) = dfipdx[Ig]*dfipdz[Ig] ;
//                     gradFgradFt(1,0) = dfipdy[Ig]*dfipdx[Ig] ;
//                     gradFgradFt(1,1) = dfipdy[Ig]*dfipdy[Ig] ; 
//                     gradFgradFt(1,2) = dfipdy[Ig]*dfipdz[Ig] ;
//                     gradFgradFt(2,0) = dfipdz[Ig]*dfipdx[Ig] ;
//                     gradFgradFt(2,1) = dfipdz[Ig]*dfipdy[Ig] ;
//                     gradFgradFt(2,2) = dfipdz[Ig]*dfipdz[Ig] ;
//                     
//                     // Element nodes 
//                     x1 = n[][0] ; y1 = n[][1] ; z1 = n[][2] ; 
//                     x2 = n[][0] ; y2 = n[][1] ; z2 = n[][2] ; 
//                     x3 = n[][0] ; y3 = n[][1] ; z3 = n[][2] ; 
//                     x4 = n[][0] ; y4 = n[][1] ; z4 = n[][2] ;    
//                     
//                     // Elements of the Jacobian of the isoparametric transformation                 
//                     dxdxi   =  dNdxi[0]*x1  + dNdxi[1]*x2 + dNdxi[2]*x3 + dNdxi[3]*x4 ; 
//                     dxdeta  =  dNdeta[0]*x1 + dNdeta[1]*x2 + dNdeta[2]*x3 + dNdeta[3]*x4 ; 
//                     dxdzeta =  dNdzeta[0]*x1 + dNdzeta[1]*x2 + dNdzeta[2]*x3 + dNdzeta[3]*x4 ; 
// 
//                     dydxi   =  dNdxi[0]*y1  + dNdxi[1]*y2 + dNdxi[2]*y3 + dNdxi[3]*y4 ; 
//                     dydeta  =  dNdeta[0]*y1 + dNdeta[1]*y2 + dNdeta[2]*y3 + dNdeta[3]*y4 ; 
//                     dydzeta =  dNdzeta[0]*y1 + dNdzeta[1]*y2 + dNdzeta[2]*y3 + dNdzeta[3]*y4 ; 
//                     
//                     dxdxi   =  dNdxi[0]*x1  + dNdxi[1]*x2 + dNdxi[2]*x3 + dNdxi[3]*x4 ; 
//                     dxdeta  =  dNdeta[0]*x1 + dNdeta[1]*x2 + dNdeta[2]*x3 + dNdeta[3]*x4 ; 
//                     dxdzeta =  dNdzeta[0]*x1 + dNdzeta[1]*x2 + dNdzeta[2]*x3 + dNdzeta[3]*x4 ;                     
//                     
//                     
//                     // Adding to the voxel summation 
//    					He = He + w[ig]*(Ne.transpose()*gradFgradFt*Ne) ; 
//    					ig++ ;    										                     
//                 }
//                 int I,J ; 
//                 int sg = 9*nbf_elem*nbf_elem* i ;
//                 // Elementary contribution to the global Hessian matrix 
//                 for (int ibf=0; ibf<nbf_elem; ibf++)
//                 {
//                     I = e(i,ibf); 
//                     for (int jbf=0; jbf<nbf_elem; jbf++)
//                         {
//            					J = e(i,jbf);
//    					
//            		            indexI[sg] = I ;
//            		            indexJ[sg] = J ;
//            		            nnz_values[sg] = He(ibf,jbf);
//            		            sg ++ ;
//            		            
//            		            indexI[sg] = I+nbf ;
//            		            indexJ[sg] = J ;
//            		            nnz_values[sg] = He(ibf+nbf_elem,jbf);
//            		            sg ++ ;
//            		            
//        		                indexI[sg] = I ;
//            		            indexJ[sg] = J +nbf;
//            		            nnz_values[sg] = He(ibf,jbf+nbf_elem);
//            		            sg ++ ;
//            		            
//            		            indexI[sg] = I+nbf ;
//            		            indexJ[sg] = J+nbf ;
//            		            nnz_values[sg] = He(ibf+nbf_elem,jbf+nbf_elem);
//            		            sg ++ ;
//            		                   
//            		            indexI[sg] = I+2*nbf ;
//            		            indexJ[sg] = J ;
//            		            nnz_values[sg] = He(ibf+2*nbf_elem,jbf);
//            		            sg ++ ;    		            
//            		            
//            		            indexI[sg] = I ;
//            		            indexJ[sg] = J+2*nbf ;
//            		            nnz_values[sg] =  He(ibf,jbf+2*nbf_elem);
//            		            sg ++ ;
//            		            
//            		            indexI[sg] = I+nbf ;
//            		            indexJ[sg] = J+2*nbf ;
//            		            nnz_values[sg] =  He(ibf+nbf_elem,jbf+2*nbf_elem);
//            		            sg ++ ;   		            
// 
//         		            indexI[sg] = I+2*nbf ;
//            		            indexJ[sg] = J+nbf ;
//            		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+nbf_elem);
//            		            sg ++ ; 
//            		               		            
//        		                indexI[sg] = I+2*nbf ;
//            		            indexJ[sg] = J+2*nbf ;
//            		            nnz_values[sg] =  He(ibf+2*nbf_elem,jbf+2*nbf_elem);
//            		            sg ++ ; 
//                        }
//                }                
//             }
//     }
// }


                                
                                                                    
