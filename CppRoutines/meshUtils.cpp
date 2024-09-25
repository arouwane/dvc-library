#include "meshUtils.h" 


SharedParameters::SharedParameters(double* image, int si, int sj, int sk, 
                     double* knotXiImage, int sknotXiImage, 
                     double* knotEtaImage, int sknotEtaImage, 
                     double* knotZetaImage, int sknotZetaImage, 	
                	 double thrsh, 
	                 int deg_xi, int deg_eta, int deg_zeta, 
                     double* knotXi, int sknotXi, 
                     double* knotEta, int sknotEta, 
                     double* knotZeta, int sknotZeta,
                     int nbg_xi, int nbg_eta, int nbg_zeta, 
                     double* xig, double* etag, double* zetag, double* wg, 
                     Eigen::MatrixXi NOELEM, 
                     int nbgTetrahedron, Eigen::MatrixXd GaussTetraBc, std::vector<double> GaussTetraW, 
                     Eigen::MatrixXd hooke)
{
    this->image = image; 
    this->si = si; 
    this->sj = sj; 
    this->sk = sk;  
    this->knotXiImage  = knotXiImage  ;
    this->sknotXiImage = sknotXiImage ;  
    this->knotEtaImage = knotEtaImage ;  
    this->sknotEtaImage = sknotEtaImage ;   
    this->knotZetaImage = knotZetaImage ;  
    this->sknotZetaImage =  sknotZetaImage ;  
    this->thrsh = thrsh ;     
    this->deg_xi  = deg_xi  ;  
    this->deg_eta = deg_eta ;  
    this->deg_zeta =  deg_zeta ;  
    this->nbf_elem_xi   = deg_xi   + 1 ; 
    this->nbf_elem_eta  = deg_eta  + 1 ;  
    this->nbf_elem_zeta = deg_zeta + 1 ; 
    this->nbf_elem  = (deg_xi   + 1)*(deg_eta  + 1)*(deg_zeta + 1) ; 
    this->knotXi   = knotXi   ;    
    this->knotEta  = knotEta  ;  
    this->knotZeta = knotZeta ;  
    this->sknotXi   = sknotXi   ; 
    this->sknotEta  = sknotEta  ;
    this->sknotZeta = sknotZeta ;  
    this->nbg_xi   = nbg_xi   ;  
    this->nbg_eta  = nbg_eta ;   
    this->nbg_zeta = nbg_zeta ;  
    this->xig    = xig   ;  
    this->etag   = etag  ;  
    this->zetag  = zetag ;  
    this->wg     = wg ;  
    this->nbgTetrahedron =  nbgTetrahedron ; 
    this->GaussTetraBc = GaussTetraBc ; 
    this->GaussTetraW  = GaussTetraW ; 
    this->NOELEM = NOELEM ; 
    this->hooke  = hooke ;  
} 

PrivateParameters::PrivateParameters( double* Cxig, double* Cetag, double* Czetag,  
                   double** bfOutputXi, double** bfOutputEta, double** bfOutputZeta, 
                   double** Nxi, double** dNxidxi, double** Neta, double** dNetadeta,  double** Nzeta,  double** dNzetadzeta,  
                   double* N, double* dNdxi, double* dNdeta, double* dNdzeta,
                   Eigen::MatrixXd* Ke, Eigen::MatrixXd* Bc )
{
    this->Cxig   = Cxig ; 
    this->Cetag  = Cetag ; 
    this->Czetag = Czetag ;  
    this->bfOutputXi   = bfOutputXi   ;  
    this->bfOutputEta  = bfOutputEta  ;
    this->bfOutputZeta = bfOutputZeta ; 
    this->Nxi = Nxi ; 
    this->dNxidxi = dNxidxi ;  
    this->Neta = Neta ;   
    this->dNetadeta = dNetadeta ;    
    this->Nzeta = Nzeta ;   
    this->dNzetadzeta = dNzetadzeta ;   
    this->N = N ;  
    this->dNdxi   = dNdxi  ;   
    this->dNdeta  = dNdeta ;   
    this->dNdzeta = dNdzeta ;  
    this->Ke = Ke ; 
    this->Bc = Bc ; 
}
                       

MeshParameters::MeshParameters(double* image, int si, int sj, int sk, 
                   double* knotXiImage, int sknotXiImage, 
                   double* knotEtaImage, int sknotEtaImage, 
                   double* knotZetaImage, int sknotZetaImage, 
                   double thrsh, 
                   int deg_xi, int deg_eta, int deg_zeta, 
                   double* knotXi, int sknotXi, 
                   double* knotEta, int sknotEta, 
                   double* knotZeta, int sknotZeta, 
                   int nbg_xi, int nbg_eta, int nbg_zeta, 
                   double* xig, double* etag, double* zetag, double* wg, 
                   double* Cxig, double* Cetag, double* Czetag,  
                   double** bfOutputXi, double** bfOutputEta, double** bfOutputZeta, 
                   double** Nxi, double** dNxidxi, double** Neta, double** dNetadeta,  double** Nzeta,  double** dNzetadzeta,  
                   double* N, double* dNdxi, double* dNdeta, double* dNdzeta, 
                   int nbgTetrahedron, Eigen::MatrixXd GaussTetraBc, std::vector<double> GaussTetraW,
                   Eigen::MatrixXd* Ke, Eigen::MatrixXd* Bc, 
                   Eigen::MatrixXi NOELEM, Eigen::MatrixXd hooke)
{
    this->image = image; 
    this->si = si; 
    this->sj = sj; 
    this->sk = sk;  
    this->knotXiImage  = knotXiImage  ;
    this->sknotXiImage = sknotXiImage ;  
    this->knotEtaImage = knotEtaImage ;  
    this->sknotEtaImage = sknotEtaImage ;   
    this->knotZetaImage = knotZetaImage ;  
    this->sknotZetaImage =  sknotZetaImage ;  
    this->thrsh = thrsh ;   
    //this->NxImage = NxImage ;  
    //this->NyImage = NyImage ;  
    //this->NzImage = NzImage ;  
    this->deg_xi  = deg_xi  ;  
    this->deg_eta = deg_eta ;  
    this->deg_zeta =  deg_zeta ;  
    this->nbf_elem_xi   = deg_xi   + 1 ; 
    this->nbf_elem_eta  = deg_eta  + 1 ;  
    this->nbf_elem_zeta = deg_zeta + 1 ; 
    this->nbf_elem  = (deg_xi   + 1)*(deg_eta  + 1)*(deg_zeta + 1) ; 
    this->knotXi   = knotXi   ;    
    this->knotEta  = knotEta  ;  
    this->knotZeta = knotZeta ;  
    this->sknotXi   = sknotXi   ; 
    this->sknotEta  = sknotEta  ;
    this->sknotZeta = sknotZeta ;  
    this->nbg_xi   = nbg_xi   ;  
    this->nbg_eta  = nbg_eta ;   
    this->nbg_zeta = nbg_zeta ;   
    this->xig    = xig   ;  
    this->etag   = etag  ;  
    this->zetag  = zetag ;  
    this->wg     = wg ;  
    this->Cxig   = Cxig ; 
    this->Cetag  = Cetag ; 
    this->Czetag = Czetag ;  
    this->bfOutputXi   = bfOutputXi   ;  
    this->bfOutputEta  = bfOutputEta  ;
    this->bfOutputZeta = bfOutputZeta ; 
    this->Nxi = Nxi ; 
    this->dNxidxi = dNxidxi ;  
    this->Neta = Neta ;   
    this->dNetadeta = dNetadeta ;    
    this->Nzeta = Nzeta ;   
    this->dNzetadzeta = dNzetadzeta ;   
    this->N = N ;  
    this->dNdxi   = dNdxi  ;   
    this->dNdeta  = dNdeta ;   
    this->dNdzeta = dNdzeta ;  
    this->nbgTetrahedron =  nbgTetrahedron ; 
    this->GaussTetraBc = GaussTetraBc ; 
    this->GaussTetraW  = GaussTetraW ; 
    this->Ke = Ke ; 
    this->Bc = Bc ; 
    this->NOELEM = NOELEM ; 
    this->hooke  = hooke ;  
}
