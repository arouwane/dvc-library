#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include "Eigen/Core"
#include <vector>

class SharedParameters
{
    public: 
    // Image parameters 
    //***************************
    double* image ; // Voxels vector 
    int si ; // Image size 
    int sj ;
    int sk ; 
    double* knotXiImage   ; 
    double* knotEtaImage  ; 
    double* knotZetaImage ; 
    int sknotXiImage   ; 
    int sknotEtaImage  ;
    int sknotZetaImage ; 
    double thrsh ;   
    //***************************  
    // B-spline mesh paramaters 
    int deg_xi     ; 
    int deg_eta    ; 
    int deg_zeta   ; 
    double* knotXi ; 
    double* knotEta  ;
    double* knotZeta ;     
    int sknotXi   ; 
    int sknotEta  ;
    int sknotZeta ;
    int nbg_xi   ; 
	int nbg_eta  ;
	int nbg_zeta ; 
	int nbf_elem_xi   ; 
	int nbf_elem_eta  ; 
	int nbf_elem_zeta ; 
	int nbf_elem  ; 
	double* xig   ; 
	double* etag  ; 
	double* zetag ; 
	double* wg    ;	
	Eigen::MatrixXi NOELEM; // Element dof connectivity 
	//************************
	// Integration of cut cells 
	// Integration rule on tetrahedrons 
	int nbgTetrahedron ; 
	Eigen::MatrixXd GaussTetraBc    ;   // Barycentric coordinates of Gauss points  
	std::vector<double> GaussTetraW ;  // Gauss weights (that must be multiplied by the element volume)   
	// Material parameter 
	Eigen::MatrixXd hooke; // Hooke matrix  
	
	SharedParameters(double* image, int si, int sj, int sk, 
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
                     Eigen::MatrixXd hooke);                      
}; 

class PrivateParameters
{
    public:
    // Element (thread) dependant parameters  
    double* Cxig  ;
	double* Cetag ;    
	double* Czetag ;
	double** bfOutputXi   ;
    double** bfOutputEta  ;
    double** bfOutputZeta ;	
    double** Nxi; 
    double** dNxidxi; 
    double** Neta; 
    double** dNetadeta; 
    double** Nzeta; 
    double** dNzetadzeta; 
	double* N; 
	double* dNdxi; 
	double* dNdeta; 
	double* dNdzeta;    
	Eigen::MatrixXd* Ke ;  // Elementary stiffness matrix 
    Eigen::MatrixXd* Bc ;  // Differential operator defined on a an integration cell ( An element contains multiple cells)  
    
    PrivateParameters( double* Cxig, double* Cetag, double* Czetag,  
                       double** bfOutputXi, double** bfOutputEta, double** bfOutputZeta, 
                       double** Nxi, double** dNxidxi, double** Neta, double** dNetadeta,  double** Nzeta,  double** dNzetadzeta,  
                       double* N, double* dNdxi, double* dNdeta, double* dNdzeta,
                       Eigen::MatrixXd* Ke, Eigen::MatrixXd* Bc );    
}; 


class MeshParameters
{
    public: 
    // Image parameters 
    //***************************
    double* image; // voxel vector  
    int si;  // Image size 
    int sj; 
    int sk; 
    double* knotXiImage   ; 
    double* knotEtaImage  ; 
    double* knotZetaImage ; 
    int sknotXiImage   ; 
    int sknotEtaImage  ;
    int sknotZetaImage ; 
    double thrsh ; 
    //double* NxImage; 
    //double* NyImage; 
    //double* NzImage; 
    //***************************
    // B-spline mesh paramaters 
    int deg_xi   ; 
    int deg_eta  ; 
    int deg_zeta ;   
    double* knotXi   ; 
    double* knotEta  ;
    double* knotZeta ; 
    int sknotXi   ; 
    int sknotEta  ;
    int sknotZeta ;
    int nbg_xi   ; 
	int nbg_eta  ;
	int nbg_zeta ; 
	int nbf_elem_xi   ; 
	int nbf_elem_eta  ; 
	int nbf_elem_zeta ;  
	int nbf_elem ; 
	double* xig   ; 
	double* etag  ; 
	double* zetag ; 
	double* wg    ;
	double* Cxig  ;
	double* Cetag ;    
	double* Czetag ;  
	double** bfOutputXi   ;
    double** bfOutputEta  ;
    double** bfOutputZeta ;
    double** Nxi; 
    double** dNxidxi; 
    double** Neta; 
    double** dNetadeta; 
    double** Nzeta; 
    double** dNzetadzeta; 
	double* N; 
	double* dNdxi; 
	double* dNdeta; 
	double* dNdzeta;   
	int nbgTetrahedron ; 
	Eigen::MatrixXd GaussTetraBc    ;   // Barycentric coordinates of Gauss points  
	std::vector<double> GaussTetraW ;  // Gauss weights (that must be multiplied by the element volume)
	Eigen::MatrixXd* Ke ;  // Elementary stiffness matrix 
    Eigen::MatrixXd* Bc ;  // Differential operator defined on a an integration cell
    Eigen::MatrixXi NOELEM; // Element dof connectivity 
    Eigen::MatrixXd hooke; // Hooke matrix 


	MeshParameters(double* image, int si, int sj, int sk, 
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
                   Eigen::MatrixXi NOELEM, Eigen::MatrixXd hooke);    
	
};




#endif
