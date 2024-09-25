#include "octreeDecomposition.h"


Cell::Cell(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    this->xmin = xmin; 
    this->xmax = xmax; 
    this->ymin = ymin; 
    this->ymax = ymax; 
    this->zmin = zmin; 
    this->zmax = zmax; 
    this->lvl = 0; 
    this->mesx = xmax-xmin; 
    this->mesy = ymax-ymin; 
    this->mesz = zmax-zmin;  
}


void Cell::DecomposeTesselationTrilinearInterp(int lvlmax, MeshParameters* meshParam, std::vector<double>* integrationPoints ) 
{
    int t;
    if (this->lvl == lvlmax)
    {   
         // First checking if the cuboid corners are > than threshold value
         // in this case, we perform an integration on a the square cell
         double eps = 1.e-8;      
         std::array<double,8> xc =  { this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps, this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps }; 
         std::array<double,8> yc =  { this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps, this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps };
         std::array<double,8> zc =  { this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmax-eps, this->zmax-eps, this->zmax-eps, this->zmax-eps };
         std::array<double,8> gl; 
         for (int i=0; i<8; i++)
         {
            gl[i] = EvaluateTrilinearInterpolationOnOnePoint(meshParam->image,meshParam->si,meshParam->sj,meshParam->sk, 
           					        meshParam->knotXiImage, meshParam->sknotXiImage,   
           					        meshParam->knotEtaImage, meshParam->sknotEtaImage, 
           					        meshParam->knotZetaImage, meshParam->sknotZetaImage,  
       					                xc[i], yc[i], zc[i]); 		                
         } 
         double thrsh = meshParam->thrsh ; 
         if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g>thrsh ;}))
         {
             // Integration the cell 
             
             // Setting Gauss integration coordinates 
             for (int ii =0; ii< meshParam->nbg_xi; ii++)
             {
               meshParam->Cxig[ii] = this->xmin + 0.5*(meshParam->xig[ii]+1)*this->mesx;
             }    
             for (int ii =0; ii< meshParam->nbg_eta; ii++)
             {
               meshParam->Cetag[ii] = this->ymin + 0.5*(meshParam->etag[ii]+1)*this->mesy;
             }  
             for (int ii =0; ii< meshParam->nbg_zeta; ii++)
             {
               meshParam->Czetag[ii] = this->zmin + 0.5*(meshParam->zetag[ii]+1)*this->mesz;
             }              
             t = 0; 
             for (int kk=0; kk<meshParam->nbg_zeta; kk++)
             {
                 for (int jj=0; jj<meshParam->nbg_eta; jj++)
                 {
                     for (int ii=0; ii<meshParam->nbg_xi; ii++)
                     {
                         // Adding the integration point to the total list of integration points 
                         integrationPoints->push_back(meshParam->Cxig[ii]); 
                         integrationPoints->push_back(meshParam->Cetag[jj]); 
                         integrationPoints->push_back(meshParam->Czetag[kk]); 
                         integrationPoints->push_back( meshParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ); 
                         t++; 
                     }
                 }
             }
         
         }
         else if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g< thrsh;})==false)
         {
             // Perform Delaunay triangulation on the last cut-cell 
             double a,b ; // Level set linear interpolation coefficients for the tesselation approximation
             double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4; // tetrahedron nodes 
             double volT;  // Volume of the terahedron   
             std::vector<Point> tesselationPoints;
             for (int i=0; i<8; i++)
             {
                 if (gl[i] > meshParam->thrsh)
                 {
                     tesselationPoints.push_back(Point(xc[i],yc[i],zc[i]));
                 }
             }
             // Distance linearization of the interface 
             // Edge 1  0-1 - x    
             if (( gl[0] > meshParam->thrsh && gl[1] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[1] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[0],xc[1],gl[0],gl[1], a,b); 
                 tesselationPoints.push_back( Point( (meshParam->thrsh-b)/a, this->ymin, this->zmin  )); 
             }
             // Edge 2  1-2 - y
             if (( gl[1] > meshParam->thrsh && gl[2] < meshParam->thrsh ) || ( gl[1] < meshParam->thrsh && gl[2] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[1],yc[2],gl[1],gl[2],a,b);
                 tesselationPoints.push_back( Point(this->xmax,(meshParam->thrsh-b)/a,this->zmin ));  
             }         
             // Edge 3  2-3 - x
             if (( gl[2] > meshParam->thrsh && gl[3] < meshParam->thrsh ) || ( gl[2] < meshParam->thrsh && gl[3] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[2],xc[3],gl[2],gl[3],a,b);
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a,this->ymax,this->zmin ));  
             }
             // Edge 4  0-3 - y
             if (( gl[0] > meshParam->thrsh && gl[3] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[3] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[0],yc[3],gl[0],gl[3],a,b);  
                 tesselationPoints.push_back( Point(this->xmin,(meshParam->thrsh-b)/a,this->zmin));  
             }
             // Edge 5  4-5 - x
             if (( gl[4] > meshParam->thrsh && gl[5] < meshParam->thrsh ) || ( gl[4] < meshParam->thrsh && gl[5] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[4],xc[5],gl[4],gl[5],a,b);  
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a,this->ymin,this->zmax)); 
             }         
             // Edge 6  5-6 - y
             if (( gl[5] > meshParam->thrsh && gl[6] < meshParam->thrsh ) || ( gl[5] < meshParam->thrsh && gl[6] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[5],yc[6],gl[5],gl[6],a,b); 
                 tesselationPoints.push_back( Point(this->xmax,(meshParam->thrsh-b)/a,this->zmax));
             }         
             // Edge 7  6-7 - x
             if (( gl[6] > meshParam->thrsh && gl[7] < meshParam->thrsh ) || ( gl[6] < meshParam->thrsh && gl[7] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[6],xc[7],gl[6],gl[7],a,b); 
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a, this->ymax, this->zmax )); 
             }
             // Edge 8  7-4 - y
             if (( gl[7] > meshParam->thrsh && gl[4] < meshParam->thrsh ) || ( gl[7] < meshParam->thrsh && gl[4] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[7],yc[4],gl[7],gl[4],a,b);
                 tesselationPoints.push_back( Point(this->xmin,(meshParam->thrsh-b)/a,this->zmax));  
             }
             // Edge 9  0-4 - z
             if (( gl[0] > meshParam->thrsh && gl[4] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[4] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[0],zc[4],gl[0],gl[4],a,b); 
                 tesselationPoints.push_back( Point(this->xmin, this->ymin , (meshParam->thrsh-b)/a  ));  
             
             }
             // Edge 10 1-5 - z
             if (( gl[1] > meshParam->thrsh && gl[5] < meshParam->thrsh ) || ( gl[1] < meshParam->thrsh && gl[5] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[1],zc[5],gl[1],gl[5],a,b); 
                 tesselationPoints.push_back( Point(this->xmax, this->ymin, (meshParam->thrsh-b)/a ));  
             }
             // Edge 11  2-6 -z
             if (( gl[2] > meshParam->thrsh && gl[6] < meshParam->thrsh ) || ( gl[2] < meshParam->thrsh && gl[6] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[2],zc[6],gl[2],gl[6],a,b); 
                 tesselationPoints.push_back( Point( this->xmax, this->ymax, (meshParam->thrsh-b)/a)); 
             }
             // Edge 12  3-7 -z 
             if (( gl[3] > meshParam->thrsh && gl[7] < meshParam->thrsh ) || ( gl[3] < meshParam->thrsh && gl[7] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[3],zc[7],gl[3],gl[7],a,b);
                 tesselationPoints.push_back( Point( this->xmin, this->ymax, (meshParam->thrsh-b)/a ));  
             }
             Delaunay dt;
             dt.insert( tesselationPoints.begin(), tesselationPoints.end() );
             // Integration over the tetrahedrons   
             for (Delaunay::Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it)
             {
                 x1 = dt.tetrahedron(it)[0].x(); y1= dt.tetrahedron(it)[0].y(); z1 = dt.tetrahedron(it)[0].z(); 
                 x2 = dt.tetrahedron(it)[1].x(); y2= dt.tetrahedron(it)[1].y(); z2 = dt.tetrahedron(it)[1].z(); 
                 x3 = dt.tetrahedron(it)[2].x(); y3= dt.tetrahedron(it)[2].y(); z3 = dt.tetrahedron(it)[2].z(); 
                 x4 = dt.tetrahedron(it)[3].x(); y4= dt.tetrahedron(it)[3].y(); z4 = dt.tetrahedron(it)[3].z(); 
 
                 //Getting the volume of the tetrahedron 
                 volT = (1./6.)*( std::abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                            (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                            (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) ) ;
                            
                 // Performing the integration over the tetrahedron   
                 // Multiply the barycentric coordinates matrix by the tetrahedron node coordinates    
                 Eigen::MatrixXd elemCoord(4,3); 
                 elemCoord(0,0)=x1; elemCoord(0,1)=y1; elemCoord(0,2)=z1;
                 elemCoord(1,0)=x2; elemCoord(1,1)=y2; elemCoord(1,2)=z2;
                 elemCoord(2,0)=x3; elemCoord(2,1)=y3; elemCoord(2,2)=z3;
                 elemCoord(3,0)=x4; elemCoord(3,1)=y4; elemCoord(3,2)=z4;
        
                 Eigen::MatrixXd pg = (meshParam->GaussTetraBc)*elemCoord ; 
                 for (int ig=0; ig<meshParam->nbgTetrahedron; ig++)
                 {
                     // Appending the integration points and weights 
                     integrationPoints->push_back(pg(ig,0)); 
                     integrationPoints->push_back(pg(ig,1)); 
                     integrationPoints->push_back(pg(ig,2)); 
                     integrationPoints->push_back( meshParam->GaussTetraW[ig]*volT );  
                 }                     
             }
         }
     }    
     else 
     {
         int nbp = std::pow(2,lvlmax-this->lvl)+1; // For the subdivision of the the boundary of a subcell
         double eps = 1.e-8; 
 		 Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nbp, this->xmin + eps, this->xmax - eps);
 		 Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(nbp, this->ymin + eps, this->ymax - eps);     
 		 Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(nbp, this->zmin + eps, this->zmax - eps);    
         int countZero = CheckCutCuboidTrilinearInterpolation(meshParam->image,meshParam->si,meshParam->sj,meshParam->sk, 
            					        meshParam->knotXiImage, meshParam->sknotXiImage,   
            					        meshParam->knotEtaImage, meshParam->sknotEtaImage, 
            					        meshParam->knotZetaImage, meshParam->sknotZetaImage, 
            					        meshParam->thrsh, x, y, z); 
         int nborder = 6*nbp*nbp; // Total number of sampling points on the boundary of the cube 
         if (countZero==0)
         {
             // Setting Gauss integration coordinates 
             for (int ii =0; ii< meshParam->nbg_xi; ii++)
             {
               meshParam->Cxig[ii] = this->xmin + 0.5*(meshParam->xig[ii]+1)*this->mesx;
             }    
             for (int ii =0; ii< meshParam->nbg_eta; ii++)
             {
               meshParam->Cetag[ii] = this->ymin + 0.5*(meshParam->etag[ii]+1)*this->mesy;
             }  
             for (int ii =0; ii< meshParam->nbg_zeta; ii++)
             {
               meshParam->Czetag[ii] = this->zmin + 0.5*(meshParam->zetag[ii]+1)*this->mesz;
             }              
             t = 0; 
             for (int kk=0; kk<meshParam->nbg_zeta; kk++)
             {
                 for (int jj=0; jj<meshParam->nbg_eta; jj++)
                 {
                     for (int ii=0; ii<meshParam->nbg_xi; ii++)
                     {
                         // Adding the integration point to the total list of integration points 
                         integrationPoints->push_back(meshParam->Cxig[ii]); 
                         integrationPoints->push_back(meshParam->Cetag[jj]); 
                         integrationPoints->push_back(meshParam->Czetag[kk]); 
                         integrationPoints->push_back( meshParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ); 
                         t++; 
                     }
                 }
             }               
         }
         else if (countZero != nborder)
         {
             // Decompose the cut cell 
            this->bottomBackLeft   = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomBackRight  = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontLeft  = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontRight = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->topBackLeft      = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topBackRight     = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontLeft     = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontRight    = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            
            this->bottomBackLeft->lvl = this->lvl + 1 ;
            this->bottomBackRight->lvl = this->lvl + 1 ;
            this->bottomFrontLeft->lvl = this->lvl + 1 ;
            this->bottomFrontRight->lvl = this->lvl + 1 ;
            this->topBackLeft->lvl = this->lvl + 1 ;
            this->topBackRight->lvl = this->lvl + 1 ;
            this->topFrontLeft->lvl = this->lvl + 1 ;
            this->topFrontRight->lvl = this->lvl + 1  ;
            
            

            this->bottomBackLeft->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->bottomBackRight->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->bottomFrontLeft->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->bottomFrontRight->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->topBackLeft->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->topBackRight->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->topFrontLeft->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints ); 
            this->topFrontRight->DecomposeTesselationTrilinearInterp(lvlmax, meshParam, integrationPoints );                
         }		        
     }
     
 
}    




void Cell::DecomposeTesselationTrilinearInterpAssembly(int ie, int spanx, int spany, int spanz, 
                                                     int lvlmax, MeshParameters* meshParam, 
                                                     double* intBf) 
{
    int t   ;
    int bfc ; 
    double wgMes ; 
    if (this->lvl == lvlmax)
    {   
         // First checking if the cuboid corners are > than threshold value
         // in this case, we perform an integration on a the square cell
         double eps = 1.e-8;      
         std::array<double,8> xc =  { this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps, this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps }; 
         std::array<double,8> yc =  { this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps, this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps };
         std::array<double,8> zc =  { this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmax-eps, this->zmax-eps, this->zmax-eps, this->zmax-eps };
         std::array<double,8> gl; 
         for (int i=0; i<8; i++)
         {
            gl[i] = EvaluateTrilinearInterpolationOnOnePoint(meshParam->image,meshParam->si,meshParam->sj,meshParam->sk, 
           					        meshParam->knotXiImage, meshParam->sknotXiImage,   
           					        meshParam->knotEtaImage, meshParam->sknotEtaImage, 
           					        meshParam->knotZetaImage, meshParam->sknotZetaImage,  
       					                xc[i], yc[i], zc[i]); 		                
         } 
         double thrsh = meshParam->thrsh ; 
         if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g>thrsh ;}))
         {
             // Integration the cell 
             
             // Setting Gauss integration coordinates
             // And evaluating the univariate basis functions  
             for (int ii =0; ii< meshParam->nbg_xi; ii++)
             {
               meshParam->Cxig[ii] = this->xmin + 0.5*(meshParam->xig[ii]+1)*this->mesx;
               dersbasisfuns(meshParam->deg_xi, meshParam->knotXi, meshParam->Cxig[ii], spanx, 1, meshParam->bfOutputXi ); 
               for (int jbf=0; jbf<meshParam->nbf_elem_xi; jbf++)
               {
                    meshParam->Nxi[ii][jbf] 	= meshParam->bfOutputXi[0][jbf];
        			meshParam->dNxidxi[ii][jbf] = meshParam->bfOutputXi[1][jbf];   
               }
             }    
             for (int ii =0; ii< meshParam->nbg_eta; ii++)
             {
               meshParam->Cetag[ii] = this->ymin + 0.5*(meshParam->etag[ii]+1)*this->mesy;
               dersbasisfuns(meshParam->deg_eta, meshParam->knotEta,meshParam->Cetag[ii], spany, 1, meshParam->bfOutputEta);
          	   for ( int jbf =0; jbf<meshParam->nbf_elem_eta; jbf++)
               {
              	   meshParam->Neta[ii][jbf] 	   = meshParam->bfOutputEta[0][jbf];
              	   meshParam->dNetadeta[ii][jbf]   = meshParam->bfOutputEta[1][jbf];
          	   }
             }  
             for (int ii =0; ii< meshParam->nbg_zeta; ii++)
             {
               meshParam->Czetag[ii] = this->zmin + 0.5*(meshParam->zetag[ii]+1)*this->mesz;
               dersbasisfuns(meshParam->deg_zeta, meshParam->knotZeta,meshParam->Czetag[ii], spanz, 1, meshParam->bfOutputZeta);
               for ( int jbf =0; jbf<meshParam->nbf_elem_zeta; jbf++)
               {
              	   meshParam->Nzeta[ii][jbf] 	    = meshParam->bfOutputZeta[0][jbf];
              	   meshParam->dNzetadzeta[ii][jbf]  = meshParam->bfOutputZeta[1][jbf];
               } 
             }        
             // Loop over the integration points of the cell
             // Adding the contribution of the cell       
             t = 0; 
             for (int kk=0; kk<meshParam->nbg_zeta; kk++)
             {
                 for (int jj=0; jj<meshParam->nbg_eta; jj++)
                 {
                     for (int ii=0; ii<meshParam->nbg_xi; ii++)
                     {  
                         // Getting the weighting of the integration point 
                         wgMes = meshParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ; 
                         
                         // Getting the trivariate basis functions 
                         bfc = 0 ; 
                         for (int kbf=0; kbf< meshParam->nbf_elem_zeta; kbf++)
                         {
                             for (int jbf=0; jbf < meshParam->nbf_elem_eta; jbf++)
                             {
                                 for (int ibf=0; ibf < meshParam->nbf_elem_xi; ibf++)
                                 {
                                     meshParam->N[bfc]       = meshParam->Nxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->Nzeta[kk][kbf]       ; 
                                     meshParam->dNdxi[bfc]   = meshParam->dNxidxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->Nzeta[kk][kbf]   ; 
                                     meshParam->dNdeta[bfc]  = meshParam->Nxi[ii][ibf]*meshParam->dNetadeta[jj][jbf]*meshParam->Nzeta[kk][kbf]  ; 
                                     meshParam->dNdzeta[bfc] = meshParam->Nxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->dNzetadzeta[kk][kbf] ; 
                                     
                                     // Add the integral value of the basis function 
                                     intBf[meshParam->NOELEM(ie,bfc)] +=  wgMes * meshParam->N[bfc] ; 
                                     bfc++; 
                                 }
                             }
                         }
                         
                         // Set the cell differential matrix 
                         for (int ibf=0; ibf< meshParam->nbf_elem; ibf++)
                         {
            				(*meshParam->Bc)( 0, ibf )                      = meshParam->dNdxi[ibf];
            				(*meshParam->Bc)( 1, ibf+meshParam->nbf_elem)   = meshParam->dNdeta[ibf] ;
            				(*meshParam->Bc)( 2, ibf+2*meshParam->nbf_elem) = meshParam->dNdzeta[ibf] ;
            				
            				(*meshParam->Bc)( 3, ibf)                       = meshParam->dNdeta[ibf] ;
            				(*meshParam->Bc)( 3, ibf+meshParam->nbf_elem)   = meshParam->dNdxi[ibf];
            				
            				(*meshParam->Bc)( 4, ibf)                       = meshParam->dNdzeta[ibf] ;
            				(*meshParam->Bc)( 4, ibf+2*meshParam->nbf_elem) = meshParam->dNdxi[ibf];  
            				
            				(*meshParam->Bc)( 5, ibf+meshParam->nbf_elem)   = meshParam->dNdzeta[ibf] ;
            				(*meshParam->Bc)( 5, ibf+2*meshParam->nbf_elem) = meshParam->dNdeta[ibf] ;  
                          
                         }
                         // Add contribution to the stiffness 
                         *(meshParam->Ke) = *(meshParam->Ke) + meshParam->Bc->transpose() * meshParam->hooke * (*meshParam->Bc) * wgMes ;
                         t++; 
                     }
                 }
             }
         }
         else if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g< thrsh;})==false)
         {
             // Perform Delaunay triangulation on the last cut-cell 
             double a,b ; // Level set linear interpolation coefficients for the tesselation approximation
             double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4; // tetrahedron nodes 
             double volT;  // Volume of the terahedron   
             std::vector<Point> tesselationPoints;
             for (int i=0; i<8; i++)
             {
                 if (gl[i] > meshParam->thrsh)
                 {
                     tesselationPoints.push_back(Point(xc[i],yc[i],zc[i]));
                 }
             }
             // Distance linearization of the interface 
             // Edge 1  0-1 - x    
             if (( gl[0] > meshParam->thrsh && gl[1] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[1] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[0],xc[1],gl[0],gl[1], a,b); 
                 tesselationPoints.push_back( Point( (meshParam->thrsh-b)/a, this->ymin, this->zmin  )); 
             }
             // Edge 2  1-2 - y
             if (( gl[1] > meshParam->thrsh && gl[2] < meshParam->thrsh ) || ( gl[1] < meshParam->thrsh && gl[2] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[1],yc[2],gl[1],gl[2],a,b);
                 tesselationPoints.push_back( Point(this->xmax,(meshParam->thrsh-b)/a,this->zmin ));  
             }         
             // Edge 3  2-3 - x
             if (( gl[2] > meshParam->thrsh && gl[3] < meshParam->thrsh ) || ( gl[2] < meshParam->thrsh && gl[3] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[2],xc[3],gl[2],gl[3],a,b);
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a,this->ymax,this->zmin ));  
             }
             // Edge 4  0-3 - y
             if (( gl[0] > meshParam->thrsh && gl[3] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[3] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[0],yc[3],gl[0],gl[3],a,b);  
                 tesselationPoints.push_back( Point(this->xmin,(meshParam->thrsh-b)/a,this->zmin));  
             }
             // Edge 5  4-5 - x
             if (( gl[4] > meshParam->thrsh && gl[5] < meshParam->thrsh ) || ( gl[4] < meshParam->thrsh && gl[5] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[4],xc[5],gl[4],gl[5],a,b);  
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a,this->ymin,this->zmax)); 
             }         
             // Edge 6  5-6 - y
             if (( gl[5] > meshParam->thrsh && gl[6] < meshParam->thrsh ) || ( gl[5] < meshParam->thrsh && gl[6] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[5],yc[6],gl[5],gl[6],a,b); 
                 tesselationPoints.push_back( Point(this->xmax,(meshParam->thrsh-b)/a,this->zmax));
             }         
             // Edge 7  6-7 - x
             if (( gl[6] > meshParam->thrsh && gl[7] < meshParam->thrsh ) || ( gl[6] < meshParam->thrsh && gl[7] > meshParam->thrsh))
             {
                 interpolateLinearly(xc[6],xc[7],gl[6],gl[7],a,b); 
                 tesselationPoints.push_back( Point((meshParam->thrsh-b)/a, this->ymax, this->zmax )); 
             }
             // Edge 8  7-4 - y
             if (( gl[7] > meshParam->thrsh && gl[4] < meshParam->thrsh ) || ( gl[7] < meshParam->thrsh && gl[4] > meshParam->thrsh))
             {
                 interpolateLinearly(yc[7],yc[4],gl[7],gl[4],a,b);
                 tesselationPoints.push_back( Point(this->xmin,(meshParam->thrsh-b)/a,this->zmax));  
             }
             // Edge 9  0-4 - z
             if (( gl[0] > meshParam->thrsh && gl[4] < meshParam->thrsh ) || ( gl[0] < meshParam->thrsh && gl[4] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[0],zc[4],gl[0],gl[4],a,b); 
                 tesselationPoints.push_back( Point(this->xmin, this->ymin , (meshParam->thrsh-b)/a  ));  
             
             }
             // Edge 10 1-5 - z
             if (( gl[1] > meshParam->thrsh && gl[5] < meshParam->thrsh ) || ( gl[1] < meshParam->thrsh && gl[5] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[1],zc[5],gl[1],gl[5],a,b); 
                 tesselationPoints.push_back( Point(this->xmax, this->ymin, (meshParam->thrsh-b)/a ));  
             }
             // Edge 11  2-6 -z
             if (( gl[2] > meshParam->thrsh && gl[6] < meshParam->thrsh ) || ( gl[2] < meshParam->thrsh && gl[6] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[2],zc[6],gl[2],gl[6],a,b); 
                 tesselationPoints.push_back( Point( this->xmax, this->ymax, (meshParam->thrsh-b)/a)); 
             }
             // Edge 12  3-7 -z 
             if (( gl[3] > meshParam->thrsh && gl[7] < meshParam->thrsh ) || ( gl[3] < meshParam->thrsh && gl[7] > meshParam->thrsh))
             {
                 interpolateLinearly(zc[3],zc[7],gl[3],gl[7],a,b);
                 tesselationPoints.push_back( Point( this->xmin, this->ymax, (meshParam->thrsh-b)/a ));  
             }
             Delaunay dt;
             dt.insert( tesselationPoints.begin(), tesselationPoints.end() );
             // Integration over the tetrahedrons   
             for (Delaunay::Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it)
             {
                 x1 = dt.tetrahedron(it)[0].x(); y1= dt.tetrahedron(it)[0].y(); z1 = dt.tetrahedron(it)[0].z(); 
                 x2 = dt.tetrahedron(it)[1].x(); y2= dt.tetrahedron(it)[1].y(); z2 = dt.tetrahedron(it)[1].z(); 
                 x3 = dt.tetrahedron(it)[2].x(); y3= dt.tetrahedron(it)[2].y(); z3 = dt.tetrahedron(it)[2].z(); 
                 x4 = dt.tetrahedron(it)[3].x(); y4= dt.tetrahedron(it)[3].y(); z4 = dt.tetrahedron(it)[3].z(); 
 
                 //Getting the volume of the tetrahedron 
                 volT = (1./6.)*( std::abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                            (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                            (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) ) ;
                            
                 // Performing the integration over the tetrahedron   
                 // Multiply the barycentric coordinates matrix by the tetrahedron node coordinates    
                 Eigen::MatrixXd elemCoord(4,3); 
                 elemCoord(0,0)=x1; elemCoord(0,1)=y1; elemCoord(0,2)=z1;
                 elemCoord(1,0)=x2; elemCoord(1,1)=y2; elemCoord(1,2)=z2;
                 elemCoord(2,0)=x3; elemCoord(2,1)=y3; elemCoord(2,2)=z3;
                 elemCoord(3,0)=x4; elemCoord(3,1)=y4; elemCoord(3,2)=z4;
        
                 Eigen::MatrixXd pg = (meshParam->GaussTetraBc)*elemCoord ; 
                 for (int ig=0; ig<meshParam->nbgTetrahedron; ig++)
                 {
 
                     // Gauss weight 
                     wgMes = meshParam->GaussTetraW[ig]*volT ; 
                     
                     // Evaluating the univariate basis functions at the gauss point 
                     dersbasisfuns(meshParam->deg_xi, meshParam->knotXi, pg(ig,0), spanx, 1, meshParam->bfOutputXi ) ; 
                     dersbasisfuns(meshParam->deg_eta, meshParam->knotEta, pg(ig,1), spany, 1, meshParam->bfOutputEta ) ;
                     dersbasisfuns(meshParam->deg_zeta, meshParam->knotZeta, pg(ig,2), spanz, 1, meshParam->bfOutputZeta ) ;
                     // Getting the trivariate basis functions 
                     
                     bfc = 0 ; 
                     for (int kbf=0; kbf< meshParam->nbf_elem_zeta; kbf++)
                     {
                         for (int jbf=0; jbf < meshParam->nbf_elem_eta; jbf++)
                         {
                             for (int ibf=0; ibf < meshParam->nbf_elem_xi; ibf++)
                             {
                                 meshParam->N[bfc]       = meshParam->bfOutputXi[0][ibf]*meshParam->bfOutputEta[0][jbf]*meshParam->bfOutputZeta[0][kbf] ; 
                                 meshParam->dNdxi[bfc]   = meshParam->bfOutputXi[1][ibf]*meshParam->bfOutputEta[0][jbf]*meshParam->bfOutputZeta[0][kbf] ; 
                                 meshParam->dNdeta[bfc]  = meshParam->bfOutputXi[0][ibf]*meshParam->bfOutputEta[1][jbf]*meshParam->bfOutputZeta[0][kbf] ; 
                                 meshParam->dNdzeta[bfc] = meshParam->bfOutputXi[0][ibf]*meshParam->bfOutputEta[0][jbf]*meshParam->bfOutputZeta[1][kbf] ; 
                                 
                                 // Add the integral value of the basis function 
                                 intBf[meshParam->NOELEM(ie,bfc)] +=  wgMes * meshParam->N[bfc] ; 
                                 bfc++; 
                             }
                         }
                     }
                     // Set the cell differential matrix 
                     for (int ibf=0; ibf< meshParam->nbf_elem; ibf++)
                     {
        				(*meshParam->Bc)( 0, ibf )                      = meshParam->dNdxi[ibf];
        				(*meshParam->Bc)( 1, ibf+meshParam->nbf_elem)   = meshParam->dNdeta[ibf] ;
        				(*meshParam->Bc)( 2, ibf+2*meshParam->nbf_elem) = meshParam->dNdzeta[ibf] ;
        				
        				(*meshParam->Bc)( 3, ibf)                       = meshParam->dNdeta[ibf] ;
        				(*meshParam->Bc)( 3, ibf+meshParam->nbf_elem)   = meshParam->dNdxi[ibf];
        				
        				(*meshParam->Bc)( 4, ibf)                       = meshParam->dNdzeta[ibf] ;
        				(*meshParam->Bc)( 4, ibf+2*meshParam->nbf_elem) = meshParam->dNdxi[ibf];  
        				
        				(*meshParam->Bc)( 5, ibf+meshParam->nbf_elem)   = meshParam->dNdzeta[ibf] ;
        				(*meshParam->Bc)( 5, ibf+2*meshParam->nbf_elem) = meshParam->dNdeta[ibf] ;  
                      
                     }
                     // Add contribution to the stiffness 
                     *(meshParam->Ke) = *(meshParam->Ke) + meshParam->Bc->transpose() * meshParam->hooke * (*meshParam->Bc) * wgMes ;
                 }                     
             }
         }
     }    
     else 
     {
         int nbp = std::pow(2,lvlmax-this->lvl)+1; // For the subdivision of the the boundary of a subcell
         double eps = 1.e-8; 
 		 Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nbp, this->xmin + eps, this->xmax - eps);
 		 Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(nbp, this->ymin + eps, this->ymax - eps);     
 		 Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(nbp, this->zmin + eps, this->zmax - eps);    
         int countZero = CheckCutCuboidTrilinearInterpolation(meshParam->image,meshParam->si,meshParam->sj,meshParam->sk, 
            					        meshParam->knotXiImage, meshParam->sknotXiImage,   
            					        meshParam->knotEtaImage, meshParam->sknotEtaImage, 
            					        meshParam->knotZetaImage, meshParam->sknotZetaImage, 
            					        meshParam->thrsh, x, y, z); 
         int nborder = 6*nbp*nbp; // Total number of sampling points on the boundary of the cube 
         if (countZero==0)
         {
              // Integration the cell 
              
              // Setting Gauss integration coordinates
              // And evaluating the univariate basis functions  
              for (int ii =0; ii< meshParam->nbg_xi; ii++)
              {
                meshParam->Cxig[ii] = this->xmin + 0.5*(meshParam->xig[ii]+1)*this->mesx;
                dersbasisfuns(meshParam->deg_xi, meshParam->knotXi, meshParam->Cxig[ii], spanx, 1, meshParam->bfOutputXi ); 
                for (int jbf=0; jbf<meshParam->nbf_elem_xi; jbf++)
                {
                    meshParam->Nxi[ii][jbf] 	= meshParam->bfOutputXi[0][jbf];
         			meshParam->dNxidxi[ii][jbf] = meshParam->bfOutputXi[1][jbf];   
                }
              }    
              for (int ii =0; ii< meshParam->nbg_eta; ii++)
              {
                meshParam->Cetag[ii] = this->ymin + 0.5*(meshParam->etag[ii]+1)*this->mesy;
                dersbasisfuns(meshParam->deg_eta, meshParam->knotEta,meshParam->Cetag[ii], spany, 1, meshParam->bfOutputEta);
           	   for ( int jbf =0; jbf<meshParam->nbf_elem_eta; jbf++)
                {
               	   meshParam->Neta[ii][jbf] 	   = meshParam->bfOutputEta[0][jbf];
               	   meshParam->dNetadeta[ii][jbf]   = meshParam->bfOutputEta[1][jbf];
           	   }
              }  
              for (int ii =0; ii< meshParam->nbg_zeta; ii++)
              {
                meshParam->Czetag[ii] = this->zmin + 0.5*(meshParam->zetag[ii]+1)*this->mesz;
                dersbasisfuns(meshParam->deg_zeta, meshParam->knotZeta,meshParam->Czetag[ii], spanz, 1, meshParam->bfOutputZeta);
                for ( int jbf =0; jbf<meshParam->nbf_elem_zeta; jbf++)
                {
               	   meshParam->Nzeta[ii][jbf] 	   = meshParam->bfOutputZeta[0][jbf];
               	   meshParam->dNzetadzeta[ii][jbf]  = meshParam->bfOutputZeta[1][jbf];
                } 
              }        
              // Loop over the integration points of the cell
              // Adding the contribution of the cell       
              t = 0; 
              for (int kk=0; kk<meshParam->nbg_zeta; kk++)
              {
                  for (int jj=0; jj<meshParam->nbg_eta; jj++)
                  {
                      for (int ii=0; ii<meshParam->nbg_xi; ii++)
                      {  
                          // Getting the weighting of the integration point 
                          wgMes = meshParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ; 
                          
                          // Getting the trivariate basis functions 
                          bfc = 0 ; 
                          for (int kbf=0; kbf< meshParam->nbf_elem_zeta; kbf++)
                          {
                              for (int jbf=0; jbf < meshParam->nbf_elem_eta; jbf++)
                              {
                                  for (int ibf=0; ibf < meshParam->nbf_elem_xi; ibf++)
                                  {
                                      meshParam->N[bfc]       = meshParam->Nxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->Nzeta[kk][kbf]  ; 
                                      meshParam->dNdxi[bfc]   = meshParam->dNxidxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->Nzeta[kk][kbf]  ; 
                                      meshParam->dNdeta[bfc]  = meshParam->Nxi[ii][ibf]*meshParam->dNetadeta[jj][jbf]*meshParam->Nzeta[kk][kbf]  ; 
                                      meshParam->dNdzeta[bfc] = meshParam->Nxi[ii][ibf]*meshParam->Neta[jj][jbf]*meshParam->dNzetadzeta[kk][kbf]  ; 
                                      
                                      // Add the integral value of the basis function 
                                      intBf[meshParam->NOELEM(ie,bfc)] +=  wgMes * meshParam->N[bfc] ; 
                                      bfc++; 
                                  }
                              }
                          }
                          // Set the cell differential matrix 
                          for (int ibf=0; ibf< meshParam->nbf_elem; ibf++)
                          {
             				(*meshParam->Bc)( 0, ibf )                      = meshParam->dNdxi[ibf];
             				(*meshParam->Bc)( 1, ibf+meshParam->nbf_elem)   = meshParam->dNdeta[ibf] ;
             				(*meshParam->Bc)( 2, ibf+2*meshParam->nbf_elem) = meshParam->dNdzeta[ibf] ;
             				
             				(*meshParam->Bc)( 3, ibf)                       = meshParam->dNdeta[ibf] ;
             				(*meshParam->Bc)( 3, ibf+meshParam->nbf_elem)   = meshParam->dNdxi[ibf];
             				
             				(*meshParam->Bc)( 4, ibf)                       = meshParam->dNdzeta[ibf] ;
             				(*meshParam->Bc)( 4, ibf+2*meshParam->nbf_elem) = meshParam->dNdxi[ibf];  
             				
             				(*meshParam->Bc)( 5, ibf+meshParam->nbf_elem)   = meshParam->dNdzeta[ibf] ;
             				(*meshParam->Bc)( 5, ibf+2*meshParam->nbf_elem) = meshParam->dNdeta[ibf] ;  
                           
                          }
                          // Add contribution to the stiffness 
                          *(meshParam->Ke) = *(meshParam->Ke) + meshParam->Bc->transpose() * meshParam->hooke * (*meshParam->Bc) * wgMes ;
                          t++; 
                      }
                  }
              }
                          
         }
         else if (countZero != nborder)
         {
             // Decompose the cut cell 
            this->bottomBackLeft   = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomBackRight  = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontLeft  = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontRight = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->topBackLeft      = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topBackRight     = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontLeft     = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontRight    = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            
            this->bottomBackLeft->lvl = this->lvl + 1 ;
            this->bottomBackRight->lvl = this->lvl + 1 ;
            this->bottomFrontLeft->lvl = this->lvl + 1 ;
            this->bottomFrontRight->lvl = this->lvl + 1 ;
            this->topBackLeft->lvl = this->lvl + 1 ;
            this->topBackRight->lvl = this->lvl + 1 ;
            this->topFrontLeft->lvl = this->lvl + 1 ;
            this->topFrontRight->lvl = this->lvl + 1  ;
            
            

            this->bottomBackLeft->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->bottomBackRight->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->bottomFrontLeft->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->bottomFrontRight->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->topBackLeft->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->topBackRight->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->topFrontLeft->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
            this->topFrontRight->DecomposeTesselationTrilinearInterpAssembly(ie, spanx, spany, spanz, lvlmax, meshParam, intBf ); 
                                                                               
         }		        
     }
     
 
} 


void Cell::DecomposeTesselationTrilinearInterpAssemblyParallel(int ie, int spanx, int spany, int spanz, 
                                                 int lvlmax, 
                                                 SharedParameters*  globalParam, 
                                                 PrivateParameters* elemParam, 
                                                 double* intBf)
{
    int t   ;
    int bfc ; 
    double wgMes ; 
    if (this->lvl == lvlmax)
    {   
         // First checking if the cuboid corners are > than threshold value
         // in this case, we perform an integration on a the square cell
         double eps = 1.e-8;      
         std::array<double,8> xc =  { this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps, this->xmin+eps, this->xmax-eps, this->xmax-eps, this->xmin+eps }; 
         std::array<double,8> yc =  { this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps, this->ymin+eps, this->ymin+eps, this->ymax-eps, this->ymax-eps };
         std::array<double,8> zc =  { this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmin+eps, this->zmax-eps, this->zmax-eps, this->zmax-eps, this->zmax-eps };
         std::array<double,8> gl; 
         for (int i=0; i<8; i++)
         {
            gl[i] = EvaluateTrilinearInterpolationOnOnePoint(globalParam->image,globalParam->si,globalParam->sj,globalParam->sk, 
           					        globalParam->knotXiImage, globalParam->sknotXiImage,   
           					        globalParam->knotEtaImage, globalParam->sknotEtaImage, 
           					        globalParam->knotZetaImage, globalParam->sknotZetaImage,  
       					                xc[i], yc[i], zc[i]); 		                
         } 
         double thrsh = globalParam->thrsh ; 
         if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g>thrsh ;}))
         {
             // Integration the cell 
             
             // Setting Gauss integration coordinates
             // And evaluating the univariate basis functions  
             for (int ii =0; ii< globalParam->nbg_xi; ii++)
             {
               elemParam->Cxig[ii] = this->xmin + 0.5*(globalParam->xig[ii]+1)*this->mesx;
               dersbasisfuns(globalParam->deg_xi, globalParam->knotXi, elemParam->Cxig[ii], spanx, 1, elemParam->bfOutputXi ); 
               for (int jbf=0; jbf<globalParam->nbf_elem_xi; jbf++)
               {
                    elemParam->Nxi[ii][jbf] 	= elemParam->bfOutputXi[0][jbf];
        			elemParam->dNxidxi[ii][jbf] = elemParam->bfOutputXi[1][jbf];   
               }
             }    
             for (int ii =0; ii< globalParam->nbg_eta; ii++)
             {
               elemParam->Cetag[ii] = this->ymin + 0.5*(globalParam->etag[ii]+1)*this->mesy;
               dersbasisfuns(globalParam->deg_eta, globalParam->knotEta,elemParam->Cetag[ii], spany, 1, elemParam->bfOutputEta);
          	   for ( int jbf =0; jbf<globalParam->nbf_elem_eta; jbf++)
               {
              	   elemParam->Neta[ii][jbf] 	   = elemParam->bfOutputEta[0][jbf];
              	   elemParam->dNetadeta[ii][jbf]   = elemParam->bfOutputEta[1][jbf];
          	   }
             }  
             for (int ii =0; ii< globalParam->nbg_zeta; ii++)
             {
               elemParam->Czetag[ii] = this->zmin + 0.5*(globalParam->zetag[ii]+1)*this->mesz;
               dersbasisfuns(globalParam->deg_zeta, globalParam->knotZeta,elemParam->Czetag[ii], spanz, 1, elemParam->bfOutputZeta);
               for ( int jbf =0; jbf<globalParam->nbf_elem_zeta; jbf++)
               {
              	   elemParam->Nzeta[ii][jbf] 	    = elemParam->bfOutputZeta[0][jbf];
              	   elemParam->dNzetadzeta[ii][jbf]  = elemParam->bfOutputZeta[1][jbf];
               } 
             }        
             // Loop over the integration points of the cell
             // Adding the contribution of the cell       
             t = 0; 
             for (int kk=0; kk<globalParam->nbg_zeta; kk++)
             {
                 for (int jj=0; jj<globalParam->nbg_eta; jj++)
                 {
                     for (int ii=0; ii<globalParam->nbg_xi; ii++)
                     {  
                         // Getting the weighting of the integration point 
                         wgMes = globalParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ; 
                         
                         // Getting the trivariate basis functions 
                         bfc = 0 ; 
                         for (int kbf=0; kbf< globalParam->nbf_elem_zeta; kbf++)
                         {
                             for (int jbf=0; jbf < globalParam->nbf_elem_eta; jbf++)
                             {
                                 for (int ibf=0; ibf < globalParam->nbf_elem_xi; ibf++)
                                 {
                                     elemParam->N[bfc]       = elemParam->Nxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->Nzeta[kk][kbf]       ; 
                                     elemParam->dNdxi[bfc]   = elemParam->dNxidxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->Nzeta[kk][kbf]   ; 
                                     elemParam->dNdeta[bfc]  = elemParam->Nxi[ii][ibf]*elemParam->dNetadeta[jj][jbf]*elemParam->Nzeta[kk][kbf]  ; 
                                     elemParam->dNdzeta[bfc] = elemParam->Nxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->dNzetadzeta[kk][kbf] ; 
                                     
                                     // Add the integral value of the basis function 
                                     intBf[globalParam->NOELEM(ie,bfc)] +=  wgMes * elemParam->N[bfc] ; 
                                     bfc++; 
                                 }
                             }
                         }
                         
                         // Set the cell differential matrix 
                         for (int ibf=0; ibf< globalParam->nbf_elem; ibf++)
                         {
            				(*elemParam->Bc)( 0, ibf )                        = elemParam->dNdxi[ibf];
            				(*elemParam->Bc)( 1, ibf+globalParam->nbf_elem)   = elemParam->dNdeta[ibf] ;
            				(*elemParam->Bc)( 2, ibf+2*globalParam->nbf_elem) = elemParam->dNdzeta[ibf] ;
            				
            				(*elemParam->Bc)( 3, ibf)                         = elemParam->dNdeta[ibf] ;
            				(*elemParam->Bc)( 3, ibf+globalParam->nbf_elem)   = elemParam->dNdxi[ibf];
            				
            				(*elemParam->Bc)( 4, ibf)                         = elemParam->dNdzeta[ibf] ;
            				(*elemParam->Bc)( 4, ibf+2*globalParam->nbf_elem) = elemParam->dNdxi[ibf];  
            				
            				(*elemParam->Bc)( 5, ibf+globalParam->nbf_elem)   = elemParam->dNdzeta[ibf] ;
            				(*elemParam->Bc)( 5, ibf+2*globalParam->nbf_elem) = elemParam->dNdeta[ibf] ;  
                          
                         }
                         // Add contribution to the stiffness 
                         *(elemParam->Ke) = *(elemParam->Ke) + elemParam->Bc->transpose() * globalParam->hooke * (*elemParam->Bc) * wgMes ;
                         t++; 
                     }
                 }
             }
         }
         else if (std::all_of(gl.begin(), gl.end(), [thrsh](double g){return g< thrsh;})==false)
         {
             // Perform Delaunay triangulation on the last cut-cell 
             double a,b ; // Level set linear interpolation coefficients for the tesselation approximation
             double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4; // tetrahedron nodes 
             double volT;  // Volume of the terahedron   
             std::vector<Point> tesselationPoints;
             for (int i=0; i<8; i++)
             {
                 if (gl[i] > globalParam->thrsh)
                 {
                     tesselationPoints.push_back(Point(xc[i],yc[i],zc[i]));
                 }
             }
             // Distance linearization of the interface 
             // Edge 1  0-1 - x    
             if (( gl[0] > globalParam->thrsh && gl[1] < globalParam->thrsh ) || ( gl[0] < globalParam->thrsh && gl[1] > globalParam->thrsh))
             {
                 interpolateLinearly(xc[0],xc[1],gl[0],gl[1], a,b); 
                 tesselationPoints.push_back( Point( (globalParam->thrsh-b)/a, this->ymin, this->zmin  )); 
             }
             // Edge 2  1-2 - y
             if (( gl[1] > globalParam->thrsh && gl[2] < globalParam->thrsh ) || ( gl[1] < globalParam->thrsh && gl[2] > globalParam->thrsh))
             {
                 interpolateLinearly(yc[1],yc[2],gl[1],gl[2],a,b);
                 tesselationPoints.push_back( Point(this->xmax,(globalParam->thrsh-b)/a,this->zmin ));  
             }         
             // Edge 3  2-3 - x
             if (( gl[2] > globalParam->thrsh && gl[3] < globalParam->thrsh ) || ( gl[2] < globalParam->thrsh && gl[3] > globalParam->thrsh))
             {
                 interpolateLinearly(xc[2],xc[3],gl[2],gl[3],a,b);
                 tesselationPoints.push_back( Point((globalParam->thrsh-b)/a,this->ymax,this->zmin ));  
             }
             // Edge 4  0-3 - y
             if (( gl[0] > globalParam->thrsh && gl[3] < globalParam->thrsh ) || ( gl[0] < globalParam->thrsh && gl[3] > globalParam->thrsh))
             {
                 interpolateLinearly(yc[0],yc[3],gl[0],gl[3],a,b);  
                 tesselationPoints.push_back( Point(this->xmin,(globalParam->thrsh-b)/a,this->zmin));  
             }
             // Edge 5  4-5 - x
             if (( gl[4] > globalParam->thrsh && gl[5] < globalParam->thrsh ) || ( gl[4] < globalParam->thrsh && gl[5] > globalParam->thrsh))
             {
                 interpolateLinearly(xc[4],xc[5],gl[4],gl[5],a,b);  
                 tesselationPoints.push_back( Point((globalParam->thrsh-b)/a,this->ymin,this->zmax)); 
             }         
             // Edge 6  5-6 - y
             if (( gl[5] > globalParam->thrsh && gl[6] < globalParam->thrsh ) || ( gl[5] < globalParam->thrsh && gl[6] > globalParam->thrsh))
             {
                 interpolateLinearly(yc[5],yc[6],gl[5],gl[6],a,b); 
                 tesselationPoints.push_back( Point(this->xmax,(globalParam->thrsh-b)/a,this->zmax));
             }         
             // Edge 7  6-7 - x
             if (( gl[6] > globalParam->thrsh && gl[7] < globalParam->thrsh ) || ( gl[6] < globalParam->thrsh && gl[7] > globalParam->thrsh))
             {
                 interpolateLinearly(xc[6],xc[7],gl[6],gl[7],a,b); 
                 tesselationPoints.push_back( Point((globalParam->thrsh-b)/a, this->ymax, this->zmax )); 
             }
             // Edge 8  7-4 - y
             if (( gl[7] > globalParam->thrsh && gl[4] < globalParam->thrsh ) || ( gl[7] < globalParam->thrsh && gl[4] > globalParam->thrsh))
             {
                 interpolateLinearly(yc[7],yc[4],gl[7],gl[4],a,b);
                 tesselationPoints.push_back( Point(this->xmin,(globalParam->thrsh-b)/a,this->zmax));  
             }
             // Edge 9  0-4 - z
             if (( gl[0] > globalParam->thrsh && gl[4] < globalParam->thrsh ) || ( gl[0] < globalParam->thrsh && gl[4] > globalParam->thrsh))
             {
                 interpolateLinearly(zc[0],zc[4],gl[0],gl[4],a,b); 
                 tesselationPoints.push_back( Point(this->xmin, this->ymin , (globalParam->thrsh-b)/a  ));  
             
             }
             // Edge 10 1-5 - z
             if (( gl[1] > globalParam->thrsh && gl[5] < globalParam->thrsh ) || ( gl[1] < globalParam->thrsh && gl[5] > globalParam->thrsh))
             {
                 interpolateLinearly(zc[1],zc[5],gl[1],gl[5],a,b); 
                 tesselationPoints.push_back( Point(this->xmax, this->ymin, (globalParam->thrsh-b)/a ));  
             }
             // Edge 11  2-6 -z
             if (( gl[2] > globalParam->thrsh && gl[6] < globalParam->thrsh ) || ( gl[2] < globalParam->thrsh && gl[6] > globalParam->thrsh))
             {
                 interpolateLinearly(zc[2],zc[6],gl[2],gl[6],a,b); 
                 tesselationPoints.push_back( Point( this->xmax, this->ymax, (globalParam->thrsh-b)/a)); 
             }
             // Edge 12  3-7 -z 
             if (( gl[3] > globalParam->thrsh && gl[7] < globalParam->thrsh ) || ( gl[3] < globalParam->thrsh && gl[7] > globalParam->thrsh))
             {
                 interpolateLinearly(zc[3],zc[7],gl[3],gl[7],a,b);
                 tesselationPoints.push_back( Point( this->xmin, this->ymax, (globalParam->thrsh-b)/a ));  
             }
             Delaunay dt;
             dt.insert( tesselationPoints.begin(), tesselationPoints.end() );
             // Integration over the tetrahedrons   
             for (Delaunay::Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it)
             {
                 x1 = dt.tetrahedron(it)[0].x(); y1= dt.tetrahedron(it)[0].y(); z1 = dt.tetrahedron(it)[0].z(); 
                 x2 = dt.tetrahedron(it)[1].x(); y2= dt.tetrahedron(it)[1].y(); z2 = dt.tetrahedron(it)[1].z(); 
                 x3 = dt.tetrahedron(it)[2].x(); y3= dt.tetrahedron(it)[2].y(); z3 = dt.tetrahedron(it)[2].z(); 
                 x4 = dt.tetrahedron(it)[3].x(); y4= dt.tetrahedron(it)[3].y(); z4 = dt.tetrahedron(it)[3].z(); 
 
                 //Getting the volume of the tetrahedron 
                 volT = (1./6.)*( std::abs( (x2-x1)*( (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1) )-
                            (y2-y1)*( (x3-x1)*(z4-z1) - (x4-x1)*(z3-z1) )+
                            (z2-z1)*( (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)) ) ) ;
                            
                 // Performing the integration over the tetrahedron   
                 // Multiply the barycentric coordinates matrix by the tetrahedron node coordinates    
                 Eigen::MatrixXd elemCoord(4,3); 
                 elemCoord(0,0)=x1; elemCoord(0,1)=y1; elemCoord(0,2)=z1;
                 elemCoord(1,0)=x2; elemCoord(1,1)=y2; elemCoord(1,2)=z2;
                 elemCoord(2,0)=x3; elemCoord(2,1)=y3; elemCoord(2,2)=z3;
                 elemCoord(3,0)=x4; elemCoord(3,1)=y4; elemCoord(3,2)=z4;
        
                 Eigen::MatrixXd pg = (globalParam->GaussTetraBc)*elemCoord ; 
                 for (int ig=0; ig<globalParam->nbgTetrahedron; ig++)
                 {
 
                     // Gauss weight 
                     wgMes = globalParam->GaussTetraW[ig]*volT ; 
                     
                     // Evaluating the univariate basis functions at the gauss point 
                     dersbasisfuns(globalParam->deg_xi, globalParam->knotXi, pg(ig,0), spanx, 1, elemParam->bfOutputXi ) ; 
                     dersbasisfuns(globalParam->deg_eta, globalParam->knotEta, pg(ig,1), spany, 1, elemParam->bfOutputEta ) ;
                     dersbasisfuns(globalParam->deg_zeta, globalParam->knotZeta, pg(ig,2), spanz, 1, elemParam->bfOutputZeta ) ;
                     // Getting the trivariate basis functions 
                     
                     bfc = 0 ; 
                     for (int kbf=0; kbf< globalParam->nbf_elem_zeta; kbf++)
                     {
                         for (int jbf=0; jbf < globalParam->nbf_elem_eta; jbf++)
                         {
                             for (int ibf=0; ibf < globalParam->nbf_elem_xi; ibf++)
                             {
                                 elemParam->N[bfc]       = elemParam->bfOutputXi[0][ibf]*elemParam->bfOutputEta[0][jbf]*elemParam->bfOutputZeta[0][kbf] ; 
                                 elemParam->dNdxi[bfc]   = elemParam->bfOutputXi[1][ibf]*elemParam->bfOutputEta[0][jbf]*elemParam->bfOutputZeta[0][kbf] ; 
                                 elemParam->dNdeta[bfc]  = elemParam->bfOutputXi[0][ibf]*elemParam->bfOutputEta[1][jbf]*elemParam->bfOutputZeta[0][kbf] ; 
                                 elemParam->dNdzeta[bfc] = elemParam->bfOutputXi[0][ibf]*elemParam->bfOutputEta[0][jbf]*elemParam->bfOutputZeta[1][kbf] ; 
                                 
                                 // Add the integral value of the basis function 
                                 intBf[globalParam->NOELEM(ie,bfc)] +=  wgMes * elemParam->N[bfc] ; 
                                 bfc++; 
                             }
                         }
                     }
                     // Set the cell differential matrix 
                     for (int ibf=0; ibf< globalParam->nbf_elem; ibf++)
                     {
        				(*elemParam->Bc)( 0, ibf )                      = elemParam->dNdxi[ibf];
        				(*elemParam->Bc)( 1, ibf+globalParam->nbf_elem)   = elemParam->dNdeta[ibf] ;
        				(*elemParam->Bc)( 2, ibf+2*globalParam->nbf_elem) = elemParam->dNdzeta[ibf] ;
        				
        				(*elemParam->Bc)( 3, ibf)                       = elemParam->dNdeta[ibf] ;
        				(*elemParam->Bc)( 3, ibf+globalParam->nbf_elem)   = elemParam->dNdxi[ibf];
        				
        				(*elemParam->Bc)( 4, ibf)                       = elemParam->dNdzeta[ibf] ;
        				(*elemParam->Bc)( 4, ibf+2*globalParam->nbf_elem) = elemParam->dNdxi[ibf];  
        				
        				(*elemParam->Bc)( 5, ibf+globalParam->nbf_elem)   = elemParam->dNdzeta[ibf] ;
        				(*elemParam->Bc)( 5, ibf+2*globalParam->nbf_elem) = elemParam->dNdeta[ibf] ;  
                      
                     }
                     // Add contribution to the stiffness 
                     *(elemParam->Ke) = *(elemParam->Ke) + elemParam->Bc->transpose() * globalParam->hooke * (*elemParam->Bc) * wgMes ;
                 }                     
             }
         }
     }    
     else 
     {
         int nbp = std::pow(2,lvlmax-this->lvl)+1; // For the subdivision of the the boundary of a subcell
         double eps = 1.e-8; 
 		 Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(nbp, this->xmin + eps, this->xmax - eps);
 		 Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(nbp, this->ymin + eps, this->ymax - eps);     
 		 Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(nbp, this->zmin + eps, this->zmax - eps);    
         int countZero = CheckCutCuboidTrilinearInterpolation(globalParam->image,globalParam->si,globalParam->sj,globalParam->sk, 
            					        globalParam->knotXiImage, globalParam->sknotXiImage,   
            					        globalParam->knotEtaImage, globalParam->sknotEtaImage, 
            					        globalParam->knotZetaImage, globalParam->sknotZetaImage, 
            					        globalParam->thrsh, x, y, z); 
         int nborder = 6*nbp*nbp; // Total number of sampling points on the boundary of the cube 
         if (countZero==0)
         {
              // Integration the cell 
              
              // Setting Gauss integration coordinates
              // And evaluating the univariate basis functions  
              for (int ii =0; ii< globalParam->nbg_xi; ii++)
              {
                elemParam->Cxig[ii] = this->xmin + 0.5*(globalParam->xig[ii]+1)*this->mesx;
                dersbasisfuns(globalParam->deg_xi, globalParam->knotXi, elemParam->Cxig[ii], spanx, 1, elemParam->bfOutputXi ); 
                for (int jbf=0; jbf<globalParam->nbf_elem_xi; jbf++)
                {
                    elemParam->Nxi[ii][jbf] 	= elemParam->bfOutputXi[0][jbf];
         			elemParam->dNxidxi[ii][jbf] = elemParam->bfOutputXi[1][jbf];   
                }
              }    
              for (int ii =0; ii< globalParam->nbg_eta; ii++)
              {
                elemParam->Cetag[ii] = this->ymin + 0.5*(globalParam->etag[ii]+1)*this->mesy;
                dersbasisfuns(globalParam->deg_eta, globalParam->knotEta,elemParam->Cetag[ii], spany, 1, elemParam->bfOutputEta);
           	   for ( int jbf =0; jbf<globalParam->nbf_elem_eta; jbf++)
                {
               	   elemParam->Neta[ii][jbf] 	   = elemParam->bfOutputEta[0][jbf];
               	   elemParam->dNetadeta[ii][jbf]   = elemParam->bfOutputEta[1][jbf];
           	   }
              }  
              for (int ii =0; ii< globalParam->nbg_zeta; ii++)
              {
                elemParam->Czetag[ii] = this->zmin + 0.5*(globalParam->zetag[ii]+1)*this->mesz;
                dersbasisfuns(globalParam->deg_zeta, globalParam->knotZeta,elemParam->Czetag[ii], spanz, 1, elemParam->bfOutputZeta);
                for ( int jbf =0; jbf<globalParam->nbf_elem_zeta; jbf++)
                {
               	   elemParam->Nzeta[ii][jbf] 	   = elemParam->bfOutputZeta[0][jbf];
               	   elemParam->dNzetadzeta[ii][jbf]  = elemParam->bfOutputZeta[1][jbf];
                } 
              }        
              // Loop over the integration points of the cell
              // Adding the contribution of the cell       
              t = 0; 
              for (int kk=0; kk<globalParam->nbg_zeta; kk++)
              {
                  for (int jj=0; jj<globalParam->nbg_eta; jj++)
                  {
                      for (int ii=0; ii<globalParam->nbg_xi; ii++)
                      {  
                          // Getting the weighting of the integration point 
                          wgMes = globalParam->wg[t]*(this->mesx)*(this->mesy)*(this->mesz)/8 ; 
                          
                          // Getting the trivariate basis functions 
                          bfc = 0 ; 
                          for (int kbf=0; kbf< globalParam->nbf_elem_zeta; kbf++)
                          {
                              for (int jbf=0; jbf < globalParam->nbf_elem_eta; jbf++)
                              {
                                  for (int ibf=0; ibf < globalParam->nbf_elem_xi; ibf++)
                                  {
                                      elemParam->N[bfc]       = elemParam->Nxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->Nzeta[kk][kbf]  ; 
                                      elemParam->dNdxi[bfc]   = elemParam->dNxidxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->Nzeta[kk][kbf]  ; 
                                      elemParam->dNdeta[bfc]  = elemParam->Nxi[ii][ibf]*elemParam->dNetadeta[jj][jbf]*elemParam->Nzeta[kk][kbf]  ; 
                                      elemParam->dNdzeta[bfc] = elemParam->Nxi[ii][ibf]*elemParam->Neta[jj][jbf]*elemParam->dNzetadzeta[kk][kbf]  ; 
                                      
                                      // Add the integral value of the basis function 
                                      intBf[globalParam->NOELEM(ie,bfc)] +=  wgMes * elemParam->N[bfc] ; 
                                      bfc++; 
                                  }
                              }
                          }
                          // Set the cell differential matrix 
                          for (int ibf=0; ibf< globalParam->nbf_elem; ibf++)
                          {
             				(*elemParam->Bc)( 0, ibf )                        = elemParam->dNdxi[ibf];
             				(*elemParam->Bc)( 1, ibf+globalParam->nbf_elem)   = elemParam->dNdeta[ibf] ;
             				(*elemParam->Bc)( 2, ibf+2*globalParam->nbf_elem) = elemParam->dNdzeta[ibf] ;
             				
             				(*elemParam->Bc)( 3, ibf)                       = elemParam->dNdeta[ibf] ;
             				(*elemParam->Bc)( 3, ibf+globalParam->nbf_elem)   = elemParam->dNdxi[ibf];
             				
             				(*elemParam->Bc)( 4, ibf)                       = elemParam->dNdzeta[ibf] ;
             				(*elemParam->Bc)( 4, ibf+2*globalParam->nbf_elem) = elemParam->dNdxi[ibf];  
             				
             				(*elemParam->Bc)( 5, ibf+globalParam->nbf_elem)   = elemParam->dNdzeta[ibf] ;
             				(*elemParam->Bc)( 5, ibf+2*globalParam->nbf_elem) = elemParam->dNdeta[ibf] ;  
                           
                          }
                          // Add contribution to the stiffness 
                          *(elemParam->Ke) = *(elemParam->Ke) + elemParam->Bc->transpose() * globalParam->hooke * (*elemParam->Bc) * wgMes ;
                          t++; 
                      }
                  }
              }
                          
         }
         else if (countZero != nborder)
         {
             // Decompose the cut cell 
            this->bottomBackLeft   = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomBackRight  = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontLeft  = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->bottomFrontRight = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, this->zmin, (this->zmin+this->zmax)/2);
            this->topBackLeft      = new Cell(this->xmin, (this->xmin+this->xmax)/2, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topBackRight     = new Cell((this->xmin+this->xmax)/2, this->xmax, this->ymin, (this->ymin+this->ymax)/2, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontLeft     = new Cell(this->xmin, (this->xmin+this->xmax)/2, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            this->topFrontRight    = new Cell((this->xmin+this->xmax)/2, this->xmax, (this->ymin+this->ymax)/2, this->ymax, (this->zmin+this->zmax)/2, this->zmax);
            
            this->bottomBackLeft->lvl = this->lvl + 1 ;
            this->bottomBackRight->lvl = this->lvl + 1 ;
            this->bottomFrontLeft->lvl = this->lvl + 1 ;
            this->bottomFrontRight->lvl = this->lvl + 1 ;
            this->topBackLeft->lvl = this->lvl + 1 ;
            this->topBackRight->lvl = this->lvl + 1 ;
            this->topFrontLeft->lvl = this->lvl + 1 ;
            this->topFrontRight->lvl = this->lvl + 1  ;
            
 
 
            this->bottomBackLeft->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->bottomBackRight->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->bottomFrontLeft->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->bottomFrontRight->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->topBackLeft->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->topBackRight->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->topFrontLeft->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
            this->topFrontRight->DecomposeTesselationTrilinearInterpAssemblyParallel(ie, spanx, spany, spanz, lvlmax, globalParam, elemParam, intBf ); 
                                                                               
         }		        
     }
     
 
}

                  
 
                   
 
                   
 
                   
                   
