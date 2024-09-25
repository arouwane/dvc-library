#ifndef QUADTREEDECOMP_H
#define QUADTREEDECOMP_H

#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
#include <cmath>
#include "tools.h"
#include "meshUtils.h"
#include "Eigen/Core"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel1;
typedef CGAL::Delaunay_triangulation_3<Kernel1> Delaunay;
typedef Delaunay::Point Point;
//typedef CGAL::Triangulation_3<Kernel1> Triangulation;
//typedef Triangulation::Point Point;


class Cell
{
    public: 
    double xmin,xmax,ymin,ymax,zmin,zmax ;
    double mesx,mesy,mesz ; 
    int lvl; 
    Cell* bottomBackLeft   ;     
    Cell* bottomBackRight  ;  
    Cell* bottomFrontLeft  ; 
    Cell* bottomFrontRight ;
    Cell* topBackLeft   ;
    Cell* topBackRight  ;
    Cell* topFrontLeft  ;
    Cell* topFrontRight ;
    Cell(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);     
    void DecomposeTesselationTrilinearInterp(int lvlmax, MeshParameters* meshParam, std::vector<double>* integrationPoints ); 
    // Assembles stiffness based on the FCM method and computes the basis function indicator 
    void DecomposeTesselationTrilinearInterpAssembly(int ie, int spanx, int spany, int spanz, 
                                                     int lvlmax, MeshParameters* meshParam, 
                                                     double* intBf); 
    void DecomposeTesselationTrilinearInterpAssemblyParallel(int ie, int spanx, int spany, int spanz, 
                                                     int lvlmax, 
                                                     SharedParameters*  globalParam, 
                                                     PrivateParameters* elemParam, 
                                                     double* intBf);                                                      
}; 

 
 
 #endif