#ifndef MEMORY_H
#define MEMORY_H

#include <iostream>


// ============================================================================
// MEMORY ALLOCATION, COPY AND DELETE FUNCTIONS 
// ============================================================================


double* array_1d(int n); 
double** array_2d(int nrows, int ncols);
double*** array_3d(int x, int y, int z, double* allElements  );

void delete_array_2d(double** arr); 
void delete_array_3d(double*** array3D, int x, double* allElements); 

void copy_array_1d(double *arr, int s1, double *copy); 
void copy_array_2d(double **arr, int s1, int s2, double *copy); 
void copy_array_3d(double ***arr, int s1, int s2, int s3, double *copy); 

void zero_1d_array(double* arr, int size); 
void identity(double** array, int size);  

void print_1d_array(double* arr, int n);  
void print_2d_array(double **arr, int n, int m);  
void print_3d_array(double ***arr, int n_arrays, int nrows, int ncols); 

//void extractSubMatrix(double* arr,int si, int imin, int imax, int jmin, int jmax, double* sub);  

//#include "memory.cpp"
#endif





 


