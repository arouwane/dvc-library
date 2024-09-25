#include "memory.h"


double* array_1d(int n)
{
    double *a = new double[n]; 
    return a; 
    // delete [] a; 
}

double** array_2d(int nrows, int ncols)
{
  int row;
  double **mat;
  mat    = new double*[nrows];
  mat[0] = new double[nrows*ncols];
  for (row = 1; row < nrows; row++)
    mat[row] = mat[row-1] + ncols;  
  return mat;
}
void delete_array_2d(double** arr)
{
    delete[] arr[0];
    delete[] arr;
}

//double *** array_3d( int ind1, int ind2, int ind3)
//{
//    double *** array = (double***)malloc((ind1*sizeof(double*))+(ind1*ind2*sizeof(double**))+(ind1*ind2*ind3*sizeof(double)));
//    for (int i=0; i< ind1; ++i){
//        array[i] = (double**)(array + ind1) + i*ind2;
//        for (int j=0; j< ind2; ++j){
//            array[i][j] = (double*)(array+ind1+ind1*ind2)+i*ind2*ind3+j*ind3; 
//            
//        }
//    }
//    return array; 
//}


double *** array_3d(int x, int y, int z, double* allElements  ) 
{
  int i, j;
  double ***array3D = new double**[x]; 
  for(i = 0; i < x; i++)
  {
    array3D[i] = new double*[y]; 
    for(j = 0; j < y; j++)
    {
        array3D[i][j] = allElements + (i * y * z) + (j * z);
    }
  }
  return array3D; 
}

void delete_array_3d(double*** array3D, int x, double* allElements)
{
    int i;  
    delete[] allElements;
    for(i = 0; i < x; i++)
    {
        delete[] array3D[i];
    }
    delete[] array3D;
}

void copy_array_1d(double *arr, int s1, double *copy)
{
    for(int i=0; i<s1; i++){
        copy[i] = arr[i]; 
    }
}

void copy_array_2d(double **arr, int s1, int s2, double *copy)
{
    int c=0; 
    for (int i=0;i<s1;i++){
        for (int j=0;j<s2;j++){
            copy[c] = arr[i][j]; 
            c++; 
            }
        }
}

void copy_array_3d(double ***arr, int s1, int s2, int s3, double *copy)
{
    int c=0;
    for (int i=0;i<s1;i++){
        for (int j=0;j<s2;j++){
            for(int k=0; k<s3; k++){
                copy[c] = arr[i][j][k]; 
                c++; 
            }
        }
    } 
}

void zero_1d_array(double* arr, int size)
{
    int i; 
    for (i=0; i< size; i++){
        arr[i] = 0 ;
    }
}

void identity(double** array, int size)
{
     int i,j;
     for (i=0; i< size; i++)
     {
         for (j=0; j<size; j++)
         {
             if (i==j)
                 array[i][j] = 1; 
             else 
                 array[i][j] = 0; 
         }
     }
}

void print_1d_array(double* arr, int n)
{
	int i;
	std::cout << "\n[ \n" ; 
	for (i =0; i <n; i++){
        std::cout << arr[i]<<"\t" ; 
		}
    std::cout << "\n ] \n" ; 
}

void print_2d_array(double **arr, int n, int m)
{
	std::cout << "\n[ \n" ; 
    for (int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            std::cout << arr[i][j] << "\t" ; 
        }
        std::cout << "\n" ; 
    }
    std::cout << "]\n" ; 
}
void print_3d_array(double ***arr, int n_arrays, int nrows, int ncols)
{
    for (int i =0; i< n_arrays; i++)
    {
        std::cout << "Array "<< i+1 ; 
        std::cout << "\n"; 
        for (int j=0; j<nrows; j++)
        {
            for(int k=0; k<ncols; k++)
            {
                std::cout << arr[i][j][k] <<"\t"; 
            }
        std::cout<<"\n";
        }
    }
}

//void extractSubMatrix(double* arr,int si, int imin, int imax, int jmin, int jmax, double* sub)
//{
//    int k =0; 	
//	for (int i=imin; i<=imax; i++){
//		for (int j=jmin; j<=jmax; j++){
//			sub[k] = arr[j+i*si]; 
//	        k++; 		
//		}
//	}
//}
