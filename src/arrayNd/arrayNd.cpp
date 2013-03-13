/*
Copyright (C) 2013 Bruno Golosio

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include "xrmc_arrayNd.h"
// Definition of functions of the namespace arrayNd for 2d and 3d array
// allocation and deallocation
//

// allocate a double 2d array with dimensions nx, ny
double **arrayNd::double_array2d(long nx, long ny)
{
  long i;
  double **arr;
  
  arr = (double**)calloc(nx, sizeof(double*));
  arr[0] = (double*)calloc(nx*ny, sizeof(double)); 
  for (i=1; i<nx; i++) {    
    arr[i] = arr[i-1] + ny; 
  }

  // return pointer
  return arr;
}


// free a double array 2d allocated by double_array2d()
int arrayNd::free_double_array2d(double **arr)
{
  free(arr[0]);
  free(arr);

  return 0;
}

// allocate a float 3d array with dimensions nx, ny, nz  
float ***arrayNd::float_array3d(long nx, long ny, long nz)
{
  long i,j;
  float ***arr;
  
  // allocate pointers to slices
  arr=(float ***) malloc(sizeof(float**)*nx);
  if (!arr) {
    printf("allocation error in float_array3d()");
    exit(EXIT_FAILURE);
  }
  // allocate rows and set pointers to them
  arr[0]=(float **) malloc(sizeof(float*)*nx*ny);
  if (!arr[0]) {
    printf("allocation error in float_array3d()");
    exit(EXIT_FAILURE);
  }
  arr[0][0]=(float *) calloc(nx*ny*nz,sizeof(float));
  if (!arr[0][0]) {
    printf("allocation error in float_array3d()");
    exit(EXIT_FAILURE);
  }
  for(j=1;j<ny;j++) arr[0][j] = arr[0][j-1] + nz;

  for(i=1; i<nx; i++) {
    arr[i] = arr[i-1] + ny;
    arr[i][0] = arr[i-1][0] + ny*nz;
    for(j=0;j<ny;j++) arr[i][j] = arr[i][j-1] + nz;
  }

  // return pointer
  return arr;
}

// allocate a double 3d array with dimensions nx, ny, nz  
double ***arrayNd::double_array3d(long nx, long ny, long nz)
{
  long i,j;
  double ***arr;
  
  // allocate pointers to slices
  arr=(double ***) malloc(sizeof(double**)*nx);
  if (!arr) {
    printf("allocation error in double_array3d()");
    exit(EXIT_FAILURE);
  }
  // allocate rows and set pointers to them
  arr[0]=(double **) malloc(sizeof(double*)*nx*ny);
  if (!arr[0]) {
    printf("allocation error in double_array3d()");
    exit(EXIT_FAILURE);
  }
  arr[0][0]=(double *) calloc(nx*ny*nz, sizeof(double));
  if (!arr[0][0]) {
    printf("allocation error in double_array3d()");
    exit(EXIT_FAILURE);
  }
  for(j=1;j<ny;j++) arr[0][j] = arr[0][j-1] + nz;

  for(i=1; i<nx; i++) {
    arr[i] = arr[i-1] + ny;
    arr[i][0] = arr[i-1][0] + ny*nz;
    for(j=0;j<ny;j++) arr[i][j] = arr[i][j-1] + nz;
  }

  // return pointer
  return arr;
}

// allocate an int 3d array with dimensions nx, ny, nz  
int ***arrayNd::int_array3d(long nx, long ny, long nz)
{
  long i,j;
  int ***arr;
  
  // allocate pointers to slices
  arr=(int ***) malloc(sizeof(int**)*nx);
  if (!arr) {
    printf("allocation error in int_array3d()");
    exit(EXIT_FAILURE);
  }
  // allocate rows and set pointers to them
  arr[0]=(int **) malloc(sizeof(int*)*nx*ny);
  if (!arr[0]) {
    printf("allocation error in int_array3d()");
    exit(EXIT_FAILURE);
  }
  arr[0][0]=(int *) calloc(nx*ny*nz, sizeof(int));
  if (!arr[0][0]) {
    printf("allocation error in int_array3d()");
    exit(EXIT_FAILURE);
  }
  for(j=1;j<ny;j++) arr[0][j] = arr[0][j-1] + nz;

  for(i=1; i<nx; i++) {
    arr[i] = arr[i-1] + ny;
    arr[i][0] = arr[i-1][0] + ny*nz;
    for(j=0;j<ny;j++) arr[i][j] = arr[i][j-1] + nz;
  }

  // return pointer
  return arr;
}

// free a float array 3d allocated by float_array3d()
int arrayNd::free_float_array3d(float ***arr)
{
        free(arr[0][0]);
        free(arr[0]);
        free(arr);

	return 0;
}

// free a double array 3d allocated by double_array3d()
int arrayNd::free_double_array3d(double ***arr)
{
        free(arr[0][0]);
        free(arr[0]);
        free(arr);

	return 0;
}

// free a int array 3d allocated by int_array3d()
int arrayNd::free_int_array3d(int ***arr)
{
        free(arr[0][0]);
        free(arr[0]);
        free(arr);

	return 0;
}
