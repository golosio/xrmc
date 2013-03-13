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
/////////////////////////////////////////
//          xrmc_arrayNd.h             //
//        31/01/2013                   //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Namespace for 2d and 3d array allocation and deallocation
//

#ifndef ARRAYNDH
#define ARRAYNDH

namespace arrayNd
{
// allocate a 2d double array with dimensions nx, ny
double **double_array2d(long nx, long ny);

// free a 2d double array allocated by double_array2d()
int free_double_array2d(double **arr);

// allocate a float, double or int 3d array with dimensions nx, ny, nz
float ***float_array3d(long nx, long ny, long nz);
double ***double_array3d(long nx, long ny, long nz);
int ***int_array3d(long nx, long ny, long nz);

// free a float, double or int 3d array
int free_float_array3d(float ***arr);
int free_double_array3d(double ***arr);
int free_int_array3d(int ***arr);

}

#endif
