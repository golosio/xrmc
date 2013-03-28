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
///////////////////////////////////
//     xrmc_math.h               //
//        07/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
//  vect3 and matr4 classes definitions
//

#ifndef XRMCMATHH
#define XRMCMATHH

#ifndef PI
#define PI 3.14159265359
#endif

//////////////////////////////////////////////////////////////////////
// vect3 class definition, derived from vect
// 3d real vector with related operations
//////////////////////////////////////////////////////////////////////
class vect3
{
 public:
  double Elem[3]; // vector components
  vect3() {} // constructor
  vect3(double v[]); // constructor from double array
  vect3(double x, double y, double z); // constructor from 3 doubles
  int Set(double  v[]); // set elements from double array
  int Set(double x, double y, double z); // set elements from 3 doubles
  vect3 operator*(double d); // product of the vector by a scalar
  vect3 operator/(double d); // division of the vector by a scalar
  vect3 operator+(vect3 v1); // sum of two vectors
  vect3 operator-(vect3 v1); // difference between two vectors
  double operator*(vect3 v1); // scalar product of two vectors
  int operator+=(vect3 v1); // sum a vector
  int operator-=(vect3 v1); // subtract a vector
  vect3 operator^(vect3 v1); // vector product
  int Normalize(); // normalize vector module to 1
  double Mod(); // module
  // insertion (output) operator
  friend std::ostream &operator<<(std::ostream &stream, vect3 v);

  // find 3d polar coordinates from cartesian coordinates of position vector
  int CartesianToPolar(double *r, double *theta, double *phi);

};

//////////////////////////////////////////////////////////////////////
// matr3 class definition
// 3 x 3 real matrix with related operations
//////////////////////////////////////////////////////////////////////
class matr3
{
 public:
  double Elem[3][3]; // matrix elements
  matr3() {} // constructor
  matr3(double **m1); // constructor from pointer to pointers
  int Set(double **m1); // set elements from pointer to pointers
  matr3 operator*(const matr3& m1); // 3 x 3 matrix multiplication
  matr3 operator-(); // minus operator
  matr3 Transpose(); // matrix transpose operator
  vect3 operator*(const vect3 &v1); // multiplication by a 3d vector
  // rotation matrix of an angle theta around axis with direction u
  static matr3 RotMatr(vect3 u, double theta); 
  static matr3 Identity(); // 3 x 3 identity matrix
  // 3 x 3 matrix insertion (output) operator
  friend std::ostream &operator<<(std::ostream &stream, matr3 m);
};

//////////////////////////////////////////////////////////////////////
// matr4 class definition
// 4 x 4 real matrix with related operations
//////////////////////////////////////////////////////////////////////
class matr4
{
 public:
  double Elem[4][4]; // matrix elements
  matr4() {} // constructor
  matr4(double **m1); // constructor from pointer to pointers
  int Set(double **m1); // set elements from pointer to pointers
  matr4 operator*(const matr4& m1); // 4 x 4 matrix multiplication
  matr4 operator-(); // minus operator
  matr4 Transpose(); // matrix transpose operator
  // rotation matrix of an angle theta around axis with direction u
  static matr4 RotMatr(vect3 u, double theta); 
  static matr4 Identity(); // 4 x 4 identity matrix
  // 4 x 4 matrix insertion (output) operator
  friend std::ostream &operator<<(std::ostream &stream, matr4 m);
  friend bool operator==(matr4 &Matr1, matr4 &Matr2);
};

// build a orthonormal basis from the two vectors vk, vi
int OrthoNormal(vect3 &vi, vect3 &vj, vect3 &vk);

#endif

