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
//     xrmc_math.cpp             //
//        07/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
//  methods of the classes vect3 and matr4
//

#include <iostream>
#include <cstring>
#include <cmath>
#include "xraylib.h"
#include "xrmc_math.h"
#include "xrmc_exception.h"

////////////////////////////////////////////////////////////////////////////
// return angle phi in 2d polar coordinates from cartesian coordinates x, y
////////////////////////////////////////////////////////////////////////////
double Angle(double x, double y)
{
  double Phi;

  if (y+x == y) {
    if (y > 0) Phi = PI/2;
    else Phi = 3*PI/2;
  }
  else {
    Phi = atan(y/x);
    if (x < 0) Phi += PI;
    if (Phi < 0) Phi += 2*PI;
  }
  return Phi;
}

//////////////////////////////////////////////////////////////////////
// build a othonormal basis from the two vectors vk, vi
//////////////////////////////////////////////////////////////////////
int OrthoNormal(vect3 &vi, vect3 &vj, vect3 &vk)
{
  vk.Normalize();
  vj=vk^vi;
  vj.Normalize();
  vi=vj^vk;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector constructor from 3 doubles
//////////////////////////////////////////////////////////////////////
vect3::vect3(double x, double y, double z)
{
  Elem[0]=x;
  Elem[1]=y;
  Elem[2]=z;
}


//////////////////////////////////////////////////////////////////////
// 3d real vector constructor from double array
//////////////////////////////////////////////////////////////////////
vect3::vect3(double v[])
{
  Elem[0]=v[0];
  Elem[1]=v[1];
  Elem[2]=v[2];
}

//////////////////////////////////////////////////////////////////////
// 3d real vector set elements from 3 doubles
//////////////////////////////////////////////////////////////////////
int vect3::Set(double x, double y, double z)
{
  Elem[0]=x;
  Elem[1]=y;
  Elem[2]=z;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector set elements from double array
//////////////////////////////////////////////////////////////////////
int vect3::Set(double v[])
{
  Elem[0]=v[0];
  Elem[1]=v[1];
  Elem[2]=v[2];
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector insertion (output) operator
//////////////////////////////////////////////////////////////////////
std::ostream &operator<<(std::ostream &stream, vect3 v)
{
  stream<<"("<<v.Elem[0]<<";"<<v.Elem[1]<<";"<<v.Elem[2]<<")";

  return stream;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector multiplication by a scalar
//////////////////////////////////////////////////////////////////////
vect3 vect3::operator*(double d)
{
  vect3 v1(Elem[0]*d,Elem[1]*d,Elem[2]*d);
  
  return v1;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector division by a scalar
//////////////////////////////////////////////////////////////////////
vect3 vect3::operator/(double d)
{
  vect3 v1(Elem[0]/d,Elem[1]/d,Elem[2]/d);
  
  return v1;
}

//////////////////////////////////////////////////////////////////////
// 3d real vectors sum
//////////////////////////////////////////////////////////////////////
vect3 vect3::operator+(vect3 v1)
{
  vect3 v2(Elem[0]+v1.Elem[0], Elem[1]+v1.Elem[1], Elem[2]+v1.Elem[2]);
  
  return v2;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector difference
//////////////////////////////////////////////////////////////////////
vect3 vect3::operator-(vect3 v1)
{
  vect3 v2(Elem[0]-v1.Elem[0], Elem[1]-v1.Elem[1], Elem[2]-v1.Elem[2]);
  
  return v2;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector sum
//////////////////////////////////////////////////////////////////////
int vect3::operator+=(vect3 v1)
{
  Elem[0] += v1.Elem[0];
  Elem[1] += v1.Elem[1];
  Elem[2] += v1.Elem[2];

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector difference
//////////////////////////////////////////////////////////////////////
int vect3::operator-=(vect3 v1)
{
  Elem[0] -= v1.Elem[0];
  Elem[1] -= v1.Elem[1];
  Elem[2] -= v1.Elem[2];
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vectors scalar product
//////////////////////////////////////////////////////////////////////
double vect3::operator*(vect3 v1)
{
  double d=Elem[0]*v1.Elem[0]+Elem[1]*v1.Elem[1]+Elem[2]*v1.Elem[2];

  return d;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector module
//////////////////////////////////////////////////////////////////////
double vect3::Mod()
{
  double mod = Elem[0]*Elem[0]+Elem[1]*Elem[1]+Elem[2]*Elem[2];

  return sqrt(mod);
}

//////////////////////////////////////////////////////////////////////
// 3d real vector normalize module to 1
//////////////////////////////////////////////////////////////////////
int vect3::Normalize()
{
  double vnorm = Mod(); 

  if (vnorm==0) return 1;
  Elem[0] /= vnorm;
  Elem[1] /= vnorm;
  Elem[2] /= vnorm;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3d real vector: vector product
//////////////////////////////////////////////////////////////////////
vect3 vect3::operator^(vect3 v1)
{
  vect3 v2(Elem[1]*v1.Elem[2]-Elem[2]*v1.Elem[1],
	   Elem[2]*v1.Elem[0]-Elem[0]*v1.Elem[2],
	   Elem[0]*v1.Elem[1]-Elem[1]*v1.Elem[0]);

  return v2;
}

//////////////////////////////////////////////////////////////////////
// evaluate 3d polar coordinates r, theta, phi of 3d position vector
//////////////////////////////////////////////////////////////////////
int vect3::CartesianToPolar(double *r, double *theta, double *phi)
{
  *theta = *phi = 0;
  *r = sqrt(Elem[0]*Elem[0] + Elem[1]*Elem[1] + Elem[2]*Elem[2]);
  if ((Elem[2] + *r) == Elem[2]) return 0;
  *phi = Angle(Elem[0], Elem[1]);
  *theta = acos(Elem[2] / (*r));

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix constructor from pointer to pointers
//////////////////////////////////////////////////////////////////////
matr3::matr3(double **m1)
{
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Elem[i][j]=m1[i][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix set element values from pointer to pointers
//////////////////////////////////////////////////////////////////////
int matr3::Set(double **m1)
{
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      Elem[i][j]=m1[i][j];
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix insertion operator (output)
//////////////////////////////////////////////////////////////////////
std::ostream &operator<<(std::ostream &stream, matr3 m)
{
  for(int i=0; i<3; i++) {
    stream<<"(";
    for(int j=0; j<2; j++) {
      stream<<m.Elem[i][j]<<";";
    }
    stream<<m.Elem[2]<<")";
  }

  return stream;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix multiplication
//////////////////////////////////////////////////////////////////////
matr3 matr3::operator*(const matr3 &m1)
{
  matr3 m2;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      m2.Elem[i][j] = Elem[i][0]*m1.Elem[0][j]+
	Elem[i][1]*m1.Elem[1][j]+
	Elem[i][2]*m1.Elem[2][j];
    }
  }
  
  return m2;
}


//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix minus operator
//////////////////////////////////////////////////////////////////////
matr3 matr3::operator-()
{
  matr3 m1;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      m1.Elem[i][j] = -Elem[i][j];
    }
  }
  
  return m1;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix transpose operator
//////////////////////////////////////////////////////////////////////
matr3 matr3::Transpose()
{
  matr3 m1;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      m1.Elem[i][j] = Elem[j][i];
    }
  }
  
  return m1;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 real matrix multiplication by a 3d vector
//////////////////////////////////////////////////////////////////////
vect3 matr3::operator*(const vect3 &v1)
{
  vect3 v2;

  for (int i=0; i<3; i++) {
    v2.Elem[i] = Elem[i][0]*v1.Elem[0]+
      Elem[i][1]*v1.Elem[1]+
      Elem[i][2]*v1.Elem[2];
  }
  
  return v2;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 matrix of rotation of an angle theta around axis with direction u
//////////////////////////////////////////////////////////////////////
matr3 matr3::RotMatr(vect3 u, double theta)
{
  matr3 m;

  double ux=u.Elem[0], uy=u.Elem[1], uz=u.Elem[2];
  double cth = cos(theta*PI/180);
  double sth = sin(theta*PI/180);

  m.Elem[0][0] = cth + ux*ux*(1-cth);
  m.Elem[0][1] = ux*uy*(1-cth) - uz*sth;
  m.Elem[0][2] = ux*uz*(1-cth) + uy*sth;

  m.Elem[1][0] = uy*ux*(1-cth) + uz*sth;
  m.Elem[1][1] = cth + uy*uy*(1-cth);
  m.Elem[1][2] = uy*uz*(1-cth) - ux*sth;

  m.Elem[2][0] = uz*ux*(1-cth) - uy*sth;
  m.Elem[2][1] = uz*uy*(1-cth) + ux*sth;
  m.Elem[2][2] = cth + uz*uz*(1-cth);

  return m;
}

//////////////////////////////////////////////////////////////////////
// 3 x 3 identity matrix
//////////////////////////////////////////////////////////////////////
matr3 matr3::Identity()
{
  matr3 m;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      m.Elem[i][j] = (i==j) ? 1 : 0;
    }
  }

  return m;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix constructor from pointer to pointers
//////////////////////////////////////////////////////////////////////
matr4::matr4(double **m1)
{
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      Elem[i][j]=m1[i][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix set element values from pointer to pointers
//////////////////////////////////////////////////////////////////////
int matr4::Set(double **m1)
{
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      Elem[i][j]=m1[i][j];
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix insertion operator (output)
//////////////////////////////////////////////////////////////////////
std::ostream &operator<<(std::ostream &stream, matr4 m)
{
  for(int i=0; i<4; i++) {
    stream<<"(";
    for(int j=0; j<3; j++) {
      stream<<m.Elem[i][j]<<";";
    }
    stream<<m.Elem[3]<<")";
  }

  return stream;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix multiplication
//////////////////////////////////////////////////////////////////////
matr4 matr4::operator*(const matr4 &m1)
{
  matr4 m2;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      m2.Elem[i][j] = Elem[i][0]*m1.Elem[0][j]+
	Elem[i][1]*m1.Elem[1][j]+
	Elem[i][2]*m1.Elem[2][j]+
	Elem[i][3]*m1.Elem[3][j];
    }
  }
  
  return m2;
}


//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix minus operator
//////////////////////////////////////////////////////////////////////
matr4 matr4::operator-()
{
  matr4 m1;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      m1.Elem[i][j] = -Elem[i][j];
    }
  }
  
  return m1;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 real matrix transpose operator
//////////////////////////////////////////////////////////////////////
matr4 matr4::Transpose()
{
  matr4 m1;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      m1.Elem[i][j] = Elem[j][i];
    }
  }
  
  return m1;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 matrix of rotation of an angle theta around axis with direction u
//////////////////////////////////////////////////////////////////////
matr4 matr4::RotMatr(vect3 u, double theta)
{
  matr4 m;

  double ux=u.Elem[0], uy=u.Elem[1], uz=u.Elem[2];
  double cth = cos(theta*PI/180);
  double sth = sin(theta*PI/180);

  m.Elem[0][0] = cth + ux*ux*(1-cth);
  m.Elem[0][1] = ux*uy*(1-cth) - uz*sth;
  m.Elem[0][2] = ux*uz*(1-cth) + uy*sth;

  m.Elem[1][0] = uy*ux*(1-cth) + uz*sth;
  m.Elem[1][1] = cth + uy*uy*(1-cth);
  m.Elem[1][2] = uy*uz*(1-cth) - ux*sth;

  m.Elem[2][0] = uz*ux*(1-cth) - uy*sth;
  m.Elem[2][1] = uz*uy*(1-cth) + ux*sth;
  m.Elem[2][2] = cth + uz*uz*(1-cth);

  for (int i=0; i<4; i++) {
    m.Elem[3][i] = m.Elem[i][3] = 0;
  }
  m.Elem[3][3] = 1;

  return m;
}

//////////////////////////////////////////////////////////////////////
// 4 x 4 identity matrix
//////////////////////////////////////////////////////////////////////
matr4 matr4::Identity()
{
  matr4 m;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      m.Elem[i][j] = (i==j) ? 1 : 0;
    }
  }

  return m;
}

bool operator==(matr4 &Matr1, matr4 &Matr2) {
	if (memcmp(Matr1.Elem, Matr2.Elem, sizeof(double)*4*4) == 0)
		return true;
	else
		return false;
}
