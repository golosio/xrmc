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
//     quadric.cpp               //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class quadric
//

#include <cmath>
#include <iostream>
#include "xrmc_geom3d.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Default constructor
//////////////////////////////////////////////////////////////////////
quadric::quadric()
{
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      Matr.Elem[i][j] = 0;
    }
  }
  NInters=0;
  tInters[0]=tInters[1]=0;
  Enter[0]=Enter[1]=0;
}

//////////////////////////////////////////////////////////////////////
// Set elements i,j and j,i to the value elem
//////////////////////////////////////////////////////////////////////
int quadric::SetElem(int i, int j, double elem)
{
  Matr.Elem[i][j] = Matr.Elem[j][i] = elem;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Congruence transform of the quadric with matrix M
//////////////////////////////////////////////////////////////////////
int quadric::Transform(matr4 M)
{
  Matr=M.Transpose()*Matr*M;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Display the content of all elements
//////////////////////////////////////////////////////////////////////
int quadric::Print()
{
  int i, j;
  
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      cout << Matr.Elem[i][j] << "\t";
    }
    cout << "\n";
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// change the sign
//////////////////////////////////////////////////////////////////////
int quadric::ChangeSign()
{
  Matr = -Matr;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// evaluates the quadratic form Mij * xi * yj
// on 4d vectors
//////////////////////////////////////////////////////////////////////
double quadric::Prod4(double *x, double *y) {
  int i, j;
  double prod = 0;

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      prod += Matr.Elem[i][j]*x[i]*y[j];
    }
  }
  return prod;
}

//////////////////////////////////////////////////////////////////////
// evaluates the quadratic form Mij * xi * yj on 3d vectors
// (x0, x1, x2) and (y0, y1, y2) 
// setting x3=1 and y3=1
//////////////////////////////////////////////////////////////////////
double quadric::Prod3(double *x, double *y) {
  int i;
  double x4[4], y4[4];
  for (i=0; i<3; i++) {
    x4[i] = x[i];    
    y4[i] = y[i];
  }
  x4[3] = y4[3] = 1; // set the missing element to 1
  return Prod4(x4, y4);
}

//////////////////////////////////////////////////////////////////////
// check if point x is inside (0) or outside (1) the quadric
//////////////////////////////////////////////////////////////////////
int quadric::Inside(vect3 x)
{
  return (Prod3(x.Elem, x.Elem)<0);
}

//////////////////////////////////////////////////////////////////////
// build ellipsoid quadric
// x0: center coordinates
// a: semi-axes
//////////////////////////////////////////////////////////////////////
int quadric::Ellipsoid(double *x0, double *a) {
  int i, j;
  double inva2;
  
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Matr.Elem[i][j] = 0;
    }
  }
  Matr.Elem[3][3] = -1;
  for (i=0; i<3; i++) {
    inva2 = 1. / (a[i]*a[i]);
    Matr.Elem[i][i] = inva2;
    Matr.Elem[3][i] = Matr.Elem[i][3] = -x0[i]*inva2;
    Matr.Elem[3][3] += x0[i]*x0[i]*inva2;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// build plane quadric
// x0: point in the plane
// u: normal to the plane
//////////////////////////////////////////////////////////////////////
int quadric::Plane(double *x0, double *u) {
  int i, j;
  double d = 0;

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      Matr.Elem[i][j] = 0;
      
    }
  }
  
  for (i=0; i<3; i++) {
    Matr.Elem[i][3] = Matr.Elem[3][i] = 0.5*u[i];
    d -= u[i]*x0[i];
  }
  Matr.Elem[3][3] = d;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Cylinder parallel to x axis
// x0[]: y, z coordinates of the center
// a[]: Ry, Rz (semi-axes of elliptical section)
////////////////////////////////////////////////////////////////////// 
int quadric::CilinderX(double *x0, double *a) {
  int i, j;
  double inva2;
  
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Matr.Elem[i][j] = 0;
    }
  }
  Matr.Elem[3][3] = -1;
  inva2 = 1. / (a[0]*a[0]);
  Matr.Elem[1][1] = inva2;
  Matr.Elem[3][1] = Matr.Elem[1][3] = -x0[0]*inva2;
  Matr.Elem[3][3] += x0[0]*x0[0]*inva2;

  inva2 = 1. / (a[1]*a[1]);
  Matr.Elem[2][2] = inva2;
  Matr.Elem[3][2] = Matr.Elem[2][3] = -x0[1]*inva2;
  Matr.Elem[3][3] += x0[1]*x0[1]*inva2;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Cylinder parallel to y axis
// x0[]: x, z coordinates of the center
// a[]: Rx, Rz (semi-axes of elliptical section)
////////////////////////////////////////////////////////////////////// 
int quadric::CilinderY(double *x0, double *a) {
  int i, j;
  double inva2;
  
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Matr.Elem[i][j] = 0;
    }
  }
  Matr.Elem[3][3] = -1;
  inva2 = 1. / (a[0]*a[0]);
  Matr.Elem[0][0] = inva2;
  Matr.Elem[3][0] = Matr.Elem[0][3] = -x0[0]*inva2;
  Matr.Elem[3][3] += x0[0]*x0[0]*inva2;

  inva2 = 1. / (a[1]*a[1]);
  Matr.Elem[2][2] = inva2;
  Matr.Elem[3][2] = Matr.Elem[2][3] = -x0[1]*inva2;
  Matr.Elem[3][3] += x0[1]*x0[1]*inva2;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Cylinder parallel to z axis
// x0[]: x, y coordinates of the center
// a[]: Rx, Ry (semi-axes of elliptical section)
////////////////////////////////////////////////////////////////////// 
int quadric::CilinderZ(double *x0, double *a) {
  int i, j;
  double inva2;
  
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Matr.Elem[i][j] = 0;
    }
  }
  Matr.Elem[3][3] = -1;
  inva2 = 1. / (a[0]*a[0]);
  Matr.Elem[0][0] = inva2;
  Matr.Elem[3][0] = Matr.Elem[0][3] = -x0[0]*inva2;
  Matr.Elem[3][3] += x0[0]*x0[0]*inva2;

  inva2 = 1. / (a[1]*a[1]);
  Matr.Elem[1][1] = inva2;
  Matr.Elem[3][1] = Matr.Elem[1][3] = -x0[1]*inva2;
  Matr.Elem[3][3] += x0[1]*x0[1]*inva2;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Find intersections of the straight line x0 + u*t with the quadric
//////////////////////////////////////////////////////////////////////
int quadric::Intersect(vect3 x0, vect3 u) {
  double x4[4], u4[4];
  double det;

  // transform x0 and u in the proper 4d vectors
  for (int i=0; i<3; i++) {
    x4[i]=x0.Elem[i];
    u4[i] = u.Elem[i];
  }
  x4[3] = 1;
  u4[3] = 0;
  
  // find coefficients of reduced formula for 2nd degree equation
  double a = Prod4(u4, u4);
  double beta = Prod4(x4, u4);
  double c = Prod4(x4, x4);

  NInters = 0;

  if ((beta+a)==beta) {  // a is 0 or comparable to 0
    if (c+beta != c) { // check that beta is not too small
      tInters[0] = -c/(2.*beta); // only one solution in this case
      NInters++; // 1 solution
    }
  }
  else { // a is not zero
    det = beta*beta - a*c; // discriminant of 2nd degree equation
    if (det > 0) { // two solutions
      det = sqrt(det);
      double sol1 = -beta + det;
      double sol2 = -beta - det;
      // if (sol1 == sol2) {
      //   tInters[0] = sol1 / a;
      //   NInters++; Reconsider in case of nonconvex quadrics
      // } 
      if (sol1 != sol2) {
	if ((sol1+a) != sol1) { // check that denominator is not too small
	  tInters[0] = sol1 / a; // first intersection
	  NInters++;
	} 
	if ((sol2+a) != sol2) { // check that denominator is not too small
	  tInters[NInters] = sol2 / a; // second intersection
	  NInters++; // two solutions
	} 
      }
    }
  }
  for (int i=0; i<NInters; i++) {
    // check the crossing direction of each intersection:
    // from inside to outside (0) or from outside to inside (1)
    Enter[i] = (a*tInters[i]+beta)<0;
  }

  cout << "tInters[0]" << tInters[0] << "\n";
  cout << "Enter[0]" << Enter[0] << "\n";

  return 0;
}

quadricarray *quadricarray::Clone(string dev_name) {
	cout << "Entering quadricarray::Clone\n";
	quadricarray *clone = new quadricarray(dev_name);
	clone->NQuadr = NQuadr;
	clone->MaxNQuadr = MaxNQuadr;
	clone->Quadr = new quadric[NQuadr];
	for (int i = 0 ; i < NQuadr ; i++)
		clone->Quadr[i] = Quadr[i];
	
	//typedef map<string, quadric*> quadric_map;
	//typedef pair<string, quadric*> quadric_map_pair;
	//typedef pair<quadric_map::iterator, bool> quadric_map_insert_pair;
	quadric_map::iterator it;

	int j = 0;
	//for (it = QuadricMap.begin() ; it != QuadricMap.end() ; it++) {
	//	clone->QuadricMap.insert(quadric_map_pair(it->first, &clone->Quadr[j++]));
	//}
	clone->QuadricMap = QuadricMap;
	return clone;
}


