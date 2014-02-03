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
//     loadquadric.cpp           //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load quadrics
//

#include <iostream>
#include <string>
#include <cmath>
#include "xrmc.h"
#include "xrmc_math.h"
#include "xrmc_geom3d.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
// Method for loading quadric array
//////////////////////////////////////////////////////////////////////
int quadricarray::Load(istream &fs)
{
  string comm="";
  double xc[3], a[3], theta, elem;
  vect3 u;
  matr4 T, R, T1;
  int i;

  cout << "Quadric array file\n";
 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    if (comm=="End") break;
    else if(comm=="MaxNQuadr") { // set the maximum number of quadrics    
      GetIntToken(fs, &MaxNQuadr);
      cout << "Maximum number of quadrics: " << MaxNQuadr << "\n";
      delete[] Quadr;
      Quadr = new quadric[MaxNQuadr+1]; // allocate quadric array
      NQuadr = 0;
    }
    else if (comm=="Ellipsoid") { // define ellipsoid quadric
      cout << "Ellipsoid ";
      MapQuadric(fs);
      for (i=0; i<3; i++) GetDoubleToken(fs, &xc[i]); // center coordinates
      printf("\txc: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      for (i=0; i<3; i++) GetDoubleToken(fs, &a[i]); // semi-axes
      printf("\ta: %g\t%g\t%g\n", a[0], a[1], a[2]);
      Quadr[NQuadr].Ellipsoid(xc, a); // build ellipsoid quadric
      // Quadr[NQuadr].Print();
      NQuadr++;
    }
    else if (comm=="Plane") { // define plane quadric 
      cout << "Plane ";
      MapQuadric(fs);
      for (i=0; i<3; i++) GetDoubleToken(fs, &xc[i]); // point coordinates
      printf("\txc: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      for (i=0; i<3; i++) GetDoubleToken(fs, &a[i]);  // normal vector
      printf("\ta: %g\t%g\t%g\n", a[0], a[1], a[2]);
      Quadr[NQuadr].Plane(xc, a); // build plane
      // Quadr[NQuadr].Print();
      NQuadr++;
    }
    else if (comm=="CylinderX") { // cylinder parallel to x
      cout << "CylinderX ";
      MapQuadric(fs);
      for (i=0; i<2; i++) GetDoubleToken(fs, &xc[i]); // y,z coordinates
      printf("\txc: %g\t%g\n", xc[0], xc[1]);
      for (i=0; i<2; i++) GetDoubleToken(fs, &a[i]); // Ry, Rz
      printf("\ta: %g\t%g\n", a[0], a[1]);
      Quadr[NQuadr].CilinderX(xc, a); // build cylinder
      // Quadr[NQuadr].Print();
      NQuadr++;
    }
    else if (comm=="CylinderY") { // cylinder parallel to y
      cout << "CylinderY ";
      MapQuadric(fs);
     for (i=0; i<2; i++) GetDoubleToken(fs, &xc[i]); // x,z coordinates
     printf("\txc: %g\t%g\n", xc[0], xc[1]);
     for (i=0; i<2; i++) GetDoubleToken(fs, &a[i]); // Rx, Rz
     printf("\ta: %g\t%g\n", a[0], a[1]);
     Quadr[NQuadr].CilinderY(xc, a); // build cylinder
     // Quadr[NQuadr].Print();
     NQuadr++;
    }
    else if (comm=="CylinderZ") { // cylinder parallel to z
      cout << "CylinderZ ";
      MapQuadric(fs);
      for (i=0; i<2; i++) GetDoubleToken(fs, &xc[i]); // x,y coordinates
      printf("\txc: %g\t%g\n", xc[0], xc[1]);
      for (i=0; i<2; i++) GetDoubleToken(fs, &a[i]); // Rx, Ry
      printf("\ta: %g\t%g\n", a[0], a[1]);
      Quadr[NQuadr].CilinderZ(xc, a); // build cylinder
      // Quadr[NQuadr].Print();
      NQuadr++;
    }
    else if (comm=="Quadric") { // generic quadric
      cout << "Quadric ";
      MapQuadric(fs);
      for (int i=0; i<4; i++) { // loop on elements
	for (int j=i; j<4; j++) { // j>=i since quadric is symmetric
	  GetDoubleToken(fs, &elem);
	  Quadr[NQuadr].SetElem(i,j, elem); // fill quadric elements 
	}
      }
      Quadr[NQuadr].Print();
      NQuadr++;
    }
    // translation of last defined quadric
    else if (comm=="Translate" && NQuadr>0) {
      cout << "Translate\n";
      //Quadr[NQuadr-1].Print();
      T = matr4::Identity(); // start from identity matrix
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &xc[i]); // translation components
	T.Elem[i][3] = -xc[i]; // set elements of translation matrix
      }
      printf("\tdx: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      Quadr[NQuadr-1].Transform(T); // congruence transform
      // Quadr[NQuadr-1].Print();
    }
    else if (comm=="Rotate" && NQuadr>0) { // rotation of last defined quadric
      cout << "Rotate\n";
      //Quadr[NQuadr-1].Print();
      T1 = T = matr4::Identity(); // start from identity matrix
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &xc[i]); // translation components
	T.Elem[i][3] = -xc[i]; // set elements of translation matrix
	T1.Elem[i][3] = xc[i]; // set elements of translation matrix
      }
      printf("\txc: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &u.Elem[i]); // get rotation axis direction
      }
      printf("\tu: %g\t%g\t%g\n", u.Elem[0], u.Elem[1], u.Elem[2]);
      u.Normalize(); // normalize axis direction
      GetDoubleToken(fs, &theta); // get angle in degrees
      cout << "\ttheta: " << theta <<endl;
      R = matr4::RotMatr(u, theta); // build rotation matrix
      Quadr[NQuadr-1].Transform(T1); // translate: rotation axis-> origin
      Quadr[NQuadr-1].Transform(R); // congruence transform: rotate
      Quadr[NQuadr-1].Transform(T); // translate back
     // Quadr[NQuadr-1].Print();
    }
    else if (comm=="ChangeSign" && NQuadr>0) { // change sign of all elements
      cout << "ChangeSign\n";
      // Quadr[NQuadr-1].Print();
      Quadr[NQuadr-1].ChangeSign();
      // Quadr[NQuadr-1].Print();
    }
    //translate all previous quadrics
    else if (comm=="TranslateAll" && NQuadr>0) {
      cout << "TranslateAll\n";
      T = matr4::Identity(); // start from identity matrix
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &xc[i]);
	T.Elem[i][3] = -xc[i];  // translation components
      }
      printf("\tdx: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      for (int iq=0; iq<NQuadr; iq++) { // loop on previous quadrics
	Quadr[iq].Transform(T); // congruence transform
      }
    }
    else if (comm=="RotateAll" && NQuadr>0) { //rotate all previous quadrics
      cout << "RotateAll\n";
      T1 = T = matr4::Identity(); // start from identity matrix
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &xc[i]); // translation components
	T.Elem[i][3] = -xc[i]; // set elements of translation matrix
	T1.Elem[i][3] = xc[i]; // set elements of translation matrix
      }
      printf("\txc: %g\t%g\t%g\n", xc[0], xc[1], xc[2]);
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &u.Elem[i]); // get rotation axis direction
      }
      printf("\tu: %g\t%g\t%g\n", u.Elem[0], u.Elem[1], u.Elem[2]);
      u.Normalize(); // normalize axis direction
      GetDoubleToken(fs, &theta); // get angle in degrees
      cout << "\ttheta: " << theta <<endl;
      R = matr4::RotMatr(u, theta); // build rotation matrix

      for (int iq=0; iq<NQuadr; iq++) { // loop on previous quadrics
	Quadr[iq].Transform(T1); // translate: rotation axis-> origin
	Quadr[iq].Transform(R); // congruence transform: rotate
	Quadr[iq].Transform(T); // translate back
      }
    }
    else {
      throw xrmc_exception("Syntax error in quadrics file\n");
      // unrecognized command
    }
    if (NQuadr>MaxNQuadr) {
      char i2ch[MAXSTRLEN];
      sprintf(i2ch, "%d", MaxNQuadr);
      string s_err="Number of quadrics greater than maximum: ";
      s_err = s_err + i2ch + "\nUse the command MaxNQuadr to "
	"increase the maximum number of quadrics";
      throw xrmc_exception(s_err);
    }
  }
  
  return 0;
}

// insert name and pointer to the quadric in the quadric map
int quadricarray::MapQuadric(istream &fs)
{
  string qname;
  quadric_map_insert_pair insert_pair;

  GetToken(fs, qname);
  cout << qname << endl;
  
  // insert name and pointer to the quadric in the quadric map
  insert_pair = QuadricMap.insert(quadric_map_pair(qname, &Quadr[NQuadr]));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Quadric ") + qname + 
			 " already inserted in quadric map\n");
  return 0;
}

//////////////////////////////////////////////////////////////////////
// set default values for quadricarray parameters
int quadricarray::SetDefault()
{
  MaxNQuadr = 10000; // maximum num. of quadrics
  Quadr = new quadric[MaxNQuadr+1];
  NQuadr = 0;

  return 0;
}
