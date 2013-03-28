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
//     loadsource.cpp            //
//        14/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load source parameters, position and orientation
//

#include <cmath>
#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_source.h"
#include "xrmc_spectrum.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "xrmc_device.h"
#include "xrmc_algo.h"

using namespace std;
using namespace gettoken;
using namespace xrmc_algo;

//////////////////////////////////////////////////////////////////////
int source::Load(FILE *fp)
  // Loads source position/orientation
{
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  char comm_ch[MAXSTRLEN], s[MAXSTRLEN];
  string comm="";

  cout << "Source position/orientation file\n";

  while (!feof(fp) && comm!="End") {
    GetToken(fp, comm_ch); // get a command/variable name from input file
    comm = comm_ch;
    // parse the command and decide what to do
    if(comm=="SpectrumName") { // set the spectrum input device name
      GetToken(fp, s);
      InputDeviceName[0] = s;
      cout << "Spectrum input device name: " << InputDeviceName[0] << "\n";
    } 
    else if(comm=="X") { // set the source coordinates
      cout << "Source position :\t";
      for (i=0; i<3; i++) {
	GetDoubleToken(fp, &X.Elem[i]);
	//cout << X[i] << "\t";
      }
      cout << X << endl;
      //cout << "\n";
    }
    else if(comm=="uk") { // uk vector, local z axis (main source direction)
      cout << "Source orientation :\n"; 
      cout << "\tuk vector (local z axis, main source direction):\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fp, &uk.Elem[i]);
	//cout << uk[i] << "\t";
      }
      cout << uk << endl;
      //cout << "\n";
    }	
    else if(comm=="ui") { // ui vector, local x axis
      cout << "Source orientation :\n"; 
      cout << "\tui vector (local x axis direction):\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fp, &ui.Elem[i]);
	//cout << ui[i] << "\t";
      }
      cout << ui << endl;
      //cout << "\n";
    }	
    else if(comm=="Divergence") { // Beam divergence (thetax, thetay)
      GetDoubleToken(fp, &Thx);
      GetDoubleToken(fp, &Thy);
      cout << "Beam divergence (thetax, thetay): "
	   << Thx << "\t" << Thy << "\n";
    }
    else if(comm=="Size") { // Source size (sigmax, sigmay, sigmaz)
                            // in local coordinate system
      SizeFlag = 1;
      GetDoubleToken(fp, &Sigmax);
      GetDoubleToken(fp, &Sigmay);
      GetDoubleToken(fp, &Sigmaz);
      cout << "Source size (sigmax, sigmay, sigmaz in local source coordinate "
	"system):\n" << Sigmax << ", " << Sigmay << ", " << Sigmaz << "\n";
    }
    else if(comm=="Rotate") { // source rotation
      cout << "Source rotation :\n"; 
      cout << "\tPoint on rotation axis x0:\t";
      for (int i=0; i<3; i++) {
	GetDoubleToken(fp, &x0.Elem[i]); // translation components
      }
      cout << x0 << endl;
      cout << "\tRotation axis direction u:\t"; 
      for (int i=0; i<3; i++) {
	GetDoubleToken(fp, &u.Elem[i]); // get rotation axis direction
      }
      u.Normalize(); // normalize axis direction
      cout << u << endl;
      cout << "\tRotation angle theta (degrees): ";
      GetDoubleToken(fp, &theta); // get angle in degrees
      cout << theta << endl;
      R = matr3::RotMatr(u, -theta); // build rotation matrix
      X -= x0; // translate source position: rotation axis-> origin
      X = R*X; // rotate source position
      uk = R*uk; // rotate source uk direction
      ui = R*ui; // rotate source ui direction
      X += x0; // translate back source position
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in source input file"); 
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
int source::SetDefault()
// set default values for source parameters
{
  InputDeviceName[0] = "Spectrum"; // input spectrum device name

  X.Set(0, -50, 0); // Source position coordinates 
  // Source orientation :
  uk.Set(0,1,0); // uk is directed as the y axis
  ui.Set(-1,0,0);// ui has direction opposite the x axis
  Thx = 0.1; // divergence (thetax)
  Thy = 0.1; // divergence (thetay)
  SizeFlag = 0;                 // Source size (sigmax, sigmay, sigmaz)      
  Sigmax = Sigmay = Sigmaz = 0; // in local coordinate system
  OrthoNormal(ui, uj, uk); // evaluates uj to form a orthonormal basis

  return 0;
}

