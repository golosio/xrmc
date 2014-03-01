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
//     loadgeom3d.cpp            //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load 3d geometric shapes
//

#include <iostream>
#include <string>
#include <algorithm>
#include "xrmc.h"
#include "xrmc_geom3d.h"
#include "xrmc_exception.h"
#include "xrmc_gettoken.h"
#include "xrmc_composition.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
// Method for loading 3d geometric shapes
//////////////////////////////////////////////////////////////////////
int geom3d::Load(istream &fs)
{
  string comm="", s, qname, s_ph_in, s_ph_out;
  int n_quadr, i;

  cout << "Geometric shapes file\n";

 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    // check if it's a command for setting an input device name
    if (ParseInputDeviceCommand(fs, comm)) continue;
    else if(comm=="X") { // set the sample region center coordinates
      cout << "Sample region center coordinates: \t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &X.Elem[i]);
	//cout << X[i] << "\t";
      }
      cout << X << endl;
      //cout << "\n";
    }
    else if(comm=="HW") { // set the sample region half-sides
      cout << "Sample region half-sides: \t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &HW[i]);
	cout << HW[i] << "\t";
      }
      cout << "\n";
    }
    else if(comm=="MaxNQVol") { // set the maximum number of 3d objects    
      GetIntToken(fs, &MaxNQVol);
      cout << "Maximum number of 3d objects: " << MaxNQVol << "\n";
      delete[] QVol;
      delete[] QVolMap;
      QVol = new qvolume[MaxNQVol+1]; // initialize 3d object array
      QVolMap = new string*[MaxNQVol+1]; // initialize 3d object array map
      NQVol = 0;
    }
    else if(comm=="Object") { // read parameters for a new 3d object
      GetToken(fs, s); // object name
      //ObjName = s;
      cout << "Object: " << s << endl;
      GetToken(fs, s_ph_in); // set internal phase
      GetToken(fs, s_ph_out);  // set external phase
      GetIntToken(fs, &n_quadr); // num. of quadrics delimiting the 3d object
      QVol[NQVol].Init(s_ph_in, s_ph_out, n_quadr); // initialize 3d object
      QVolMap[NQVol] = new string[n_quadr]; // initialize map of quadrics
      // delimiting the 3d object
      cout << "\tPhase Inside: " << QVol[NQVol].PhaseInName << endl;
      cout << "\tPhase Outside: " << QVol[NQVol].PhaseOutName << endl;
      cout << "\tNum. of quadrics: " << QVol[NQVol].NQuadr << endl;
     
      if (find(used_phases.begin(), used_phases.end(), s_ph_in) == used_phases.end()) {
      	used_phases.push_back(s_ph_in);
      }
      if (find(used_phases.begin(), used_phases.end(), s_ph_out) == used_phases.end()) {
      	used_phases.push_back(s_ph_out);
      }

      for(i=0; i<n_quadr; i++) { // loop on quadrics delimiting the 3d object
	GetToken(fs, s);
	qname = s; // quadric name
	QVolMap[NQVol][i] = qname; // put the quadric in the map
	cout << "\tQuadric: " << qname << "\n";
      }
      NQVol++;
      if (NQVol>MaxNQVol) {
	char i2ch[MAXSTRLEN];
	sprintf(i2ch, "%d", MaxNQVol);
	string s_err="Number of 3d objects greater than maximum: ";
	s_err = s_err + i2ch + "\nUse the command MaxNQVol to "
	  "increase the maximum number of objects";
	throw xrmc_exception(s_err);
      }
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in geometric shapes input file"); 
    }
  }
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// set default values for geom3d parameters
int geom3d::SetDefault()
{
  InputDeviceName[0] = "QuadricArray"; // Quadric array input device name
  InputDeviceName[1] = "Composition"; // composition device name
  X.Set(0,0,0); // Sample region center coordinates
  HW[0] = HW[1] = HW[2] = 20; // Sample region half-sides
  MaxNQVol = 10000; // maximum num. of 3d objects
  QVol = new qvolume[MaxNQVol+1]; // initialize 3d object array
  QVolMap = new string*[MaxNQVol+1]; // initialize 3d object array map
  NQVol = 0; // num. of 3d objects

  return 0;
}
