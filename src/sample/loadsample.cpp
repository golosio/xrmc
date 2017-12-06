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
//     loadsample.cpp          //
//        08/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load sample interaction parameters
//
#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "xrmc_sample.h"
#include "xrmc_photon.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
int sample::Load(istream &fs)
  // Loads sample parameters file
{
  string comm="";

  cout << "Sample parameters file\n";

  // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    // check if it's a command for setting an input device name
    if (ParseInputDeviceCommand(fs, comm)) continue;
    // flag for weighted step length
    else if(comm=="WeightedStepLength") { 
      GetIntToken(fs, &WeightedStepLength);
      cout << "Weighted step length (0/1): " << WeightedStepLength << "\n";
    }
    // flag for activating/deactivating fluorescent emission
    else if(comm=="FluorFlag") { 
      GetIntToken(fs, &photon::FluorFlag);
      cout << "Activate Fluorescence (0/1): " << photon::FluorFlag << "\n";
    }
    else if(comm=="ScattOrderNum") { // Maximum scattering order
      GetIntToken(fs, &ScattOrderNum);
      cout << "Maximum scattering order: " << ScattOrderNum << "\n";
      ScattOrderNum++;
      // remove previously allocated array
      if (PhotonNum!=NULL) delete[] PhotonNum;
      PhotonNum = new int[ScattOrderNum]; // allocate new array
      for (int i=0; i<ScattOrderNum; i++) {
	GetIntToken(fs, &PhotonNum[i]); // Multiplicity of simulated events
	cout << "Multiplicity of simulated events for scattering order "
	     << i << ": " << PhotonNum[i] << "\n";
      }
    }
    else if(comm=="Rlim") { // limit on x0-x1 distance in forced detection
                // below which a distance-weighted integration method is used
      GetDoubleToken(fs, &Rlim);
      cout << "Rlim: " << Rlim << "\n";
    }
    else if(comm=="RWFract") { // fraction of the events in forced detection
                // on which a distance-weighted integration method is used
      GetDoubleToken(fs, &RWFract);
      cout << "RWFract: " << RWFract << "\n";
    }
    else if(comm=="RlimDet") { // limit on x0-x1 distance in forced detection
             // for the last step before detection
             // below which a distance-weighted integration method is used
      GetDoubleToken(fs, &RlimDet);
      cout << "RlimDet: " << RlimDet << "\n";
    }
    else if(comm=="RWFractDet") { // fraction of the events in forced detection
              // for the last step before detection
              // on which a distance-weighted integration method is used
      GetDoubleToken(fs, &RWFractDet);
      cout << "RWFractDet: " << RWFractDet << "\n";
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empty string\n";
    }
    else {
      throw xrmc_exception("syntax error in sample input file"); 
    }
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
int sample::SetDefault()
// set default values for sample parameters
{
  InputDeviceName[0] = "Source"; // source input device name
  InputDeviceName[1] = "Geom3D"; // geom3d input device name
  InputDeviceName[2] = "Composition"; // composition device name
  WeightedStepLength = 0; // flag for weighted step length
  photon::FluorFlag = 1; // flag for activating fluorescent emission
  ScattOrderNum = 1; // Num. of scattering orders (including transmission)
  PhotonNum = new int[1];
  PhotonNum[0] = 1; // Multiplicity of simulated events for transmission
  Rlim = 0.1; // 1 mm
  RWFract = 0.1; // 10% of the events
  RlimDet = 0.1; // 1 mm
  RWFractDet = 0.1; // 10% of the events

  return 0;
}
