/*
Copyright (C) 2014 Bruno Golosio and Tom Schoonjans

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
//          loadradionuclide.cpp       //
//            25/06/2014               //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
// Load radionuclide parameters 
//

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_radionuclide.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include <xraylib.h>

using namespace std;
using namespace gettoken;

int radionuclide::Load(istream &fs)
{
  int i;
  string comm="";
  string source="", unit;

  cout << "Radionuclide Parameters\n";

  // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    cout << "x" << comm << "x" << endl;
    // parse the command and decide what to do
    if(comm=="Unit") { 
	GetToken(fs, unit);
	if (unit=="mCi") {
		cout << "mCi unit set" << endl;
		Unit = 1;
	}	
	else if (unit=="Ci") {
		cout << "Ci unit set" << endl;
		Unit = 2;
	}
	else if (unit=="GBq") {
		cout << "GBq unit set" << endl;
		Unit = 3;
	}
	else if (unit=="Bq") {
		cout << "Bq unit set" << endl;
		Unit = 4;
	}
        else throw xrmc_exception("Invalid Unit detected.\n");
    }
    else if (comm == "Activity") {
    	GetDoubleToken(fs, &Activity);
	if (Activity <= 0.0)
		throw xrmc_exception("Activity must be a strictly positive real number\n");
	else cout << "Activity" << "\t" << Activity << endl;
    }
    else if (comm == "RadioNuclideSource") {
    	GetToken(fs, RadioNuclide);
	rnd = GetRadioNuclideDataByName(RadioNuclide.c_str(), NULL);
	if (rnd != NULL) {
		cout << "RadioNuclide set to " << RadioNuclide << endl;
	}
	else {
		throw xrmc_exception("Unknown radionuclide detected\n");
	}
    }
    else if(comm=="LoopFlag") { // flag for modality of energy extraction
      GetIntToken(fs, &LoopFlag);
      if (LoopFlag==0)
	cout << "Extract random energies on the whole spectrum\n";
      else if (LoopFlag==1)
	cout << "Loop on all lines and sampling points\n";
      else throw xrmc_exception("Wrong loop flag.\n");
    }
    else if(comm=="LinePhotonNum") {
      // Number of samples for each line in the spectrum 
      GetIntToken(fs, &LinePhotonNum);
      cout << "Number of samples for each line in the spectrum: " 
	   << LinePhotonNum << "\n";
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empty string\n";
    }
    else {
      throw xrmc_exception("syntax error in radionuclide input file"); 
    }
  }
  return 0;
}

int radionuclide::SetDefault() 
{
	cout << "Setting radionuclide defaults\n";
	Unit = 1;
	Activity = 100.0;
	RadioNuclide = "55Fe";


	spectrum::SetDefault();
	return 0;
}

