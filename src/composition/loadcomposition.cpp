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
//     loadcomposition.cpp       //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load sample phases composition and density
//

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_composition.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"

using namespace std;
using namespace gettoken;

int composition::Load(istream &fs)
{
  /* Loads sample phases composition and density */
  int n_elem;
  int z_elem;
  double w, rho;
  string comm;

  cout << "Composition file\n";

 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    if (comm=="End") break;
    else if(comm=="MaxNPhases") { // set the maximum number of phases    
      GetIntToken(fs, &MaxNPhases);
      cout << "Maximum number of phases: " << MaxNPhases << "\n";
      delete[] Ph;
      Ph = new phase[MaxNPhases+1]; // allocate phase array
      NPhases = 1;
      // initializes phase 0 as vacuum
      Ph[0].NElem = 0;
      Ph[0].Rho = 0;
    }
    else if(comm=="Phase") { // read parameters for a new phase
      MapPhase(fs);
      GetToken(fs, comm);
      if (comm != "NElem")
	throw xrmc_exception("NElem variable initialization not found"); 
      GetIntToken(fs, &n_elem); // read number of elements in the phase
      Ph[NPhases].NElem = n_elem;
      cout << "Num. of elements: " << Ph[NPhases].NElem << endl;
      cout << "\tZ\tweight fract.\n";
      Ph[NPhases].W = new double[n_elem];  // initializes weight fraction,
      Ph[NPhases].Z = new int[n_elem];     // atomic numbers and absorption
      Ph[NPhases].MuAtom = new double[n_elem]; // coefficient arrays
      
      for(int elem_idx=0; elem_idx<n_elem; elem_idx++) { //read element Z and w
	GetIntToken(fs, &z_elem);
	GetDoubleToken(fs, &w);
	Ph[NPhases].Z[elem_idx] = z_elem;
	Ph[NPhases].W[elem_idx] = w/100; // convert from percent to fraction
	cout << "\t" << Ph[NPhases].Z[elem_idx] << "\t"
	     << Ph[NPhases].W[elem_idx] << endl;
      }
      GetToken(fs, comm);
      if (comm != "Rho")
	throw xrmc_exception("Rho variable initialization not found"); 
      GetDoubleToken(fs, &rho); // read mass density of the phase in g/c3
      Ph[NPhases].Rho = rho;
      cout << "\tRho: " << Ph[NPhases].Rho << endl;    
      NPhases++;
    }
    else {
      throw xrmc_exception("Syntax error in composition file\n");
      // unrecognized command
    }
    if (NPhases>MaxNPhases) {
      char i2ch[MAXSTRLEN];
      sprintf(i2ch, "%d", MaxNPhases);
      string s_err="Number of phases greater than maximum: ";
      s_err = s_err + i2ch + "\nUse the command MaxNPhases to "
	"increase the maximum number of phases";
      throw xrmc_exception(s_err);
    }
  }

  return 0;
}

// insert name and pointer to the phase in the phase map
int composition::MapPhase(istream &fs)
{
  string ph_name;
  phase_map_insert_pair insert_pair;

  if (!GetToken(fs, ph_name))
    throw xrmc_exception(string("Syntax error reading phase name\n"));
  // insert name and pointer to the phase in the phase map
  insert_pair = PhaseMap.insert(phase_map_pair(ph_name, NPhases));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Phase ") + ph_name + 
			 " already inserted in phase map\n");
  cout << "Phase: " << ph_name << endl; 

  return 0;
}

//////////////////////////////////////////////////////////////////////
// set default values for composition parameters
int composition::SetDefault()
{
  phase_map_insert_pair insert_pair;

  MaxNPhases = 10000; // maximum num. of phases
  Ph = new phase[MaxNPhases+1];
  NPhases = 1;
  // initializes phase 0 as vacuum
  Ph[0].NElem = 0;
  Ph[0].Rho = 0;
  string ph_name="Vacuum";
  insert_pair = PhaseMap.insert(phase_map_pair(ph_name, 0));

  return 0;
}
