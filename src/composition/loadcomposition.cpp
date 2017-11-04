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
//        31/10/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load sample phases composition and density
//

#include <sstream>
#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_composition.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include <xraylib.h>
#include <cstring>
#include <algorithm>

using namespace std;
using namespace gettoken;

int composition::Load(istream &fs)
{
  /* Loads sample phases composition and density */
  int n_comp;
  double w, rho, fact;
  string comm, comp;
  struct compoundData *cd;

  cout << "Composition file\n";

  // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    if (comm=="End") break;
    else if(comm=="Phase") { // read parameters for a new phase
      int i_phase = Ph.size();
      Ph.push_back(phase());
      string ph_name = MapPhase(fs, i_phase);
      GetToken(fs, comm);
      if (comm != "NElem")
	throw xrmc_exception("NElem variable initialization not found"); 
      GetIntToken(fs, &n_comp); // read number of elements in the phase
      
      int Zelem;
      for(int i_comp=0; i_comp<n_comp; i_comp++) { //read element Z and w
	GetToken(fs, comp);
	GetDoubleToken(fs, &w);
	istringstream buffer(comp);
	if ((buffer >> Zelem)) {
	  //classic mode: integer found
	  vector<int>::iterator it=find(Ph[i_phase].Z.begin(),
					Ph[i_phase].Z.end(), Zelem);
	  int i_elem = distance(Ph[i_phase].Z.begin(), it);
	  if (i_elem<Ph[i_phase].NElem()) {
	    Ph[i_phase].Z[i_elem] = Zelem;
	    Ph[i_phase].W[i_elem] += w/100.0;
	  }
	  else {
	    Ph[i_phase].Z.push_back(Zelem);
	    Ph[i_phase].W.push_back(w/100.0);
	  }
	}
	else {
	  //modern mode: chemical formula
	  if ((cd = CompoundParser(comp.c_str())) == NULL){
	    throw xrmc_exception("Cannot parse compound\n");
	  }
	  int Zcd;
	  for (Zcd=0 ; Zcd<cd->nElements; Zcd++) {
	    vector<int>::iterator it=find(Ph[i_phase].Z.begin(),
					  Ph[i_phase].Z.end(), Zcd);
	    int i_elem = distance(Ph[i_phase].Z.begin(), it);
	    if (i_elem<Ph[i_phase].NElem()) {
	      Ph[i_phase].W[i_elem] += w*cd->massFractions[Zcd]/100.0;
	    }
	    else {
	      Ph[i_phase].Z.push_back(cd->Elements[Zcd]);
	      Ph[i_phase].W.push_back(w*cd->massFractions[Zcd]/100.0);
	    }
	  }
	  FreeCompoundData(cd);
	}	
      }
      Ph[i_phase].MuAtom.assign (Ph[i_phase].NElem(),0);
      cout << "Num. of elements: " << Ph[i_phase].NElem() << endl;
      cout << "\tZ\tweight fract.\n";
      for (int i_elem=0 ; i_elem<Ph[i_phase].NElem(); i_elem++) {
	cout << "\t" << Ph[i_phase].Z[i_elem] << "\t"
	     << Ph[i_phase].W[i_elem] << endl;
      }
      GetToken(fs, comm);
      if (comm != "Rho")
	throw xrmc_exception("Rho variable initialization not found"); 
      GetDoubleToken(fs, &rho); // read mass density of the phase in g/c3
      Ph[i_phase].Rho = rho;
      cout << "\tRho: " << Ph[i_phase].Rho << endl;

       // create a material with only one phase
       int i_mat = Mater.size();
       Mater.push_back(material());
       Mater[i_mat].iPh.push_back(i_phase);
       Mater[i_mat].Fact.push_back(1);
       MapMater(ph_name, i_mat);
    }
    else if(comm=="Mixture") { // read parameters for a new mixture
      int i_mat = Mater.size();
      Mater.push_back(material());
      MapMater(fs, i_mat);
      GetToken(fs, comm);
      if (comm != "NComp")
	throw xrmc_exception("NComp variable initialization not found"); 
      GetIntToken(fs, &n_comp); // read number of compounds in the mixture
       for(int i_comp=0; i_comp<n_comp; i_comp++) { //read compounds
	string comp_name;
	GetToken(fs, comp_name);
	phase_map::iterator it = PhaseMap.find(comp_name);
	if (it==PhaseMap.end()) {
	  throw xrmc_exception(string("Compound ") + comp_name + 
			       " not found in composition map\n");
	}
	// get phase index from the phase map
	int i_ph = (*it).second;
	GetDoubleToken(fs, &w);
	Mater[i_mat].iPh.push_back(i_ph);
	Mater[i_mat].Fact.push_back(w);
      }
      cout << "Num. of compounds: " << Mater[i_mat].NPh() << endl;
      cout << "\tcomp.\tweight fract.\n";
      for (int i_comp=0 ; i_comp<Mater[i_mat].NPh(); i_comp++) {
	int i_ph = Mater[i_mat].iPh[i_comp];
	string comp_name="";
	for (phase_map::iterator it=PhaseMap.begin();
	     it!=PhaseMap.end(); ++it) {
	  if (it->second==i_ph) {
	    comp_name = it->first;
	    break;
	  }
	}
	if (comp_name=="") {
	  throw xrmc_exception("Compound index not found in composition map\n");
	}
	cout << "\t" << comp_name << "\t"
	     << Mater[i_mat].Fact[i_comp] << endl;
      }
      GetToken(fs, comm);
      if (comm != "Rho")
	throw xrmc_exception("Rho variable initialization not found"); 
      GetDoubleToken(fs, &rho); // read mass density of the mixture in g/c3
      cout << "\tRho: " << rho << endl;

      for (int i_comp=0 ; i_comp<Mater[i_mat].NPh(); i_comp++) {
	int i_ph = Mater[i_mat].iPh[i_comp];
	double rho_ph = Ph[i_ph].Rho;
	Mater[i_mat].Fact[i_comp] *= rho/rho_ph;
      }
    }
    else {
      throw xrmc_exception("Syntax error in composition file\n");
      // unrecognized command
    }
  }

  return 0;
}

  // insert name and pointer to the phase in the phase map
string composition::MapPhase(istream &fs, int i_ph)
{
  string ph_name;
  phase_map_insert_pair insert_pair;
  
  if (!GetToken(fs, ph_name))
    throw xrmc_exception(string("Syntax error reading phase name\n"));
  // insert name and pointer to the phase in the phase map
  insert_pair = PhaseMap.insert(phase_map_pair(ph_name, i_ph));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Phase ") + ph_name + 
			 " already inserted in phase map\n");
  cout << "Phase: " << ph_name << endl; 

  return ph_name;
}

// insert name and pointer to the material in the material map
int composition::MapMater(string mater_name, int i_mat)
{
  phase_map_insert_pair insert_pair;

  // insert name and pointer to the material in the material map
  insert_pair = MaterMap.insert(phase_map_pair(mater_name, i_mat));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Material ") + mater_name + 
			 " already inserted in material map\n");
  //cout << "Material: " << mater_name << endl; 

  return 0;
}

// same function as above but reading material name from file
string composition::MapMater(istream &fs, int i_mat)
{
  string mat_name;
  phase_map_insert_pair insert_pair;
  
  if (!GetToken(fs, mat_name))
    throw xrmc_exception(string("Syntax error reading material name\n"));
  // insert name and pointer to the material in the material map
  insert_pair = MaterMap.insert(phase_map_pair(mat_name, i_mat));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Material ") + mat_name + 
			 " already inserted in material map\n");
  cout << "Material: " << mat_name << endl; 

  return mat_name;
}

//////////////////////////////////////////////////////////////////////
// set default values for composition parameters
int composition::SetDefault()
{
  phase_map_insert_pair insert_pair;
  
  Ph.clear();
  Mater.clear();

  // initializes phase 0 as vacuum
  Ph.push_back(phase());
  Ph[0].Rho = 0;
  string ph_name="Vacuum";
  insert_pair = PhaseMap.insert(phase_map_pair(ph_name, 0));

   // initializes material 0 as vacuum
  Mater.push_back(material());
  Mater[0].iPh.push_back(0);
  Mater[0].Fact.push_back(0.0);
  MapMater(ph_name, 0);

  //add NIST compound database
  struct compoundDataNIST *cdn;
  char **list = GetCompoundDataNISTList(NULL);

  for (int i = 0 ; list[i] != NULL ; i++) {
	cdn = GetCompoundDataNISTByIndex(i);
	//cout << cdn->name << endl;
  	xrlFree(list[i]);
	int i_ph = Ph.size();
	
	Ph.push_back(phase());
	Ph[i_ph].Z = vector<int>(cdn->Elements,
				    cdn->Elements+cdn->nElements);
	Ph[i_ph].W = vector<double>(cdn->massFractions,
				    cdn->massFractions+cdn->nElements);
	Ph[i_ph].Rho = cdn->density;
	Ph[i_ph].MuAtom.assign (Ph[i_ph].NElem(),0);

  	insert_pair = PhaseMap.insert(phase_map_pair(cdn->name, i_ph));

	// create a material with only one phase
	int i_mat = Mater.size();
	Mater.push_back(material());
	Mater[i_mat].iPh.push_back(i_ph);
	Mater[i_mat].Fact.push_back(1);
	MapMater(cdn->name, i_mat);

	FreeCompoundDataNIST(cdn);	
  }
  xrlFree(list);

  return 0;
}
