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
//     composition.cpp           //
//        31/10/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the classes composition and phase
//
#include "xrmc_composition.h"
#include "xraylib.h"
#include "xrmc_algo.h"
#include <iostream>
#include <cstring>

using namespace std;
using namespace xrmc_algo;

const double phase::KD=4.15179082788e-4; // constant for computing Delta

 // destructor
composition::~composition() {
  for (unsigned int i=0; i<Ph.size(); i++) {
    if (Ph[i].W!=NULL) free(Ph[i].W);
    if (Ph[i].Z!=NULL) free(Ph[i].Z);
    if (Ph[i].MuAtom!=NULL) free(Ph[i].MuAtom);
  }
}

// constructor
composition::composition(string dev_name) {
  Runnable = false;
  NInputDevices=0;
  SetDevice(dev_name, std::string("composition"));
}

// Evaluates the absorption coefficient Mu of each phase at energy E
int composition::Mu(double E)
{
  for (unsigned int i=0; i<Ph.size(); i++) {
    Ph[i].Mu(E);
  }

  return 0;
}

// Evaluates the coefficient delta of each phase at energy E
// 1-delta is the real part of the refractive index
int composition::Delta(double E)
{
  for (unsigned int i=0; i<Ph.size(); i++) {
    Ph[i].Delta(E);
  }

  return 0;
}

// Evaluates the absorption coefficient Mu of the phase at energy E
int phase::Mu(double E)
{
  double mu = 0;
  for (int i=0; i<NElem; i++) { // loop on all elements of the phase
    MuAtom[i] = CS_Total(Z[i], E); // total cross section  for atomic num. Z
    mu += W[i]*MuAtom[i]; // sum of C.S. weighted by the weight fractions
  }
  LastMu = mu*Rho; // linear absorption coefficient (Rho is the mass density)

  return 0;
}

// Evaluates the delta coefficient of the phase at energy E
int phase::Delta(double E)
{
  double delta = 0;
  for (int i=0; i<NElem; i++) {
    // for details about the calculation see references in the documentation
    double num = W[i]*KD*(Z[i]+Fi(Z[i],E));
    double denom = AtomicWeight(Z[i])*E*E;
    // check that the denominator is not much smaller than numerator
    if (num+denom == num) delta = 0;
    else delta += num/denom;
  }
  LastDelta = delta*Rho; // delta coefficient (Rho is the mass density)

  return 0;
}

// extract the atomic species with which the interaction will occur
// Zelem: atomic number; mu_atom: total cross section for this element
int phase::AtomType(int *Zelem, double *mu_atom)
{
  // check that Rho is not 0 and that it is not much smaller than numerator
  if (Rho==0 || LastMu+Rho==LastMu) return -1;

  double mu_compound = LastMu / Rho; // absorption coefficient of the phase
  double R = Rnd_r(Rng)*mu_compound; // random num. between 0 and mu_compound

  double sum = 0;
  int i = 0;
  // extract the index of the element with which the interaction will occur
  while (sum <= R) { // stop when the cumulative distribution is >= R
    sum += W[i]*MuAtom[i];
    i++;
  }
  if (i > NElem) i = NElem;
  i--;
  *Zelem = Z[i]; // atomic number of the element
  *mu_atom = MuAtom[i]; // total cross section for the element

  return 0;
}


composition *composition::Clone(string dev_name) {
	//cout << "Entering composition::Clone\n";
	composition *clone = new composition(dev_name);

	cout << "Clone " << dev_name << endl;
	for (unsigned int i = 0 ; i < Ph.size(); i++) {
	  clone->Ph.push_back(Ph[i]);
	  cout << "ph " << i << endl; 
	  cout << "N " << clone->Ph[i].NElem << endl;
	  for (int j=0; j<clone->Ph[i].NElem; j++) {
	    cout << "j " << j << " " << clone->Ph[i].Z[j] << " "
		 << clone->Ph[i].W[j] << endl;
	  }

	}
	//PhaseMap
	//typedef map<string, int> phase_map;
	//typedef pair<string, int> phase_map_pair;
	//typedef pair<phase_map::iterator, bool> phase_map_insert_pair;
	phase_map::iterator it;

	for (it = PhaseMap.begin(); it != PhaseMap.end() ; it++) {
		clone->PhaseMap.insert(phase_map_pair(it->first, it->second));
		cout << "pm " << it->first << " " << it->second << endl;
	}



	return clone;
}


phase& phase::operator= (const phase &Phase) {
	//cout << "Entering phase assignment operator\n";

	if (this == &Phase)
		return *this;

	NElem = Phase.NElem;
	Rho = Phase.Rho;
	Z = (int *) malloc(sizeof(int)*NElem); 
	memcpy(Z, Phase.Z, sizeof(int)*NElem);
	W = (double *) malloc(sizeof(double)*NElem); 
	memcpy(W, Phase.W, sizeof(double)*NElem);
	MuAtom = (double *) malloc(sizeof(double)*NElem);
	memcpy(MuAtom, Phase.MuAtom, sizeof(double)*NElem);
	LastMu = Phase.LastMu;
	LastDelta = Phase.LastDelta;
	
	return *this;
}

int composition::SetRng(randmt_t *rng)
{
  Rng = rng;

  for (unsigned int i=0; i<Ph.size(); i++) {
    Ph[i].Rng = rng;
  }

  return 0;
}
int composition::ReduceMap(vector<string> used_phases) {

  vector <phase> Ph_new;
  phase_map PhaseMap_new;
  phase_map::iterator it;

  cout << "Begin inside reduce Ph " << endl;
  for (unsigned int i = 0 ; i < Ph.size(); i++) {
    cout << "ph " << i << endl; 
    cout << "N " << Ph[i].NElem << endl;
    for (int j=0; j<Ph[i].NElem; j++) {
      cout << "j " << j << " " << Ph[i].Z[j] << " "
	   << Ph[i].W[j] << endl;
    }  
  }

  for (int i=0 ; i<used_phases.size() ; i++) {
    it = PhaseMap.find(used_phases.at(i));
    if (it==PhaseMap.end()) {
      throw xrmc_exception(string("Phase ") + used_phases.at(i)
                           + " not found in phase map\n");
    }

    cout << "ok ph " << it->second << endl; 
    cout << "N " << Ph[it->second].NElem << endl;
    for (int j=0; j<Ph[it->second].NElem; j++) {
      cout << "j " << j << " " << Ph[it->second].Z[j] << " "
	   << Ph[it->second].W[j] << endl;
    }  

    
    Ph_new.push_back(Ph[it->second]);
    cout << "ok ph_new " << i << endl; 
    cout << used_phases.at(i) << endl;
    cout << "N " << Ph_new[i].NElem << endl;
    for (int j=0; j<Ph_new[i].NElem; j++) {
      cout << "j " << j << " " << Ph_new[i].Z[j] << " "
	   << Ph_new[i].W[j] << endl;
    }  

    PhaseMap_new.insert(phase_map_pair(used_phases.at(i), i));
  }

  cout << "Reduce " << endl;
  for (unsigned int i = 0 ; i < Ph_new.size(); i++) {
    cout << "ph " << i << endl; 
    cout << used_phases.at(i) << endl;
    cout << "N " << Ph_new[i].NElem << endl;
    for (int j=0; j<Ph_new[i].NElem; j++) {
      cout << "j " << j << " " << Ph_new[i].Z[j] << " "
	   << Ph_new[i].W[j] << endl;
    }  
  }
  cout << "Ph " << endl;
  for (unsigned int i = 0 ; i < Ph.size(); i++) {
    cout << "ph " << i << endl; 
    cout << "N " << Ph[i].NElem << endl;
    for (int j=0; j<Ph[i].NElem; j++) {
      cout << "j " << j << " " << Ph[i].Z[j] << " "
	   << Ph[i].W[j] << endl;
    }  
  }

  //delete []Ph;
  Ph = Ph_new;
  PhaseMap.clear();
	
  for (it = PhaseMap_new.begin(); it != PhaseMap_new.end() ; it++) {
    PhaseMap.insert(phase_map_pair(it->first, it->second));
  }

  return 0;
}
