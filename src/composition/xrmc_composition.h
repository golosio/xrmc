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
/////////////////////////////////////////
//          xrmc_composition.h         //
//        31/01/2013                   //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Definition of the classes composition and phase
//

#ifndef COMPOSITIONH
#define COMPOSITIONH
#include <string>
#include "xrmc_device.h"

// phase class definition, member variables and functions
class phase
{
 public:
  int NElem; // number of elements in the phase
  double Rho; // mass density of the phase in g/cm3
  int *Z; // atomic number array
  double *W; // weight fraction array
  double *MuAtom; // atomic interactions total cross section array
  double LastMu; // linear absorption coefficient
  double LastDelta; // delta coefficient 
  static const double KD; // constant used for computing Delta

  ~phase() { // destructor
    if (W!=NULL) delete[] W;
    if (Z!=NULL) delete[] Z;
    if (MuAtom!=NULL) delete[] MuAtom;
  }
  phase() { // constructor
    W = NULL;
    Z = NULL;
    MuAtom = NULL;
    NElem = 0;
    Rho = LastMu = LastDelta = 0;
  }
  int Mu(double E);// Evaluates the absorption coefficient at energy E
  int Delta(double E); // Evaluates the delta coefficient at energy E
  int AtomType(int *Z, double *mu_atom);// extract the atomic species with
                                        // which the interaction will occur


  phase& operator= (const phase &Phase);
};

typedef map<string, int> phase_map;
typedef pair<string, int> phase_map_pair;
typedef pair<phase_map::iterator, bool> phase_map_insert_pair;

// composition class definition, member variables and functions
class composition : public xrmc_device
{
 public:
  int NPhases; // number of phases (materials)
  int MaxNPhases; // maximum number of phases (materials)
  phase *Ph;   // phase array
  phase_map PhaseMap; // map of phases with their names

  ~composition(); // destructor
  composition(std::string dev_name);  // constructor
  int Load(istream &fs); // Loads sample phases composition and density
  int SetDefault(); // Set default values for composition parameters
  // insert name and index of the phase in the phase map
  int MapPhase(istream &fs);

  int Mu(double E); // Evaluates the absorption coefficient of each phase
  int Delta(double E); // Evaluates the delta coefficient of each phase
  composition *Clone(string dev_name);
};

#endif

