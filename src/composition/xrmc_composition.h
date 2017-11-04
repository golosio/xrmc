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
//        31/10/2017                   //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Definition of the classes composition and phase
//

#ifndef COMPOSITIONH
#define COMPOSITIONH
#include <string>
#include <cstdlib>
#include "xrmc_device.h"

// phase class definition, member variables and functions
class phase
{
 public:
  //int NElem; // number of elements in the phase
  double Rho; // mass density of the phase in g/cm3
  std::vector<int> Z; // atomic number array
  std::vector<double> W; // weight fraction array
  std::vector<double> MuAtom; // atomic interactions total cross section array
  double LastMu; // linear absorption coefficient
  double LastDelta; // delta coefficient 
  static const double KD; // constant used for computing Delta
  randmt_t *Rng;

  phase() { // constructor
    Rho = LastMu = LastDelta = 0;
  }

  inline int NElem() { // number of elements in the phase
    return Z.size();
  }
  int Mu(double E);// Evaluates the absorption coefficient at energy E
  int Delta(double E); // Evaluates the delta coefficient at energy E
  int AtomType(int *Z, double *mu_atom);// extract the atomic species with
                                        // which the interaction will occur
};


// structure defining a material (compound or mixture) 
class material
{
 public:
  vector<int> iPh;
  vector<double> Fact;

  inline int NPh() { // number of compounds in the material
    return iPh.size();
  }

};

typedef map<string, int> phase_map;
typedef pair<string, int> phase_map_pair;
typedef pair<phase_map::iterator, bool> phase_map_insert_pair;

// composition class definition, member variables and functions
class composition : public xrmc_device
{
 public:
  //int MaxNPhases; // maximum number of phases (materials)
  vector<phase> Ph;   // phase array
  phase_map PhaseMap; // map of phases with their names

  vector<material> Mater; // material array
  phase_map MaterMap; // map of materials with their names

  virtual ~composition(); // destructor
  composition(std::string dev_name);  // constructor
  virtual int Load(istream &fs); // Loads sample phases composition and density
  virtual int SetDefault(); // Set default values for composition parameters
  // insert name and index of the phase in the phase map
  std::string MapPhase(istream &fs, int i_ph);
  int MapMater(std::string mater_name, int i_mat);
  std::string MapMater(istream &fs, int i_mat);

  int Mu(double E); // Evaluates the absorption coefficient of each phase
  int Delta(double E); // Evaluates the delta coefficient of each phase
  virtual composition *Clone(string dev_name);
  virtual int SetRng(randmt_t *rng);
  int ReduceMap(vector<string> used_phases);
  int ReduceMaterMap(vector<string> used_mater);

};

#endif

