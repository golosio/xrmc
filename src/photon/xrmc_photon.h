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
//             photon.h                //
//           08/02/2013                //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Definition of the class photon
// This is the main class for x-ray photon transport and
// interaction with matter simulation
//

#ifndef XRMCPHOTONH
#define XRMCPHOTONH
#include "xrmc_math.h"
#include "randmt.h"

#define FLUORESCENCE 0
#define COHERENT 1
#define INCOHERENT 2

#define PHOTO_NO_FLUOR 1

class sample;

//////////////////////////////////////////////////////////////////////
// Definition of the class photon with member variables and methods
//////////////////////////////////////////////////////////////////////
class photon
{
 public:
  static int FluorFlag; // flag for fluorescent emission activation
  double w; // event weight
  double E; // photon energy
  vect3 x; // photon position
  vect3 ui, uj, uk; // photon local coord. system
  // uk: direction; ui: polarization vector
  randmt_t *Rng;
  sample *mySample;

  photon() {Rng=NULL;} // constructor
// move the photon in the direction uk by a distance step_length
  int MoveForward(double step_length); 

  int ChangeDirection(double theta, double phi); // Update the photon direction

  // Evaluates the photon next interaction type and position
  int MonteCarloStep(sample *Sample, int *iZ, int *iType);

  // Evaluates energy deposition and mu at photon end point
  int EnergyDeposition(sample *Sample, int *iZ, int *iType, double *mu_x1,
		       double *Edep);

  // Evaluates absorption coefficient at photon end point
  int Mu_x1(sample *Sample, int *iZ, int *iType, double *mu_x1);

  // cross sections of the three interaction types with the extracted element
  int CSInteractions(int Z, double *mu_interaction, double *cs_tot);

  // cross sections of the interaction types with energy deposition
  int CSInteractionsEdep(int Z, double *cs_interaction,
			 double *cs_tot_Edep);

  // extract interaction type (elastic/inelastic scattering or fluorescence) 
  int InteractionType(double *cs_interaction, double cs_tot);

  // Extract the line of fluorescent emission
  int SetFluorescenceEnergy(int Z);

  // Check the interaction type and launch the corresponding method
  int Scatter(int Z, int interaction_type);

  // As before in case of forced direction
  int Scatter(int Z, int interaction_type, vect3 v_r);

  int Fluorescence(); // Fluorescent emission
  int Coherent(int Z); // Coherent (elastic) scattering
  int Incoherent(int Z); // Incoherent (inelastic) scattering
  int ComptonEnergyDoppler(int Z, double theta);
};
#endif

