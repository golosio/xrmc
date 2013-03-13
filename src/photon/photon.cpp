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
//             photon.cpp              //
//            08/02/1013               //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Methods of the class photon
// These are the main method for x-ray photon transport and
// interaction with matter simulation
//

#include <cmath>
#include <string>
#include <iostream>
#include <stdio.h>
#include "xrmc_photon.h"
#include "xrmc_sample.h"
#include "xrmc_composition.h"
#include "xraylib.h"
#include "xrmc_math.h"
#include "xrmc_algo.h"
#include "xrmc_fluo_lines.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

int photon::FluorFlag; // flag for fluorescent emission activation 

//////////////////////////////////////////////////////////////////////
// Evaluates the photon next interaction type and position
// The photon starting coordinates, direction and polarization
// are the vectors x, uk and ui. The energy is E
// The interaction occurs with the device "Sample"
// At the end of the procedure, x, uk, ui and E are updated to
// the values at the interaction position
// iType is the interaction type (fluorescence emission, elastic
// scattering or inelastic scattering), iZ the atomic number.
//////////////////////////////////////////////////////////////////////
int photon::MonteCarloStep(sample *Sample, int *iZ, int *iType)
{
  double step_length, weight;
  phase *ph_compound;
  int step_idx, interaction_type, Z;
  double mu_atom, cs_interaction[3], cs_tot; 

  // evaluates the intersections of the photon trajectory with
  // the phases inside the sample and the absorption coefficient
  // at each step
  Sample->LinearAbsorption(x, uk);
  if (Sample->Path->NSteps==0) { // check if it does not intersect any object
    w = 0;
    return 0; 
  }

  // Check if the flag for weighted step length is activated. In this case,
  // extract the next interaction position using the weighted event approach
  if (Sample->WeightedStepLength == 1) 
    step_length = Sample->Path->WeightedStepLength(&step_idx, &weight);
  else // otherwise use the standard MC approach 
    step_length = Sample->Path->StepLength(&step_idx, &weight);

  w *= weight; // update the event weight
  if (w==0) return 0;

  // move the photon in the direction uk by a distance step_length
  MoveForward(step_length);

  // phase where the next interaction will occur
  ph_compound = &Sample->Comp->Ph[Sample->Path->iPh0[step_idx]];
  // extract the atomic species that the photon will interact with
  ph_compound->AtomType(&Z, &mu_atom);
  // cross sections of the three interaction types with the extracted element
  CSInteractions(Z, cs_interaction, &cs_tot);
  if (cs_tot+mu_atom == cs_tot) { // check for division overflow
    w=0;
    return 0;
  }
  // cs_tot includes all interactions excluding photoelectric absorption
  // without fluorescent emission. Basically, the photon is forced to be
  // scattered or emitted, it cannot "die" in the interaction point
  w *= cs_tot / mu_atom; // update the event weight
  // extract interaction type (elastic/inelastic scattering or fluorescence) 
  interaction_type = InteractionType(cs_interaction, cs_tot);

  // evaluate the energy for fluorescence emission
  if (interaction_type == FLUORESCENCE) SetFluorescenceEnergy(Z);
  *iZ = Z; // atomic species that the photon interact with
  *iType = interaction_type;

  return 0;
}

// move the photon in the direction uk by a distance step_length
int photon::MoveForward(double step_length)
{
  x += uk*step_length;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// cross sections of the three interaction types with the extracted element
//////////////////////////////////////////////////////////////////////
int photon::CSInteractions(int Z, double *cs_interaction, double *cs_tot)
{
  int i;
  
  cs_interaction[FLUORESCENCE] = 0;
  if (FluorFlag==1) { // check if fluorescence emission is activated
    for (i=0; i<NLines[Z]; i++) { // loop on possible fluorescence lines
      // cumulative sum of the line cross sections
      cs_interaction[FLUORESCENCE] += CS_FluorLine(Z, Line[Z][i], E);
    }
  }
  // coherent and incoherent cross sections
  cs_interaction[COHERENT] = CS_Rayl(Z, E);
  cs_interaction[INCOHERENT] = CS_Compt(Z, E);
  // total interaction cross section
  *cs_tot = cs_interaction[FLUORESCENCE] + cs_interaction[COHERENT] +
    cs_interaction[INCOHERENT];

  return 0;
}

//////////////////////////////////////////////////////////////////////
// extract interaction type (elastic/inelastic scattering or fluorescence) 
//////////////////////////////////////////////////////////////////////
int photon::InteractionType(double *cs_interaction, double cs_tot)
{
  double R;

  R = Rnd()*cs_tot; // random number between 0 and total cross section

  // determine in which of the three intervals the random number is comprised
  // and return the index of the interval
  if (R < cs_interaction[0]) return 0;
  if ((cs_interaction[2]==0) || 
      (R < (cs_interaction[0] + cs_interaction[1]))) {
    return 1;
  }
  return 2;
}

//////////////////////////////////////////////////////////////////////
// Update the photon direction. The two angles theta and phi are
// the polar and the azimuthal angle of the new direction, respectively,
// in the local coordinate system associated to the photon ui, uj, uk 
//////////////////////////////////////////////////////////////////////
int photon::ChangeDirection(double theta, double phi)
{
  double x, y, z, ui_uk;
  int i;

  x = sin(theta)*cos(phi); // new direction components 
  y = sin(theta)*sin(phi); // in the local coordinate system
  z = cos(theta);
  uj = uk^ui;
  uk = ui*x + uj*y + uk*z; //new direction vector in the absolute coord system
  uk.Normalize(); // normalize the new direction vector to 1
  for(;;) {
    ui_uk = ui*uk;
    ui -= uk*ui_uk; // component of ui perpendicular to the new direction

    // normalize ui. If it is not zero, its\'s done 
    if (ui.Normalize() == 0)  break;
    for (i=0; i<3; i++) { // otherwise extract a random direction
      ui.Elem[i] = 2.0*Rnd() - 1.0; //  and repeat the loop
    }
  }
  uj=uk^ui;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Extract the line of fluorescent emission
// and update the photon energy accordingly
//////////////////////////////////////////////////////////////////////
int photon::SetFluorescenceEnergy(int Z)
{
  int i_line;
  double sum, R, cs_line_sum[MAXNLINES+1];

  sum = 0;
  cs_line_sum[0] = 0;
  for (i_line=0; i_line<NLines[Z]; i_line++) { // loop on fluorescent lines
    sum += CS_FluorLine(Z, Line[Z][i_line], E);
    cs_line_sum[i_line+1] = sum; // cumulative sum of their cross sections
  }
  R = Rnd()*sum; // random number between 0 and total fluor. cross section
  Locate(R, cs_line_sum, NLines[Z], &i_line); // extract a line

  E = LineEnergy(Z, Line[Z][i_line]); // fluorescent line energy

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Check the interaction type and launch the corresponding method
// to update the photon direction, polarization and energy
//////////////////////////////////////////////////////////////////////
int photon::Scatter(int Z, int interaction_type)
{
  char s[3];

  switch (interaction_type)
    {
    case FLUORESCENCE:
      Fluorescence();
      break;
    case COHERENT : 
      Coherent(Z);
      break;
    case INCOHERENT : 
      Incoherent(Z);
      break;
    default :
      sprintf(s, "%d", interaction_type);
      throw xrmc_exception(string("Wrong Interaction Type Index: ")
			   + s + "\n");
    }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Check the interaction type and launch the corresponding method
// In this case, the photon is forced to go in the direction v_r
// Update the photon direction, polarization and energy
//////////////////////////////////////////////////////////////////////
int photon::Scatter(int Z, int interaction_type, vect3 v_r)
{
  const double th_min=1e-10;
  vect3 rel_r;
  double r, theta, phi, weight, num, denom;
  char s[3];

  // evaluates the components of the new direction v_r
  // in the local coordinate system of the photon ui, uj, uk
  rel_r.Set(ui*v_r, uj*v_r, uk*v_r);
  // evaluates the polar coordinates of the new direction
  // in the local coordinate system of the photon
  rel_r.CartesianToPolar(&r, &theta, &phi); 

  switch (interaction_type)
    {
    case FLUORESCENCE :
      weight = 1 / (4*PI); // isotropic angular distribution
      break;
    case COHERENT :
      // the weight is the ratio between differential and total cross section
      num = DCSP_Rayl(Z, E, theta, phi);
      denom = CS_Rayl(Z, E);
      if (num+denom==num) weight = 0; // check for division overflow
      else weight = num/denom;
      break;
    case INCOHERENT :
      if (theta<th_min) theta=th_min;
      // the weight is the ratio between differential and total cross section
      num = DCSP_Compt(Z, E, theta, phi);
      denom = CS_Compt(Z, E);
      if (num+denom==num) weight = 0; // check for division overflow
      else weight = num/denom;
      E = ComptonEnergy(E, theta); // update the photon energy
      break;
    default :
      sprintf(s, "%d", interaction_type);
      throw xrmc_exception(string("Wrong Interaction Type Index: ")
			   + s + "\n");
    }
  w *= weight; // update the weight

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Fluorescent emission
// The photon is emitted with random direction and random polarization
//////////////////////////////////////////////////////////////////////
int photon::Fluorescence()
{
  double r;

  do {
      uk.Elem[0] = Rnd()-0.5; // random direction
      uk.Elem[1] = Rnd()-0.5; // elements between -0.5
      uk.Elem[2] = Rnd()-0.5; // and +0.5
      r = uk.Mod();
  } while (r>0.5 || r==0); // uniform distribution inside a r=1/2 sphere
  uk = uk/r; // normalize uk to 1

  do {
    do {
      ui.Elem[0] = Rnd()-0.5; // random polarization vector
      ui.Elem[1] = Rnd()-0.5; // elements between -0.5
      ui.Elem[2] = Rnd()-0.5; // and +0.5
      r = ui.Mod();
    } while (r>0.5 || r==0); // uniform distribution inside a r=1/2 sphere
    uj = uk^ui; // vector perpendicular to both ui and uk 
    r = uj.Mod();
  } while (r==0);
  uj = uj/r; // normalize uj to 1
  ui = uj^uk; // build a orthonormal basis

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Coherent scattering
// Weighted event method
// The photon is emitted with random direction weighted with a weight
// proportional to the differential cross section
//////////////////////////////////////////////////////////////////////
int photon::Coherent(int Z)
{
  double theta, phi, cos_theta, weight;

  cos_theta = 2*Rnd()-1;   // random polar angles theta and phi with uniform
  theta = acos(cos_theta); //  distribution on the whole 4*PI solid angle
  phi = 2*PI*Rnd();
  double num = 4.*PI*DCSP_Rayl(Z, E, theta, phi); // differential cross section
  double denom = CS_Rayl(Z, E); // total cross section
  if (num+denom==num) weight = 0; // check for division overflow
  else weight = num/denom;
  w *= weight; // update the event weight
  ChangeDirection(theta, phi);  // update the photon direction

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Coherent scattering
// Weighted event method
// The photon is emitted with random direction weighted with a weight
// proportional to the differential cross section
//////////////////////////////////////////////////////////////////////
int photon::Incoherent(int Z)
{
  const double th_min=1e-10;
  double theta, phi, cos_theta, weight;

  cos_theta = 2*Rnd()-1;   // random polar angles theta and phi with uniform
  theta = acos(cos_theta); //  distribution on the whole 4*PI solid angle
  if (theta<th_min) theta=th_min;
  phi = 2*PI*Rnd();
  double num = 4.*PI*DCSP_Compt(Z, E, theta, phi); // differential cross section
  double denom = CS_Compt(Z, E); // total cross section
  if (num+denom==num) weight = 0; // check for division overflow
  else weight = num/denom;
  w *= weight; // update the event weight
  E = ComptonEnergy(E, cos(theta)); // update the photon energy

  ChangeDirection(theta, phi); // update the photon direction

  return 0;
}

