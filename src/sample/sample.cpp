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
//     sample.cpp                //
//        12/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the classes sample and path
//

#include <iostream>
#include <cmath>
#include "xrmc_sample.h"
#include "xrmc_composition.h"
#include "xrmc_algo.h"
#include "xrmc_photon.h"
#include "xrmc_math.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

// destructor
sample::~sample() { // destructor
  if (PhotonNum!=NULL) delete[] PhotonNum;
  if (Path!=NULL) delete Path;
}

// constructor
sample::sample(string dev_name) {
  Runnable = false;
  NInputDevices = 3;
  PhotonNum = NULL;
  Path = NULL;
  SetDevice(dev_name, "sample");
}

//////////////////////////////////////////////////////////////////////
// method for casting input devices
// to basesource, geom3d and composition types
/////////////////////////////////////////////////////////////////////
int sample::CastInputDevices()
{
  // cast InputDevice[0] to type basesource*
  Source = dynamic_cast<basesource*>(InputDevice[0]);
  if (Source==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0] +
			 " cannot be casted to type basesource\n");

  // cast InputDevice[1] to type geom3d*
  Geom3D = dynamic_cast<geom3d*>(InputDevice[1]);
  if (Geom3D==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[1] +
              " cannot be casted to type geom3d\n");

  // cast InputDevice[2] to type composition*
  Comp = dynamic_cast<composition*>(InputDevice[2]);
  if (Comp==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[2] +
              " cannot be casted to type composition\n");

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for linking an input device
/////////////////////////////////////////////////////////////////////
/*
int sample::LinkInputDevice(string command, xrmc_device *dev_pt)
{
  if (command=="Geom3DName") {
    // cast it to type geom3d*
    Geom3D = dynamic_cast<geom3d*>(dev_pt);
    if (Geom3D==0)
      throw xrmc_exception(string("Device cannot be casted to type "
				  "geom3d\n"));
  }
  else if (command=="SourceName") {
    // cast it to type basesource*
    Source = dynamic_cast<basesource*>(dev_pt);
    if (Source==0)
      throw xrmc_exception(string("Device cannot be casted to type "
				  "basesource\n"));
  }
  else if (command=="CompName") {
    // cast it to type composition*
    Comp = dynamic_cast<composition*>(dev_pt);
    if (Comp==0)
      throw xrmc_exception(string("Device cannot be casted to type "
				  "composition\n"));
  }
  else
    throw xrmc_exception(string("Unrecognized command: ") + command + "\n");

  return 0;
}
*/

//////////////////////////////////////////////////////////////////////
// sample run initialization method
//////////////////////////////////////////////////////////////////////
int sample::RunInit()
{
  int mns;
 
  if (Path!=NULL) delete Path;
  Path = new path; // allocate the path member variable
  mns = 2*(Geom3D->QArr->NQuadr+4)+2; // maximum number of intersections
  Path->MaxNSteps = mns;            // (steplengths)
  Path->t = new double[mns]; // intersections (distances from starting point)
  Path->Step = new double[mns]; // steplengths between adjacent intersections
  Path->iPh0 = new int[mns]; // index of phase before intersection
  Path->iPh1 = new int[mns]; // index of phase after intersection
  Path->Mu = new double[mns]; // absorption coeff. in each step
  Path->Delta = new double[mns]; // delta coeff. in each step
  Path->SumMuS = new double[mns]; // cumulative sum of Mu * steplength
  Path->SumS = new double[mns]; // cumulative sum of steplengths

  return 0;
}

//////////////////////////////////////////////////////////////////////
// sample cleaning after run method
//////////////////////////////////////////////////////////////////////
//int sample::RunFree()
//{
//  // deallocate arrays related to intersection steps
//  delete[] Path->t;
//  delete[] Path->Step;
//  delete[] Path->iPh0;
//  delete[] Path->iPh1;
//  delete[] Path->Mu;
//  delete[] Path->Delta;
//  delete[] Path->SumMuS;
//  delete[] Path->SumS;
//  delete Path; // deallocate the path member variable

//  return 0;
//}

// initialize loop on events
int sample::Begin()
{
  Source->Begin();// initialize source loop on events
  ScattOrderIdx = PhotonIdx = 0; // set scatt. order index and event index to 0

  return 0;
}

// next step of the event loop
int sample::Next()
{
  // check if the loop is finished
  if (ScattOrderIdx >= ScattOrderNum) return 1;

  Source->Next();// next step of the source event loop
  if (Source->End()) { // check if the loop on source events is finished
    Source->Begin(); // if yes, reinitialize source loop on events, ...
    PhotonIdx++;     // increase the event index, ...
    // check if the event index reached the event multiplicity
    if (PhotonIdx >= PhotonNum[ScattOrderIdx]) {
      PhotonIdx = 0;   // if yes, reset the event index
      ScattOrderIdx++; // and start with the next scattering order
    }
  }
   
  return 0;
}

// check if the end of the loop is reached
bool sample::End()
{
  if (ScattOrderIdx >= ScattOrderNum) return true;
  else return false;
}

// event multiplicity
int sample::EventMulti()
{
  int em = 0;
  for (int is=0; is<ScattOrderNum; is++) {
    em += PhotonNum[is];
  }

  return em*Source->EventMulti();
}

// number of scattering orders
int sample::ModeNum()
{
  return ScattOrderNum;
}

//////////////////////////////////////////////////////////////////////
// generate an event with a photon forced to end on the point x1
// and using ModeIdx orders of scattering
//////////////////////////////////////////////////////////////////////
int sample::Out_Photon_x1(photon *Photon, vect3 x1, int *ModeIdx)
{
  *ModeIdx = ScattOrderIdx; // set the scattering order
  return Out_Photon_x1(Photon, x1); // generate the event
}

//////////////////////////////////////////////////////////////////////
// method for evaluating the intersections of a straight trajectory
//  x0 + u*t with the sample
//////////////////////////////////////////////////////////////////////
int sample::Intersect(vect3 x0, vect3 u)
{
  int i;

  // find the intersections of the trajectory with the 3d objects
  Geom3D->Intersect(x0, u, Path->t, Path->iPh0, Path->iPh1, 
		    &Path->NSteps);
  for (i=0; i<Path->NSteps; i++) { // loop on the intersections
    // evaluate the steplengths
    Path->Step[i] = (i!=0) ? (Path->t[i] - Path->t[i-1]) : Path->t[0];
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for evaluating the absorption coefficient in each step
// of the intersections of a straight trajectory
//  x0 + u*t with the sample
//////////////////////////////////////////////////////////////////////
double sample::LinearAbsorption(vect3 x0, vect3 u)
{
  // find the intersections of the trajectory with the 3d objects
  Geom3D->Intersect(x0, u, Path->t, Path->iPh0, Path->iPh1, 
		    &Path->NSteps);
  Path->StepMu(Comp); // evaluate the absorption coefficient in each step

  return Path->MuL; // return the cumulative sum of Mu * steplength
}

//////////////////////////////////////////////////////////////////////
// Analogous to the previous function for t<tmax; tmax must be >= 0
//////////////////////////////////////////////////////////////////////
double sample::LinearAbsorption(vect3 x0, vect3 u, double tmax)
{
  int imax;
  // find the intersections of the trajectory with the 3d objects
  Geom3D->Intersect(x0, u, Path->t, Path->iPh0, Path->iPh1, 
		    &Path->NSteps);

  if (Path->NSteps>0 && tmax>=0) {
    if (tmax<Path->t[0]) {
      Path->NSteps = 1;
      Path->t[0] = tmax;
      Path->iPh1[0] = Path->iPh0[0];
    }
    else {
      Locate(tmax, Path->t, Path->NSteps, &imax); // locate tmax in the inters.
      if (imax < Path->NSteps-1) { // and eventually reduce the number of
	Path->NSteps = imax + 2;   // intersections and limit them to tmax
	Path->t[imax+1] = tmax;
	Path->iPh1[imax+1] = Path->iPh0[imax+1];
      }
    }
  }

  Path->StepMu(Comp); // evaluate the absorption coefficient in each step

  return Path->MuL; // return the cumulative sum of Mu * steplength
}

int sample::LinearMuDelta(vect3 x0, vect3 u)
{
  // find the intersections of the trajectory with the 3d objects
  Geom3D->Intersect(x0, u, Path->t, Path->iPh0, Path->iPh1, 
		    &Path->NSteps);
  Path->StepMuDelta(Comp); // evaluate the delta coefficient in each step

  return 0;
}

// generate a random position in the sample region
vect3 sample::RandomPoint()
{
  vect3 v = Geom3D->X; // center of the sample region
  v.Elem[0] += (-1. + 2.*Rnd())*Geom3D->HW[0]; // uniform probability
  v.Elem[1] += (-1. + 2.*Rnd())*Geom3D->HW[1]; // distribution
  v.Elem[2] += (-1. + 2.*Rnd())*Geom3D->HW[2]; // in a parallelepiped

  return 0;
}

//////////////////////////////////////////////////////////////////////
// evaluate the absorption coefficient at each step of the intersection
//////////////////////////////////////////////////////////////////////
int path::StepMu(composition *comp)
{
  int i;

  MuL = 0;

  for (i=0; i<NSteps; i++) { // loop on intersections
    Step[i] = (i!=0) ? (t[i] - t[i-1]) : t[0]; // steplength of intersection
    Mu[i] = comp->Ph[iPh0[i]].LastMu; // absorption coefficient of the phase
    MuL += Mu[i]*Step[i]; // cumulative sum of Mu * steplength

  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// evaluate absorption and delta coefficients at each step of the intersection
//////////////////////////////////////////////////////////////////////
int path::StepMuDelta(composition *comp)
{
  int i;

  MuL = 0;
  DeltaL = 0;

  for (i=0; i<NSteps; i++) { // loop on intersections
    Step[i] = (i!=0) ? (t[i] - t[i-1]) : t[0]; // steplength of intersection
    Mu[i] = comp->Ph[iPh0[i]].LastMu; // absorption coefficient of the phase
    Delta[i] = comp->Ph[iPh0[i]].LastDelta; // delta coefficient of the phase
    MuL += Mu[i]*Step[i]; // cumulative sum of Mu * steplength
    DeltaL += Delta[i]*Step[i]; // cumulative sum of Delta * steplength
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Extract the next interaction position using the standard MC approach
//////////////////////////////////////////////////////////////////////
double path::StepLength(int *step_idx, double *p_abs)
{
  double logargmin = 1e-10;
  double sum_mu_s, sum_s, step_length, R, mu_s_length, logarg;
  int i, m;

  SumMuS[0] = SumS[0] = sum_mu_s = sum_s = 0;
  for (i=0; i<NSteps; i++) { // loop on intersections
    sum_mu_s += Mu[i]*Step[i];
    sum_s += Step[i];
    SumMuS[i+1] = sum_mu_s; // cumulative sum of mu * steplength
    SumS[i+1] = sum_s; // cumulative sum of steplengths
  }
  *p_abs = 1. - exp(-sum_mu_s); // total absorption probability
  
  do {
    R = Rnd(); // random number 0-1
    logarg = 1 - R*(*p_abs); // use the inverse cumulative distribution
  } while (logarg<logargmin); // check for logarith underflow
  
  mu_s_length = -log(logarg); // sum of mu*L up to the interaction position
  Locate(mu_s_length, SumMuS, NSteps, &m); // locate the step index
  *step_idx = m;
  double num = mu_s_length - SumMuS[m];
  double denom = Mu[m];
  if (num+denom==num) step_length=SumS[m+1]; // check for division overflow
  else {
    double s=num/denom; // path length to the interaction point in the last step
    step_length = SumS[m] + s; // total path length to the interaction point
  }
  
  return step_length; // distance from starting point to interaction point
}

///////////////////////////////////////////////////////////////////////////////
// Extract the next interaction position using the weighted steplength approach
///////////////////////////////////////////////////////////////////////////////
double path::WeightedStepLength(int *step_idx, double *weight)
{
  double sum, step_length, l, MuS, w;
  int i, m;

  m = (int)floor(Rnd()*NSteps);  // choose a step at random
  if (m >= NSteps) m = NSteps - 1;

  sum = 0;
  step_length = 0;
  for (i=0; i<m; i++) { // loop on all steps before the one extracted
    step_length += Step[i]; // sum of steplengths
    sum += Mu[i]*Step[i]; // sum of mu *steplength
  }
  MuS = Mu[m]*Step[m]; // mu * steplength on the step extracted
  w = exp(-sum)*NSteps;
  l = Rnd()*Step[m]; // random position on last step
  step_length += l; // total steplength up to the interaction position
  *weight = w*exp(-Mu[m]*l) * MuS; // weight of the event
  *step_idx = m;

  return step_length;
}


//////////////////////////////////////////////////////////////////////
// generate an event with a photon forced to end on the point x1
//////////////////////////////////////////////////////////////////////
int sample::Out_Photon_x1(photon *Photon, vect3 x1)
{
  int Z, interaction_type;
  vect3 vr;
  double tmax=0;

  if (PhotonNum[ScattOrderIdx]==0) {
    Photon->w = 0;
    return 0;
  }
  if (ScattOrderIdx == 0) { // transmission
    // ask the source to generate photon directed torward the point x1
    Source->Out_Photon_x1(Photon, x1);
    vr = x1 - Photon->x; // end position relative to photon starting position
    tmax = vr.Mod(); // maximum intersection distance
    vr.Normalize(); // normalized direction
    Comp->Mu(Photon->E); // absorption coefficients at photon energy
  }
  else {
    Source->Out_Photon(Photon); // at least one scattering process
    Comp->Mu(Photon->E); // absorption coefficients at photon energy
    // loop on scattering interactions up to the scattering order
    for (int is=1; is<=ScattOrderIdx; is++) {
      // Evaluates the photon next interaction type and position
      Photon->MonteCarloStep(this, &Z, &interaction_type);
      if (Photon->w == 0) break;
      if (is<ScattOrderIdx) { // check that it is not the last interaction
	// update the photon direction, polarization and energy
	Photon->Scatter(Z, interaction_type);
	// update absorption coefficients using the new energy value
	if (interaction_type != COHERENT) Comp->Mu(Photon->E);
      }   
    }
    if (Photon->w != 0) {
      vr = x1 - Photon->x; //end position relative to last interaction position
      tmax = vr.Mod(); // maximum intersection distance
      vr.Normalize(); // normalized direction
      // update the photon direction, polarization and energy
      // the photon is forced to go in the direction v_r
      Photon->Scatter(Z, interaction_type, vr);
      // update absorption coefficients using the new energy value
      if (interaction_type != COHERENT) Comp->Mu(Photon->E);
    }
  }

  // weight the event with the survival probability and divide by multiplicity
  if (Photon->w>0) {
    Photon->w *= exp(-LinearAbsorption(Photon->x, vr, tmax))
      /PhotonNum[ScattOrderIdx];
  }

  return 0;
}

