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
//     xrmc_sample.h             //
//        12/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// sample and path classes definition
//

#ifndef SAMPLEH
#define SAMPLEH

#include <string>
#include "xrmc_geom3d.h"
#include "xrmc_composition.h"
#include "xrmc_device.h"
#include "xrmc_source.h"
#include "phcdevice.h"
#include "randmt.h"

using namespace std;

class photon;

// path class definition, member variables and functions
class path
{
 public:
  int MaxNSteps; // maximum number of intersections
  int NSteps; // number of intersections of a trajectory with the sample
  double *t; // array of the intersections
             // (distances from the starting coordinate)
  double *Step; // steplength between intersections
  int *iPh0; // array of the entrance phase indexes
  int *iPh1; // array of the exit phase indexes
  double *Mu; // absorption coefficient at each step
  double MuL; // sum of mu * steplength
  double *Delta; // delta coefficient
  double DeltaL; // sum of delta * steplength
  double *SumMuS; // cumulative sum of mu * steplength
  double *SumS; // cumulative sum of steplengths
  randmt_t *Rng;

  ~path() { // destructor
    if(t!=NULL) delete[] t;
    if(Step!=NULL) delete[] Step;
    if(iPh0!=NULL) delete[] iPh0;
    if(iPh1!=NULL) delete[] iPh1;
    if(Mu!=NULL) delete[] Mu;
    if(Delta!=NULL) delete[] Delta;
    if(SumMuS!=NULL) delete[] SumMuS;
    if(SumS!=NULL) delete[] SumS;
  }
  path() { // constructor
    t = Step = Mu = Delta = SumMuS = SumS = NULL;
    iPh0 = iPh1 = NULL;
    NSteps = 0;
    Rng = NULL;
  }
  // evaluate the absorption coefficient at each step of the intersection
  int StepMu(composition *comp);
  // evaluate absorption and delta coefficients at each step of the intersection
  int StepMuDelta(composition *comp);
  // Extract the next interaction position using the standard MC approach
  double StepLength(int *step_idx, double *weight);
  //Extract the next interaction position using the weighted steplength approach
  double WeightedStepLength(int *step_idx, double *weight);
  path *Clone();
};

// sample class definition, member variables and functions
class sample : public basesource
{
 public:
  path *Path; //object storing all info about intersections
  composition *Comp; // input composition device
  int ScattOrderNum; // number of scattering orders
  int *PhotonNum; // event multiplicity for each scattering order 
  int ScattOrderIdx; // scattering order index
  int PhotonIdx; // event index
  int WeightedStepLength; // flag for weighted steplength extraction method
  basesource *Source; // input source device

  virtual ~sample(); // destructor
  sample(std::string dev_name); // constructor
  // method for casting input devices
  virtual int CastInputDevices();
  virtual int Load(istream &fs); // load sample parameters from file
  virtual int SetDefault(); // set default values for sample parameters
  virtual int Begin(); // begin event loop
  virtual int Next();  // next step in event loop
  virtual bool End();  // check for end of event loop
  virtual long long EventMulti(); // event multiplicity
  virtual int ModeNum(); // return the number of modes (scattering orders)
  virtual int Intersect(vect3 x0, vect3 u); // evaluate intersection of a trajectory
                                    // with sample objects
  //evaluate the absorption coefficient at each step of the intersections
  double LinearAbsorption(vect3 x0, vect3 u);
  // Analogous to the previous function for t<tmax
  double LinearAbsorption(vect3 x0, vect3 u, double tmax);
  //evaluate absorption and delta coefficient at each step of the intersections
  int LinearMuDelta(vect3 x0, vect3 u);

// weight the event with the survival probability
  int PhotonSurvivalWeight(photon *Photon, double tmax);
  // simulates the photon history up to the last interaction point
  int PhotonHistory(photon *Photon, int &Z, int &interaction_type);
  // generate an event with a photon forced to end on the point x1
  virtual int Out_Photon_x1(photon *Photon, vect3 x1);
  virtual int Out_Photon_x1(photon *Photon, vect3 x1, int *ModeIdx);
  // simulates an event up to the last interaction point
  virtual int Out_Photon(photon *Photon);
  virtual int Out_Photon(photon *Photon, int *ModeIdx);
  virtual int Out_Phase_Photon_x1(photon *Photon, vect3 x1);
  virtual double GetPhC_E0();
  virtual int PhCOn();
  virtual int PhCOff();

  virtual int RunInit(); // sample run initialization method
 // set the random number generator structure
  virtual int SetRng(randmt_t *rng);
  virtual basesource *Clone(string dev_name);
 private:
  geom3d *Geom3D; // input geom3d device

  vect3 RandomPoint(); // generate a random position in the sample region

};

sample *LoadGeom3DSample(char *QuadricFileName, char *GeomFileName,
                         char *CompFileName);

#endif
