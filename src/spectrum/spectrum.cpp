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
//          spectrum.cpp               //
//           14/02/2013                //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Methods of the class spectrum
//

#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include "xrmc_spectrum.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"

using namespace std;
using namespace xrmc_algo;

 // destructor
spectrum::~spectrum() {
  if (ContinuousEne!=NULL) delete[] ContinuousEne;
  if (ContSIntensity[0]!=NULL) delete[] ContSIntensity[0];
  if (ContSIntensity[1]!=NULL) delete[] ContSIntensity[1];
  if (IntervalIntensity[0]!=NULL) delete[] IntervalIntensity[0];
  if (IntervalIntensity[1]!=NULL) delete[] IntervalIntensity[1];
  if (IntervalWeight[0]!=NULL) delete[] IntervalWeight[0];
  if (IntervalWeight[1]!=NULL) delete[] IntervalWeight[1];
  if (IntervalCumul!=NULL) delete[] IntervalCumul;
  if (LineEne!=NULL) delete[] LineEne;
  if (LineSigma!=NULL) delete[] LineSigma;
  if (LineIntensity[0]!=NULL) delete[] LineIntensity[0];
  if (LineIntensity[1]!=NULL) delete[] LineIntensity[1];
  if (LineWeight[0]!=NULL) delete[] LineWeight[0];
  if (LineWeight[1]!=NULL) delete[] LineWeight[1];
  if (LineCumul!=NULL) delete[] LineCumul;
}
// constructor
spectrum::spectrum(string dev_name) {
  Runnable = false;
  NInputDevices = 0;
  ContinuousEne = ContSIntensity[0] = ContSIntensity[1] =
    IntervalIntensity[0] = IntervalIntensity[1] = IntervalWeight[0] =
    IntervalWeight[1] = IntervalCumul = LineEne = LineSigma =
    LineIntensity[0] = LineIntensity[1] = LineWeight[0] = 
    LineWeight[1] = LineCumul = NULL;
  EneContinuousNum = EneLineNum = 0;
  rng = NULL;
  SetDevice(dev_name, "spectrum");
}

// initialize loop on events
int spectrum::Begin()
{
  // check if flag for loop on all lines and on all intervals is enabled
  if (LoopFlag==0) LoopIdx = 0; // if not, set loop index to zero
  else { // if yes, initialize all nested loop indexes
    PolIdx = ModeIdx = ContinuousPhotonIdx = IntervalIdx
      = LinePhotonIdx = LineIdx = 0;
    // check if continuous part of the spectrum is defined
    // if not, set mode index to 1 (discrete lines mode)
    if (ContinuousPhotonNum==0 || EneContinuousNum<2) ModeIdx = 1;
  }

  return 0;
}

// next step of the event loop
int spectrum::Next()
{
  // check if flag for loop on all lines and on all intervals is enabled
  if (LoopFlag==0) LoopIdx = 1;  // if not, set loop index to 1
  else { // if yes, update all nested loop indexes
    if (ModeIdx == 0) { // continuous spectrum mode
      IntervalIdx++; // increase index of events within the interval 
      if (IntervalIdx == EneContinuousNum-1) { // check for last event
	IntervalIdx = 0; // if yes, reset index to zero ...
	ContinuousPhotonIdx++; // go to the next interval
	if (ContinuousPhotonIdx
	    == ContinuousPhotonNum) { // check if the last interval is reached
	  // if yes and if discrete lines are present
	  // initialize loop on discrete lines
	  if (LinePhotonNum>0 && EneLineNum>0) {
	    LinePhotonIdx = LineIdx = 0;
	    ModeIdx = 1; // discrete lines mode
	  }
	  else { // if there are no discrete lines reinitialize loop on
	         // intervals of continuous spectrum ...
	    ContinuousPhotonIdx = IntervalIdx = 0;
	    PolIdx++; // with the other polarization type
	  }
	}
      }
    }
    else if (ModeIdx == 1) { // discrete lines mode
      LineIdx++; // increase index of events within the line 
      if (LineIdx == EneLineNum) { // check for last event
	LineIdx = 0; // if yes, reset index to zero ...
	LinePhotonIdx++; // go to the next line
	if (LinePhotonIdx==LinePhotonNum) { //check if the last line is reached
	  // if yes and if continuous spectrum is present
	  // initialize loop on intervals of continuous spectrum
	  if (ContinuousPhotonNum>0 && EneContinuousNum>1) {
	    ContinuousPhotonIdx = IntervalIdx = 0;
	    ModeIdx = 0; // continuous spectrum mode
	  }
	  // if there is no continuous spectrum reinitialize loop on
	  // discrete lines
	  else LinePhotonIdx = LineIdx = 0;
	  PolIdx++; // switch polarization type
	}
      }
    } 
    else return 1;
  }

  return 0;
}

// check if the end of the event loop is reached
bool spectrum::End()
{
  if (LoopFlag==0) return (LoopIdx==1);
  else return (PolIdx > 1); // end when all lines and intervals have been
                            // completed with both polarization types 
}

// event multiplicity
int spectrum::EventMulti()
{
  if (LoopFlag==0) return 1;
  else {
    return 2*(EneContinuousNum*ContinuousPhotonNum +
	      EneLineNum*LinePhotonNum);
  }
}

//////////////////////////////////////////////////////////////////////
// Method for resampling the continuous spectrum
//////////////////////////////////////////////////////////////////////
int spectrum::Resample()
{
  // initialize arrays of energies and intensities of new sampling points
  double *Ene = new double[ResampleNum];
  double *xI = new double[ResampleNum];
  double *yI = new double[ResampleNum];

  double dE = (Emax-Emin)/(ResampleNum-1); // resampling energy step
  double E = Emin; // start with resampling minimum energy
  for (int i=0; i<ResampleNum; i++) { // loop on resampling points
    Ene[i] = E; // resampling point energy
    ExtractSpectrum(E, &xI[i], &yI[i]); // evaluate x, y intensities
                                        // at energy E
    E += dE; // increase energy to the next resampling point
  }
  delete[] ContinuousEne; // deallocate previous arrays
  delete[] ContSIntensity[0];
  delete[] ContSIntensity[1];
  ContinuousEne = Ene; // reinitialize arrays to the resampled ones
  ContSIntensity[0] = xI;
  ContSIntensity[1] = yI;
  EneContinuousNum = ResampleNum;

  return 0;
}
//////////////////////////////////////////////////////////////////////
// method for spectrum initialization before run
//////////////////////////////////////////////////////////////////////
int spectrum::RunInit()
{
  int i;
  double intensity;

  if ((EneLineNum + EneContinuousNum) == 0)
    throw xrmc_exception("Number of discrete lines and number of sampling "
			 "points in the continuous spectrum are both null.\n");

  if (IntervalIntensity[0]!=NULL) delete[] IntervalIntensity[0];
  if (IntervalIntensity[1]!=NULL) delete[] IntervalIntensity[1];
  if (IntervalWeight[0]!=NULL) delete[] IntervalWeight[0];
  if (IntervalWeight[1]!=NULL) delete[] IntervalWeight[1];
  if (IntervalCumul!=NULL) delete[] IntervalCumul;
  // if continuous spectrum is present allocate its arrays
  if (EneContinuousNum > 1) {
    IntervalIntensity[0] = new double[EneContinuousNum-1];
    IntervalIntensity[1] = new double[EneContinuousNum-1];
    IntervalWeight[0] = new double[EneContinuousNum-1];
    IntervalWeight[1] = new double[EneContinuousNum-1];
    IntervalCumul = new double[2*EneContinuousNum-1];
  }
  
  if (LineWeight[0]!=NULL) delete[] LineWeight[0];
  if (LineWeight[1]!=NULL) delete[] LineWeight[1];
  if (LineCumul!=NULL) delete[] LineCumul;
  
  // if discrete lines are present allocate their arrays
  if (EneLineNum > 0) {
    LineWeight[0] = new double[EneLineNum];
    LineWeight[1] = new double[EneLineNum];
    LineCumul = new double[2*EneLineNum+1];
  }
  
  // initialize maximum intensity, total continuous intensity
  // and range of continuous energies
  MaxIntensity = ContinuousIntensity = ContinuousEnergyRange = 0;
  
  // if continuous spectrum is present fill its arrays  
  if (EneContinuousNum > 1) {
    ContinuousEnergyRange = ContinuousEne[EneContinuousNum-1]-ContinuousEne[0];
    cout << "Intervals \n";
    for (i=0; i<EneContinuousNum; i++) { // loop on sampling points
      if (i<EneContinuousNum-1) {
	// the intensity of each interval between adjacent sampling points
	// is evaluated as the area of the trapezium
	IntervalIntensity[0][i] = // x polarization
	  (ContSIntensity[0][i+1] + ContSIntensity[0][i])
	  * (ContinuousEne[i+1] - ContinuousEne[i]) / 2;
	IntervalIntensity[1][i] =  // y polarization
	  (ContSIntensity[1][i+1] + ContSIntensity[1][i])
	  * (ContinuousEne[i+1] - ContinuousEne[i]) / 2;
	// the integral of the continuous spectrum is approximated using
	// trapezoidal areas
	ContinuousIntensity += IntervalIntensity[0][i]+IntervalIntensity[1][i]; 
      }
      intensity = ContSIntensity[0][i] + ContSIntensity[1][i];
      if (i==0 || intensity>MaxIntensity) MaxIntensity = intensity;
    }
  }
  cout << "Total Continuous Intensity: " << ContinuousIntensity << "\n";
  cout << "Maximum Continuous Intensity: " << MaxIntensity << "\n";
  
  DiscreteIntensity = 0;
  for (i=0; i<EneLineNum; i++) { // loop on discrete lines
    // sum contributions of both polarization types
    DiscreteIntensity += LineIntensity[0][i] + LineIntensity[1][i];
  }
  cout << "Total Discrete Intensity: " << DiscreteIntensity << "\n";
  TotalIntensity = ContinuousIntensity + DiscreteIntensity; 

  // evaluate interval weights and cumulative distribution function
  // for continuous spectrum
  if (EneContinuousNum>1) {
    IntervalCumul[0] = 0;
    for (i=0; i<EneContinuousNum-1; i++) { // loop on intervals
      // interval weights for both polarization types
      IntervalWeight[0][i] = IntervalIntensity[0][i] / TotalIntensity
	/ ContinuousPhotonNum;
      IntervalWeight[1][i] = IntervalIntensity[1][i] / TotalIntensity
	/ ContinuousPhotonNum;
      // cumulative sum for interval i, polarization type 0
      IntervalCumul[2*i+1] = IntervalCumul[2*i]
	+ IntervalIntensity[0][i]/ContinuousIntensity;
      // cumulative sum for interval i, polarization type 1
      IntervalCumul[2*i+2] = IntervalCumul[2*i+1] 
      + IntervalIntensity[1][i]/ContinuousIntensity;    
    }
    if (IntervalCumul[2*EneContinuousNum-2] > 1)
      IntervalCumul[2*EneContinuousNum-2] = 1;
  }

  // evaluate weights and cumulative distribution function
  // for discrete lines
  if (EneLineNum>0) {
    LineCumul[0] = 0;
    for (i=0; i<EneLineNum; i++) { // loop on discrete lines
      // line weights for both polarization types
      LineWeight[0][i] = LineIntensity[0][i] / TotalIntensity
	/ LinePhotonNum;
      LineWeight[1][i] = LineIntensity[1][i] / TotalIntensity
	/ LinePhotonNum;
      // cumulative sum for line i, polarization type 0
      LineCumul[2*i+1] = LineCumul[2*i]
	+ LineIntensity[0][i]/DiscreteIntensity;
      // cumulative sum for line i, polarization type 1
      LineCumul[2*i+2] = LineCumul[2*i+1] 
	+ LineIntensity[1][i]/DiscreteIntensity;    
    }
    if (LineCumul[2*EneLineNum] > 1)
      LineCumul[2*EneLineNum] = 1;
  }

  return 0;
}


//////////////////////////////////////////////////////////////////////
// Generates a random energy value and polarization type according 
// to the spectrum distribution
//////////////////////////////////////////////////////////////////////
int spectrum::ExtractEnergy(double *weight, double *Energy, int *polarization)
{
  double R;
 
  // check if flag for loop on all lines and all intervals is activated
  if (LoopFlag == 0) { // if not, extract energy and polarization
                       // using the cumulative distribution approach
    *weight = 1;
    // decide if a discrete line or continuous spectrum will be used
    R = (rng == NULL ? Rnd() : Rnd_r(rng))*(ContinuousIntensity + DiscreteIntensity);
    if (R <= ContinuousIntensity && EneContinuousNum>1)
      ContinuousRandomEnergy(Energy, polarization);
    else DiscreteRandomEnergy(Energy, polarization);
  }
  else { // if yes, extract energy from current interval or line
    *polarization = PolIdx;
    if (ModeIdx == 0) {
      IntervalRandomEnergy(Energy, IntervalIdx, PolIdx);
      *weight = IntervalWeight[PolIdx][IntervalIdx];
    }
    else {
      *Energy = LineEne[LineIdx];
      *weight = LineWeight[PolIdx][LineIdx];
    } 
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Extract the energy value and polarization type from a specified 
// interval of the continuous spectrum distribution
//////////////////////////////////////////////////////////////////////
int spectrum::IntervalRandomEnergy(double *E, int interval_idx, int pol_idx)
{
  // check if a random energy in the interval must be extracted
  if (RandomEneFlag == 1) {
    double y1 = ContSIntensity[pol_idx][interval_idx]; //left side height
    double y2 = ContSIntensity[pol_idx][interval_idx+1]; // right side height
    double R = (rng == NULL ? Rnd() : Rnd_r(rng));
    if (y2 != y1) {
      // extract the energy value in the interval using a linear probability
      // distribution
      double x = (sqrt(y1*y1*(1-R) + y2*y2*R) - y1);
      if ((x + y2 - y1) != x) R = x / (y2 - y1);
    }
    *E = ContinuousEne[interval_idx]*(1.-R) + ContinuousEne[interval_idx+1]*R;
  }
  // otherwise return the middle energy of the interval
  else *E = (ContinuousEne[interval_idx] + ContinuousEne[interval_idx+1])/2; 

  return 0;
}

/////////////////////////////////////////////////////////////////////////////
// returns the intensity of the continuous spectrum at the given energy value
/////////////////////////////////////////////////////////////////////////////
int spectrum::ExtractSpectrum(double trial_energy, double *x_intensity,
		     double *y_intensity)  
{
  int j;
  double t, u;

  // find the index of the interval where the specified energy is located
  Locate(trial_energy, ContinuousEne, EneContinuousNum, &j);
  if (j>=0 && j<EneContinuousNum-1) {
    // use a linear interpolation to evaluate the x, y intensities
    // at the specified energy value
    t = (trial_energy-ContinuousEne[j])/(ContinuousEne[j+1]-ContinuousEne[j]);
    u = 1 - t;
    *x_intensity = ContSIntensity[0][j]*u + ContSIntensity[0][j+1]*t;
    *y_intensity = ContSIntensity[1][j]*u + ContSIntensity[1][j+1]*t;
  }
  else {
    *x_intensity = *y_intensity = 0;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Generates a random energy value according to the discrete part
// of the spectrum
//////////////////////////////////////////////////////////////////////
int spectrum::DiscreteRandomEnergy(double *Energy, int *polarization)
{
  int i, j;
  double R, E, E0;

  // generate a random number and locate it in the cumulative distribution
  R = (rng == NULL ? Rnd() : Rnd_r(rng));
  Locate(R, LineCumul, 2*EneLineNum, &j);
  *polarization = j % 2; // polarization type depends on j being even or odd
  i = j / 2; // line index
  E = LineEne[i]; // line energy
  if (LineSigma[i] != 0) { // line width (standard deviation)
    E0 = E;
    do {
      // if sigma !=0 generate the energy using a Gaussian distribution
      E = E0 + LineSigma[i]*(rng == NULL ? GaussRnd() : GaussRnd_r(rng));
    } while (E <= 0); // E must be positive
  }
  *Energy = E;

  return 0;
}
//////////////////////////////////////////////////////////////////////
// Generates a random energy value according to the continuous part
// of the spectrum
//////////////////////////////////////////////////////////////////////
int spectrum::ContinuousRandomEnergy(double *Energy, int *polarization)
{
  int i, j;
  double R;

  // generate a random number and locate it in the cumulative distribution
  R = (rng == NULL ? Rnd() : Rnd_r(rng));
  Locate(R, IntervalCumul, 2*(EneContinuousNum-1), &j);
  *polarization = j % 2; // polarization type depends on j being even or odd
  i = j / 2; // interval index
  // generate a random energy within the interval i
  IntervalRandomEnergy(Energy, i, *polarization);

  return 0;
}

spectrum *spectrum::Clone(string dev_name) {
	//cout << "Entering spectrum::Clone\n";
	spectrum *clone = new spectrum(dev_name);
	clone->PolarizedFlag = PolarizedFlag;
	clone->LoopFlag = LoopFlag;
	clone->RandomEneFlag = RandomEneFlag;
	clone->ResampleFlag = ResampleFlag;
	clone->EneContinuousNum = EneContinuousNum;
	if (ResampleFlag && ContinuousEne) {
		clone->ContinuousEne =new double[ResampleNum];
		memcpy(clone->ContinuousEne, ContinuousEne, sizeof(double)*ResampleNum);
	}
	else if (ContinuousEne) {
		clone->ContinuousEne =new double[EneContinuousNum];
		memcpy(clone->ContinuousEne, ContinuousEne, sizeof(double)*EneContinuousNum);
	}
	else {
		clone->ContinuousEne = NULL;
	}
	
	if (ResampleFlag && ContSIntensity[0]) {
		clone->ContSIntensity[0] =new double[ResampleNum];
		memcpy(clone->ContSIntensity[0], ContSIntensity[0], sizeof(double)*ResampleNum);
	}
	else if (ContSIntensity[0]) {
		clone->ContSIntensity[0] =new double[EneContinuousNum];
		memcpy(clone->ContSIntensity[0], ContSIntensity[0], sizeof(double)*EneContinuousNum);
	}
	else {
		clone->ContSIntensity[0] = NULL;
	}
	if (ResampleFlag && ContSIntensity[1]) {
		clone->ContSIntensity[1] =new double[ResampleNum];
		memcpy(clone->ContSIntensity[1], ContSIntensity[1], sizeof(double)*ResampleNum);
	}
	else if (ContSIntensity[1]) {
		clone->ContSIntensity[1] =new double[EneContinuousNum];
		memcpy(clone->ContSIntensity[1], ContSIntensity[1], sizeof(double)*EneContinuousNum);
	}
	else {
		clone->ContSIntensity[1] = NULL;
	}

	clone->ContinuousEnergyRange = ContinuousEnergyRange;
	clone->ContinuousIntensity = ContinuousIntensity;
	clone->MaxIntensity = MaxIntensity;
	if (EneContinuousNum > 1) {
    		clone->IntervalIntensity[0] = new double[EneContinuousNum-1];
    		clone->IntervalIntensity[1] = new double[EneContinuousNum-1];
    		clone->IntervalWeight[0] = new double[EneContinuousNum-1];
    		clone->IntervalWeight[1] = new double[EneContinuousNum-1];
    		clone->IntervalCumul = new double[2*EneContinuousNum-1];
		memcpy(clone->IntervalIntensity[0], IntervalIntensity[0], sizeof(double)*(EneContinuousNum-1));
		memcpy(clone->IntervalIntensity[1], IntervalIntensity[1], sizeof(double)*(EneContinuousNum-1));
		memcpy(clone->IntervalWeight[0], IntervalWeight[0], sizeof(double)*(EneContinuousNum-1));
		memcpy(clone->IntervalWeight[1], IntervalWeight[1], sizeof(double)*(EneContinuousNum-1));
		memcpy(clone->IntervalCumul, IntervalCumul, sizeof(double)*(2*EneContinuousNum-1));
	}
	else {
    		clone->IntervalIntensity[0] = 
    		clone->IntervalIntensity[1] = 
    		clone->IntervalWeight[0] = 
    		clone->IntervalWeight[1] = 
    		clone->IntervalCumul = NULL; 
	}

	clone->EneLineNum = EneLineNum;
  	if (EneLineNum > 0) {
    		clone->LineWeight[0] = new double[EneLineNum];
    		clone->LineWeight[1] = new double[EneLineNum];
    		clone->LineCumul = new double[2*EneLineNum+1];
		clone->LineEne = new double[EneLineNum];
		clone->LineSigma = new double[EneLineNum];
		clone->LineIntensity[0] = new double[EneLineNum];
		clone->LineIntensity[1] = new double[EneLineNum];
		memcpy(clone->LineWeight[0], LineWeight[0], sizeof(double)*EneLineNum);
		memcpy(clone->LineWeight[1], LineWeight[1], sizeof(double)*EneLineNum);
		memcpy(clone->LineCumul, LineCumul, sizeof(double)*(2*EneLineNum+1));
		memcpy(clone->LineEne, LineEne, sizeof(double)*EneLineNum);
		memcpy(clone->LineSigma, LineSigma, sizeof(double)*EneLineNum);
		memcpy(clone->LineIntensity[0], LineIntensity[0], sizeof(double)*EneLineNum);
		memcpy(clone->LineIntensity[1], LineIntensity[1], sizeof(double)*EneLineNum);
  	}
	else {
    		clone->LineWeight[0] = 
    		clone->LineWeight[1] = 
    		clone->LineIntensity[0] = 
    		clone->LineIntensity[1] = 
    		clone->LineEne = 
    		clone->LineSigma = 
    		clone->LineCumul = NULL;
	}
	clone->DiscreteIntensity = DiscreteIntensity;
	clone->ResampleNum = ResampleNum;
	clone->Emin = Emin;
	clone->Emax = Emax;
	clone->ContinuousPhotonNum = ContinuousPhotonNum;
	clone->ContinuousPhotonIdx = ContinuousPhotonIdx;
	clone->LinePhotonNum = LinePhotonNum;
	clone->LinePhotonIdx = LinePhotonIdx;
	clone->PolIdx = PolIdx;
	clone->ModeIdx = ModeIdx;
	clone->IntervalIdx = IntervalIdx;
	clone->LineIdx = LineIdx;
	clone->TotalIntensity = TotalIntensity;

	return clone;
}

