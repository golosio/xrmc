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
//          spectrum.h                 //
//            14/02/2013               //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// spectrum class definition
//

#ifndef SpectrumH
#define SpectrumH
#include <string>
#include "xrmc_device.h"
#include "randmt.h"
using namespace std;

// spectrum class definition, member variables and functions
class spectrum : public xrmc_device
{
 private:
  int PolarizedFlag; // flag for polarized(1) / unpolarized(0) beam
  int LoopFlag; // flag for loop on all lines and all intervals of the spectrum
  int RandomEneFlag; // flag for extracting random energies in the intervals
  int ResampleFlag; // flag for resampling the continuous spectrum

  int EneContinuousNum; // Number of sampling points in the continuous spectrum
  double *ContinuousEne; // energies of sampling points
  double *ContSIntensity[2]; // intensity at sampling points
  double ContinuousEnergyRange; // energy range of the continuous spectrum
  double ContinuousIntensity; // total continuous intensity
  double MaxIntensity;        // maximum continuous intensity
  double *IntervalIntensity[2]; // intensities of continuous spectrum intervals
  double *IntervalWeight[2]; // weights of continuous spectrum intervals
  double *IntervalCumul; // cumulative distribution  of the intervals

  int EneLineNum; // number of discrete lines of the spectrum
  double *LineEne; // energies of the discrete lines
  double *LineSigma; // widths (sigma) of the discrete lines
  double *LineIntensity[2]; // intensities of the discrete lines
  double DiscreteIntensity; // total intensity of the discrete spectrum
  double *LineWeight[2]; // weights of the discrete lines
  double *LineCumul; // cumulative distribution  of the lines

  int ResampleNum; // number of resampling points of the continuous spectrum
  double Emin, Emax; // minimum and maximum resampling energy

  int ContinuousPhotonNum; // num. of events to be extracted for each interval
  int ContinuousPhotonIdx; // index of the event extracted for the interval
  int LinePhotonNum; // num. of events to be extracted for each discrete line
  int LinePhotonIdx; // index of the event extracted for the line
  int PolIdx; // index of polarization type (0: x, 1: y)
  int ModeIdx; // mode index: continuous spectrum (0) or discrete lines(1)
  int IntervalIdx; // index of the interval of the continuous spectrum
  int LineIdx; // index of the discrete line

  int Resample(); // method for resampling the continuous spectrum
 public:
  double TotalIntensity;
  randmt_t *rng;

  ~spectrum() { // destructor
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
  spectrum(string dev_name) {
    ContinuousEne = ContSIntensity[0] = ContSIntensity[1] =
      IntervalIntensity[0] = IntervalIntensity[1] = IntervalWeight[0] =
      IntervalWeight[1] = IntervalCumul = LineEne = LineSigma =
      LineIntensity[0] = LineIntensity[1] = LineWeight[0] = 
      LineWeight[1] = LineCumul = NULL;
    EneContinuousNum = EneLineNum = 0;
    rng = NULL;
    xrmc_device(dev_name, "spectrum");
  }
  int Load(FILE *fp); // load spectrum parameters from file
  int Begin(); // begin the event loop
  int Next();  // next step of the event loop  
  bool End(); // check if the end of the event loop is reached
  int EventMulti(); // event multiplicity
  int Init(); // initialize the spectrum
  int SetDefault(); // set default values for spectrum parameters

  // Generates a random energy value and polarization type
  int ExtractEnergy(double *weight, double *Energy, int *polarization);

  // Generates a random energy value from the continuous  spectrum
  int ContinuousRandomEnergy(double *Energy, int *polarization);

  // returns the intensity of the continuous spectrum at the given energy value
  int ExtractSpectrum(double trial_energy, double *x_intensity,
		     double *y_intensity);

  // Generates a random energy value from the discrete part of the spectrum
  int DiscreteRandomEnergy(double *Energy, int *polarization);

  // Extract the energy value and polarization type from a specified 
  // interval of the continuous spectrum distribution
  int IntervalRandomEnergy(double *E, int interval_idx, int pol_idx);

  spectrum *Clone(string dev_name);
}; 

#endif
