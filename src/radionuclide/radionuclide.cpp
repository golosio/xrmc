/*
Copyright (C) 2014 Bruno Golosio and Tom Schoonjans

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
//          radionuclide.cpp           //
//           25/06/2014                //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
// Methods of the class radionuclide 
//

#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include "xrmc_radionuclide.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include <cstdlib>
#include <xraylib.h>

using namespace std;
using namespace xrmc_algo;


radionuclide::radionuclide(string dev_name) : spectrum(dev_name){
  Runnable = false;
  NInputDevices = 0;
  ContinuousEne = ContSIntensity[0] = ContSIntensity[1] =
    IntervalIntensity[0] = IntervalIntensity[1] = IntervalWeight[0] =
    IntervalWeight[1] = IntervalCumul = LineEne = LineSigma =
    LineIntensity[0] = LineIntensity[1] = LineWeight[0] = 
    LineWeight[1] = LineCumul = NULL;
  EneContinuousNum = EneLineNum = 0;
  Rng = NULL;
  rnd = NULL;
  SetDevice(dev_name, "radionuclide");
} 

radionuclide::~radionuclide() {
	FreeRadioNuclideData(rnd);
	rnd = NULL;
	//now the spectrum destructor should be called...
}

int radionuclide::RunInit() {
	double realActivity;

	if (Unit == 1) {
		//mCi
		realActivity = Activity*3.7E7;
	}
	else if (Unit == 2) {
		//Ci
		realActivity = Activity*3.7E10;
	}
	else if (Unit == 3) {
		//GBq
		realActivity = Activity*1E9;
	}
	else if (Unit == 4) {
		//Bq -> SI!
		realActivity = Activity;
	}
                                                                                   

	//copy everything to spectrum members
	LineEne = new double[rnd->nXrays+rnd->nGammas];
	LineSigma = new	double[rnd->nXrays+rnd->nGammas];
	LineIntensity[0] = new double[rnd->nXrays+rnd->nGammas];
	LineIntensity[1] = new double[rnd->nXrays+rnd->nGammas];
	EneLineNum = rnd->nXrays+rnd->nGammas;
	int i;
	for (i = 0 ; i < rnd->nXrays ; i++) {
		LineEne[i] = LineEnergy(rnd->Z_xray, rnd->XrayLines[i]);
		LineSigma[i] = 0.0;
		LineIntensity[0][i] = LineIntensity[1][i] = realActivity*rnd->XrayIntensities[i]/2.0;
	}
	for (i = 0 ; i < rnd->nGammas; i++) {
		LineEne[i+rnd->nXrays] = rnd->GammaEnergies[i];
		LineSigma[i+rnd->nXrays] = 0.0;
		LineIntensity[0][i+rnd->nXrays] = LineIntensity[1][i+rnd->nXrays] = realActivity*rnd->GammaIntensities[i]/2.0;
	}

	for (i = 0 ; i < EneLineNum ; i++)
		cout << LineEne[i] << "\t" << LineIntensity[0][i] << endl;

	spectrum::RunInit();


	return 1;
}
