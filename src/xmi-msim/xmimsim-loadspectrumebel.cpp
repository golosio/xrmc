/*
Copyright (C) 2013 Bruno Golosio and Tom Schoonjans

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
//          loadspectrum_ebel.cpp      //
//            24/05/2013               //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
// Load spectrum Ebel parameters
//

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_spectrum_ebel.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include <xraylib.h>
#include <unistd.h>

using namespace std;
using namespace gettoken;

int spectrum_ebel::Load(istream &fs)
{
  int i;
  string comm="";
  string material="";

  cout << "Ebel Spectrum Parameters\n";

  // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    if(comm=="TransmissionFlag") { 
	GetIntToken(fs, &TransmissionFlag);
	if (TransmissionFlag == 0) cout << "No transmission tube" << endl;
	else if (TransmissionFlag == 1) cout << "Transmission tube" << endl;
        else throw xrmc_exception("Wrong transmission flag.\n");
    }
    else if (comm == "TubeCurrent") {
    	GetDoubleToken(fs, &TubeCurrent);
	if (TubeCurrent <= 0.0)
		throw xrmc_exception("TubeCurrent must be a strictly positive real number\n");
	else cout << "TubeCurrent" << "\t" << TubeCurrent << endl;
    }
    else if (comm == "TubeVoltage") {
    	GetDoubleToken(fs, &TubeVoltage);
	if (TubeVoltage <= 0.0)
		throw xrmc_exception("TubeVoltage must be a strictly positive real number\n");
	else cout << "TubeVoltage" << "\t" << TubeVoltage << endl;
    }
    else if (comm == "ElectronAngle") {
    	GetDoubleToken(fs, &ElectronAngle);
	if (ElectronAngle <= 0.0)
		throw xrmc_exception("ElectronAngle must be a strictly positive real number\n");
	else cout << "ElectronAngle" << "\t" << ElectronAngle << endl;
    }
    else if (comm == "XrayAngle") {
    	GetDoubleToken(fs, &XrayAngle);
	if (XrayAngle <= 0.0)
		throw xrmc_exception("XrayAngle must be a strictly positive real number\n");
	else cout << "XrayAngle" << "\t" << XrayAngle << endl;
    }
    else if (comm == "IntervalWidth") {
    	GetDoubleToken(fs, &IntervalWidth);
	if (IntervalWidth <= 0.0)
		throw xrmc_exception("IntervalWidth must be a strictly positive real number\n");
	else cout << "IntervalWidth" << "\t" << IntervalWidth << endl;
    }
    else if (comm == "AnodeMaterial") {
    	GetToken(fs, material);
	AnodeMaterial = SymbolToAtomicNumber(material.c_str(), NULL);
	if (AnodeMaterial == 0) 
		throw xrmc_exception("AnodeMaterial must be a chemical symbol\n");
	else cout << "AnodeMaterial" << "\t" << AnodeMaterial << endl; 
    }
    else if (comm == "FilterMaterial") {
    	GetToken(fs, material);
	FilterMaterial = SymbolToAtomicNumber(material.c_str(), NULL);
	if (FilterMaterial == 0) 
		throw xrmc_exception("FilterMaterial must be a chemical symbol\n");
	else cout << "FilterMaterial" << "\t" << FilterMaterial << endl; 
    }
    else if (comm == "WindowMaterial") {
    	GetToken(fs, material);
	WindowMaterial = SymbolToAtomicNumber(material.c_str(), NULL);
	if (WindowMaterial == 0) 
		throw xrmc_exception("WindowMaterial must be a chemical symbol\n");
	else cout << "WindowMaterial" << "\t" << WindowMaterial << endl; 
    }
    else if (comm == "AnodeDensity") {
  	GetDoubleToken(fs, &AnodeDensity);
	if (AnodeDensity <= 0.0)
		throw xrmc_exception("AnodeDensity must be a strictly positive real number\n");
	else cout << "AnodeDensity" << "\t" << AnodeDensity << endl;
    }
    else if (comm == "FilterDensity") {
  	GetDoubleToken(fs, &FilterDensity);
	if (FilterDensity <= 0.0)
		throw xrmc_exception("FilterDensity must be a strictly positive real number\n");
	else cout << "FilterDensity" << "\t" << FilterDensity << endl;
    }
    else if (comm == "WindowDensity") {
  	GetDoubleToken(fs, &WindowDensity);
	if (WindowDensity <= 0.0)
		throw xrmc_exception("WindowDensity must be a strictly positive real number\n");
	else cout << "WindowDensity" << "\t" << WindowDensity << endl;
    }
    else if (comm == "AnodeThickness") {
  	GetDoubleToken(fs, &AnodeThickness);
	if (AnodeThickness <= 0.0)
		throw xrmc_exception("AnodeThickness must be a strictly positive real number\n");
	else cout << "AnodeThickness" << "\t" << AnodeThickness << endl;
    }
    else if (comm == "FilterThickness") {
  	GetDoubleToken(fs, &FilterThickness);
	if (FilterThickness <= 0.0)
		throw xrmc_exception("FilterThickness must be a strictly positive real number\n");
	else cout << "FilterThickness" << "\t" << FilterThickness << endl;
    }
    else if (comm == "WindowThickness") {
  	GetDoubleToken(fs, &WindowThickness);
	if (WindowThickness <= 0.0)
		throw xrmc_exception("WindowThickness must be a strictly positive real number\n");
	else cout << "WindowThickness" << "\t" << WindowThickness << endl;
    }
    else if (comm == "TransmissionEfficiencyFile") {
	GetToken(fs, TransmissionEfficiencyFile);
	if (access(TransmissionEfficiencyFile.c_str(), R_OK) != 0)
		throw xrmc_exception("TransmissionEfficiencyFile is not readable\n");
	else cout << "TransmissionEfficiencyFile" << "\t" << TransmissionEfficiencyFile << endl;
    }
    else if(comm=="Resample") { // resample continuous part of spectrum
      	GetIntToken(fs, &ResampleFlag);
      	if (ResampleFlag==0) {
      		cout << "Resample continuous spectrum disabled\n";
		continue;
	}
      	else if (ResampleFlag==1) {
      	    cout << "Resample continuous spectrum\n";
            // if (RandomEneFlag != 1) { // check
      	    GetDoubleToken(fs, &Emin);
            cout << "Emin: " << Emin << "\n";
            GetDoubleToken(fs, &Emax);
            cout << "Emax: " << Emax << "\n";
            GetIntToken(fs, &ResampleNum);
            cout << "Number of resampling points: " << ResampleNum << "\n";
            //Resample(); -> will have to be done later...
       }
       else throw xrmc_exception("Wrong value for Resample flag\n");
    }
    else if(comm=="RandomEneFlag") {
      // flag for enabling(1)/disabling(0) random energy in each interval
      GetIntToken(fs, &RandomEneFlag);
      if (RandomEneFlag==0) cout << "Random energy in intervals disabled\n";
      else if (RandomEneFlag==1) cout << "Random energy in intervals enabled\n";
      else throw xrmc_exception("Wrong random energy in intervals flag.\n");
    }
    else if(comm=="LoopFlag") { // flag for modality of energy extraction
      GetIntToken(fs, &LoopFlag);
      if (LoopFlag==0)
	cout << "Extract random energies on the whole spectrum\n";
      else if (LoopFlag==1)
	cout << "Loop on all lines and sampling points\n";
      else throw xrmc_exception("Wrong loop flag.\n");
    }
    else if(comm=="ContinuousPhotonNum") {
      // Number of samples for each interval in the continuous spectrum 
      GetIntToken(fs, &ContinuousPhotonNum);
      cout <<"Number of samples for each interval in the continuous spectrum: " 
	   << ContinuousPhotonNum << "\n";
    }
    else if(comm=="LinePhotonNum") {
      // Number of samples for each line in the spectrum 
      GetIntToken(fs, &LinePhotonNum);
      cout << "Number of samples for each line in the spectrum: " 
	   << LinePhotonNum << "\n";
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empty string\n";
    }
    else {
      throw xrmc_exception("syntax error in spectrum input file"); 
    }
  }
  return 0;
}

int spectrum_ebel::SetDefault() 
{
	cout << "Setting spectrum_ebel defaults\n";
	TubeCurrent = 1.0;
	TubeVoltage = 30.0;
	AnodeMaterial = 47; //Ag
	AnodeDensity = 10.5;
	AnodeThickness = 0.0002;	
	FilterMaterial = 0;
	FilterDensity = 0.0;
	FilterThickness = 0.0;
	WindowMaterial = 0;
	WindowDensity = 0.0;
	WindowThickness = 0.0;
	ElectronAngle = 60.0;
	XrayAngle = 60.0;
	IntervalWidth = 0.1;
	TransmissionFlag = 0;
	spectrum::SetDefault();

	UnitSolidAngleFlag = 1;
	return 0;
}

