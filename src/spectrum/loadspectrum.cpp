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
//          loadspectrum.cpp           //
//            14/02/2013               //
//     Author : Bruno Golosio          //
/////////////////////////////////////////
// Load spectrum parameters
//

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_spectrum.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
// method for loading spectrum parameters from a file
//////////////////////////////////////////////////////////////////////
int spectrum::Load(FILE *fp)
{
  int i;
  char comm_ch[MAXSTRLEN], spectrum_file[MAXSTRLEN];
  string comm="";
  FILE *sp_fp;

  cout << "Spectrum File\n";

  while (!feof(fp) && comm!="End") {
    GetToken(fp, comm_ch); // get a command/variable name from input file
    comm = comm_ch;
    // parse the command and decide what to do
    if(comm=="PolarizedFlag") { // flag for unpolarized(0) or polarized(1) beam
      GetIntToken(fp, &PolarizedFlag);
      if (PolarizedFlag == 0) cout << "Unpolarized beam\n";
      else if (PolarizedFlag == 1) cout << "Polarized beam\n";
      else throw xrmc_exception("Wrong polarization flag.\n");
    }
    else if(comm=="LoopFlag") { // flag for modality of energy extraction
      GetIntToken(fp, &LoopFlag);
      if (LoopFlag==0)
	cout << "Extract random energies on the whole spectrum\n";
      else if (LoopFlag==1)
	cout << "Loop on all lines and sampling points\n";
      else throw xrmc_exception("Wrong loop flag.\n");
    }
    else if(comm=="ContinuousPhotonNum") {
      // Number of samples for each interval in the continuous spectrum 
      GetIntToken(fp, &ContinuousPhotonNum);
      cout <<"Number of samples for each interval in the continuous spectrum: " 
	   << ContinuousPhotonNum << "\n";
    }
    else if(comm=="LinePhotonNum") {
      // Number of samples for each line in the spectrum 
      GetIntToken(fp, &LinePhotonNum);
      cout << "Number of samples for each line in the spectrum: " 
	   << LinePhotonNum << "\n";
    }
    else if(comm=="RandomEneFlag") {
      // flag for enabling(1)/disabling(0) random energy in each interval
      GetIntToken(fp, &RandomEneFlag);
      if (RandomEneFlag==0) cout << "Random energy in intervals disabled\n";
      else if (RandomEneFlag==1) cout << "Random energy in intervals enabled\n";
      else throw xrmc_exception("Wrong random energy in intervals flag.\n");
    }
    else if(comm=="Lines") { // read discrete lines
      GetIntToken(fp, &EneLineNum); // number of discrete lines
      cout << "Number of lines in the spectrum: " << EneLineNum << "\n";
      if (EneLineNum > 0) {
	LineEne = new double[EneLineNum];
	LineSigma = new double[EneLineNum];
	LineIntensity[0] = new double[EneLineNum];
	LineIntensity[1] = new double[EneLineNum];
	cout << "Energy Lines :\n";
	for (i=0; i<EneLineNum; i++) {
	  GetDoubleToken(fp, &LineEne[i]);
	  GetDoubleToken(fp, &LineSigma[i]);
	  GetDoubleToken(fp, &LineIntensity[0][i]);
	  if (PolarizedFlag == 0) {
	    cout << LineEne[i] << "\t" << LineSigma[i] << "\t"
		 <<  LineIntensity[0][i] << "\n";
	    LineIntensity[0][i] /= 2;
	    LineIntensity[1][i] = LineIntensity[0][i];
	  }
	  else {
	    GetDoubleToken(fp, &LineIntensity[1][i]);
	    cout << LineEne[i] << "\t" << LineSigma[i] << "\t"
		 <<  LineIntensity[0][i] << "\t" << LineIntensity[1][i] << "\n";
	  }
	}
      }
    }
    else if(comm=="ContinuousSpectrum") { // read continuous part of spectrum
      GetIntToken(fp, &EneContinuousNum);
      cout << "Number of sampling points in the continuous spectrum: "
	   << EneContinuousNum << "\n";
      if (EneContinuousNum < 0 || EneContinuousNum==1)
	throw xrmc_exception("Number of of sampling points in the continuous "
			     "spectrum must be >1 or 0.\n");
      if (EneContinuousNum > 0) {
	ContinuousEne = new double[EneContinuousNum];
	ContSIntensity[0] = new double[EneContinuousNum];
	ContSIntensity[1] = new double[EneContinuousNum];
	GetToken(fp, spectrum_file); // continuous spectrum file name
	cout << "Continuous spectrum file name: " << spectrum_file << "\n";
	if ((sp_fp = fopen(spectrum_file,"r")) == NULL)
	  throw xrmc_exception("Continuous spectrum file can not be opened.");
	cout << "Continuos spectrum :\n";
	for (i=0; i<EneContinuousNum; i++) {
	  GetDoubleToken(sp_fp, &ContinuousEne[i]);
	  GetDoubleToken(sp_fp, &ContSIntensity[0][i]);
	  if (PolarizedFlag == 0) {
	    cout << ContinuousEne[i] << "\t" << ContSIntensity[0][i] << "\n";
	    ContSIntensity[0][i] /= 2;
	    ContSIntensity[1][i] = ContSIntensity[0][i];
	  }
	  else {
	    GetDoubleToken(sp_fp, &ContSIntensity[1][i]);
	    cout << ContinuousEne[i] << "\t" << ContSIntensity[0][i] << "\t"
		 << ContSIntensity[1][i] << "\n";
	  }
	}
	fclose(sp_fp);
      }
    }
    else if(comm=="Resample") { // resample continuous part of spectrum
      cout << "Resample continuous spectrum\n";
      ResampleFlag = 1;
      if (EneContinuousNum<=1)
	throw xrmc_exception("Number of of sampling points in the continuous "
			     "spectrum must be >1 for resampling.\n");
      // if (RandomEneFlag != 1) { // check
      cout << "Resample continuous spectrum\n";
      GetDoubleToken(fp, &Emin);
      cout << "Emin: " << Emin << "\n";
      GetDoubleToken(fp, &Emax);
      cout << "Emax: " << Emax << "\n";
      GetIntToken(fp, &ResampleNum);
      cout << "Number of resampling points: " << ResampleNum << "\n";
      Resample();
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in spectrum input file"); 
    }
  }   
  if ((EneLineNum + EneContinuousNum) == 0)
    throw xrmc_exception("Number of discrete lines and number of sampling "
			 "points in the continuous spectrum are both null.\n");


  Init();
  return 0;
}

//////////////////////////////////////////////////////////////////////
int spectrum::SetDefault()
// set default values for spectrum parameters
{
  PolarizedFlag=0; // flag for unpolarized(0) or polarized(1) beam
  LoopFlag=0; // flag for modality of energy extraction
  // 0: Extract random energies on the whole spectrum

  // Number of samples for each interval in the continuous spectrum 
  ContinuousPhotonNum=1;

  // Number of samples for each line in the spectrum 
  LinePhotonNum=1;

  // flag for enabling(1)/disabling(0) random energy in each interval
  RandomEneFlag=1;

  EneLineNum=0; // number of discrete lines in the spectrum

  // Number of sampling points in the continuous spectrum
  EneContinuousNum=0;

  ResampleFlag=0; // Do not Resample continuous spectrum

  return 0;
}
