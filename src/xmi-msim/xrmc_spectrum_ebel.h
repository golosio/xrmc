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
//          xrmc_spectrum_ebel.h       //
//            24/05/2013               //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
// spectrum class definition
//

#ifndef XRMC_SPECTRUM_EBEL_H
#define XRMC_SPECTRUM_EBEL_H
#include "xrmc_spectrum.h"
#include <xmi_msim.h>

typedef int (*XmiMsimTubeEbel)(xmi_layer *tube_anode, xmi_layer *tube_window,
                  xmi_layer *tube_filter, double tube_voltage,
                  double tube_current, double tube_angle_electron,
                  double tube_angle_xray, double tube_delta_energy,
                  double tube_solid_angle, int tube_transmission,
		  size_t tube_nefficiencies, double *tube_energies, double *tube_efficiencies,
                  xmi_excitation **ebel_spectrum
                  );

typedef gboolean (*XmiMsimTransmissionEfficiencyRead)(
		const char *filename,
		size_t *nefficiencies,
		double **energies,
		double **efficiencies,
		GError **error
		);

class spectrum_ebel : public spectrum 
{
 public:
  double TubeCurrent;
  double TubeVoltage;
  int AnodeMaterial;
  double AnodeDensity;
  double AnodeThickness;
  int TransmissionFlag;
  int WindowMaterial;
  double WindowDensity;
  double WindowThickness;
  int FilterMaterial;
  double FilterDensity;
  double FilterThickness;
  double ElectronAngle;
  double XrayAngle;
  double IntervalWidth;
  std::string TransmissionEfficiencyFile;

  spectrum_ebel(std::string dev_name);
  virtual int RunInit(); // initialize the spectrum before run
  virtual int SetDefault(); // set default values for spectrum parameters
  virtual int Load(std::istream &fs); // load spectrum parameters from file

 private:


};

#endif
