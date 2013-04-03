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

#ifndef DETECTOR_CONVOLUTE_H
#define DETECTOR_CONVOLUTE_H


#include "xrmc_device.h"
#include "xrmc_composition.h"
#include <xmi_msim.h>
#include <gmodule.h>
#include <stdlib.h>

typedef int (*XmiCheckXrmcXmimsimPlugin) (void);
typedef int (*XmiMsimDetectorConvolute) (double ***Image, double ***convolutedImage, struct xmi_layer *det_absorber, struct xmi_detector *xd, int ModeNum, int N, int NBins);


class detectorconvolute : public detectorarray
{
 public:
  struct xmi_detector *xd;
  double ***convolutedImage;
  struct xmi_layer *det_absorber;
  
  ~detectorconvolute(); //destructor
  detectorconvolute(string dev_name); //constructor

  int CastInputDevices(); // cast input device method
  int Load(istream &fs); // load detector parameters, position, orientation
  int Save(string file_name); // save the acquired image in a file
  int SaveUnconvoluted(string file_name); // save the acquired, unconvoluted image in a file
  int SetDefault(); // set the default values for detector parameters
  int Run(); // calculate detector convolution
  int Clear(); // clear 
  int RunInit(); // detectorarray initialization before run

 private:
  detectorarray *DetectorArray;
  composition *Composition;
  phase *CrystalPhase;
  phase *WindowPhase;
  //string CompositionName;
  string CrystalPhaseName;
  string WindowPhaseName;
  double CrystalThickness;
  double WindowThickness;
  GModule *xrmc_xmimsim;
  
};
#endif
