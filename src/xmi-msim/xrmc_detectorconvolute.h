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


class detectorconvolute : public detectorarray
{
 public:
  struct xmi_detector *xd;
  double ***convolutedImage;
  struct xmi_layer *exc_absorber;


  //destructor
  ~detectorconvolute() {
 	if (xd->n_crystal_layers > 0 && xd->crystal_layers != NULL) {
		for (int i = 0 ; i < xd->n_crystal_layers ; i++) {
			free(xd->crystal_layers[i].Z);
			free(xd->crystal_layers[i].weight);
		}
		free(xd->crystal_layers);
	}
	if (exc_absorber != NULL) {
		free(exc_absorber->Z);
		free(exc_absorber->weight);
		free(exc_absorber);
	}

	if (xd != NULL)
		free(xd);
	if (convolutedImage != NULL)
		free(convolutedImage);
  }
  //constructor
  detectorconvolute(string dev_name) : detectorarray(dev_name) {
  	xd = (struct xmi_detector *) malloc(sizeof(struct xmi_detector));
 	xd->detector_type = -1;
	xd->live_time = 0.0;
	xd->pulse_width = 0.0;
	xd->gain = 0.0;
	xd->zero = 0.0;
	xd->fano = 0.0;
	xd->noise = 0.0;
	xd->max_convolution_energy = 0.0;
	xd->n_crystal_layers = 0;
	xd->crystal_layers = NULL;
	convolutedImage = NULL;
	CrystalThickness = 0.0;
	WindowThickness = 0.0;
	exc_absorber = NULL;
  	Composition = NULL; 
  	CrystalPhase = NULL;
  	WindowPhase = NULL;
	//xrmc_device(dev_name, "detectorconvolute");
  }
  int ImportDevice(xrmc_device_map *dev_map); // import device method
  int Load(FILE *fp); // load detector parameters, position, orientation
  int Save(string file_name); // save the acquired image in a file
  int SaveUnconvoluted(string file_name); // save the acquired, unconvoluted image in a file
  int SetDefault(); // set the default values for detector parameters
  int Run(); // calculate detector convolution
  int Clear(); // clear 
  int Init(); // detectorarray initialization

 private:
  detectorarray *DetectorArray;
  composition *Composition;
  phase *CrystalPhase;
  phase *WindowPhase;
  string DetectorArrayName;
  string CompositionName;
  string CrystalPhaseName;
  string WindowPhaseName;
  double CrystalThickness;
  double WindowThickness;
  GModule *xrmc_xmimsim;
  
};
#endif
