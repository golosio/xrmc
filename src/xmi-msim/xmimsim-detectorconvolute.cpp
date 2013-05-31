/*
Copyright (C) 2013 Tom Schoonjans and Bruno Golosio 

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

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xrmc_math.h"
#include "xrmc_detector.h"
#include "xrmc_photon.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include "xrmc_arrayNd.h"
#include "xrmc_detectorconvolute.h"
#include "xrmc_loadxmimsim.h"
#include <glib/gprintf.h> 

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;

  //destructor
detectorconvolute::~detectorconvolute() {
  if (xd->n_crystal_layers > 0 && xd->crystal_layers != NULL) {
    for (int i = 0 ; i < xd->n_crystal_layers ; i++) {
      free(xd->crystal_layers[i].Z);
      free(xd->crystal_layers[i].weight);
    }
    free(xd->crystal_layers);
  }
  if (det_absorber != NULL) {
    free(det_absorber->Z);
    free(det_absorber->weight);
    free(det_absorber);
  }
  
  if (xd != NULL)
    free(xd);
  //if (ConvolutedImage != NULL)
  //  free(ConvolutedImage);
}

//constructor
detectorconvolute::detectorconvolute(string dev_name)
  : detectorarray(dev_name) {
  //SaveDataName[0]="UnConvolutedImage";
  //SaveDataName.push_back("ConvolutedImage");
  NInputDevices = 2;
  InputDeviceCommand.push_back("CompositionName");
  InputDeviceDescription.push_back("Composition input device name");

  SetDevice(dev_name, "detectorconvolute");

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
  ConvolutedImage = NULL;
  CrystalThickness = 0.0;
  WindowThickness = 0.0;
  det_absorber = NULL;
  Composition = NULL; 
  CrystalPhase = NULL;
  WindowPhase = NULL;
  //xrmc_device(dev_name, "detectorconvolute");

  LoadXrmcXmimsimPlugin();
}

int detectorconvolute::RunInit() {
        ConvolveFlag = 1;
        detectorarray::RunInit();

	xd->live_time = ExpTime;
	xd->zero = 0.0;
	xd->gain = Emax/NBins;

	//initialize XMI-MSIM detector and det_absorber structs
	if (xd->crystal_layers != NULL) {
		free(xd->crystal_layers->Z);
		free(xd->crystal_layers->weight);
		free(xd->crystal_layers);
	}
	xd->n_crystal_layers = 1;
	xd->crystal_layers = (struct xmi_layer *) malloc(sizeof(struct xmi_layer)); 
	xd->crystal_layers[0].n_elements = CrystalPhase->NElem;
	xd->crystal_layers[0].Z = (int *) malloc(sizeof(int)*CrystalPhase->NElem); 
	xd->crystal_layers[0].weight = (double *) malloc(sizeof(double)*CrystalPhase->NElem); 
	memcpy(xd->crystal_layers[0].Z, CrystalPhase->Z, sizeof(int)*CrystalPhase->NElem);
	memcpy(xd->crystal_layers[0].weight, CrystalPhase->W, sizeof(double)*CrystalPhase->NElem);
	xd->crystal_layers[0].thickness = CrystalThickness;
	xd->crystal_layers[0].density = CrystalPhase->Rho;

	if (det_absorber != NULL) {
		free(det_absorber->Z);
		free(det_absorber->weight);
		free(det_absorber);
	}
	if (WindowPhase != NULL) {
		det_absorber = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
		det_absorber->n_elements = WindowPhase->NElem;
		det_absorber->Z = (int *) malloc(sizeof(int)*WindowPhase->NElem);
		det_absorber->weight = (double *) malloc(sizeof(double)*WindowPhase->NElem);
		memcpy(det_absorber->Z, WindowPhase->Z,sizeof(int)*WindowPhase->NElem);
		memcpy(det_absorber->weight, WindowPhase->W,sizeof(double)*WindowPhase->NElem);
		det_absorber->thickness = WindowThickness;
		det_absorber->density = WindowPhase->Rho;
	}
	else {
		det_absorber = NULL;
	}

	return 0;
}

int detectorconvolute::Clear()
{
  detectorarray::Clear();

  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iy=0; iy<NY; iy++) {
      for (int ix=0; ix<NX; ix++) {
	for (int ibin=0; ibin<NBins; ibin++) {
	  ConvolutedImage[mode_idx*NBins+ibin][iy][ix] = 0;
	}
      }
    }
  }

  return 0;
}

int detectorconvolute::CastInputDevices() {
	detectorarray::CastInputDevices();
	// cast it to type composition*
	Composition = dynamic_cast<composition*>(InputDevice[1]);
	if (Composition==0)
		throw xrmc_exception(string("Device ") + InputDeviceName[1]
			 + " cannot be casted to type composition\n");

	//now get the crystal phase
	phase_map::iterator it2 = Composition->PhaseMap.find(CrystalPhaseName);
	if (it2 == Composition->PhaseMap.end())
		throw xrmc_exception(string("Phase ") + CrystalPhaseName
			+ " not found in composition map\n");

	int i_phase = (*it2).second;
	CrystalPhase = &(Composition->Ph[i_phase]);

	//and then the window phase
	if (WindowPhaseName != "Vacuum") {
		it2 = Composition->PhaseMap.find(WindowPhaseName);
		if (it2 == Composition->PhaseMap.end())
			throw xrmc_exception(string("Phase ") + WindowPhaseName
				+ " not found in composition map\n");

		i_phase = (*it2).second;
		WindowPhase = &(Composition->Ph[i_phase]);
	}
	else {
		WindowPhase = NULL;
	}

  	//Init();  // initialize detectorconvolute

	return 0;
}

int detectorconvolute::Run() {
	detectorarray::Run();
	XmiMsimDetectorConvolute xmi_msim_detector_convolute;

    	if (!g_module_symbol(xrmc_xmimsim, "xmi_msim_detector_convolute", (gpointer *) &xmi_msim_detector_convolute)) {
        	g_module_close(xrmc_xmimsim);
		g_fprintf(stderr, "GModule error message: %s\n", g_module_error());
    		throw xrmc_exception("Could not get symbol xmi_msim_detector_convolute from module xrmc-xmimsim\n");
    	}

    	if (xmi_msim_detector_convolute == NULL) {
        	g_module_close(xrmc_xmimsim);
    		throw xrmc_exception("Symbol xmi_msim_detector_convolute from module xrmc-xmimsim is NULL\n");
    	}
	
	if (xmi_msim_detector_convolute(Image, ConvolutedImage, det_absorber, xd, ModeNum, NBins, NY, NX) == 0)
    		throw xrmc_exception("Error in xmi_msim_detector_convolute\n");


	return 0;
}


