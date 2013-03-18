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
#include <glib/gprintf.h> 

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;


int detectorconvolute::Init() {
	

	convolutedImage = double_array3d(ModeNum, N, NBins);	
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
    for (int ipix=0; ipix<N; ipix++) {
      for (int ibin=0; ibin<NBins; ibin++) {
	convolutedImage[mode_idx][ipix][ibin] = 0;
      }
    }
  }

  return 0;
}
int detectorconvolute::ImportDevice(xrmc_device_map *dev_map) {
	detectorarray::ImportDevice(dev_map);

	// check if the input device is defined
	xrmc_device_map::iterator it;
	xrmc_device *dev_pt;

	it = dev_map->find(CompositionName);

	// if not display error and exit
	if (it==dev_map->end())
		throw xrmc_exception(string("Device ") + CompositionName
			 + " not found in device map\n");

	// get device pointer from the device map
	dev_pt = (*it).second;
	// cast it to type basesource*
	Composition = dynamic_cast<composition*>(dev_pt);
	if (Composition==0)
		throw xrmc_exception(string("Device ") + CompositionName
			 + " cannot be casted to type composition\n");

	//now get the crystal phase
	phase_map::iterator it2 = Composition->PhaseMap.find(CrystalPhaseName);
	if (it2 == Composition->PhaseMap.end())
		throw xrmc_exception(string("Device") + CrystalPhaseName
			+ " not found in composition map\n");

	int i_phase = (*it2).second;
	CrystalPhase = &(Composition->Ph[i_phase]);

	//and then the window phase
	if (WindowPhaseName != "Vacuum") {
		it2 = Composition->PhaseMap.find(WindowPhaseName);
		if (it2 == Composition->PhaseMap.end())
			throw xrmc_exception(string("Device") + WindowPhaseName
				+ " not found in composition map\n");

		i_phase = (*it2).second;
		WindowPhase = &(Composition->Ph[i_phase]);
	}
	else {
		WindowPhase = NULL;
	}

	//Initialize xrmc_xmimsim
    	XmiCheckXrmcXmimsimPlugin xmi_check_xrmc_xmimsim_plugin;
    	gchar *module_path;
    	gchar *plugin_dir;

    	if (!g_module_supported()) {
    		throw xrmc_exception("GModule not supported\n");
    	}
   	//environment variable could also be useful here
	if ((plugin_dir = (gchar *) g_getenv("XRMC_XMIMSIM_MODULE")) == NULL)
    		plugin_dir = g_strdup(XRMC_XMIMSIM_LIB);
	
    	module_path = g_strdup_printf("%s" G_DIR_SEPARATOR_S "%s.%s", plugin_dir, "xrmc-xmimsim", G_MODULE_SUFFIX);
    	xrmc_xmimsim = g_module_open(module_path, (GModuleFlags) 0);
    	g_free(module_path);

    	if (!xrmc_xmimsim) {
		g_fprintf(stderr,"GModule error message: %s\n", g_module_error());
    	throw xrmc_exception("Could not load module xrmc-xmimsim\n");
    }

    if (!g_module_symbol(xrmc_xmimsim, "xmi_check_xrmc_xmimsim_plugin", (gpointer *) &xmi_check_xrmc_xmimsim_plugin)) {
        g_module_close(xrmc_xmimsim);
	g_fprintf(stderr, "GModule error message: %s\n", g_module_error());
    	throw xrmc_exception("Could not get symbol xmi_check_xrmc_xmimsim_plugin from module xrmc-xmimsim\n");
    }

    if (xmi_check_xrmc_xmimsim_plugin == NULL) {
        g_module_close(xrmc_xmimsim);
    	throw xrmc_exception("Symbol xmi_check_xrmc_xmimsim_plugin from module xrmc-xmimsim is NULL\n");
    }

    if (xmi_check_xrmc_xmimsim_plugin() == 1) {
        g_fprintf(stdout,"xrmc_xmimsim plugin is functional\n");
    }
	


  	Init();  // initialize detectorconvolute

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
	
	if (xmi_msim_detector_convolute(Image, convolutedImage, det_absorber, xd, ModeNum, N, NBins) == 0)
    		throw xrmc_exception("Error in xmi_msim_detector_convolute\n");


	return 0;
}


