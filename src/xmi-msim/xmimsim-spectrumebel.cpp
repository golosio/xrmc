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
//          spectrum_ebel.cpp          //
//           24/05/2013                //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
// Methods of the class spectrum
//

#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include "xrmc_spectrum_ebel.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include <xmi_msim.h>
#include <cstdlib>
#include "xrmc_loadxmimsim.h"
#include <glib/gstdio.h>

using namespace std;
using namespace xrmc_algo;


spectrum_ebel::spectrum_ebel(string dev_name) : spectrum(dev_name){
  Runnable = false;
  NInputDevices = 0;
  ContinuousEne = ContSIntensity[0] = ContSIntensity[1] =
    IntervalIntensity[0] = IntervalIntensity[1] = IntervalWeight[0] =
    IntervalWeight[1] = IntervalCumul = LineEne = LineSigma =
    LineIntensity[0] = LineIntensity[1] = LineWeight[0] = 
    LineWeight[1] = LineCumul = NULL;
  EneContinuousNum = EneLineNum = 0;
  Rng = NULL;
  SetDevice(dev_name, "spectrum_ebel");
} 


int spectrum_ebel::RunInit() {

	LoadXrmcXmimsimPlugin();

  	struct xmi_excitation *exc;
  	struct xmi_layer *anode;
  	struct xmi_layer *window = NULL;
  	struct xmi_layer *filter = NULL;
	//anode etc
	anode = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
	anode->Z = (int *) malloc(sizeof(int));
	anode->weight = (double *) malloc(sizeof(double));
	anode->n_elements = 1;
	anode->Z[0] = AnodeMaterial;
	anode->weight[0] = 1.0;
	anode->density = AnodeDensity;
	anode->thickness = AnodeThickness;

	if (FilterDensity > 0.0 && FilterThickness > 0.0 && FilterMaterial > 0) {
		cout << "Filter will be used for X-ray tube" << endl;
		filter = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
		filter->Z = (int *) malloc(sizeof(int));
		filter->weight = (double *) malloc(sizeof(double));
		filter->n_elements = 1;
		filter->Z[0] = FilterMaterial;
		filter->weight[0] = 1.0;
		filter->density = FilterDensity;
		filter->thickness = FilterThickness;
	}
	
	if (WindowDensity > 0.0 && WindowThickness > 0.0 && WindowMaterial > 0) {
		cout << "Window will be used for X-ray tube" << endl;
		window = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
		window->Z = (int *) malloc(sizeof(int));
		window->weight = (double *) malloc(sizeof(double));
		window->n_elements = 1;
		window->Z[0] = WindowMaterial;
		window->weight[0] = 1.0;
		window->density = WindowDensity;
		window->thickness = WindowThickness;
	}

	XmiMsimTubeEbel xmi_msim_tube_ebel;
        if (!g_module_symbol(xrmc_xmimsim, "xmi_msim_tube_ebel", (gpointer *) &xmi_msim_tube_ebel)) {
                g_module_close(xrmc_xmimsim);                                      
                g_fprintf(stderr, "GModule error message: %s\n", g_module_error());
                throw xrmc_exception("Could not get symbol xmi_msim_tube_ebel from module xrmc-xmimsim\n");
        }
    
        if (xmi_msim_tube_ebel == NULL) {                                 
                g_module_close(xrmc_xmimsim);
                throw xrmc_exception("Symbol xmi_msim_tube_ebel from module xrmc-xmimsim is NULL\n");
        }
                                                                                   
        if (xmi_msim_tube_ebel(anode, window, filter, TubeVoltage,
                  TubeCurrent, ElectronAngle,
                  XrayAngle, IntervalWidth,
                  1.0, TransmissionFlag,
                  &exc) == 0)
                throw xrmc_exception("Error in xmi_msim_tube_ebel\n");

	//copy everything to spectrum members
	LineEne = new double[exc->n_discrete];
	LineSigma = new	 double[exc->n_discrete];
	LineIntensity[0] = new double[exc->n_discrete];
	LineIntensity[1] = new double[exc->n_discrete];
	EneLineNum = exc->n_discrete;
	int i;
	for (i = 0 ; i < exc->n_discrete ; i++) {
		LineEne[i] = exc->discrete[i].energy;
		LineSigma[i] = 0.0;
		LineIntensity[0][i] = exc->discrete[i].horizontal_intensity;
		LineIntensity[1][i] = exc->discrete[i].vertical_intensity;
	}
	
	ContinuousEne = new double[exc->n_continuous];
	ContSIntensity[0] =  new double[exc->n_continuous];
	ContSIntensity[1] =  new double[exc->n_continuous];
	EneContinuousNum = exc->n_continuous;
	for (i = 0 ; i < exc->n_continuous ; i++) {
		ContinuousEne[i] = exc->continuous[i].energy;
		ContSIntensity[0][i] = exc->continuous[i].horizontal_intensity;
		ContSIntensity[1][i] = exc->continuous[i].vertical_intensity;
	}

	for (i = 0 ; i < EneContinuousNum ; i++)
		cout << ContinuousEne[i] << "\t" << ContSIntensity[0][i] << endl;

	for (i = 0 ; i < EneLineNum ; i++)
		cout << LineEne[i] << "\t" << LineIntensity[0][i] << endl;




	if (ResampleFlag)
		Resample();	

	spectrum::RunInit();



}
