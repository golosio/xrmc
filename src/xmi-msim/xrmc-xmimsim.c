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

#include <xmi_msim.h>
#include <gmodule.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <xraylib.h>
#include <math.h>
#include <glib/gstdio.h>
#include <omp.h>

G_MODULE_EXPORT int xmi_check_xrmc_xmimsim_plugin(void) {
	//very simple function to check if the xrmc_xmimsim plugin works
	//and to start random number acquisition
	if (xmi_start_random_acquisition() == 0) {
		return 0;
	}
	//load xml catalog
	if (xmi_xmlLoadCatalog() == 0) {
		return 0;
	}
	g_fprintf(stdout,"XML catalog loaded\n");

	SetErrorMessages(0);
	return 1;
}


G_MODULE_EXPORT int xmi_msim_detector_convolute(double ***Image, double ***convolutedImage, struct xmi_layer *det_absorber, struct xmi_detector *xd, int ModeNum, int NBins, int NY, int NX) {


	struct xmi_main_options options;
	struct xmi_escape_ratios *escape_ratios_def=NULL;
	char *xmi_input_string;
	char *xmimsim_hdf5_escape_ratios = NULL;
	struct xmi_input *input;
	xmi_inputFPtr inputFPtr;
	int i, j, k;
	double *channels_conv_temp;
	gchar *hdf5_file=NULL;

	g_fprintf(stdout,"ModeNum: %i\n", ModeNum);
	g_fprintf(stdout,"NBins: %i\n", NBins);
	g_fprintf(stdout,"NY: %i\n", NY);
	g_fprintf(stdout,"NX: %i\n", NY);


	options.use_M_lines = 1;
	options.use_self_enhancement = 0;
	options.use_cascade_auger = 1;
	options.use_cascade_radiative = 1;
	options.use_variance_reduction = 1;
	options.use_optimizations = 1;
	options.use_sum_peaks = 0;
	options.use_poisson = 0;
	options.verbose = 1;
#if XMI_MSIM_VERSION_MAJOR >= && XMI_MSIM_VERSION_MINOR >= 1
	options.extra_verbose = 1;
	options.omp_num_threads = omp_get_max_threads();
#endif

	if (xd->pulse_width > 0.0)
		options.use_sum_peaks = 1;

	input = xmi_init_empty_input();

	//modify input
	//put in one layer into composition
	input->composition->n_layers = 1;
	input->composition->layers = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
	input->composition->layers[0].n_elements = 1;
	input->composition->layers[0].Z = (int *) malloc(sizeof(int));
	input->composition->layers[0].weight = (double *) malloc(sizeof(double));
	input->composition->layers[0].Z[0] = 14;
	input->composition->layers[0].weight[0] = 1.0;
	input->composition->layers[0].density = 1;
	input->composition->layers[0].thickness = 1;

	xmi_free_absorbers(input->absorbers);
	input->absorbers = (struct xmi_absorbers *) malloc(sizeof(struct xmi_absorbers));
	input->absorbers->n_exc_layers = 0;
	input->absorbers->exc_layers = NULL;
	if (det_absorber != NULL) {
		input->absorbers->det_layers = det_absorber;
		input->absorbers->n_det_layers = 1;
	}

	//detector
	if (input->detector->n_crystal_layers > 0) {
		for (i = 0 ; i < input->detector->n_crystal_layers ; i++) 
			xmi_free_layer(input->detector->crystal_layers+i);
		free(input->detector->crystal_layers);
	}

	free(input->detector);
	input->detector = xd;


	
	xmi_print_input(stdout, input);

	//get full path to XMI-MSIM HDF5 data file
	if (xmi_get_hdf5_data_file(&hdf5_file) == 0) {
		return 0;
	}



	//read escape ratios
	if (xmi_get_escape_ratios_file(&xmimsim_hdf5_escape_ratios, 1) == 0)
		return 0;

	if (options.verbose)
		g_fprintf(stdout,"Querying %s for escape peak ratios\n",xmimsim_hdf5_escape_ratios);

	//check if escape ratios are already precalculated
#if XMI_MSIM_VERSION_MAJOR >= && XMI_MSIM_VERSION_MINOR >= 1
	if (xmi_find_escape_ratios_match(xmimsim_hdf5_escape_ratios , input, &escape_ratios_def, options) == 0)
#else
	if (xmi_find_escape_ratios_match(xmimsim_hdf5_escape_ratios , input, &escape_ratios_def) == 0)
#endif
		return 0;
	if (escape_ratios_def == NULL) {
		if (options.verbose)
			g_fprintf(stdout,"Precalculating escape peak ratios\n");
		//doesn't exist yet
		//convert input to string
		if (xmi_write_input_xml_to_string(&xmi_input_string,input) == 0) {
			return 0;
		}
		xmi_escape_ratios_calculation(input, &escape_ratios_def, xmi_input_string,hdf5_file,options);
		//update hdf5 file
		if( xmi_update_escape_ratios_hdf5_file(xmimsim_hdf5_escape_ratios , escape_ratios_def) == 0)
			return 0;
		else if (options.verbose)
			g_fprintf(stdout,"%s was successfully updated with new escape peak ratios\n",xmimsim_hdf5_escape_ratios);
	}
	else if (options.verbose)
		g_fprintf(stdout,"Escape peak ratios already present in %s\n",xmimsim_hdf5_escape_ratios);
	xmi_input_C2F(input,&inputFPtr);

	//initialization
	if (xmi_init_input(&inputFPtr) == 0) {
		return 0;
	}


	//apply detector window absorption if necessary	
	double *blbs = (double * ) malloc(sizeof(double)*NBins);
	for (i = 0 ; i < NBins ; i++) 
		blbs[i] = 1.0;

	double sum;
	if (det_absorber != NULL) {
		for (i = 0 ; i < NBins ; i++) {
			sum = 0.0;
			for (j = 0 ; j < det_absorber[0].n_elements ; j++)
				sum += CS_Total_Kissel(det_absorber[0].Z[j], i*xd->gain)*det_absorber[0].weight[j];
			blbs[i] *= exp(-1.0 * det_absorber[0].density * det_absorber[0].thickness * sum);

		}
	}


	for (i = 0 ; i < NBins ; i++) { 
		sum = 0.0;
		for (j = 0 ; j < xd->crystal_layers[0].n_elements ; j++)
			sum += CS_Total_Kissel(xd->crystal_layers[0].Z[j], i*xd->gain)*xd->crystal_layers[0].weight[j];
		blbs[i] *= (1.0 - exp(-1.0*xd->crystal_layers[0].thickness*xd->crystal_layers[0].density*sum));
	}

	double *abscorrImage = (double *) malloc(sizeof(double)*NBins);

	for (i = 0 ; i < ModeNum ; i++) {
	  for (iy = 0 ; iy < NY ; iy++) {
	    for (ix = 0 ; ix < NX ; ix++) {
	      for (k = 0 ; k < NBins ; k++) 
		abscorrImage[k] = Image[i*NBins+k][iy][ix] * blbs[k];
	      xmi_detector_convolute(inputFPtr, abscorrImage, &channels_conv_temp, NBins, options, escape_ratios_def);
	      for (k = 0 ; k < NBins ; k++) 
		convolutedImage[i*NBins+k][iy][ix] = channels_conv_temp[k];
	      xmi_deallocate(channels_conv_temp);
	    }
	  }
	}

	free(abscorrImage);
	free(blbs);
	
	//xmi_free_input_F(&inputFPtr);
	//xmi_free_input(input);
	if (xmi_end_random_acquisition() == 0) {
	  return 0;
	}
	return 1;
}
