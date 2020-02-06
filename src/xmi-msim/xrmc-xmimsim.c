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
	if (xmi_xmlLoadCatalog(NULL) == 0) {
		return 0;
	}
	g_fprintf(stdout,"XMI-MSIM XML catalog loaded\n");
	return 1;
}


G_MODULE_EXPORT int xmi_msim_detector_convolute(double ***Image, double ***ConvolutedImage, xmi_layer *det_absorber, xmi_detector *xd, int ModeNum, int NBins, int NY, int NX) {


	xmi_main_options *options = xmi_main_options_new();
	xmi_escape_ratios *escape_ratios_def=NULL;
	xmi_input *input = NULL;
	char *xmi_input_string;
	char *xmimsim_hdf5_escape_ratios = NULL;
	xmi_inputFPtr inputFPtr;
	int i, j, k;
	double *channels_conv_temp;
	gchar *hdf5_file=NULL;
	GError *error = NULL;

	g_fprintf(stdout,"ModeNum: %i\n", ModeNum);
	g_fprintf(stdout,"NBins: %i\n", NBins);
	g_fprintf(stdout,"NY: %i\n", NY);
	g_fprintf(stdout,"NX: %i\n", NY);


	xmi_init_hdf5();

	options->use_M_lines = 1;
	options->use_cascade_auger = 1;
	options->use_cascade_radiative = 1;
	options->use_variance_reduction = 1;
	options->use_sum_peaks = 0;
	options->use_poisson = 0;
	options->verbose = 1;
	options->extra_verbose = 0;
	options->omp_num_threads = omp_get_max_threads();
	options->use_escape_peaks = 1;

	if (xd->pulse_width > 0.0)
		options->use_sum_peaks = 1;

	xd->nchannels = NBins;

	input = xmi_input_init_empty();
	input->general->outputfile = strdup("non-existent-file.xmso");

	//modify input
	//put in one layer into composition
	input->composition->n_layers = 1;
	input->composition->layers = (xmi_layer *) malloc(sizeof(xmi_layer));
	input->composition->layers[0].n_elements = 1;
	input->composition->layers[0].Z = (int *) malloc(sizeof(int));
	input->composition->layers[0].weight = (double *) malloc(sizeof(double));
	input->composition->layers[0].Z[0] = 14;
	input->composition->layers[0].weight[0] = 1.0;
	input->composition->layers[0].density = 1;
	input->composition->layers[0].thickness = 1;
	input->composition->reference_layer = 1;

	xmi_absorbers_free(input->absorbers);
	input->absorbers = (xmi_absorbers *) malloc(sizeof(xmi_absorbers));
	input->absorbers->n_exc_layers = 0;
	input->absorbers->exc_layers = NULL;
	if (det_absorber != NULL) {
		input->absorbers->det_layers = det_absorber;
		input->absorbers->n_det_layers = 1;
	}

	//detector
	xmi_detector_free(input->detector);
	input->detector = xd;


	
	xmi_input_print(input, stdout);

	//get full path to XMI-MSIM HDF5 data file
	if (xmi_get_hdf5_data_file(&hdf5_file) == 0) {
		return 0;
	}



	//read escape ratios
	if (xmi_get_escape_ratios_file(&xmimsim_hdf5_escape_ratios, 1) == 0)
		return 0;

	if (options->verbose)
		g_fprintf(stdout,"Querying %s for escape peak ratios\n",xmimsim_hdf5_escape_ratios);

	//check if escape ratios are already precalculated
	if (xmi_find_escape_ratios_match(xmimsim_hdf5_escape_ratios , input, &escape_ratios_def, options) == 0)
		return 0;
	if (escape_ratios_def == NULL) {
		if (options->verbose)
			g_fprintf(stdout,"Precalculating escape peak ratios\n");
		//doesn't exist yet
		//convert input to string
		if (!xmi_input_write_to_xml_string(input, &xmi_input_string, &error)) {
			g_fprintf(stderr, "xmi_write_input_xml_to_string error: %s", error->message);
			return 0;
		}
		xmi_escape_ratios_calculation(input, &escape_ratios_def, xmi_input_string,hdf5_file,options, xmi_get_default_escape_ratios_options());
		//update hdf5 file
		if (xmi_update_escape_ratios_hdf5_file(xmimsim_hdf5_escape_ratios , escape_ratios_def) == 0)
			return 0;
		else if (options->verbose)
			g_fprintf(stdout,"%s was successfully updated with new escape peak ratios\n",xmimsim_hdf5_escape_ratios);
	}
	else if (options->verbose)
		g_fprintf(stdout,"Escape peak ratios already present in %s\n",xmimsim_hdf5_escape_ratios);
	xmi_input_C2F(input,&inputFPtr);

	//initialization
	if (xmi_init_input(&inputFPtr) == 0) {
		return 0;
	}

	double *abscorrImage = (double *) malloc(sizeof(double)*NBins);
	int ix, iy;

	for (iy = 0 ; iy < NY ; iy++) {
	  for (ix = 0 ; ix < NX ; ix++) {
	    for (k = 0 ; k < NBins ; k++) {
	      abscorrImage[k] = 0.0;
	    }
	    for (i = 0 ; i < ModeNum ; i++) {
	      for (k = 0 ; k < NBins ; k++) {
		abscorrImage[k] += Image[i*NBins+k][iy][ix];
	      }
	      xmi_detector_convolute_spectrum(inputFPtr, abscorrImage, &channels_conv_temp, options, escape_ratios_def, i);
	      for (k = 0 ; k < NBins ; k++) {
	        ConvolutedImage[i*NBins+k][iy][ix] = channels_conv_temp[k];
	      }
	      g_free(channels_conv_temp);
	    }
	  }
	}

	for (iy = 0 ; iy < NY ; iy++) {
	  for (ix = 0 ; ix < NX ; ix++) {
	    for (i = ModeNum-1 ; i > 0 ; i--) {
	      for (k = 0 ; k < NBins ; k++) {
	        ConvolutedImage[i*NBins+k][iy][ix] = MAX(ConvolutedImage[i*NBins+k][iy][ix]-ConvolutedImage[(i-1)*NBins+k][iy][ix], 0.0);
	      }
	    }
	  }
	}


	free(abscorrImage);
	xmi_free_escape_ratios(escape_ratios_def);
	xmi_main_options_free(options);
	
	//these next two lines need to see fixes first in XMI-MSIM
	//before they can be used again
	//xmi_free_input_F(&inputFPtr);
	//xmi_free_input(input);
	//if (xmi_end_random_acquisition() == 0) {
	//  return 0;
	//}
	return 1;
}

G_MODULE_EXPORT int xmi_msim_tube_ebel(xmi_layer *tube_anode, xmi_layer *tube_window,
                  xmi_layer *tube_filter, double tube_voltage,
                  double tube_current, double tube_angle_electron,
                  double tube_angle_xray, double tube_delta_energy,
                  double tube_solid_angle, int tube_transmission,
		  size_t tube_nefficiencies, double *tube_energies, double *tube_efficiencies,
                  xmi_excitation **ebel_spectrum
                  ) {
	return xmi_tube_ebel(tube_anode, tube_window,
				tube_filter, tube_voltage,
				tube_current, tube_angle_electron,
				tube_angle_xray, tube_delta_energy,
				tube_solid_angle, tube_transmission,
				tube_nefficiencies, tube_energies, tube_efficiencies,
				ebel_spectrum);
}

G_MODULE_EXPORT gboolean xmi_msim_transmission_efficiency_read(const char *filename, size_t *nefficiencies, double **energies, double **efficiencies, GError **error) {
	return xmi_transmission_efficiency_read(
			filename,
			nefficiencies,
			energies,
			efficiencies,
			error
			);
}
