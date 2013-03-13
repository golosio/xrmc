/*
Copyright (C) 2013 Antonio Brunetti

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
//---------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "detector_routines.h"
#include "xraylib.h"

//---------------------------------------------------------------------------

int detector_convolution_SDD(double *spectrum_eff, char *spectrum_out,
			     double Emin, double Emax, int spectrum_size)
{
  FILE *file_out;
  int i,j;
  double sigma_var;
  double *conv_nucleus;
  double *spectrum_conv;
  double factor;
  double E;
  int conv_lenght=60;
  conv_nucleus=(double *)calloc(conv_lenght,sizeof(double));
  spectrum_conv=(double *)calloc(spectrum_size,sizeof(double));

  for (i=0;i<spectrum_size;i++) spectrum_conv[i]=0.0;

  for (i=0;i<spectrum_size;i++) {
    sigma_var=sqrt(125.0*125.0-120*120.0+2440.0*i*0.03)/2.3548/30;
    factor=1.0/sigma_var/sqrt(2.0*M_PI);
    for (j=0;j<conv_lenght;j++)
      conv_nucleus[j]=factor*exp(-1.0*j*j/(2.0*sigma_var*sigma_var));

    for (j=-conv_lenght;j<conv_lenght;j++)
      if ( ((i+j)>=0) && ((i+j)<spectrum_size))
	spectrum_conv[i+j]+=spectrum_eff[i]*conv_nucleus[abs(j)];
  }

  file_out=fopen(spectrum_out,"w");
  for (i=0;i<spectrum_size;i++) {
    E = Emin +(Emax-Emin)*i/spectrum_size;
    fprintf(file_out, "%e\t%e\n", E, spectrum_conv[i]);
  }
  fclose(file_out);
  free(conv_nucleus);
  free (spectrum_conv);

  return 0;
}


double *detector_efficiency(char *spectrum_in, int Z,double thickness,
			    int Z_window,double window_thick,
			    double Emin, double Emax, int spectrum_size)
{
  FILE *file_in;
  int i;
  double *spectrum_temp,*spectrum_eff;
  double density,density1,density_air=0.0012,attenuation,
    attenuation_window,attenuation_air,energy,
    air_thickness;
  
  XRayInit();
  
  air_thickness=0;
  spectrum_temp=(double *)calloc(spectrum_size,sizeof(double));
  spectrum_eff=(double *)calloc(spectrum_size,sizeof(double));
  
  for (i=0;i<spectrum_size;i++) spectrum_eff[i]=0.0;
  file_in=fopen(spectrum_in,"rb");  
  fseek(file_in, sizeof(double)*spectrum_size, SEEK_SET);  
  if (fread(spectrum_temp,sizeof(double),spectrum_size,file_in)
      !=spectrum_size) {
    printf("Error reading input file\n");
    exit(EXIT_FAILURE);
  }
  fclose(file_in);
  
  density=2.33;
  density1=1.848;
  for (i=0;i<spectrum_size;i++) {
    energy=Emin +(Emax-Emin)*i/spectrum_size;
    if (energy>0.0) {
      attenuation=CS_Total(Z, energy);
      attenuation_window=CS_Total(Z_window, energy);
      attenuation_air=0.7804*CS_Total(7, energy)+0.20946*CS_Total(8, energy)+
	0.00934*CS_Total(18, energy);
      spectrum_eff[i]=spectrum_temp[i]
	*(1.0-exp(-attenuation*thickness*density))
	*exp(-attenuation_window*window_thick*density1
	     -attenuation_air*air_thickness*density_air);
      
    }
    else {
      spectrum_eff[i]=0.0;
    }   
  }

  free (spectrum_temp);
  
  return spectrum_eff;
}
