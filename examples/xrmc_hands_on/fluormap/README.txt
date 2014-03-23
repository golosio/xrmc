//////////////////////////////////////////////////////////////////////
// Fluorescence mapping example
//////////////////////////////////////////////////////////////////////

This example simulates a fluorescence mapping experiment
using an ideal detector.
The sample is an iron parallelepiped containing 4 cilinders made of a
Cu-Sn alloy with different relative concentrations.

//////////////////////////////////////////////////////////////////////
// Visualizing the sample
//////////////////////////////////////////////////////////////////////
This is just for visualization purpose. Since XRMC does not include a tool
for 3d visualization, a trick to visualize the sample is to simulate a
radiographic image.
Run:
xrmc input_image.dat
the output will be stored in the file image.dat
To visualize it using ImageJ select:
 File->Import->Raw
choose image.dat
then in the form set:
Image type: 64-bit Real
Width 400, Height 400, Little-endian byte order
The sample is sligthly rotated just for visualization purpose.

//////////////////////////////////////////////////////////////////////
// Executing the simulation
//////////////////////////////////////////////////////////////////////
The input file translate.dat performs a 20 micron translation of all quadrics.
If you have a look at the main input file input.dat, you can see that
after loading all devices, it runs the simulation, saves the output,
translates the sample, runs another simulation and so on, until the scan
is complete.
To execute the simulation, run:
 xrmc input.dat
The output files will be stored in the directory output

//////////////////////////////////////////////////////////////////////
// Visualizing the map of the Sn signal using ImageJ
//////////////////////////////////////////////////////////////////////
It's a bit tricky, but it is possible to visualize the map of the Sn signal
using ImageJ. Select:
File->Import->Raw
in the directory "output", select output_0.dat
then in the form set:
Image type: 64-bit Real
Width 800, Height 1, Little-endian byte order, open all files in folder
The first 400 pixels are not used as they refer to transmission.
The pixels from 400 to 799 are the 1st scattering order signal.
There are 400 energy bins from 0 to 40 keV.
To see the measured spectrum on one step of the scan, choose
Edit->Selection->Select all
Analyze->Plot profile
To see all spectra in an image, select the imported stack,
then on ImageJ
Images->Stack->Make montage
Columns 1
Rows 41
Scale factor 1
First slice 1
Last slice 41

The Ka line of Sn is at  25.192 keV, which corresponds to
the column n. 400+252=652
Choose a region around the Sn Ka line
Edit->Selection->Specify
Width 20
Height 41
x coord 640
y coord 0

then plot a vertical profile:
Edit->Options->Profile plot options
Check "vertical profile"
Analyze->plot profile

