/*
Copyright (C) 2013 Bruno Golosio

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
///////////////////////////////////
//     loadnewdevice.cpp         //
//        29/03/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class device
//
#include <config.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_device.h"
#include "xrmc_spectrum.h"
#include "xrmc_source.h"
#include "xrmc_composition.h"
#include "xrmc_geom3d.h"
#include "xrmc_sample.h"
#include "xrmc_detector.h"
#include "beamsource.h"
#include "beamscreen.h"
#include "anisotropicsource.h"
#include "intensityscreen.h"
#ifdef HAVE_XMIMSIM
#include "xrmc_detectorconvolute.h"
#endif

using namespace std;
using namespace gettoken;

int xrmc_device::LoadNewDevice(istream &dev_fs,  xrmc_device*& dev_pt)
{
  string dev_type, dev_name;
  xrmc_device_map_insert_pair insert_pair;

  GetToken(dev_fs, dev_type); // read the device type
  cout << "Device type: " << dev_type << "\n";
  GetToken(dev_fs, dev_name); // read the device name
  cout << "Device name: " << dev_name << "\n";

  // compare the device type against known types
  // and creates the device with the appropriate derived class constructor
  if (dev_type=="spectrum") {
    dev_pt = new spectrum(dev_name);
  }
  else if (dev_type=="source") {
    dev_pt = new source(dev_name);
  }
  else if (dev_type=="composition") {
    dev_pt = new composition(dev_name);
  }
  else if (dev_type=="quadricarray") {
    dev_pt = new quadricarray(dev_name);
  }
  else if (dev_type=="geom3d") {
    dev_pt = new geom3d(dev_name);
  }
  else if (dev_type=="sample") {
    dev_pt = new sample(dev_name);
  }
  else if (dev_type=="detectorarray") {
    dev_pt = new detectorarray(dev_name);
  }
  else if (dev_type=="beamsource") {
    dev_pt = new beamsource(dev_name);
  }
  else if (dev_type=="beamscreen") {
    dev_pt = new beamscreen(dev_name);
  }
  else if (dev_type=="anisotropicsource") {
    dev_pt = new anisotropicsource(dev_name);
  }
  else if (dev_type=="intensityscreen") {
    dev_pt = new intensityscreen(dev_name);
  }
  else if (dev_type=="detectorconvolute") {
#ifdef HAVE_XMIMSIM
    dev_pt = new detectorconvolute(dev_name);
#else
    throw xrmc_exception(string("Device type ") + dev_type + " not supported.\nRecompile XRMC with the XMI-MSIM plug-in\n");
#endif
  }
  else // unknown device type
    throw xrmc_exception(string("Device type ") + dev_type + " not known.\n");
   
  dev_pt->SetDefault(); // set default values for device parameters
  dev_pt->Load(dev_fs); // load device parameters from device file
                        // and initialize the device

  return 0;
}
