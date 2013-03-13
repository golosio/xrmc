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
//     loaddevice.cpp            //
//        06/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// xrmc method for loading device
//

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_exception.h"
#include "xrmc_device.h"
#include "xrmc_gettoken.h"
#include "xrmc_spectrum.h"
#include "xrmc_source.h"
#include "xrmc_composition.h"
#include "xrmc_geom3d.h"
#include "xrmc_sample.h"
#include "xrmc_detector.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
// xrmc method for loading device
// fp is a opinter to the main input file
//////////////////////////////////////////////////////////////////////
int xrmc::LoadDevice(FILE *fp)
{
  char file_name[MAXSTRLEN], dev_type[MAXSTRLEN], dev_name[MAXSTRLEN];
  string type_str;
  FILE *dev_fp;
  xrmc_device *dev_pt;
  xrmc_device_map_insert_pair insert_pair;

  GetToken(fp, file_name); // read from fp the name of the device input file
  cout << "Device file: " << file_name << "\n";
  if ((dev_fp = fopen(file_name,"r")) == NULL)
    throw xrmc_exception("Device file can not be opened.");
  GetToken(dev_fp, dev_type); // read the device type
  cout << "Device type: " << dev_type << "\n";
  GetToken(dev_fp, dev_name); // read the device name
  cout << "Device name: " << dev_name << "\n";
  type_str = dev_type;

  // compare the device type against known types
  // and creates the device with the appropriate derived class constructor
  if (type_str=="spectrum") {
    dev_pt = new spectrum(dev_name);
  }
  else if (type_str=="source") {
    dev_pt = new source(dev_name);
  }
  else if (type_str=="composition") {
    dev_pt = new composition(dev_name);
  }
  else if (type_str=="quadricarray") {
    dev_pt = new quadricarray(dev_name);
  }
  else if (type_str=="geom3d") {
    dev_pt = new geom3d(dev_name);
  }
  else if (type_str=="sample") {
    dev_pt = new sample(dev_name);
  }
  else if (type_str=="detectorarray") {
    dev_pt = new detectorarray(dev_name);
  }
  else // unknown device type
    throw xrmc_exception(string("Device type ") + type_str + " not known.\n");
  // insert name and pointer to the new device in the device map
  insert_pair = DeviceMap.insert(xrmc_device_map_pair(dev_name,
						       dev_pt));
  if(insert_pair.second == false) // check that it was not already inserted
    throw xrmc_exception(string("Device ") + dev_name + 
	      " already inserted in device map\n");
   
  dev_pt->SetDefault(); // set default values for device parameters
  dev_pt->Load(dev_fp); // load device parameters from device file
                        // and initialize the device
  fclose(dev_fp);
  
  return 0;
}
