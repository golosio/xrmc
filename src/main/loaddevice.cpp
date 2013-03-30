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

#include <config.h>
#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_exception.h"
#include "xrmc_device.h"
#include "xrmc_gettoken.h"
				    /*
#include "xrmc_spectrum.h"
#include "xrmc_source.h"
#include "xrmc_composition.h"
#include "xrmc_geom3d.h"
#include "xrmc_sample.h"
#include "xrmc_detector.h"
#ifdef HAVE_XMIMSIM
#include "xrmc_detectorconvolute.h"
#endif
				    */
using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
// xrmc method for loading device
// fp is a opinter to the main input file
//////////////////////////////////////////////////////////////////////
int xrmc::LoadDevice(FILE *fp)
{
  char file_name[MAXSTRLEN], comm[MAXSTRLEN], dev_name[MAXSTRLEN];
  string comm_str;
  FILE *dev_fp;
  xrmc_device_map_insert_pair insert_pair;
  xrmc_device *dev_pt;

  GetToken(fp, file_name); // read from fp the name of the device input file
  cout << "Device file: " << file_name << "\n";
  if ((dev_fp = fopen(file_name,"r")) == NULL)
    throw xrmc_exception("Device file can not be opened.");
  GetToken(dev_fp, comm); // read the command
  comm_str = comm;
  if (comm_str=="Newdevice") {
    if (xrmc_device::LoadNewDevice(dev_fp, dev_pt)==0) {
      // insert name and pointer to the new device in the device map
      insert_pair = DeviceMap.insert(xrmc_device_map_pair(dev_pt->Name,
							  dev_pt));
      if(insert_pair.second == false) // check that it was not already inserted
	throw xrmc_exception(string("Device ") + dev_pt->Name + 
			     " already inserted in device map\n");
    }
  }
  else if (comm_str=="Device") {
    GetToken(dev_fp, dev_name); // read the device name
    cout << "Device name: " << dev_name << "\n";
    xrmc_device_map::iterator it = DeviceMap.find(dev_name);
     // check that the device was defined in device map
    if (it==DeviceMap.end()) {
      throw xrmc_exception(string("Device ") + dev_name + 
			   " not found in device map!");
    }
    dev_pt = (*it).second;
    dev_pt->Load(dev_fp); // load device parameters from device file
  }
  else // unknown command
    throw xrmc_exception(string("Unrecognized command ") + comm_str +
			 " in device file.\n");
  
  return 0;
}

