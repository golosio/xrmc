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
//     device.cpp              //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class device
//
#include <iostream>
#include <string>
#include <vector>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_device.h"

using namespace std;
using namespace gettoken;

// initialization of the class device
int xrmc_device::SetDevice(string dev_name, string dev_type)
{
  Name = dev_name; // device name
  DeviceType = dev_type; // device type
  if (NInputDevices>0) {
    InputDevice = vector<xrmc_device*>(NInputDevices,NULL);
    InputDeviceName = vector<string>(NInputDevices,"");
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for linking input devices
/////////////////////////////////////////////////////////////////////
int xrmc_device::LinkInputDevices(xrmc_device_map *dev_map)
{
  // loop on input devices
  for (int id=0; id<NInputDevices; id++) {
    // check if the input device is defined
    xrmc_device_map::iterator it = dev_map->find(InputDeviceName[id]);

  // if not display error and exit
    if (it==dev_map->end())
      throw xrmc_exception(string("Device ") + InputDeviceName[id]
			 + " not found in device map\n");

    // get device pointer from the device map
    InputDevice[id] = (*it).second;
  }
  CastInputDevices();

  return 0;
}

//////////////////////////////////////////////////////////////////////
// recursive link method
//////////////////////////////////////////////////////////////////////
int xrmc_device::RecursiveLink(xrmc_device_map *dev_map)
{
  LinkInputDevices(dev_map);

  // loop on input devices
  for (int id=0; id<NInputDevices; id++) {
    InputDevice[id]->RecursiveLink(dev_map);
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// recursive run initialization method
//////////////////////////////////////////////////////////////////////
int xrmc_device::RecursiveRunInit()
{
  // loop on input devices
  for (int id=0; id<NInputDevices; id++) {
    InputDevice[id]->RecursiveRunInit();
  }
  RunInit();

  return 0;
}

//////////////////////////////////////////////////////////////////////
// recursive cleaning after run method
//////////////////////////////////////////////////////////////////////
/*
int xrmc_device::RecursiveRunFree()
{
  RunFree();
  // loop on input devices
  for (int id=0; id<NInputDevices; id++) {
    InputDevice[id]->RecursiveRunFree();
  }

  return 0;
}
*/

//////////////////////////////////////////////////////////////////////
// parse device file for input device commands
//////////////////////////////////////////////////////////////////////
bool xrmc_device::ParseInputDeviceCommand(FILE *fp, string comm_str)
{
  char s[MAXSTRLEN];

  for(int i=0; i<NInputDevices; i++) {
    if (comm_str==InputDeviceCommand[i]) { // set the input device name
      GetToken(fp, s);
      InputDeviceName[i] = s;
      cout << InputDeviceDescription[i] << ": " << InputDeviceName[i] << "\n";
      return true;
    } 
  }

  return false;
}

