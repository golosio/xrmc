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
//        xrmc.cpp               //
//        06/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// xrmc class member functions definitions
//

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xrmc_exception.h"
#include "xrmc_device.h"
#include "xrmc_gettoken.h"
#ifdef HAVE_XMIMSIM
#include "xrmc_detectorconvolute.h"
#endif

using namespace std;
using namespace gettoken;

/*---------------------------------------------------------------------------*/
// Run method
// file_name: main input file name
// 
int xrmc::Run(string file_name)
{
  FILE *fp;
  char s[MAXSTRLEN];
  string command;
  
  cout.precision(7); // max number of digits on floating-point values cout
  cout << "Files used by the simulation: " <<  file_name << "\n";
  if ((fp = fopen(file_name.c_str(),"r")) == NULL)
    throw xrmc_exception("Input file not found.\n");

  while(1) { 
    GetToken(fp, s); // read a command from input file
    if (feof(fp)) break;
    command = s;
    // parse the command and decide what to do
    if (command=="LoadParams") { // read simulation parameter file
      LoadParams(fp);
    }
    else if (command=="Load") { // load a new device
      LoadDevice(fp);
    }
    else if (command=="Link") { // link all previously loaded devices
      LinkDevices();
    }
    else if (command=="Run") {// launch the Run method on a device
      RunDevice(fp);
    }
    else if (command=="Save") { // launch the Save method on a device
      SaveDevice(fp);
    }
    else if (command=="SaveUnconvoluted") {
#ifdef HAVE_XMIMSIM
      SaveUnconvolutedDevice(fp);
#else
      throw xrmc_exception("Command SaveUnconvoluted is not supported.\nRecompile XRMC with the XMI-MSIM plug-in\n");
#endif
    }
    else if (command!="End")
      throw xrmc_exception(string("Syntax error: ") + command
			   + "\ncommand not found.\n");
  }
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for linking all previously defined devices
//////////////////////////////////////////////////////////////////////
int xrmc::LinkDevices()
{
  xrmc_device *dev_pt;

  // loop on device map
  for (xrmc_device_map::iterator it=DeviceMap.begin(); it!=DeviceMap.end();
       it++) {
    dev_pt=it->second; // pointer to the device
    dev_pt->LinkInputDevices(&DeviceMap); //launch LinkInputDevice method on it
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// RunDevice method
// fp is the main input file
//////////////////////////////////////////////////////////////////////
int xrmc::RunDevice(FILE *fp)
{
  char dev_name[MAXSTRLEN];
  xrmc_device *dev_pt;

  GetToken(fp, dev_name); // device on which the Run method should be launched  
  cout << "Run device: " << dev_name << "\n";
  xrmc_device_map::iterator it = DeviceMap.find(dev_name);
   // check that the device was defined in device map
  if (it==DeviceMap.end())
    throw xrmc_exception("Device not found in device map.\n");

  dev_pt = (*it).second;
   // check that the run method is defined for this device
  if (dev_pt->Runnable==false)
    throw xrmc_exception("Run method is not defined for this device.\n");
  dev_pt->RecursiveLink(&DeviceMap); // recursive link
  dev_pt->RecursiveRunInit(); // recursive run initialization
  dev_pt->Run(); // launch the Run method on the device

  return 0;
}
//////////////////////////////////////////////////////////////////////
// SaveDevice method
// fp is the main input file
//////////////////////////////////////////////////////////////////////
int xrmc::SaveDevice(FILE *fp)
{
  char dev_name[MAXSTRLEN], file_name[MAXSTRLEN];
  xrmc_device *dev_pt;

  GetToken(fp, dev_name); // device on which the Save method should be launched
  GetToken(fp, file_name);
  cout << "Save device: " << dev_name << " to file " << file_name << "\n";
  xrmc_device_map::iterator it = DeviceMap.find(dev_name);
   // check that the device was defined in device map
  if (it==DeviceMap.end())
    throw xrmc_exception("Device not found in device map.\n");

  dev_pt = (*it).second;
  dev_pt->Save(file_name); // launch the Save method on the device

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for deleting all previously defined devices
//////////////////////////////////////////////////////////////////////
int xrmc::DeleteDevices()
{
  xrmc_device *dev_pt;

  // loop on device map
  for (xrmc_device_map::iterator it=DeviceMap.begin(); it!=DeviceMap.end();
       it++) {
    dev_pt=it->second; // pointer to the device
    delete dev_pt; // delete the device
  }

  return 0;
}

