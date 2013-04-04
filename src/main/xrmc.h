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
//        xrmc.h                 //
//        06/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// xrmc class definition
//
#ifndef XRMCH
#define XRMCH

#include <string>
#include <fstream>
#include "xrmc_device.h"

using namespace std;

#define MAXSTRLEN 1000

//////////////////////////////////////////////////////////////////////
// xrmc class definition
// This is the main class of the program
//////////////////////////////////////////////////////////////////////
class xrmc
{
 public:
  ~xrmc() { DeleteDevices(); }// destructor
  int Run(string file_name); // method for running the simulation
 private:
  xrmc_device_map DeviceMap; // map of devices used by the simulation

  int LoadDevice(istream &fs); //method for loading a (new or existing) device
  int LinkDevices(); // method for linking all previously defined devices
  int RunDevice(istream &fs); // launch the run method on a device
  int SaveData(istream &fs); // launching the save method on a device
  int LoadParams(istream &fs); // load the main simulation parameters
  int DeleteDevices(); // delete all previously defined devices

};

#endif

