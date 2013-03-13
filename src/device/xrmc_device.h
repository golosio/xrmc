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
//     xrmc_device.h             //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// device class definition
//
#ifndef XRMC_DEVICE_H
#define XRMC_DEVICE_H

#include <string>
#include <map>
#include "xrmc_math.h"

using namespace std;

class xrmc_device;

typedef map<string, xrmc_device*> xrmc_device_map;
typedef pair<string, xrmc_device*> xrmc_device_map_pair;
typedef pair<xrmc_device_map::iterator, bool> xrmc_device_map_insert_pair;

///////////////////////////////////////////////
// device class definition
///////////////////////////////////////////////
class xrmc_device
{
 public:
  string Name; // device name
  string DeviceType; // device type
  xrmc_device() {}; // constructors
  xrmc_device(string dev_name, string dev_type);
  //virtual methods of the class 
  virtual ~xrmc_device() {}; // destructor
  virtual int ImportDevice(xrmc_device_map *dev_map) {return 0;}
  virtual int Load(FILE *fp) {return 0;} // load device from input file
  virtual int SetDefault() {return 0;} // set default values for device params
  virtual int Save(string file_name) {return 0;} // save device output to a file
  virtual int Begin() {LoopIdx = 0; return 0;} // start loop on events
  virtual int Next() {LoopIdx = 1; return 0;}  // next step of the loop
  virtual bool End() {return (LoopIdx==1);}    // end of the loop 
  virtual int Run() {return 0;} // generic virtual run method
  virtual int EventMulti() {return 1;} // event multiplicity
 protected:
  int LoopIdx; // index of event loop
};

///////////////////////////////////////////////
// bodydevice class definition
//
// bodydevice is derived from the class device
// besides the device member variables and functions,
// a bodydevise is characterized by a position
// and a local coordinate system
///////////////////////////////////////////////
class bodydevice : public xrmc_device
{
 public:
  vect3 X; // bodydevice position
  vect3 ui, uj, uk; // local coordinate system axis directions
};

#endif
