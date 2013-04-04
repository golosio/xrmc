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

//#define UNREF( x ) (void)x

#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "xrmc_math.h"
#include "xrmc_exception.h"

class xrmc_device;

typedef std::map<std::string, xrmc_device*> xrmc_device_map;
typedef std::pair<std::string, xrmc_device*> xrmc_device_map_pair;
typedef std::pair<xrmc_device_map::iterator, bool> xrmc_device_map_insert_pair;

///////////////////////////////////////////////
// device class definition
///////////////////////////////////////////////
class xrmc_device
{
 public:
  std::string Name; // device name
  std::string DeviceType; // device type
  bool Runnable;
  std::vector<std::string> SaveDataName;
  int NInputDevices;
  std::vector<xrmc_device*> InputDevice;
  std::vector<std::string> InputDeviceName;
  std::vector<std::string> InputDeviceDescription;
  std::vector<std::string> InputDeviceCommand;
  virtual ~xrmc_device() {};
  // initialization
  int SetDevice(std::string dev_name, std::string dev_type);
  // method for linking input devices
  int LinkInputDevices(xrmc_device_map *dev_map);
  // recursive link method
  int RecursiveLink(xrmc_device_map *dev_map);
  // recursive run initialization method
  int RecursiveRunInit();
  // recursive cleaning after run method
  //int RecursiveRunFree();
  //method for loading a new device
  static int LoadNewDevice(istream &dev_fs, xrmc_device*& dev_pt);
  //////////////////////////////////////////////////////////////////////
  //virtual methods of the class 
  //////////////////////////////////////////////////////////////////////
  // method for casting input devices to derived classes
  virtual int CastInputDevices() {return 0;}
  // run initialization method
  virtual int RunInit() {return 0;}
  // cleaning after run method
  //virtual int RunFree();
  //virtual int Load(FILE*)=0;// {return 0;} // load device from input file
  virtual int Load(std::istream&) {return 0;} // load device from input file
  virtual int SetDefault() {return 0;} // set default values for device params
  // save device output to a file
  virtual int SaveData(string, string) {return 0;}
  virtual int Begin() {LoopIdx = 0; return 0;} // start loop on events
  virtual int Next() {LoopIdx = 1; return 0;}  // next step of the loop
  virtual bool End() {return (LoopIdx==1);}    // end of the loop 
  virtual int Run() {return 0;} // generic virtual run method
  virtual int EventMulti() {return 1;} // event multiplicity
  //virtual int LinkInputDevice(string command, xrmc_device *dev_pt)
  //{
  //  throw xrmc_exception(string("Device of type ") + DeviceType + 
  //			 "does not have input devices\n");
  //}
  virtual bool ParseInputDeviceCommand(std::istream &fs, std::string str);
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
