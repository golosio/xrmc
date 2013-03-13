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
#include <string>
#include "xrmc_device.h"

using namespace std;

// constructor of the class device
xrmc_device::xrmc_device(string dev_name, string dev_type)
{
  Name = dev_name; // device name
  DeviceType = dev_type; // device type
}



