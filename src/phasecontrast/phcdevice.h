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
//       phcdevice.h             //
//        19/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// basesource and source classes definition
//
#ifndef PHCDEVICEH
#define PHCDEVICEH

#include <string>
#include "phcdevice.h"

///////////////////////////////////////////////
// phcdevice class definition
//
// generic (virtual) class for phase contrast device
// phcdevice is derived from the class device
///////////////////////////////////////////////
class phcdevice
{
 public:
  bool PhCFlag;

  // generate an event with a photon directed toward the position x1
  // virtual int Out_Phase_Photon_x1(photon *, vect3) {return 0;}
  virtual double GetPhC_E0()=0;
  virtual int PhCOn() {PhCFlag=true; return 0;}
  virtual int PhCOff() {PhCFlag=false; return 0;}
  //virtual phcdevice *Clone(string) {return NULL;};
};

#endif
