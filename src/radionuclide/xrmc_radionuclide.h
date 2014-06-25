/*
Copyright (C) 2014 Bruno Golosio and Tom Schoonjans

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
/////////////////////////////////////////
//          xrmc_radionuclide.h        //
//            24/06/2014               //
//     Author : Tom Schoonjans         //
/////////////////////////////////////////
//

#ifndef XRMC_RADIONUCLIDE_H
#define XRMC_RADIONUCLIDE_H
#include "xrmc_spectrum.h"
#include <xraylib.h>


class radionuclide : public spectrum 
{
 public:
  int Unit; //1: mCi, 2: Ci, 3: GBq, 4: Bq
  double Activity;
  std::string RadioNuclide; 

  radionuclide(std::string dev_name);
  virtual int RunInit(); // initialize the spectrum before run
  virtual int SetDefault(); // set default values for spectrum parameters
  virtual int Load(std::istream &fs); // load spectrum parameters from file
  virtual ~radionuclide();

 private:
  struct radioNuclideData *rnd;


};

#endif
