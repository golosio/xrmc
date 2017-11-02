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
//     geom3d.cpp                //
//        31/10/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class geom3d
//

#include <string>
#include <iostream>
#include "xrmc_geom3d.h"
#include "xrmc_algo.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

bool geom3d::MapReduced = false;
// destructor
geom3d::~geom3d() {
  for (int i=0; i<NQVol; i++) {
    if (QVolMap[i]!=NULL) delete[] QVolMap[i];
  }
  if (QVolMap!=NULL) delete[] QVolMap;
  if (QVol!=NULL) delete[] QVol;
  if (QArr!=NULL) {
	delete QArr;
	QArr = NULL;
  }
  if (Comp!=NULL) {
	delete Comp;
	Comp = NULL;
  }
}

// constructor
geom3d::geom3d(string dev_name) {
  Runnable = false;
  NInputDevices = 2;
  InputDeviceCommand.push_back("QArrName");
  InputDeviceDescription.push_back("Quadric array input device name");
  InputDeviceCommand.push_back("CompName");
  InputDeviceDescription.push_back("Composition input device name");

  QVolMap = NULL;
  QVol = NULL;
  NQVol = 0;
  SetDevice(dev_name, "geom3d");
}

// method for casting input devices to types quadricarray and composition
int geom3d::CastInputDevices()
{
  // cast InputDevice[0] to type quadricarray*
  QArr = dynamic_cast<quadricarray*>(InputDevice[0]);
  if (QArr==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type quadricarray\n");

  // cast InputDevice[1] to type composition*
  Comp = dynamic_cast<composition*>(InputDevice[1]);
  if (Comp==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[1] +
              " cannot be casted to type composition\n");

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for linking input devices
/////////////////////////////////////////////////////////////////////
/*
int geom3d::LinkInputDevice(string command, xrmc_device *dev_pt)
{
  if (command=="QArrName") {
    // cast it to type quadricarray*
    QArr = dynamic_cast<quadricarray*>(dev_pt);
    if (QArr==0)
      throw xrmc_exception(string("Device cannot be casted to type "
				  "quadricarray\n"));
  }
  else if (command=="CompName") {
    // cast it to type composition*
    Comp = dynamic_cast<composition*>(dev_pt);
    if (Comp==0)
      throw xrmc_exception(string("Device cannot be casted to type "
				  "composition\n"));
  }
  else
    throw xrmc_exception(string("Unrecognized command: ") + command + "\n");

  return 0;
}
*/

// method for initializing the 3d geometry before a run
int geom3d::RunInit()
{

  //clean up the materials
  if (!MapReduced) {
    Comp->ReduceMaterMap(used_mater);
    MapReduced = true;
  }

  for(int iqv=0; iqv<NQVol; iqv++) { // loop on 3d objects
    for(int iq=0; iq<QVol[iqv].NQuadr; iq++) { // loop on quadrics delimiting
                                               // the 3d object
      // initialize the pointers of the quadrics delimiting the 3d object
      // using the 3d object quadric map
      string qname = QVolMap[iqv][iq];
      // check if the quadric name is present in quadric map
      quadric_map::iterator it = QArr->QuadricMap.find(qname);
      // if not display error and exit
      if (it==QArr->QuadricMap.end())
	throw xrmc_exception(string("Quadric ") + qname
			     + " not found in quadric map\n");
      // get quadric pointer from the quadric map
      quadric *quadr_pt = (*it).second;

      QVol[iqv].Quadr[iq] = quadr_pt;

    }
    string s_mat_in = QVol[iqv].MaterInName;

    // check if the material name is present in material map
    phase_map::iterator it = Comp->MaterMap.find(s_mat_in);
    // if not display error and exit
    if (it==Comp->MaterMap.end())
      throw xrmc_exception(string("Material ") + s_mat_in
			   + " not found in material map\n");
    // get material index from the material map
    int i_mat_in = (*it).second;
    QVol[iqv].iMaterIn = i_mat_in; // set the material inside

    string s_mat_out = QVol[iqv].MaterOutName;
    // check if the material name is present in material map
    it = Comp->MaterMap.find(s_mat_out);
    // if not display error and exit
    if (it==Comp->MaterMap.end())
      throw xrmc_exception(string("Material ") + s_mat_out
			   + " not found in material map\n");
    // get material index from material map
    int i_mat_out = (*it).second;
    QVol[iqv].iMaterOut = i_mat_out; // set the material inside
  }
 
  return 0;
}
//////////////////////////////////////////////////////////////////////
// method for evaluating the intersections of a straight trajectory
// with the 3d objects
//////////////////////////////////////////////////////////////////////
// input:
// x0 and u are the starting coordinates and direction of the trajectory
// output:
// t is the array of the intersections (distances from the starting coordinate)
// imat0 is the array of the entrance material indexes
// imat1 is the array of the exit material indexes
// n_inters is the number of intersections
 //////////////////////////////////////////////////////////////////////
int geom3d::Intersect(vect3 x0, vect3 u, double *t, int *imat0, int *imat1,
		       int *n_inters)
{
  int i;
  for (i=0; i<QArr->NQuadr; i++) { // loop on all quadrics
    QArr->Quadr[i].Intersect(x0, u); // evaluate intersection of the trajectory
                                     // with the quadric
  }
  *n_inters = 0;
  for (i=0; i<NQVol; i++) { // loop on 3d objects
 // evaluate intersection of the trajectory with the 3d objects
    QVol[i].Intersect(x0, u, t, imat0, imat1, n_inters);
  }
  SortInters(t, imat0, imat1, *n_inters); // sort the intersections
                                        // with growing values of t

  return 0;
}

geom3d *geom3d::Clone(string dev_name) {
	//cout << "Entering geom3d::Clone\n";
	geom3d *clone = new geom3d(dev_name);

	*clone = *this;
	clone->HW[0] = HW[0];
	clone->HW[1] = HW[1];
	clone->HW[2] = HW[2];
	clone->QArr = QArr->Clone(QArrName);
	clone->Comp = Comp->Clone(CompName);

	clone->QVol = new qvolume[NQVol];
	clone ->QVolMap = new string*[NQVol];	
	for (int i = 0 ; i < NQVol ; i++) {
	  clone->QVol[i] = QVol[i];
	  clone->QVolMap[i] = new string[QVol[i].NQuadr];
	  for (int j = 0 ; j < QVol[i].NQuadr ; j++)
	    clone->QVolMap[i][j] = QVolMap[i][j];
	}
	clone->RunInit();
	
	return clone;
}
