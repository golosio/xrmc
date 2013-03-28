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
//     qvolume.cpp               //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class qvolume
//

#include <cmath>
#include <iostream>
#include "xrmc_geom3d.h"
using namespace std;

//////////////////////////////////////////////////////////////////////
// qvolume initialization method
//////////////////////////////////////////////////////////////////////
int qvolume::Init(string s_ph_in, string s_ph_out, int n_quadr)
{
  NQuadr = n_quadr; // set the number of quadrics
  PhaseInName = s_ph_in; // set the phase inside
  PhaseOutName = s_ph_out; // set the phase outside
  Quadr = new quadric*[NQuadr]; // allocate pointer-to-quadric array

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Find intersections of the straight line x0 + u*t
// with the quadrics delimiting the object
//////////////////////////////////////////////////////////////////////
int qvolume::Intersect(vect3 x0, vect3 u, double *t, int *iph0, int *iph1,
			int *n_inters)
{
  int n_inters_q=0;
  int enter_q=0;
  double t_q=0;
  double t01max=0, t10min=0;
  int t01flag=0, t10flag=0;

  for (int iq=0; iq<NQuadr; iq++) { // loop on quadrics delimiting the object
    n_inters_q = Quadr[iq]->NInters; // num. of intersections with quadric iq
    // cout << "iq, n_inters_q: " << iq << "\t" << n_inters_q <<"\n";
    if (n_inters_q==0){ // if the trajectory does not intersect the quadric 
      if (Quadr[iq]->Inside(x0)) { // check if x0 is inside the quadric
	continue;
      }
      else {
	t01flag = 0; // if not, then the line does not intersect the object
	break;
      }
    }
    for(int ii=0; ii<n_inters_q; ii++) { // loop on intersections with quadric
      t_q = Quadr[iq]->tInters[ii]; // intersection parametric coordinate t
      // cout << "t_q: " << t_q << "\n";
      enter_q = Quadr[iq]->Enter[ii]; // check if the quadric is crossed from
                                      // outside to inside or viceversa
      //cout << "enter_q: " << enter_q << "\n";
      //cout << "t_q: " << t_q << "\n";
      // t01flag: flag for at least one crossing from outside to inside
      // t01max: maximum value of t for this type of crossing
      if (enter_q==1 && (t01flag==0 || t_q>t01max)) {
	t01flag = 1;
	t01max = t_q;
      }
      // t10flag: flag for at least one crossing from inside to outside
      // t10min: minimum value of t for this type of crossing
      if (enter_q==0 && (t10flag==0 || t_q<t10min)) {
	t10flag = 1;
	t10min = t_q;
      }
    }
    //cout << "t01max: " << t01max << "\n";
    //cout << "t10min: " << t10min << "\n";
  }

  if (t01flag && t10flag && t10min>t01max) {
    // if this condition is satisfied then the line crosses the object
    // t01max and t10min are the parametric coordinates of the intersections 
    if (t01max > 0) { // if the 1st intersection is in the forward direction
      t[*n_inters] = t01max; // store it
      iph0[*n_inters] = iPhaseOut;
      iph1[*n_inters] = iPhaseIn;
      (*n_inters)++;
    }
    if (t10min > 0) { // if the 2nd intersection is in the forward direction
      t[*n_inters] = t10min; // store it
      iph0[*n_inters] = iPhaseIn;
      iph1[*n_inters] = iPhaseOut;
      (*n_inters)++;
    }
  }
  //cout << "n_inters: " << *n_inters << "\n";

  return 0;
}


qvolume& qvolume::operator= (const qvolume &QVolume) {
	//cout << "Entering qvolume assignment operator\n";

	if (this == &QVolume)
		return *this;
	NQuadr = QVolume.NQuadr;
	iPhaseIn = QVolume.iPhaseIn;
	iPhaseOut = QVolume.iPhaseOut;
	PhaseInName = QVolume.PhaseInName;
	PhaseOutName = QVolume.PhaseOutName;
	Quadr = new quadric*[NQuadr];

	//cout << "Leaving qvolume assignment operator\n";

}
