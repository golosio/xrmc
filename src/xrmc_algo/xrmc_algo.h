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
//     xrmc_algo.h               //
//        07/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
//  xrmc generic algorithms
//  namespace xrmc_algo definition

#ifndef XRMCALGOH
#define XRMCALGOH
#include "randmt.h"

#define Rnd() (rand_unif())
#define Rnd_r(x) (rand_unif_r(x))
#define GaussRnd_r(x) (rand_normal_r(x))
#define GaussRnd() (rand_normal())

namespace xrmc_algo
{  
  // find position of a real value in a sorted array
  int Locate(double x, double x_arr[], int N, int *idx);
  //double GaussRnd(); // random number with Gaussian distribution
  // integral of the real function func from a to b
  double Integrate(double (*func)(double), double a, double b);
  int SortInters(double *t, int *iph0, int *iph1, int n_inters);
}

#endif
