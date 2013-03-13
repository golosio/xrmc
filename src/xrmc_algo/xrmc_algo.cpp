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
//     xrmc_algo.cpp             //
//        07/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
//  xrmc generic algorithms
//  functions of the namespace xrmc_algo
//

#include <cmath>
#include <iostream>
#include <stdio.h>
#include "xrmc_algo.h"
#include "xrmc_exception.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Locate
// find the position of a real (double) value x
// in a sorted array x_arr[N] using the binary search algorithm
//
//////////////////////////////////////////////////////////////////////
int xrmc_algo::Locate(double x, double x_arr[], int N, int *idx)
{

  if (x<x_arr[0]) { // x is smaller than lower limit
    *idx=0;
    return -1;
  }
  if (x>=x_arr[N-1]) { // x is greater or equal to upper limit
    *idx=N-1;
    return 1;
  }
  int id=0, iu=N-1; // lower and upper index of the subarray to search
  while (iu-id>1) { // search until the size of the subarray to search is >1
    int im = (id + iu)/2; // use the midpoint for equal partition
    // decide which subarray to search
    if (x>=x_arr[im]) id=im; // change min index to search upper subarray
    else iu=im; // change max index to search lower subarray
  }
  *idx=id;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// Generate a pseudo-random number
// according to a Gaussian probability distribution
// with standsrd deviation sigma=1
// Once every two calls, this function generates two values, X,Y
// according to a 2d Gaussian distribution, and returns the value of X
// At odd calls, it returns the value of Y generated in the previous call
//////////////////////////////////////////////////////////////////////
double xrmc_algo::GaussRnd()
{
  double x, y, r2, inv_cumul, X;
  static int PrevFlag = 0; //1 if there is a ready value of Y from previous call
  static double Y; // previously generated value of Y
  
  if  (PrevFlag==1) { // there is a ready value of Y from previous call
    PrevFlag = 0;
    return Y;
  }
  
  // generates a uniform 2d distribution in a circle
  do {
    // first generate a uniform distribution in a 2 x 2 square
    x = 2.0*Rnd()-1.0; // x = r*cos(theta) in polar coordinates
    y = 2.0*Rnd()-1.0; // y = r*sin(theta)
    r2 = x*x+y*y; // r^2: square of polar coordinate r
  } while (r2>=1.0 || r2==0.0); //repeat if x,y are not in the circle with r=1
  // inverse cumulative distr. function for R/r, with R=sqrt(X*X+Y*Y)
  // i.e. polar coordinate R of (X,Y)
  inv_cumul = sqrt(-2.0*log(r2)/r2);
  X = inv_cumul*x; // X = R*cos(theta) = (R/r)*x
  Y = inv_cumul*y; // Y = R*sin(theta) = (R/r)*y
  PrevFlag = 1;

  return X;
}

//////////////////////////////////////////////////////////////////////
// Computes the integral of the real function func from a to b
// using a sum on equally spaced samples.
// it's used only once, therefore we don't need a more efficient algo 
//////////////////////////////////////////////////////////////////////
double xrmc_algo::Integrate(double (*func)(double), double a, double b)
{
  const int NSAMPLES=1000000; // n. of samples used to evaluate the integral
  double dx=(b-a)/NSAMPLES; // sample step
  double integral=0;
  for (int i=0; i<NSAMPLES; i++) { // loop on sampling points
    double x = dx*i +dx/2;
    integral += func(x); // sum function values on sampl. pts
  }
  integral *= dx;

  return integral;
}

//////////////////////////////////////////////////////////////////////
// Sort intersection of trajectory with quadrics
// t is the unordered array of intersection distances
// iph0 is the entrance phase index at each intersection
// iph1 is the exit phase index at each intersection
// n_inters is the nymber of intersections
//////////////////////////////////////////////////////////////////////
int xrmc_algo::SortInters(double *t, int *iph0, int *iph1, int n_inters)
{
  int i, j;
  double tj;
  int iph0j, iph1j;
  
  for (j=1; j<n_inters; j++) { // loop starting from the 2nd intersection 
    tj = t[j]; // intersection distance from starting position
    iph0j = iph0[j]; // entrance phase index 
    iph1j = iph1[j]; // exit phase index 
    i = j - 1; // previous intersection index
    while (i>=0 && t[i]>tj) { // if previous intersection is more distant
      // switch the two intersections
      t[i+1] = t[i];
      iph0[i+1] = iph0[i];
      iph1[i+1] = iph1[i];
      i--; // compare the two previous intersections
    }
    t[i+1] = tj; // set with stored temporary values
    iph0[i+1] = iph0j;
    iph1[i+1] = iph1j;
  }
  return 0;
}
