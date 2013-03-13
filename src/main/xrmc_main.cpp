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
//        xrmc_main.cpp          //
//        06/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Main File
// X-Ray Monte Carlo

#include <string>
#include <iostream>
#include "xrmc.h"
#include "xraylib.h"
#include "xrmc_exception.h"

using namespace std;

int ReadArg(int argc, char *argv[]);
int Banner();
/*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
/* Main Procedure. A single command line input input argument specifies
   the name of the input data file, which in turns contains the names
   of all files used by the simulation. */
{
  xrmc XRMC; // create object of class xrmc, main class for the simulation
  time_t time0, time1;
  double time_diff;
  clock_t cpu_time;
  double cpu_time_diff;
  
  try {
    Banner();
    XRayInit(); // xraylib libraries initialization
    ReadArg(argc, argv); // check command line arguments
    cout << argv[argc-1] << "\n";
    
    cout << "Start simulation...\n";
    time0 = time(NULL); // starting time
    cpu_time = clock(); // starting cpu time

    XRMC.Run(string(argv[argc-1])); // launch the simulation

    // cpu time difference converted to seconds
    cpu_time_diff = (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
    time1 = time(NULL); // final time
    time_diff = difftime(time1, time0); // time difference
    cout << "Elapsed time: " << time_diff << " sec\n";
    cout << "Elapsed CPU time: " << cpu_time_diff << " sec\n";
  }
  catch (xrmc_exception &e){ // handle possible runtime errors
    cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  catch (bad_alloc&) {
    cerr << "Error allocating memory." << "\n";
    return 1;
  }
  catch (...) {
    cerr << "Unrecognized error\n";
    return 1;
}

  return 0;
}

///////////////////////////////////////////////////////////////////////
int ReadArg(int argc, char *argv[])
  // Check command line input arguments
{
  if (argc != 2) // there should be two arguments
    // otherwise write usage message 
   throw xrmc_exception(string("Wrong input arguments.\nUsage : ")
			+ string(argv[0]) + " file-name\n");

 return 0;
}

// display banner
int Banner()
{
  cout << "xrmc (X-ray Monte Carlo)\n";
  cout << "Software for X-ray imaging and spectroscopy experiment simulation\n";
  cout << "based on Monte Carlo method with variance reduction techniques\n";
  cout << "Authors:\n";
  cout << "B. Golosio, A. Brunetti, G. L. Masala, P. Oliva, T. Schoonjans\n";
  cout << "Università degli Studi di Sassari\n";
  cout << "Please refer to the manual for usage conditions and instructions\n";
  cout << "-----------------------------------------------------------------\n";
  cout << endl;

  return 0;
}
    
