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
//////////////////////////////////////////////////
//               loadparams.cpp                 //
//                06/02/2013                    //
//         author: Bruno Golosio                //
//////////////////////////////////////////////////
// xrmc method for loading simulation parameters
//

#include <time.h>
#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_algo.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"

using namespace std;
using namespace gettoken;
using namespace xrmc_algo;

//////////////////////////////////////////////////////////////////////
// method for loading simulation parameters
// fp is a pointer to the main input file
//////////////////////////////////////////////////////////////////////
int xrmc::LoadParams(FILE *fp)
{
  FILE *par_fp;
  char file_name[MAXSTRLEN], comm_ch[MAXSTRLEN];
  string comm="";
  time_t t1;
  long seed;

  GetToken(fp, file_name); // read the name of the parameter file
  cout << "Parameters file: " << file_name << "\n";
  if ((par_fp = fopen(file_name,"r")) == NULL)
    throw xrmc_exception("Parameters file not found.\n");

  GetToken(par_fp, comm_ch); // get a command/variable name from input file
  comm = comm_ch;
  if(comm!="Seed")
    throw xrmc_exception("syntax error in parameter file. "
			 "Seed command expected"); 
  else {
    GetLongToken(par_fp, &seed); // read the starting seed for random numbers
    if (seed == 0) { // if it is zero use the time from system clock as seed
      (void) time(&t1);
      seed=(long)t1; // convert t1 to long type
    } 
    cout << "Random seed: " <<  seed << "\n";
    init_randmt(seed);
  }

  fclose(par_fp);

  return 0;
}
