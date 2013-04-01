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
/////////////////////////////////////////
//          gettoken.h                 //
//     Author : Bruno Golosio          //
//            06/02/2013               //
/////////////////////////////////////////
// definition of namespace gettoken
// gets character arrays, integer values and real values from file
// skipping comments, blanks, ... 

namespace gettoken
{
/*---------------------------------------------------------------------------*/
  bool GetToken(istream &fs, string &str); // read string from stream
/*---------------------------------------------------------------------------*/
  bool GetIntToken(istream &fs, int *i); // read integer from stream
/*---------------------------------------------------------------------------*/
  bool GetLongToken(istream &fs, long *l); // read long integer from stream
/*---------------------------------------------------------------------------*/
  bool GetDoubleToken(istream &fs, double *d); // read double from stream
/*---------------------------------------------------------------------------*/
}

