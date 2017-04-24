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
//          gettoken.cpp               //
//     Author : Bruno Golosio          //
//           06/02/2013                //
/////////////////////////////////////////
// Functions of the namespace gettoken
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "xrmc_exception.h"
#include "xrmc_gettoken.h"

//////////////////////////////////////////////////////////////////////
// Reads string from stream
//////////////////////////////////////////////////////////////////////
bool gettoken::GetToken(istream &fs, string &str)
{
  enum {FREE, IN_TOKEN, IN_TOKEN_BS, IN_COMMENT, DONE};

  char c, state=FREE;
  string start_comment = ";#";
  string end_line = "\n\r";
  str="";
  while (state!=DONE && fs.good()) {
    fs.get(c);
    switch (state)
      {
      case FREE:
	if (start_comment.find(c)!=string::npos) {
	  state = IN_COMMENT;
	}
	else if (c=='\\') {
	  state = IN_TOKEN_BS;
	}
	else if (isgraph(c)) {
	  str = str + c;
	  state = IN_TOKEN;
	}
	break;
      case IN_TOKEN:
	if (!isgraph(c) || start_comment.find(c)!=string::npos) {
	  state = DONE;
	}
	else if (c=='\\') {
	  state = IN_TOKEN_BS;
	}
	else {
	  str = str + c;
	}
	break;
      case IN_TOKEN_BS:
	if (isprint(c)) {
	  str = str + c;
	  state = IN_TOKEN;
	}
	else
	  state = DONE;
	break;
      case IN_COMMENT:
	if (end_line.find(c)!=string::npos) {
	  state = FREE;
	}
	//default:
      }
  }
  if(str!="") return true;
  else return false;
}

//////////////////////////////////////////////////////////////////////
// Reads long integer value from file stream
//////////////////////////////////////////////////////////////////////
bool gettoken::GetLongToken(istream &fs, long *l)
{
 string s;
 GetToken(fs, s); // get string token
 istringstream buffer(s);
 if (!(buffer >> (*l))) { // convert it to long integer
   std::ostringstream msg;
   msg << "Cannot convert token " << buffer.str() << " to long" << std::endl;
   throw xrmc_exception(msg.str());
 }

 return true;
}

//////////////////////////////////////////////////////////////////////
// Reads integer value from file stream
//////////////////////////////////////////////////////////////////////
bool gettoken::GetIntToken(istream &fs, int *i)
{
 string s;
 GetToken(fs, s); // get string token
 istringstream buffer(s);
 if (!(buffer >> (*i))) { // convert it to integer
   std::ostringstream msg;
   msg << "Cannot convert token " << buffer.str() << " to int" << std::endl;
   throw xrmc_exception(msg.str());
 }

 return true;
}

//////////////////////////////////////////////////////////////////////
// Reads double value from file stream
//////////////////////////////////////////////////////////////////////
bool gettoken::GetDoubleToken(istream &fs, double *d)
{
 string s;
 GetToken(fs, s); // get string token
 istringstream buffer(s);
 if (!(buffer >> (*d))) { // convert it to double
   std::ostringstream msg;
   msg << "Cannot convert token " << buffer.str() << " to double" << std::endl;
   throw xrmc_exception(msg.str());
 }

 return true;
}

