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

#include <stdio.h>
#include <string.h>
#include <xrmc_exception.h>
#include <xrmc_gettoken.h>
/*---------------------------------------------------------------------------*/

//////////////////////////////////////////////////////////////////////
// get an array of characters (representing a string) from file fp 
//////////////////////////////////////////////////////////////////////
int gettoken::GetToken(FILE *fp, char *s)
{
  char ch=' ';

  // read characters from file fp
  // skip blanks, tab and newlines 
  while (!feof(fp) && strchr(" \n\r\t", (ch=(char)fgetc(fp))) != NULL);
  while (!feof(fp) && strchr("!#;", ch) != NULL) {
    // if ch represents the beginning of a comment
    // skip characters until the end of the line
    while (!feof(fp) && strchr("\n\r", (ch=(char)fgetc(fp))) == NULL);
    // and again skip blanks, tab and newlines 
    while (!feof(fp) && strchr(" \n\r\t", (ch=(char)fgetc(fp))) != NULL);     
  }
  // check for end of file 
  if (feof(fp) && strchr("!#; \n\r\t", ch) != NULL) return 1;
  // at this point ch should be the starting character of a string
  fseek(fp, -1, SEEK_CUR); // step back one character
  int i=fscanf(fp, "%s", s); // and read the string
  if (i != 1) { // string is empty
    printf(" ...\n");
    return 1;
  }
  return 0;
}  

/*---------------------------------------------------------------------------*/

// get integer value from file fp
int gettoken::GetIntToken(FILE *fp, int *i)
{
 char s[100];
 GetToken(fp, s); // get character array
 if (sscanf(s, "%d", i) != 1) // convert it to integer
   throw xrmc_exception("Cannot convert token to integer.\n");

 return 0;
}
/*---------------------------------------------------------------------------*/

// get float value from file fp
int gettoken::GetFloatToken(FILE *fp, float *f)
{
 char s[100];
 GetToken(fp, s); // get character array
 if (sscanf(s, "%f", f) != 1) // convert it to float
   throw xrmc_exception("Cannot convert token to float.\n"); 

 return 0;
}

/*---------------------------------------------------------------------------*/
// get double value from file fp
int gettoken::GetDoubleToken(FILE *fp, double *lf)
{
 char s[100];

 GetToken(fp, s); // get character array
 if (sscanf(s, "%lf", lf) != 1) // convert it to double
     throw xrmc_exception("Cannot convert token to double.\n"); 

 return 0;
}

/*---------------------------------------------------------------------------*/
// get long integer value from file fp
int gettoken::GetLongToken(FILE *fp, long *li)
{
 char s[100];
 GetToken(fp, s); // get character array
 if (sscanf(s, "%ld", li) != 1) // convert it to long integer
     throw xrmc_exception("Cannot convert token to long.\n");

 return 0;
}

/*---------------------------------------------------------------------------*/

