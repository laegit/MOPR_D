/*
getQAP.cc (C) Joshua Knowles 12/3/2002

This program reads in a file written by makeQAPuni.cc or makeQAPrl.cc. This code
can be embedded into your favourite QAP search algorithm.

Compile: gcc getQAP.cc -o getQAP -lm

Usage: ./getQAP qap_problem_file


** Please contact me - Joshua Knowles - if you have any comments, suggestions
or questions regarding this program or multiobjective QAP problems. My email
address is jknowles@ulb.ac.be

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version. 

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details. 

   The GNU General Public License is available at:
      http://www.gnu.org/copyleft/gpl.html
   or by writing to: 
        The Free Software Foundation, Inc., 
        675 Mass Ave, Cambridge, MA 02139, USA.  


This program also uses Professor Peter Ross' (Napier University) getdouble() and getint() functions.

*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>


// true and false used in getdouble() and getint()
#define GDFALSE 0  
#define GDTRUE 1




double getdouble(FILE * file, double *valaddr, int stopateol);
int getint(FILE * f, int *valaddr, int stopateol);
long getlong(FILE * f, long *valaddr, int stopateol);


// The following file i/o functions are written by Professor Peter Ross 
// who is at Napier University, Edinburgh

/***************************************************************/
/* Get the next number from the input: put it in the location  */
/* addressed by second argument. This function returns 0 on    */
/* EOF. If stopateol is true, it returns -1 when it hits \n    */
/* (after which some other procedure has to read past the \n), */
/* otherwise it continues looking for the next number.         */
/* A number has an optional sign, perhaps followed by digits,  */
/* perhaps followed by a decimal point, perhaps followed by    */
/* more digits. There must be a digit somewhere for it to count*/
/* as a number. So it would read any of:                       */
/*  -.5                                                        */
/*  -0.5                                                       */
/*  -.5.7                                                      */
/* as minus-a-half. In the last case, it would read .7 next    */
/* time around.                                                */
/*   There doesn't seem to be a neat and reliable way to do    */
/* all this, including stopateol, using scanf?                 */
/***************************************************************/

double getdouble(FILE * file, double *valaddr, int stopateol)

{
  int c;
  int found = GDFALSE, indecimal = GDFALSE;
  int sign = +1;
  double n = 0.0, p = 1.0;

// First find what looks like start of a number - the first digit. 
// And note any sign and whether we just passed a decimal point.   
  do {
    c = fgetc(file);
    if (c == EOF)
      return (0);
    else if (stopateol && c == '\n')
      return (-1);
    else if (c == '+' || c == '-') {
      sign = (c == '+') ? +1 : -1;
      c = fgetc(file);
      if (c == EOF)
	return (0);
      else if (stopateol && c == '\n')
	return (-1);
    }
    if (c == '.') {
      indecimal = GDTRUE;
      c = fgetc(file);
      if (c == EOF)
	return (0);
      else if (stopateol && c == '\n')
	return (-1);
    }
    if (c >= '0' && c <= '9') {
      found = GDTRUE;
    } else {
      sign = +1;
      indecimal = GDFALSE;
    }
  } while (!found);

  // Now we've got digit(s) ...
  do {
    n = 10.0 * n + c - '0';
    p = 10.0 * p;
    c = fgetc(file);

    if ((c < '0') || (c > '9')) {
      found = GDFALSE;
  // We've run out. If we already saw a decimal point, return now 
      if (indecimal) {
	if (c != EOF)
	  ungetc(c, file);
	*valaddr = sign * n / p;
	return (1);
      } else
	p = 1.0;
    }
  } while (found);

  // We ran out and we didn't see a decimal point, so is this a decimal? 
  if (c != '.') {
  // No, give it back to caller  
    if (c != EOF)
      ungetc(c, file);
    *valaddr = sign * n;
    return (1);
  } else {
// It is. Step past it, carry on hoping for more digits 
    c = fgetc(file);
    while (c >= '0' && c <= '9') {
      n = 10.0 * n + c - '0';
      p = p * 10.0;
      c = fgetc(file);
    }
// We've run out of digits but we have a number to give 
    if (c != EOF)
      ungetc(c, file);
    *valaddr = sign * n / p;
    return (1);
  }
}

// Use getdouble() above but convert result to int. 
int getint(FILE * f, int *valaddr, int stopateol)
{
  int r;
  double x;
  r = (int)getdouble(f, &x, stopateol);
  *valaddr = (int) x;
  return (r);
}

// Use getdouble() above but convert result to long. 
long getlong(FILE * f, long *valaddr, int stopateol)
{
  int r;
  double x;
  r = (long)getdouble(f, &x, stopateol);
  *valaddr = (long) x;
  return (r);
}

// end Peter Ross' file-reading functions
