/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUREMOVE:  $Date: Julyh 2000  */

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUREMOVE- remove zero traces                                        ", 
  " suremove< stdin > stdout [optional parameters]          		",
  "                                                                     ",
  " If only some interval should be touch use hmin and hmax             ",
  " hmin=0 if tr.offset < hmin	then pass even if dead    		",
  " hmax=1000000 if tr.offset > hmax then pass even if dead   		",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr;

int main(int argc, char **argv)
{
  int j;
  register int it;
  int count;
  float hmin;
  float hmax;
  float h;
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
  if (!getparfloat("hmin",&hmin)) hmin=0;
  if (!getparfloat("hmax",&hmax)) hmin=1000000;

 // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");

  count=0;
  j=0;
  do {   // Loop over traces
    j++; 
    h=fabs((float) tr.offset);
    if (tr.trid!=2){ puttr(&tr); count++;}
    else if ((h > hmax)||(h < hmin)) { puttr(&tr); count++;}
  }while (gettr(&tr));
  fprintf(stderr,"number of good traces = %d\n", count);
  return EXIT_SUCCESS;
}




















