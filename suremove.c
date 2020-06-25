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
  " 									",
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
  int nt;
  float sum;
  int count;
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  //////////////////////////////////////////////////////////////////
 // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");
  if (!tr.ns) err("ns header field must be set");
  nt = (int) tr.ns;

  /////////////////////////////////////////////////////////////////  
  count=0;
  j=0;
  do {   // Loop over traces
    j++; 
    sum=0; for (it=0;it<nt;it++) sum+=tr.data[it];
    if(fabs(sum)> 1e-7){  puttr(&tr); count++;}
  }while (gettr(&tr));
  fprintf(stderr,"number of good traces = %d\n", count);
  return EXIT_SUCCESS;
}




















