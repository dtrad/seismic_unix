/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUADDZEROS:  $Date: Julyh 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUADDZEROS - pad zeros from t=0 to t=f1                             ", 
  " suaddzeros < stdin > stdout [optional parameters]          		",
  " Example:                                                            ",
  " suwind tmin=5.4 < datain | .... | suaddzeros f1=5.4 > dataout       ",
  " 									",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr,tr2;

int main(int argc, char **argv)
{
  int j;
  register int it;
  int nt;
  float sum;
  int count;
  int if1;
  float dt;
  float f1;
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  //////////////////////////////////////////////////////////////////
 // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.dt) err("dt header field must be set");
  if (!getparfloat("f1",&f1)) f1=tr.f1;
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  //f1+=dt; 
  if1=(int) (f1/dt);

  /////////////////////////////////////////////////////////////////  
  count=0;
  j=0;
  do {   // Loop over traces
    j++; 
    for (it=0;it<if1-1;it++) tr2.data[it]=0; 
    for (it=if1-1;it<(nt+if1);it++) tr2.data[it]=tr.data[it-if1-1];
    memcpy(&tr2,&tr,HDRBYTES);
    tr2.f1=0;
    tr2.ns=nt+if1+1;
    tr2.delrt=0;
    puttr(&tr2);
  }while (gettr(&tr));

  return EXIT_SUCCESS;
}




















