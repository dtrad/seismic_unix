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
  int j, ih;
  register int it;
  int nt;
  int ntfinal;
  int ntr;
  float sum;
  int count;
  float dt;
  int nhbeg;
  int nhend;
  int ntraces;
  float dh;

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  //////////////////////////////////////////////////////////////////
 // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.dt) err("dt header field must be set");
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  ntr= (int) tr.ntr;
  if (!getparint("ntfinal",&ntfinal)) ntfinal=2*nt;
  if (!getparint("nhbeg",&nhbeg)) nhbeg=0;
  if (!getparint("nhend",&nhend)) nhend=0;
  if (!getparfloat("dh",&dh)) dh=50;

  /////////////////////////////////////////////////////////////////  
  count=0;
  j=0;
  ntraces=ntr+nhbeg+nhend;
  

  if (nhbeg){
    for (ih=0;ih<nhbeg;ih++){
      j++;
      for (it=0;it<ntfinal;it++) tr2.data[it]=0;
      memcpy(&tr2,&tr,HDRBYTES);
      tr2.ntr=ntraces;
      tr2.ns=ntfinal;
      tr2.offset=tr.offset-(nhbeg-ih)*dh;
      puttr(&tr2);     
    }
  }

  do {   // Loop over traces
    j++; 
    for (it=nt;it<ntfinal;it++) tr.data[it]=0;
    tr.ns=ntfinal;
    tr.ntr=ntraces;
    puttr(&tr);
  }while (gettr(&tr));

  if (nhend){
    for (ih=0;ih<nhend;ih++){
      j++;
      for (it=0;it<ntfinal;it++) tr2.data[it]=0;
      memcpy(&tr2,&tr,HDRBYTES);
      tr2.ntr=ntraces;
      tr2.ns=ntfinal;
      tr2.offset=tr.offset+(ih+1)*dh;
      puttr(&tr2);     
    }
  }


  return EXIT_SUCCESS;
}




















