/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUMYNMO:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUMYNMO - simple nmo test                                           ", 
  "	                                                         	",
  " 	   								",
  " sumynmo < stdin > stdout [optional parameters]          		",
  " 									",
  " vel = 2000                                                          ",
  "                          						",
  "                                                   			",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr,trf; 

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,float smute);

int main(int argc, char **argv)
{
  int j,i,k,ix;
  register int it;
  float *d;
  float *t; 
  int nt;
  float dt;	
  int invert;
  float *vnmo;
  float *tnmo;
  float *vel;
  int nvel;
  int ntime;
  float tn;
  float h;
  float smute;

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  if (!getparint("invert", &invert)) invert = 0;  
  if (!getparfloat("smute", &smute)) smute = 2.0;  
  //////////////////////////////////////////////////////////////////
  nvel = countnparval(1,"vnmo");nvel=(MAX(nvel,1));
  if ((vnmo=ealloc1float(nvel))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");
  if (!getnparfloat(1,"vnmo",vnmo)) vnmo[0]=2000;

  ntime = countnparval(1,"tnmo");ntime=(MAX(ntime,1));
  if ((tnmo=ealloc1float(ntime))==NULL)
    fprintf(stderr,"***Space for tnmo could not be allocated\n");
  if (!getnparfloat(1,"tnmo",tnmo)) tnmo[0]=0;
        
  if (ntime!=nvel){
    fprintf(stderr,"ntime=%d,nvel=%d\n",ntime,nvel);
    err("nvel and ntime must be equal\n");
  }
 // Get info from first trace 
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;

  /////////////////////////////////////////////////////////////////  
 
  // Create velocity function
  if ((vel=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");

  for (it=1; it<ntime; ++it)
    if (tnmo[it]<=tnmo[it-1]) err("tnmo values must increase monotonically");

  if ((vel=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");
  
  for (it=0,tn=0; it<nt; ++it,tn+=dt) 
    intlin(ntime,tnmo,vnmo,vnmo[0],vnmo[nvel-1],1,&tn,&vel[it]);  

  // Allocate memory for data and model

  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
        
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 

  for(it=0;it<nt;it++) t[it]=it*dt;

  fprintf(stderr,"**************************\n");  
  
  j=0;
  do {   // Loop over traces
    j++; 
    h=(float) tr.offset;
    if (invert) nmo(d,tr.data,t,h,vel,invert,nt,dt,smute);
    else nmo(tr.data,d,t,h,vel,invert,nt,dt,smute);	  
    memcpy((void *) tr.data,(const void *) d,nt*sizeof(float));
    puttr(&tr);
  }while (gettr(&tr));

  free1float(t);
  free1float(d);
  free1float(vel);
  free1float(tnmo);
  free1float(vnmo);

  return EXIT_SUCCESS;
}




















