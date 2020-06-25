/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTEMPLATE1:  $Date: June 2000  */


#include "su.h"
#include "segy.h"
#include "header.h"
#include "radoninv.h"

void radon00inv(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float **model, float *q, int nq, float smute,float nmofactor, int pseudohyp, float depth, float fmax, int rtmethod);


/*********************** self documentation **********************/
char *sdoc[] = {
  " suradoninv data RTdata > stdout [optional parameters]   		",
  "        								",
  "        								",
  " 	This program reads two files, data and RTdata. The first one is ",
  "     the original data file and the second one is the Radon file     ",
  "     The inverse RT is performed on file2 and the headers from       ",
  "     file1 are used.                                                 ",
  "     Note that file1 provides only headers, so there is no use of	",
  "     the data themselves.             				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt                          */
/**************** end self doc ***********************************/

segy tr,tr2; 
FILE *fp1;		/* file pointer for first file		*/
FILE *fp2;		/* file pointer for second file		*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

int main(int argc, char **argv)
{
  int j,ih,iq;
  register int it;
  int nt;
  int ntr;
  int nh;
  int nq;
  float dt;
  float *h;
  float *q;
  int plot;
  float **d;
  float **m;

  float *velint;/* array[nt] of vel for a particular trace */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */

  int pseudohyp;
  float depth;
  float nmofactor;
  float  *t;     // time axis for input and output 
  float smute;
  float fmax;
  int rtmethod;

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparint("pseudohyp",&pseudohyp)) pseudohyp=0;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;

  /* Open two files given as arguments for reading */
  fp1 = efopen(argv[1], "r");
  fp2 = efopen(argv[2], "r");  
  tracefp = etmpfile();
  headerfp = etmpfile();


  // Read first file and save it to a temporal file 
  if (!fgettr(fp1,&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set"); 
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;  
  if (!(ntr=(unsigned long int) tr.ntr)) err("***ntr must be set\n");
  h=ealloc1float(ntr);

  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);


  j=0;
  do{   // Loop over traces    
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);  	   
    h[j]=tr.offset;
    j++;
  }while(fgettr(fp1, &tr));
  nh=j;
  erewind(tracefp);
  erewind(headerfp);
  fprintf(stderr,"nh=%d,nt=%d\n",nh,nt);
  float cdpgather=tr.cdp;

  d=ealloc2float(nt,nh);

  // Read second file and save the data to a 2D array m
 
  if (!fgettr(fp2,&tr)) err("can't read first trace");
  if (!(nq=(unsigned long int) tr.ntr)) err("***ntr must be set\n");

  m=ealloc2float(nt,nq);
  q=ealloc1float(nq);
  t=ealloc1float(nt);

  iq=0;
  do{   // Loop over traces    
    memcpy((void *) m[iq],(const void *) tr.data,nt*sizeof(float));
    q[iq]=tr.f2;
    iq++;
  }while(fgettr(fp2, &tr));
  nq=iq;
  fprintf(stderr,"nh2=%d\n",nq);
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////
  /* compute new square slowness and anis function */
  //#include "getvelocities.h" 
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  fprintf(stderr,"ncdp=%d\n",ncdp);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  radoninv0(h,nh,d ,t,nt,dt,velint,m,q,nq,smute,nmofactor,depth,fmax,rtmethod);
  
  for (ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data, (const void *) d[ih], nt*sizeof(float));
    puttr(&tr);
  }

  free2float(d);
  free2float(m);
  free1float(q);
  free1float(t);
  free1float(velint);
  free2float(ovv);
  free1float(h);
  free1float(cdp);
  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);

  return EXIT_SUCCESS;
}





























