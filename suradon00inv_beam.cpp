/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTEMPLATE1:  $Date: June 2000  */


#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonhybridinv.h"

void radon00inv(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float **model, float *q, int nq, float smute,float nmofactor, int pseudohyp, float depth);

void radon00inv_beam(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float **model, float *q, int nq, float smute, float nmofactor, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2,  int mute1, int mute2, float fmax1, float fmax2, float *ffilter, float *amps);


/*********************** self documentation **********************/
char *sdoc[] = {
  " sumixgathers file1 file2 > stdout [optional parameters]   		",
  "        								",
  "        								",
  " 	This template is for reading a model file generated from some   ",
  "     and from the model generate the original data                   ",
  "     The original data file is read to put back the original offset  ",
  "     but no computation is done with the tr.data                     ",
  "        								",
  "                                      				",
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
  int j,k,ih,iq;
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

  /// For radon_beam
  int comb_method;
  float qmin1;
  float qmin2;
  int nq1;
  int nq2;
  int rtmethod1;
  int rtmethod2;
  float depth1;
  float depth2;
  float factor1;
  float factor2;
  int mute1;
  int mute2;
  float fmax1;
  float fmax2;
  float *ffilter;
  float *amps;
  

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparint("pseudohyp",&pseudohyp)) pseudohyp=0;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("smute",&smute)) smute=2;

  // The following are parameters use for radon_beam with two operators
  if (!getparint("comb_method",&comb_method)) comb_method=0;
  if (!getparfloat("depth1",&depth1)) depth1=1000;
  if (!getparfloat("depth2",&depth2)) depth2=1000;
  if (!getparfloat("qmin1",&qmin1)) qmin1=-1e-4;
  if (!getparfloat("qmin2",&qmin2)) qmin2=-5e-4;
  if (!getparint("nq1", &nq1))  nq1 = nq/2+20;
  if (!getparint("nq2", &nq2))  nq2 = nq-nq1;
  if (!getparfloat("factor1",&factor1)) factor1=4;
  if (!getparfloat("factor2",&factor2)) factor2=4;
  if (!getparint("rtmethod1", &rtmethod1))  rtmethod1 = 1;
  if (!getparint("rtmethod2", &rtmethod2))  rtmethod2 = 3;
  if (!getparint("mute1",&mute1)) mute1=0;
  if (!getparint("mute2",&mute2)) mute2=0;
  if (!getparfloat("fmax1",&fmax1)) fmax1=20;
  if (!getparfloat("fmax2",&fmax2)) fmax2=70;
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  /* Get frequencies that define the filter */

  int npoly;
  int namps;
  if ((npoly = countparval("ffilter"))!=0) {
    ffilter = ealloc1float(npoly);
    getparfloat("ffilter",ffilter);
  } else npoly = 0;

  if ((namps = countparval("amps"))!=0) {
    amps = ealloc1float(namps);
    getparfloat("amps",amps);
  } else  namps = 0;

  if (!(namps==npoly)) 	err("number of f values must = number of amps values");

  for (k=0;k<npoly;k++) fprintf(stderr,"ffilter[%d]=%f\n",k,ffilter[k]);
  for (k=0;k<namps;k++) fprintf(stderr,"amps[%d]=%f\n",k,amps[k]);

  ////////////////////////////////////////////////////////////////////////

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
    //fprintf(stderr,"q[%d]=%e\n",iq,q[iq]);
    iq++;
  }while(fgettr(fp2, &tr));
  nq=iq;
  fprintf(stderr,"nq=%d\n",nq);
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  
  // Here we do something with m[iq][it] and put the result to d[ih][it]
  // Just for testing we make data1 = data2 
  ///////////////////////////////////////////////////
  /* compute new square slowness and anis function */
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  radon00inv_beam(h,nh,d,t,nt,dt,velint,m,q,nq,smute,nmofactor,nq1,
		  nq2,rtmethod1,rtmethod2,depth1,depth2,mute1,mute2,fmax1,fmax2,ffilter,
		  amps);
  
  
  for (ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data, (const void *) d[ih], nt*sizeof(float));
    puttr(&tr);
  }

  
  free1float(cdp);
  free2float(ovv);
  free2float(d);
  free2float(m);
  free1float(q);
  free1float(t);
  free1float(velint);
  free1float(h);

  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);

  return EXIT_SUCCESS;
}





























