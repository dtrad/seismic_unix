/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTEMPLATE1:  $Date: June 2000  */


#include "su.h"
#include "segy.h"
#include "header.h"
#include "radoninv.h"


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
  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  int j,ih,iq;
  register int it;
  int nt;
  int ntr;
  int nh=0;
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
  TRACE; 
  if (!getparint("plot",&plot)) plot=0;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparint("pseudohyp",&pseudohyp)) pseudohyp=0;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  //if (STREQ(offsetfile,"")) offsetfile=NULL;

  TRACE;
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
  nh = ntr;

  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);

  h=ealloc1float(5*nh); // allocate larger for now
  TRACE;
  /************************** output header information *****************/
  // If offsetfile is given interpolation is requested 
  // The first data file is read to obtained headers
  // If offsetfile is not given the offset is taken from this file

  j=0;
  do{   // Loop over traces    
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);  	   
    if (!offsetfile) h[j]=tr.offset;
    j++;
  }while(fgettr(fp1, &tr));
  erewind(tracefp);
  erewind(headerfp);
  fprintf(stderr,"Using %s offsetfile \n",offsetfile);
  // nh is given by offsetfile if it is given, otherwise by datafile1
  if (offsetfile) nh=read_ascii_file(offsetfile,h);

  else  nh=j;

  h=erealloc1float(h,nh);
  /**********************************************************************/

  fprintf(stderr,"nh=%d,nt=%d\n",nh,nt);
  float cdpgather=tr.cdp;

  /************************** input data information ********************/
  // Read second file and save the data to a 2D array m
  if (!fgettr(fp2,&tr)) err("can't read first trace");
  if (!(nq=(unsigned long int) tr.ntr)) err("***ntr must be set\n");

  d=ealloc2float(nt,nh);
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
  fprintf(stderr,"nq=%d\n",nq);

  for(it=0;it<nt;it++) t[it]=0+it*dt;

  /***********************************************************************/
  /* compute new square slowness and anis function */
  //#include "getvelocities.h" 
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  fprintf(stderr,"ncdp=%d\n",ncdp);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  TRACE;
  radoninv0(h,nh,d ,t,nt,dt,velint,m,q,nq,smute,nmofactor,depth,fmax,rtmethod);
  TRACE;  
  /*********************************************************************** 
      If offsetfile is NULL the headers from datafile1 are put them back 
      Otherwise the headers will need to be edited later
  ***********************************************************************/
  for (ih=0;ih<nh;ih++){
    if (!offsetfile || ih<ntr) efread(&tr,HDRBYTES,1,headerfp);
    tr.offset=(int) h[ih];
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

  fprintf(stderr,"%s Finished successfully\n",__FILE__);

  return EXIT_SUCCESS;
}





























