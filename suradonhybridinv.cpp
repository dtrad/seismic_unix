/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTEMPLATE1:  $Date: June 2000  */


#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonhybridinv.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " SURADONHYBRIDINV - Frequency domain inverse  Radon Transform  using ",
  "                    two  different operators.                        ",
  " suradonhybridinv  file1 file2 > stdout [optional parameters]     	",
  "                                                                     ",
  " Input must be sorted by offset.                                     ",
  "                                                                     ",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " nq=100          Size of Radon EOM                                   ",
  " nq1=nq/2+20     Size of model1                                      ",
  " nq2=nq-nq1      Size of model2                                      ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1.9   nmofactor * equiv offset is used for nmo		",
  " filter=0        =1 mutes model1, =2 mutes model2,                   ",
  "      	    =3 apply polygonal filter on 1                      ",
  "		    =4 apply polygonal filter on 2                      ",
  " rtmethod1= 1         						",
  " rtmethod2= 3         						",
  " fmax1 = 20        							",
  " fmax1 = 70        							",
  " depth1 = 1       For the pseudo offset of Foster and Mosher		",
  " depth2 = 3       Only used if rtmethod1 = 3 or rtmethod =3 		",
  "                                      				",
  "        								",
  " 	This program reads two files, the first one is the original     ",
  "     data file, which provides the offset for reconstruction         ",
  "     and the headers for the output traces. The second file is the   ",
  "     Radon transform of file1 computed with program suradonhybrid.   ",
  "     The q parameters are read from the file2 in the header word f2,	",
  "     The filter parameter is used to mute the radon space before	",
  "     reconstruction.							",
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

  float depth;
  float nmofactor;
  float  *t;     // time axis for input and output 
  float smute;

  /// For radon_2op
  int nq1;
  int nq2;
  int rtmethod1;
  int rtmethod2;
  float depth1;
  float depth2;
  float fmax1;
  float fmax2;
  // Filtering
  int filter;  /* =1 mutes model1, =2 mutes model2, 
		  =3 apply polygonal filter on 1
		  =4 apply polygonal filter on 2 */

  float *ffilter=0;
  float *amps=0;
  

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("smute",&smute)) smute=2;
  // The following are parameters use for radon_beam with two operators
  if (!getparfloat("depth1",&depth1)) depth1=1000;
  if (!getparfloat("depth2",&depth2)) depth2=1000;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("nq1", &nq1))  nq1 = nq/2+20;
  if (!getparint("nq2", &nq2))  nq2 = nq-nq1;
  if (!getparint("rtmethod1", &rtmethod1))  rtmethod1 = 1;
  if (!getparint("rtmethod2", &rtmethod2))  rtmethod2 = 3;
  if (!getparfloat("fmax1",&fmax1)) fmax1=20;
  if (!getparfloat("fmax2",&fmax2)) fmax2=70;
  if (!getparint("filter",&filter)) filter=0;
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

  int cdpgather=tr.cdp;
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

  radonhybridinv(d,h,nh,t,nt,dt,m,q,nq,velint,nq1,nq2,smute,nmofactor,rtmethod1,rtmethod2,depth1,depth2,fmax1,fmax2,filter,npoly,ffilter,amps);  
  TRACE;
  for (ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data, (const void *) d[ih], nt*sizeof(float));
    puttr(&tr);
  }
  TRACE;
  
  free1float(cdp);
  free2float(ovv);
  free2float(d);
  free2float(m);
  free1float(q);
  free1float(t);
  free1float(velint);
  free1float(h);
  TRACE;
  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);
  TRACE;
  return EXIT_SUCCESS;
}





























