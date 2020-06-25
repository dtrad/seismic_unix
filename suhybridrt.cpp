/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON00FORW  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonhybrid.h"
#include <signal.h>
#include <math.h>
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONHYBRID - Frequency domain Radon Transform  using two         ",
  " 	            different operators.	       			",
  " 	   								",
  " suradfreq < stdin > stdout [optional parameters]          		",
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
  " symmetric = 0    If = 1 defines q1 axis to be symmetric without	",
  "                  actually passing for zero                          ",
  " 	This program reads a data gather sorted by offset and produces  ",
  "     in the file = modelfile a Radon transform panel with two differents",
  "     model spaces, corresponding to two different Radon operators    ",
  "     The nature of each operator is defined by the corresponding     ",
  "     parameters rtmethod 1 and 2, depth 1 and 2, nq 1 and 2, etc.  	",
  "     The model obeys the following problem d=L1 m1 + L2 m2           ",
  "     The inverse RT gather (reconstructed) is output to stdoutput    ", 
  "     The filter parameter is used to mute the radon space before	",
  "     reconstruction.							",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset            */
/**************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
FILE *offsetfile; 
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr;
  cwp_String modelfile=""; /* output sufile for the model */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of offset traces
  int nt;   // number of time samples 
  float dt; // Time sampling
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  ///////////////////////////////////////////////////////////////
  int  k, ih, iq; // General counters 
  register int it;
  float *d;      /* single trace */
  float **data;   /* Common Scattering gather */
  float  *m;    /* Final modeled trace   */
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float t0;      // First useful time
  float *velint;/* array[nt] of vel for a particular trace */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  char *tmpdir;		/* directory path for tmp files		*/

  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 

  int verbose;		/* flag for echoing info		*/
  float fmax;
  int nq;  
  float *q;

  inv_par inv; // struct defined in inversion_par.h 
  float **model;            // Temporal array for the model
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  float smute;
  float nmofactor;
  int ntr;
  /// For radon_beam
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
  float fmax1;
  float fmax2;
  int filter;
  float *ffilter=0;
  float *amps=0;
  int symmetric;


  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
   // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  
  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (!(ntr=tr.ntr)) err("***ntr must be set\n");
  fmax=0.8/(2*dt);
  nh=ntr;

  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 5e-2;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 5e-2;
  if (!getparfloat("step", &inv.step))  inv.step = 0.9;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 5;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparint("norm",&inv.norm)) inv.norm=0;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";

  // The following are parameters use for radon_beam with two operators
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
  if (!getparint("filter",&filter)) filter=0;
  if (!getparfloat("fmax1",&fmax1)) fmax1=20;
  if (!getparfloat("fmax2",&fmax2)) fmax2=70;
  if (!getparint("symmetric",&symmetric))  symmetric = 0;
  ////////////////////////////////////////////////////////////////////////
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

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);

  /* Look for user-supplied tmpdir */

  if (!getparstring("tmpdir",&tmpdir) &&
      !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
  if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
    err("you can't write in %s (or it doesn't exist)", tmpdir);  
  
    // Store traces in the disk for later use
  if (STREQ(tmpdir,"")){
    tracefp = etmpfile();
    headerfp = etmpfile();
    if (verbose) warn("using tmpfile() call");
  }else{ /* user-supplied tmpdir */
    char directory[BUFSIZ];
    strcpy(directory, tmpdir);
    strcpy(tracefile, temporary_filename(directory));
    strcpy(headerfile, temporary_filename(directory));
    /* Trap signals so can remove temp files */
    signal(SIGINT,  (void (*) (int)) closefiles);
    signal(SIGQUIT, (void (*) (int)) closefiles);
    signal(SIGHUP,  (void (*) (int)) closefiles);
    signal(SIGTERM, (void (*) (int)) closefiles);
    tracefp = efopen(tracefile, "w+");
    headerfp = efopen(headerfile, "w+");
    istmpdir=cwp_true;		
    if (verbose) warn("putting temporary files in %s", directory);
  }
  
  // Allocate memory for data and model
  d=ealloc1float(nt);
  data=ealloc2float(nt,nh);
  model=alloc2float(nt,nq);
  m=ealloc1float(nt);
  t=ealloc1float(nt);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  Wm=alloc1float(nt*nq);
  Wd=alloc1float(nt*nh);

    
  // arrays for nmo model
  // Velint will contain velocity law for one cmp
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////////////////////////
  
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ntr = 0;
  do {
    ntr++;
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);
  ih=0;
  nh=ntr;
  
  for (ih=0;ih<nh;ih++){
    efread(tr.data, FSIZE, nt, tracefp);
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    efread(&tr,HDRBYTES,1,headerfp);
    h[ih]=tr.offset;
  }
  int cdpgather=tr.cdp;
  /* compute new square slowness and anis function */
  
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);
  
  radonhybrid(data,h,nh,t,nt,dt,model,q,nq,velint,inv,qmin1,qmin2,nq1,nq2,factor1,factor2,smute,nmofactor,rtmethod1,rtmethod2,depth1,depth2,fmax1,fmax2,filter,npoly,ffilter,amps,symmetric);

  // Output trace
  erewind(tracefp);
  erewind(headerfp);  
    
  for (ih=0;ih<nh;ih++){ 
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
    puttr(&tr);
  }

  if (1){
    save_gather(model,nq,nt,0.004,"model2");
    system("sufft  < model2 | suamp | suxwigb title=\"model2\" perc=99 & ");
  }

  
  for (iq=0;iq<nq;iq++){
    memcpy((void *) tr.data,(const void *) model[iq],nt*sizeof(float));
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nq;
    tr.ns=nt;
    tr.tracl=iq+1;
    tr.tracr=iq+1;
    tr.f2=q[iq];
    fputtr(modelfilep,&tr);    
  }

  free1float(cdp);
  free2float(ovv);
  free1float(Wd);
  free1float(Wm);
  free1float(q);
  free1float(velint);
  free1float(t);
  free1float(h);
  free1float(m);
  free2float(data);
  free2float(model);
  free1float(d);
  efclose(headerfp);
  efclose(tracefp);
  efclose(modelfilep);

  if (istmpdir) eremove(headerfile);
  if (istmpdir) eremove(tracefile);  

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}






























#include "su.h"
#include "radonhybrid.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, inv_par inv, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps, int symmetric)
{
  int it, ih;
  float *dtemp;
  float *Wd;
  int testadj=0;
  float dq1=0;
  float dq2=0;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2; 
  TRACE;
  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);
  
  fprintf(stderr,"fmax1=%f,fmax2=%f\n",fmax1,fmax2);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  }

  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  
  

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);
  
  dataweigths(h,nh,Wd,TRUE);
  TRACE;
  radon_param_2op(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2, nq1, nq2, depth1, depth2,
		   rtmethod1, rtmethod2, factor1, factor2,&dq1,&dq2, symmetric);
  TRACE;

  radon_wtcgls_2op(data,ph1,ph2,nh,t,nt,dt,model,q,nq,q1,nq1,q2,nq2,inv,Wd,testadj,fmax1,fmax2);
  TRACE;
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);  

  radoninv_2op(data,ph1,ph2,nh,t,nt,model,q,nq,q1,nq1,q2,nq2,fmax1,fmax2,filtering,npoly,ffilter,amps);
  
  /////////////////////////////

  plotgather(data[0],nt,nh,"data_before_inmo");

  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }    
  }

  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(Wd);
  free1float(dtemp);

  return;

}
























#include "su.h"
#include "radonhybrid.h"


/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/


void radon_wtcgls_2op(float **data, float *ph1, float *ph2, int nh, float *t, int nt, float dt, float **model, float *q, int nq, float *q1, int nq1, float *q2, int nq2, inv_par inv, float *Wd, int testadj, float fmax1, float fmax2)
{
  int iq,iter;
  complex czero;  czero.r=czero.i=0;
  float w=0,df;
  int freq, nf,nfft;
  int maxfreq1, maxfreq2;
  complex **m2;
  complex **d2;
  complex **L;
  complex **L1;
  complex **L2;
  float *Wm;
  float  *Jtot;
  float fmin=0; // Minimum freq to compute in Hz;
  float sigmam;
  float sigmad;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;

  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)
  Jtot=ealloc1float(20);
  Wm=ealloc1float(nq);

  //for (iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);

  zero_vector(Jtot,inv.iter_end);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
  
  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);
  

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
    
     //if (testadj) test=testadj_rad_f(L,LH);
    
    for (iter=1;iter<=inv.iter_end;iter++){
      weights_inv(m2[freq],nq,inv.norm,sigmam,Wm,iter);
      if (iter==1 && freq < 10 ){
	for (iq=0;iq<nq1;iq++) Wm[iq]*=2;
	for (iq=nq1;iq<nq;iq++) Wm[iq]*=0.5;
      }
      wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,0,inv.step,inv.itercg);
    }
    
    fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
  }
  //freqweighting(m2,nf,nq,df,fmin,fmax1);

  if (nq2>0){
    fprintf(stderr,"From now on one operator only\n");
    //////////////////////////////////////////////////////////////
    // For Ground roll,
    // Between maxfreq1 and freqmax2   we use only one operator
    complex **m2w=window(m2,0,nf-1,nq1,nq-1);
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      
      for (iter=1;iter<=inv.iter_end;iter++){
	if (iter==2) 
	  deviations(m2w[freq],nq2,d2[freq],nh,inv.norm,quantil1,quantil2,&sigmam,&sigmad);
	weights_inv(m2w[freq],nq2,inv.norm,sigmam,Wm,iter);   
	wtcgls(d2[freq],L2,m2w[freq],Wm,Wd,nh,nq2,0,inv.step,inv.itercg);
      }
      fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
    }
    freqweighting(m2w,nf,nq2,df,fmin,fmax2);
   }

   fprintf(stderr,"w=%f\n",w);
   
   
   for (iq=nq1;iq<nq;iq++) m2[freq][iq]/=nq;      
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  



   free1float(Wm);
   free1float(Jtot);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}








































#include "su.h"
#include "radonhybrid.h"

void radoninv_2op(float **data,  float *ph1, float *ph2, int nh, float *t, int nt, 
		  float **model, float *q,  int nq, float *q1, int nq1, float *q2, int nq2, 
		  float fmax1, float fmax2, int filter, int npoly, float *f, float *amps)
{
  int freq;
  complex **m2;
  complex **d2;
  complex czero; czero.r=czero.i=0;
  complex **L;
  complex **L1;
  complex **L2;
  float w,df;
  float dt=t[1]-t[0];
  int nfft;
  int nf;
  int maxfreq1, maxfreq2;
  float fmin=0;
  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh);  
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);

  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  if (filter==3) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);
  //if (filter==4) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
 
    Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
    fprintf(stderr,":");
    
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
  }
  if (nq2>0){ // Between maxfreq1 and freqmax2   we use only one operator
    fprintf(stderr,"From now on one operator only\n");
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      Atimesx(d2[freq],L2,&m2[freq][nq1],nh,nq2,FALSE);
      fprintf(stderr,".");      
    }
   }

  fprintf(stderr,"freq=%f\n",freq*df);      

  freqweighting(d2,nf,nh,df,fmin,fmax2);
  
  fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf); 

  TRACE;
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  TRACE;
  return;

}






























#include "su.h"
#include "stddef.h"

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod)
/* Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    Daniel Trad- UBC- 16-2-99
*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   if (rtmethod==2) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1) { //LRT
	  dq= 1/(fmax*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod,
float factor)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   
   if (rtmethod==2 || rtmethod==5) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1 || rtmethod==3) { //LRT
	  dq= 1/(fmax*fabs(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}
	















#include "su.h"
#include "stddef.h"



void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);
void interval(float *pos, int nh, float *dx_max, float *dx_av);


void radon_param_2op(float fmax1, float fmax2, float *h, int nh, float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1,
		      float factor2,float *pdq1, float *pdq2, int symmetricq1)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{	
  float qmaxt;
  float qmax;
  float dx_max;
  float dx_av;
  float dq1;
  float dq2;
  int iq;

  fprintf(stderr,"qmin1=%e,qmin2=%e\n",qmin1,qmin2);

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
  radon_param(fmax1,h,nh,dx_av,qmin1,&qmaxt,&qmax,&dq1,nq1,rtmethod1,factor1);
    
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax1,dq1);

  if (symmetricq1){
    for (iq=0;iq<nq1/2;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
    //fprintf(stderr,"q[%d]=%f\n",nq1/2-1,q1[nq1/2-1]);
    q1[nq1/2]=-q1[nq1/2-1];
    fprintf(stderr,"q[%d]=%f\n",nq1/2,q1[nq1/2]);
    for (iq=nq1/2+1;iq<nq1;iq++){
      q1[iq]=q1[iq-1]+dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
  }
  else{
    //qmin2=q1[nq1-1]+dq1;
    for (iq=0;iq<nq1;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q1[%d]=%f\n",iq,q1[iq]);
    }
  }
  
  radon_param(fmax2,h,nh,dx_av,qmin2,&qmaxt,&qmax,&dq2,nq1,rtmethod2,factor2);

  for (iq=0;iq<nq2;iq++) q2[iq]=qmin2+iq*dq2;

  for (iq=0;iq<nq1;iq++) q[iq]=q1[iq];
  for (iq=0;iq<nq2;iq++){
    q[iq+nq1]=q2[iq];
    fprintf(stderr,"q2[%d]=%e\n",iq,q2[iq]);
  }
  *pdq1=dq1;
  *pdq2=dq2;

  return;
}










#include "su.h"

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t)
{
	static int indx=0;
	int it;
	float a1,a2;

	/* if before first cdp, constant extrapolate */
	if (cdpt<=cdp[0]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[0][it];
			oa1t[it] = oa1[0][it];
			oa2t[it] = oa2[0][it];
		      };
	
	/* else if beyond last cdp, constant extrapolate */
	} else if (cdpt>=cdp[ncdp-1]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[ncdp-1][it];
			oa1t[it] = oa1[ncdp-1][it];
			oa2t[it] = oa2[ncdp-1][it];
		      };
	
	/* else, linearly interpolate */
	} else {
		xindex(ncdp,cdp,cdpt,&indx);
		a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
		a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
		for (it=0; it<nt; ++it) {
			ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
			oa1t[it] = a1*oa1[indx][it]+a2*oa1[indx+1][it];
			oa2t[it] = a1*oa2[indx][it]+a2*oa2[indx+1][it];
		      };
	}
}

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt)
{
  static int indx=0;
  int it;
  float a1,a2;
  
  /* if before first cdp, constant extrapolate */
  if (cdpt<=cdp[0]) for (it=0; it<nt; ++it) ovvt[it] = ovv[0][it];
  /* else if beyond last cdp, constant extrapolate */
  else if (cdpt>=cdp[ncdp-1]) for (it=0; it<nt; ++it) ovvt[it] = ovv[ncdp-1][it];
  /* else, linearly interpolate */
  else {
    xindex(ncdp,cdp,cdpt,&indx);
    a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
    a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
    for (it=0; it<nt; ++it) ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
  }
}











/*
  It computes the term Cm for the model weight 
  for the Hessian (L' Cd^{-1} L+ Cm^{-1})
  This term corresponds to the probability model, so that
  the distribution parameters sigma and norm are passed.
  Wm is a vector of dimension nx, but in fact is the diagonal 
  of the Wm matrix of size nx x nx.
  Inqut 
          m: model
          nx: number of model traces
          norm: implemented 1 Huber, 0 Cauchy, else L2
  Output
          Wm 
          eps1: standard deviation of the model

  Daniel Trad- 14 March 2000. UBC- Canada
  Based in Sacchi, 1996. phD thesis. UBC. Canada
*/
#include "su.h"
#include "math.h"
#include "Complex.h"

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=MAX(fabs(m[i]),sigmam);
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}
void weights_cgfft(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=sqrt(MAX(fabs(m[i]),sigmam));
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=sqrt(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter)
{ 
      int i;
      
      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
      
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=1./sqrt(abs(m[i])*sigmam);
      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=1./sqrt(sigmam*sigmam+abs(m[i]*m[i]));
      else if(norm==2){ // Mask
	for (i=0;i<nx;i++)
	  Wm[i]=1+200./(1+exp(1*(abs(m[i])-sigmam)+0.5));
      }
      return;
}

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter)
{ 
  /* 
     It computes the model weights using a window in the model space m(f,q).
     For example, Wm[iq] is a function of m(f,iq), m(f-1,iq), ..., m(f-buffer,iq)
     Particularly useful for the dealiased RT (Hermman et al.)

  */     
      int  iq, iw;
      float maveg;

      if (iter==1){ 
	for (iq=0;iq<nq;iq++) Wm[iq]=1;
	return;
      }
      
      if (norm==1){ 
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(fabs(maveg)*sigmam);
	}
      }
      else if(norm==0){
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(sigmam*sigmam+fabs(maveg*maveg));
	}
      }
      else if(norm==2){ // Mask
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1+200./(1+exp(1*(fabs(maveg)-sigmam)+0.5));
	}
      }
      return;
}

void modelweight(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=fabs(m[isamax(nx,m,1)]);
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=fabs(m[i]);
      }
      else if(norm==0){
	if (maxm>1e-4) 
	  for (i=0;i<nx;i++){
	    // Solved !!!!!!!!!!!!!!!!!
	    // The right Wm from Cauchy is 
	    // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // But if M^-1 ATA x = M^-1 AT b is solved instead
	    // of the satndard form M=WmT *Wm 
	    Wm[i]=(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // Actually it works even better with (I don't know why)
	    //Wm[i]=Wm[i]*Wm[i];
	    //if (Wm[i]>2) Wm[i]=2; 
	  }
	else for (i=0;i<nx;i++) Wm[i]=1e-3;
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}

void modelweight_inv(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=m[isamax(nx,m,1)];
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=1./MIN(fabs(m[i]),eps1);
      }
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
 	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=1./(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}





void modelweight(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
            
      if (norm==1)
	for (i=0;i<nx;i++) Wm[i]=abs(m[i]);

      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=eps1+pow(abs(m[i]),2.0);

      else for (i=0;i<nx;i++) Wm[i]=1.;     


      return;
}

void modelweight_inv(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
  float *mabs;
  int ix;
  float sigma;

  mabs=ealloc1float(nx);
  if (norm==0 || norm ==1){ 
    for (ix=0;ix<nx;ix++) mabs[ix]=abs(m[ix]);
    sigma=quest(eps1,nx,mabs);
    //sigma=MAX(sigma,1e-3);
    //sigma=eps1;
    if (norm==1) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(MAX(mabs[ix],sigma)))+1e-7;
    else if(norm==0) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(sigma*sigma+pow(mabs[ix],2.0)))+1e-7;
  }
  else for (ix=0;ix<nx;ix++) Wm[ix]=1.;     

  free1float(mabs);
  return;
}



void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=1;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=fabs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=fabs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);
      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);
      free1float(mabs);
      free1float(dabs);

      return;
}

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=0;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=abs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=abs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);

      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);

      free1float(mabs);
      free1float(dabs);

      return;
}


void dataweigths(float *pos, int nh, float *Wd, int add)
{
  int ih;
  float *dh;
  float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
  int verbose=0;

  dh=ealloc1float(nh);

  for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   

  dh[0]=pos[1]-pos[0];
  dh[nh-1]=pos[nh-1]-pos[nh-2];
  
  for (ih=0;ih<nh;ih++){
    if (add) Wd[ih]*=MIN(dh[ih],dhav);
    else  Wd[ih]=MIN(dh[ih],dhav);
    if (verbose) fprintf(stderr," Wd[%d]=%f\n",ih,Wd[ih]);
  }

  free1float(dh);
  return;

}









































//#include <math.h>
#include "su.h"

void nmo(float *d,float *m,float *t,float h,float *vel,int invert,int nt,float dt,float smute) 
{
  register int it;
  float moveout;
  float t2;
  float vel2;
  float dt2=dt*dt;
  float *ttn;
  float *tnt;
  //  float smute=1.5; /* zero samples with NMO stretch exceeding smute */
  float osmute;	   /* 1/smute */
  int lmute=25;    /* length in samples of linear ramp for mute */
  int itmute=0;	   /* zero samples with indices less than itmute */
  float *atn;	   /* amplitude a(tn) for NMO */
  float *at;	   /* amplitude a(t) for inverse NMO */
  int sscale=0;	   /* if non-zero, apply NMO stretch scaling */
  ttn=ealloc1float(nt);
  at = ealloc1float(nt);
  atn = ealloc1float(nt);
  tnt=ealloc1float(nt);

  if (invert) for(it=0;it<nt;it++) d[it]=0;
  else for (it=0;it<nt;it++) m[it]=0;

  for (it=0;it<nt;it++){
    vel2=vel[it]*vel[it];
    t2=t[it]*t[it];
    moveout=h*h/vel2;
    ttn[it]=sqrt(t2/dt2+moveout/dt2);
  }
  /* compute inverse of stretch factor a(tn) */
  atn[0] = ttn[1]-ttn[0];
  for (it=1; it<nt; ++it)
    atn[it] = ttn[it]-ttn[it-1];
  //fprintf(stderr,"smute=%f\n",smute);
  /* determine index of first sample to survive mute */
  osmute = 1.0/smute;
  for (it=0,itmute=0; it<nt && atn[it]<osmute; ++it)
    itmute++;
  
  if (invert){ 
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    /* adjust mute time */
    itmute = (int) (1.0+ttn[itmute]);
    itmute = MIN(nt-2,itmute);
    
    /* compute a(t) */
    if (sscale) {
      for (it=itmute+1; it<nt; ++it)
	at[it] = tnt[it]-tnt[it-1];
      at[itmute] = at[itmute+1];
    }

    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      d[it] = 0.0;
			
    /* if specified, undo NMO stretch factor scaling */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	d[it] *= at[it];


  }
  else{
    ints8r(nt,1.0,0,d,0.0,0.0,nt,ttn,m);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      m[it] = 0.0;
    
    /* apply linear ramp */
    for (it=itmute; it<itmute+lmute && it<nt; ++it)
      m[it] *= (float)(it-itmute+1)/(float)lmute;
    
    /* if specified, scale by the NMO stretch factor */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	m[it] *= atn[it];
  }

  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}










#include "su.h"

void radon_matrix(complex *R, complex **l,complex **lh,float *g,float *q,int nh,int nq,float w,float *dh)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];

        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;

        for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      lh[iq][ih]=(1./Aperture)*dh[ih]*co;
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
	for (iq=0;iq<nq;iq++){
	  R[iq].r=0;
	  R[iq].i=0;
	  for (ih=0;ih<nh;ih++)
	    R[iq]+=lh[0][ih]*l[ih][iq]; //Top row of LL=LH*L
	  //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
	}

        return;
}

void radon_matrix(complex **l, float *g,float *q,int nh,int nq,float w)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];


	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
        return;
}

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth)
{
  int ih;
  
  if (rtmethod==1) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih];
  else if (rtmethod==2) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih]*h[ih];
  else if (rtmethod==3) 
    for (ih=0;ih<nh;ih++) g[ih]=sqrt(h[ih]*h[ih]+depth*depth)-depth;   

  return;

}


void radon_matrix_irrq(complex **L, float *h, float *q,int nh,int nq,float w)
{
        
        register int ih;
	register int iq;  
        complex  arg;
	arg.r=0;
        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.i=w*h[ih]*q[iq];
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}


void radon_matrix_2op(complex **L, complex **L1, complex **L2, int nh, int nq1, int nq2)
{
  int ih, iq;

  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq1;iq++)
      L[ih][iq]=L1[ih][iq];
    for (iq=0;iq<nq2;iq++) 
      L[ih][iq+nq1]=L2[ih][iq];
  }	


  return;
}





#include "su.h"
void interval(float *x,int lx,float *pmx, float *pax)
{
   int i;
   float dx, mx, ax;		
   mx=0;
   ax=0;	
     for (i=1;i<lx;i++){
	 dx=fabs(x[i]-x[i-1]);
	 if (dx>mx) 
	        mx=dx;
	 ax=ax+dx;
     }
     ax=ax/(lx-1);
     *pmx=mx;
     *pax=ax;
     return;
}
#include "su.h"
#include "clibrary.h"
#include <math.h>

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	    float *Wd,int nh,int nq, float tol, float step, int itercg)
{
  /* This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
     Notice that LH=FH WdT and L= Wd F
     so that Wd only needs to be used the first time 
  */
 
  float normb,dq,dq2,nit,beta, alpha, alphanum, alphaden;
  int j,in,num;
  complex czero;

  register int i;
  complex *q;
  complex *q1;
  complex *s;
  complex *x1;
  complex *z;
  complex *z1;
  complex *r;
  complex *Az;
  complex *xold;
  float *eta;
  float *rho;
  float *gcv;
  float J;
  complex *r2;
  complex *xtemp;

  q=ealloc1complex(nq); 
  q1=ealloc1complex(nq);
  s=ealloc1complex(nq);
  x1=ealloc1complex(nq);
  z=ealloc1complex(nq);
  z1=ealloc1complex(nq);
  r=ealloc1complex(nh);
  Az=ealloc1complex(nh);
  eta=ealloc1float(nq);
  rho=ealloc1float(nq);
  gcv=ealloc1float(nq);
  r2=ealloc1complex(nh);
  xtemp=ealloc1complex(nq);
  xold=ealloc1complex(nq);

  czero.r=czero.i=0;
  for (i=0;i<nq;i++) x[i]=czero;
  normb=sqrt(rcdot(nh,b,b));
  //xequaly(r,b,nh);
  for (i=0;i<nh;i++) r[i]=Wd[i]*b[i];
  for (i=0;i<nh;i++) r2[i]=r[i]*Wd[i];  
  Atimesx(r2,L,s,nh,nq,1);
  
  nit=MIN(itercg,nq);
  for(i=0;i<nq;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  xequaly(z,q,nq);
  dq=rcdot(nq,s,q);
  xequaly(z1,q1,nq);
  for(i=0;i<nq;i++) x1[i]=czero;       
  for (j=0;j<nit;j++){
    Atimesx(Az,L,z,nh,nq,0);            
    for (i=0;i<nh;i++) Az[i]*=Wd[i];  
    alphanum=dq;
    alphaden=rcdot(nh,Az,Az);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1.e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e\n",
	      alphanum,alphaden);
      break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
         
    //// Update model u and residuals

    for(i=0;i<nh;i++) r[i]-=alpha*Az[i];
    for(i=0;i<nq;i++) xold[i]=x[i];    
    for(i=0;i<nq;i++) x[i]+=alpha*z[i];    

    for(i=0;i<nh;i++) r2[i]=r[i]*Wd[i];
    Atimesx(r2,L,s,nh,nq,1);

    for(i=0;i<nq;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=rcdot(nq,s,q);
    beta=dq2/dq;
    dq=dq2;
    for (i=0;i<nq;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(rcdot(nq,s,s))/normb;
    //fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);
    for (i=0;i<nq;i++) xtemp[i]=x[i]/Wm[i];
    J=rcdot(nh,r2,r2)+rcdot(nq,xtemp,xtemp);
    for (i=0;i<nq;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(rcdot(nq,x1,x1));
    if ((tol==0) && (j>1)){ // GCV criteria
       in = j;
       for (i=0;i<in;i++){
       num=(nq-i)*(nq-i); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if ((gcv[j-2]<gcv[j-1])&&(gcv[j-1]<gcv[j])) { 
         if (1) fprintf(stderr,"GCV Criteria, iteration %d\n",j-1);
	 for(i=0;i<nq;i++) x[i]=xold[i];    
         nit = j-1;
         break;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        if (1) fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        break;
       }
    }          
  }
  if (0) fprintf(stderr,"j=%d\n",j);

  free1complex(xold);
  free1complex(xtemp);
  free1complex(r2);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1complex(Az);
  free1complex(r);
  free1complex(z1);
  free1complex(z);
  free1complex(x1);
  free1complex(s);
  free1complex(q1);
  free1complex(q);

 

  return(J);
}















#include "su.h"

 
void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv)
{
  /* This function reads the par file in the command line and computes the vector
  ovv. This vecotr can then be used with interpovv to get the interpolated velocities at the desired cdp
 input :
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 

The calling function requires this:

  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv);

  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  free1float(cdp);
  free2float(ovv);
  free1float(velint);

  */
  int it;
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////

  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  //fprintf(stderr,"countparname=%d\n",countparname("vnmo"));
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
  }
  
  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
  } 
}





void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1,
		   float **oa2)
{
  /* This function reads the par file in the command line and computes the vectors
  ovv, oa1, and oa2. These 3 vecotrs can then be sued with interpovv to get the interpolated velocities at the desired cdp
 input :
    float **oa1;	 array[ncdp][nt] of anis1 functions 
    float **oa2;	 array[ncdp][nt] of anis2 functions 
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float *velint; array[nt] of vel for a particular trace 
    float *oa1t;	 array[nt] of anis1 for a particular trace 
    float *oa2t;	 array[nt] of anis2 for a particular trace 

The calling function requires this:

  float **oa1;
  float **oa2;
  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;
  float *oa1t;	
  float *oa2t;	

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv,oa1,oa2);

  interpovv(nt,ncdp,cdp,ovv,oa1,oa2,cdpgather,velint,oa1t,oa2t);

  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(velint);

  */

  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  int i, it;

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////


  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
    if (countparname("anis1")!=ncdp &&
	countparname("anis1")!=0)
      err("an anis1 array must be specified for each cdp, "
	  "or omitted at all");
    if (countparname("anis2")!=ncdp &&
	countparname("anis2")!=0)
      err("an anis2 array must be specified for each cdp, "
	  "or omitted at all");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
    if (countparname("anis1")>1)
      err("only one (or no) anis1 array must be specified");
    if (countparname("anis2")>1)
      err("only one (or no) anis2 array must be specified");    
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    nanis1 = countnparval(icdp+1,"anis1");
    nanis2 = countnparval(icdp+1,"anis2");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nanis1!=nvnmo && nanis1 != 0)
      err("number of vnmo and anis1 values must be equal");
    if (nanis2!=nvnmo && nanis2 != 0)
      err("number of vnmo and anis2 values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    if (nanis1==0) nanis1 = nvnmo;
    if (nanis2==0) nanis2 = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    anis1 = ealloc1float(nvnmo);
    anis2 = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    if (!getnparfloat(icdp+1,"anis1",anis1)) 
      for (i=0; i<nvnmo; i++) anis1[i] = 0.0;
    if (!getnparfloat(icdp+1,"anis2",anis2))
      for (i=0; i<nvnmo; i++) anis2[i] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis1,anis1[0],anis1[nanis1-1],1,&tn,&oa1[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis2,anis2[0],anis2[nanis2-1],1,&tn,&oa2[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
    free1float(anis1);
    free1float(anis2);
  }

  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    aoa1 = oa1[jcdp];
    aoa2 = oa2[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
      oa1[icdp+1] = oa1[icdp];
      oa2[icdp+1] = oa2[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
    oa1[icdp+1] = aoa1;
    oa2[icdp+1] = aoa2;
  } 

}








#include "su.h"
#include "radonhybrid.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, inv_par inv, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps, int symmetric)
{
  int it, ih;
  float *dtemp;
  float *Wd;
  int testadj=0;
  float dq1=0;
  float dq2=0;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2; 
  TRACE;
  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);
  
  fprintf(stderr,"fmax1=%f,fmax2=%f\n",fmax1,fmax2);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  }

  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  
  

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);
  
  dataweigths(h,nh,Wd,TRUE);
  TRACE;
  radon_param_2op(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2, nq1, nq2, depth1, depth2,
		   rtmethod1, rtmethod2, factor1, factor2,&dq1,&dq2, symmetric);
  TRACE;

  radon_wtcgls_2op(data,ph1,ph2,nh,t,nt,dt,model,q,nq,q1,nq1,q2,nq2,inv,Wd,testadj,fmax1,fmax2);
  TRACE;
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);  

  radoninv_2op(data,ph1,ph2,nh,t,nt,model,q,nq,q1,nq1,q2,nq2,fmax1,fmax2,filtering,npoly,ffilter,amps);
  
  /////////////////////////////

  plotgather(data[0],nt,nh,"data_before_inmo");

  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }    
  }

  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(Wd);
  free1float(dtemp);

  return;

}
























#include "su.h"
#include "radonhybrid.h"


/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/


void radon_wtcgls_2op(float **data, float *ph1, float *ph2, int nh, float *t, int nt, float dt, float **model, float *q, int nq, float *q1, int nq1, float *q2, int nq2, inv_par inv, float *Wd, int testadj, float fmax1, float fmax2)
{
  int iq,iter;
  complex czero;  czero.r=czero.i=0;
  float w=0,df;
  int freq, nf,nfft;
  int maxfreq1, maxfreq2;
  complex **m2;
  complex **d2;
  complex **L;
  complex **L1;
  complex **L2;
  float *Wm;
  float  *Jtot;
  float fmin=0; // Minimum freq to compute in Hz;
  float sigmam;
  float sigmad;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;

  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)
  Jtot=ealloc1float(20);
  Wm=ealloc1float(nq);

  //for (iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);

  zero_vector(Jtot,inv.iter_end);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
  
  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);
  

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
    
     //if (testadj) test=testadj_rad_f(L,LH);
    
    for (iter=1;iter<=inv.iter_end;iter++){
      weights_inv(m2[freq],nq,inv.norm,sigmam,Wm,iter);
      if (iter==1 && freq < 10 ){
	for (iq=0;iq<nq1;iq++) Wm[iq]*=2;
	for (iq=nq1;iq<nq;iq++) Wm[iq]*=0.5;
      }
      wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,0,inv.step,inv.itercg);
    }
    
    fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
  }
  //freqweighting(m2,nf,nq,df,fmin,fmax1);

  if (nq2>0){
    fprintf(stderr,"From now on one operator only\n");
    //////////////////////////////////////////////////////////////
    // For Ground roll,
    // Between maxfreq1 and freqmax2   we use only one operator
    complex **m2w=window(m2,0,nf-1,nq1,nq-1);
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      
      for (iter=1;iter<=inv.iter_end;iter++){
	if (iter==2) 
	  deviations(m2w[freq],nq2,d2[freq],nh,inv.norm,quantil1,quantil2,&sigmam,&sigmad);
	weights_inv(m2w[freq],nq2,inv.norm,sigmam,Wm,iter);   
	wtcgls(d2[freq],L2,m2w[freq],Wm,Wd,nh,nq2,0,inv.step,inv.itercg);
      }
      fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
    }
    freqweighting(m2w,nf,nq2,df,fmin,fmax2);
   }

   fprintf(stderr,"w=%f\n",w);
   
   
   for (iq=nq1;iq<nq;iq++) m2[freq][iq]/=nq;      
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  



   free1float(Wm);
   free1float(Jtot);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}








































#include "su.h"
#include "radonhybrid.h"

void radoninv_2op(float **data,  float *ph1, float *ph2, int nh, float *t, int nt, 
		  float **model, float *q,  int nq, float *q1, int nq1, float *q2, int nq2, 
		  float fmax1, float fmax2, int filter, int npoly, float *f, float *amps)
{
  int freq;
  complex **m2;
  complex **d2;
  complex czero; czero.r=czero.i=0;
  complex **L;
  complex **L1;
  complex **L2;
  float w,df;
  float dt=t[1]-t[0];
  int nfft;
  int nf;
  int maxfreq1, maxfreq2;
  float fmin=0;
  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh);  
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);

  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  if (filter==3) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);
  //if (filter==4) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
 
    Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
    fprintf(stderr,":");
    
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
  }
  if (nq2>0){ // Between maxfreq1 and freqmax2   we use only one operator
    fprintf(stderr,"From now on one operator only\n");
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      Atimesx(d2[freq],L2,&m2[freq][nq1],nh,nq2,FALSE);
      fprintf(stderr,".");      
    }
   }

  fprintf(stderr,"freq=%f\n",freq*df);      

  freqweighting(d2,nf,nh,df,fmin,fmax2);
  
  fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf); 

  TRACE;
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  TRACE;
  return;

}






























#include "su.h"
#include "stddef.h"

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod)
/* Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    Daniel Trad- UBC- 16-2-99
*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   if (rtmethod==2) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1) { //LRT
	  dq= 1/(fmax*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod,
float factor)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   
   if (rtmethod==2 || rtmethod==5) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1 || rtmethod==3) { //LRT
	  dq= 1/(fmax*fabs(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}
	















#include "su.h"
#include "stddef.h"



void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);
void interval(float *pos, int nh, float *dx_max, float *dx_av);


void radon_param_2op(float fmax1, float fmax2, float *h, int nh, float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1,
		      float factor2,float *pdq1, float *pdq2, int symmetricq1)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{	
  float qmaxt;
  float qmax;
  float dx_max;
  float dx_av;
  float dq1;
  float dq2;
  int iq;

  fprintf(stderr,"qmin1=%e,qmin2=%e\n",qmin1,qmin2);

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
  radon_param(fmax1,h,nh,dx_av,qmin1,&qmaxt,&qmax,&dq1,nq1,rtmethod1,factor1);
    
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax1,dq1);

  if (symmetricq1){
    for (iq=0;iq<nq1/2;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
    //fprintf(stderr,"q[%d]=%f\n",nq1/2-1,q1[nq1/2-1]);
    q1[nq1/2]=-q1[nq1/2-1];
    fprintf(stderr,"q[%d]=%f\n",nq1/2,q1[nq1/2]);
    for (iq=nq1/2+1;iq<nq1;iq++){
      q1[iq]=q1[iq-1]+dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
  }
  else{
    //qmin2=q1[nq1-1]+dq1;
    for (iq=0;iq<nq1;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q1[%d]=%f\n",iq,q1[iq]);
    }
  }
  
  radon_param(fmax2,h,nh,dx_av,qmin2,&qmaxt,&qmax,&dq2,nq1,rtmethod2,factor2);

  for (iq=0;iq<nq2;iq++) q2[iq]=qmin2+iq*dq2;

  for (iq=0;iq<nq1;iq++) q[iq]=q1[iq];
  for (iq=0;iq<nq2;iq++){
    q[iq+nq1]=q2[iq];
    fprintf(stderr,"q2[%d]=%e\n",iq,q2[iq]);
  }
  *pdq1=dq1;
  *pdq2=dq2;

  return;
}










#include "su.h"

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t)
{
	static int indx=0;
	int it;
	float a1,a2;

	/* if before first cdp, constant extrapolate */
	if (cdpt<=cdp[0]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[0][it];
			oa1t[it] = oa1[0][it];
			oa2t[it] = oa2[0][it];
		      };
	
	/* else if beyond last cdp, constant extrapolate */
	} else if (cdpt>=cdp[ncdp-1]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[ncdp-1][it];
			oa1t[it] = oa1[ncdp-1][it];
			oa2t[it] = oa2[ncdp-1][it];
		      };
	
	/* else, linearly interpolate */
	} else {
		xindex(ncdp,cdp,cdpt,&indx);
		a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
		a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
		for (it=0; it<nt; ++it) {
			ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
			oa1t[it] = a1*oa1[indx][it]+a2*oa1[indx+1][it];
			oa2t[it] = a1*oa2[indx][it]+a2*oa2[indx+1][it];
		      };
	}
}

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt)
{
  static int indx=0;
  int it;
  float a1,a2;
  
  /* if before first cdp, constant extrapolate */
  if (cdpt<=cdp[0]) for (it=0; it<nt; ++it) ovvt[it] = ovv[0][it];
  /* else if beyond last cdp, constant extrapolate */
  else if (cdpt>=cdp[ncdp-1]) for (it=0; it<nt; ++it) ovvt[it] = ovv[ncdp-1][it];
  /* else, linearly interpolate */
  else {
    xindex(ncdp,cdp,cdpt,&indx);
    a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
    a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
    for (it=0; it<nt; ++it) ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
  }
}











/*
  It computes the term Cm for the model weight 
  for the Hessian (L' Cd^{-1} L+ Cm^{-1})
  This term corresponds to the probability model, so that
  the distribution parameters sigma and norm are passed.
  Wm is a vector of dimension nx, but in fact is the diagonal 
  of the Wm matrix of size nx x nx.
  Inqut 
          m: model
          nx: number of model traces
          norm: implemented 1 Huber, 0 Cauchy, else L2
  Output
          Wm 
          eps1: standard deviation of the model

  Daniel Trad- 14 March 2000. UBC- Canada
  Based in Sacchi, 1996. phD thesis. UBC. Canada
*/
#include "su.h"
#include "math.h"
#include "Complex.h"

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=MAX(fabs(m[i]),sigmam);
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}
void weights_cgfft(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=sqrt(MAX(fabs(m[i]),sigmam));
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=sqrt(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter)
{ 
      int i;
      
      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
      
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=1./sqrt(abs(m[i])*sigmam);
      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=1./sqrt(sigmam*sigmam+abs(m[i]*m[i]));
      else if(norm==2){ // Mask
	for (i=0;i<nx;i++)
	  Wm[i]=1+200./(1+exp(1*(abs(m[i])-sigmam)+0.5));
      }
      return;
}

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter)
{ 
  /* 
     It computes the model weights using a window in the model space m(f,q).
     For example, Wm[iq] is a function of m(f,iq), m(f-1,iq), ..., m(f-buffer,iq)
     Particularly useful for the dealiased RT (Hermman et al.)

  */     
      int  iq, iw;
      float maveg;

      if (iter==1){ 
	for (iq=0;iq<nq;iq++) Wm[iq]=1;
	return;
      }
      
      if (norm==1){ 
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(fabs(maveg)*sigmam);
	}
      }
      else if(norm==0){
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(sigmam*sigmam+fabs(maveg*maveg));
	}
      }
      else if(norm==2){ // Mask
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1+200./(1+exp(1*(fabs(maveg)-sigmam)+0.5));
	}
      }
      return;
}

void modelweight(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=fabs(m[isamax(nx,m,1)]);
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=fabs(m[i]);
      }
      else if(norm==0){
	if (maxm>1e-4) 
	  for (i=0;i<nx;i++){
	    // Solved !!!!!!!!!!!!!!!!!
	    // The right Wm from Cauchy is 
	    // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // But if M^-1 ATA x = M^-1 AT b is solved instead
	    // of the satndard form M=WmT *Wm 
	    Wm[i]=(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // Actually it works even better with (I don't know why)
	    //Wm[i]=Wm[i]*Wm[i];
	    //if (Wm[i]>2) Wm[i]=2; 
	  }
	else for (i=0;i<nx;i++) Wm[i]=1e-3;
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}

void modelweight_inv(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=m[isamax(nx,m,1)];
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=1./MIN(fabs(m[i]),eps1);
      }
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
 	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=1./(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}





void modelweight(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
            
      if (norm==1)
	for (i=0;i<nx;i++) Wm[i]=abs(m[i]);

      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=eps1+pow(abs(m[i]),2.0);

      else for (i=0;i<nx;i++) Wm[i]=1.;     


      return;
}

void modelweight_inv(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
  float *mabs;
  int ix;
  float sigma;

  mabs=ealloc1float(nx);
  if (norm==0 || norm ==1){ 
    for (ix=0;ix<nx;ix++) mabs[ix]=abs(m[ix]);
    sigma=quest(eps1,nx,mabs);
    //sigma=MAX(sigma,1e-3);
    //sigma=eps1;
    if (norm==1) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(MAX(mabs[ix],sigma)))+1e-7;
    else if(norm==0) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(sigma*sigma+pow(mabs[ix],2.0)))+1e-7;
  }
  else for (ix=0;ix<nx;ix++) Wm[ix]=1.;     

  free1float(mabs);
  return;
}



void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=1;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=fabs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=fabs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);
      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);
      free1float(mabs);
      free1float(dabs);

      return;
}

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=0;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=abs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=abs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);

      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);

      free1float(mabs);
      free1float(dabs);

      return;
}


void dataweigths(float *pos, int nh, float *Wd, int add)
{
  int ih;
  float *dh;
  float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
  int verbose=0;

  dh=ealloc1float(nh);

  for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   

  dh[0]=pos[1]-pos[0];
  dh[nh-1]=pos[nh-1]-pos[nh-2];
  
  for (ih=0;ih<nh;ih++){
    if (add) Wd[ih]*=MIN(dh[ih],dhav);
    else  Wd[ih]=MIN(dh[ih],dhav);
    if (verbose) fprintf(stderr," Wd[%d]=%f\n",ih,Wd[ih]);
  }

  free1float(dh);
  return;

}









































//#include <math.h>
#include "su.h"

void nmo(float *d,float *m,float *t,float h,float *vel,int invert,int nt,float dt,float smute) 
{
  register int it;
  float moveout;
  float t2;
  float vel2;
  float dt2=dt*dt;
  float *ttn;
  float *tnt;
  //  float smute=1.5; /* zero samples with NMO stretch exceeding smute */
  float osmute;	   /* 1/smute */
  int lmute=25;    /* length in samples of linear ramp for mute */
  int itmute=0;	   /* zero samples with indices less than itmute */
  float *atn;	   /* amplitude a(tn) for NMO */
  float *at;	   /* amplitude a(t) for inverse NMO */
  int sscale=0;	   /* if non-zero, apply NMO stretch scaling */
  ttn=ealloc1float(nt);
  at = ealloc1float(nt);
  atn = ealloc1float(nt);
  tnt=ealloc1float(nt);

  if (invert) for(it=0;it<nt;it++) d[it]=0;
  else for (it=0;it<nt;it++) m[it]=0;

  for (it=0;it<nt;it++){
    vel2=vel[it]*vel[it];
    t2=t[it]*t[it];
    moveout=h*h/vel2;
    ttn[it]=sqrt(t2/dt2+moveout/dt2);
  }
  /* compute inverse of stretch factor a(tn) */
  atn[0] = ttn[1]-ttn[0];
  for (it=1; it<nt; ++it)
    atn[it] = ttn[it]-ttn[it-1];
  //fprintf(stderr,"smute=%f\n",smute);
  /* determine index of first sample to survive mute */
  osmute = 1.0/smute;
  for (it=0,itmute=0; it<nt && atn[it]<osmute; ++it)
    itmute++;
  
  if (invert){ 
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    /* adjust mute time */
    itmute = (int) (1.0+ttn[itmute]);
    itmute = MIN(nt-2,itmute);
    
    /* compute a(t) */
    if (sscale) {
      for (it=itmute+1; it<nt; ++it)
	at[it] = tnt[it]-tnt[it-1];
      at[itmute] = at[itmute+1];
    }

    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      d[it] = 0.0;
			
    /* if specified, undo NMO stretch factor scaling */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	d[it] *= at[it];


  }
  else{
    ints8r(nt,1.0,0,d,0.0,0.0,nt,ttn,m);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      m[it] = 0.0;
    
    /* apply linear ramp */
    for (it=itmute; it<itmute+lmute && it<nt; ++it)
      m[it] *= (float)(it-itmute+1)/(float)lmute;
    
    /* if specified, scale by the NMO stretch factor */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	m[it] *= atn[it];
  }

  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}










#include "su.h"

void radon_matrix(complex *R, complex **l,complex **lh,float *g,float *q,int nh,int nq,float w,float *dh)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];

        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;

        for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      lh[iq][ih]=(1./Aperture)*dh[ih]*co;
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
	for (iq=0;iq<nq;iq++){
	  R[iq].r=0;
	  R[iq].i=0;
	  for (ih=0;ih<nh;ih++)
	    R[iq]+=lh[0][ih]*l[ih][iq]; //Top row of LL=LH*L
	  //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
	}

        return;
}

void radon_matrix(complex **l, float *g,float *q,int nh,int nq,float w)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];


	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
        return;
}

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth)
{
  int ih;
  
  if (rtmethod==1) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih];
  else if (rtmethod==2) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih]*h[ih];
  else if (rtmethod==3) 
    for (ih=0;ih<nh;ih++) g[ih]=sqrt(h[ih]*h[ih]+depth*depth)-depth;   

  return;

}


void radon_matrix_irrq(complex **L, float *h, float *q,int nh,int nq,float w)
{
        
        register int ih;
	register int iq;  
        complex  arg;
	arg.r=0;
        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.i=w*h[ih]*q[iq];
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}


void radon_matrix_2op(complex **L, complex **L1, complex **L2, int nh, int nq1, int nq2)
{
  int ih, iq;

  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq1;iq++)
      L[ih][iq]=L1[ih][iq];
    for (iq=0;iq<nq2;iq++) 
      L[ih][iq+nq1]=L2[ih][iq];
  }	


  return;
}





#include "su.h"
void interval(float *x,int lx,float *pmx, float *pax)
{
   int i;
   float dx, mx, ax;		
   mx=0;
   ax=0;	
     for (i=1;i<lx;i++){
	 dx=fabs(x[i]-x[i-1]);
	 if (dx>mx) 
	        mx=dx;
	 ax=ax+dx;
     }
     ax=ax/(lx-1);
     *pmx=mx;
     *pax=ax;
     return;
}
#include "su.h"
int taper(float **data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;
  if (ntaper<=0) return EXIT_SUCCESS;
  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    //taper[k] = s*s;
    taper[k]= pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 

  if(flag==4){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  if(flag==5){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih][it]*=scalemute;
      }
    } 
  
  
  return EXIT_SUCCESS;
}

int taper(float *data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;

  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    taper[k] = pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih*nt+it]*=scalemute;
      }
    } 

  if(flag==5)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=0;   
  
  return EXIT_SUCCESS;
}


int taper(float *data, int nh, int ntaper)
{
  float *taper;
  int k, ih;
  float s;
  
  taper=ealloc1float(ntaper);

  for (k=0;k<ntaper;++k){ 
    s=sin(k*PI/(2*ntaper));
    taper[k]=s;
  }  
  taper[0]=1e-5; 
  /* Taper at the left end of the data set */
  for (ih = 0; ih < ntaper; ++ih) data[ih] /= taper[ih];

  /* Taper at the right end of the data set */
  for (ih=nh - ntaper ; ih < nh; ++ih)  data[ih] /= taper[nh - ih - 1];

  return EXIT_SUCCESS;

}










#include "su.h"
#include "clibrary.h"
#include <math.h>

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	    float *Wd,int nh,int nq, float tol, float step, int itercg)
{
  /* This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
     Notice that LH=FH WdT and L= Wd F
     so that Wd only needs to be used the first time 
  */
 
  float normb,dq,dq2,nit,beta, alpha, alphanum, alphaden;
  int j,in,num;
  complex czero;

  register int i;
  complex *q;
  complex *q1;
  complex *s;
  complex *x1;
  complex *z;
  complex *z1;
  complex *r;
  complex *Az;
  complex *xold;
  float *eta;
  float *rho;
  float *gcv;
  float J;
  complex *r2;
  complex *xtemp;

  q=ealloc1complex(nq); 
  q1=ealloc1complex(nq);
  s=ealloc1complex(nq);
  x1=ealloc1complex(nq);
  z=ealloc1complex(nq);
  z1=ealloc1complex(nq);
  r=ealloc1complex(nh);
  Az=ealloc1complex(nh);
  eta=ealloc1float(nq);
  rho=ealloc1float(nq);
  gcv=ealloc1float(nq);
  r2=ealloc1complex(nh);
  xtemp=ealloc1complex(nq);
  xold=ealloc1complex(nq);

  czero.r=czero.i=0;
  for (i=0;i<nq;i++) x[i]=czero;
  normb=sqrt(rcdot(nh,b,b));
  //xequaly(r,b,nh);
  for (i=0;i<nh;i++) r[i]=Wd[i]*b[i];
  for (i=0;i<nh;i++) r2[i]=r[i]*Wd[i];  
  Atimesx(r2,L,s,nh,nq,1);
  
  nit=MIN(itercg,nq);
  for(i=0;i<nq;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  xequaly(z,q,nq);
  dq=rcdot(nq,s,q);
  xequaly(z1,q1,nq);
  for(i=0;i<nq;i++) x1[i]=czero;       
  for (j=0;j<nit;j++){
    Atimesx(Az,L,z,nh,nq,0);            
    for (i=0;i<nh;i++) Az[i]*=Wd[i];  
    alphanum=dq;
    alphaden=rcdot(nh,Az,Az);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1.e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e\n",
	      alphanum,alphaden);
      break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
         
    //// Update model u and residuals

    for(i=0;i<nh;i++) r[i]-=alpha*Az[i];
    for(i=0;i<nq;i++) xold[i]=x[i];    
    for(i=0;i<nq;i++) x[i]+=alpha*z[i];    

    for(i=0;i<nh;i++) r2[i]=r[i]*Wd[i];
    Atimesx(r2,L,s,nh,nq,1);

    for(i=0;i<nq;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=rcdot(nq,s,q);
    beta=dq2/dq;
    dq=dq2;
    for (i=0;i<nq;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(rcdot(nq,s,s))/normb;
    //fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);
    for (i=0;i<nq;i++) xtemp[i]=x[i]/Wm[i];
    J=rcdot(nh,r2,r2)+rcdot(nq,xtemp,xtemp);
    for (i=0;i<nq;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(rcdot(nq,x1,x1));
    if ((tol==0) && (j>1)){ // GCV criteria
       in = j;
       for (i=0;i<in;i++){
       num=(nq-i)*(nq-i); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if ((gcv[j-2]<gcv[j-1])&&(gcv[j-1]<gcv[j])) { 
         if (1) fprintf(stderr,"GCV Criteria, iteration %d\n",j-1);
	 for(i=0;i<nq;i++) x[i]=xold[i];    
         nit = j-1;
         break;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        if (1) fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        break;
       }
    }          
  }
  if (0) fprintf(stderr,"j=%d\n",j);

  free1complex(xold);
  free1complex(xtemp);
  free1complex(r2);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1complex(Az);
  free1complex(r);
  free1complex(z1);
  free1complex(z);
  free1complex(x1);
  free1complex(s);
  free1complex(q1);
  free1complex(q);

 

  return(J);
}















#include "su.h"

 
void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv)
{
  /* This function reads the par file in the command line and computes the vector
  ovv. This vecotr can then be used with interpovv to get the interpolated velocities at the desired cdp
 input :
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 

The calling function requires this:

  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv);

  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  free1float(cdp);
  free2float(ovv);
  free1float(velint);

  */
  int it;
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////

  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  //fprintf(stderr,"countparname=%d\n",countparname("vnmo"));
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
  }
  
  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
  } 
}





void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1,
		   float **oa2)
{
  /* This function reads the par file in the command line and computes the vectors
  ovv, oa1, and oa2. These 3 vecotrs can then be sued with interpovv to get the interpolated velocities at the desired cdp
 input :
    float **oa1;	 array[ncdp][nt] of anis1 functions 
    float **oa2;	 array[ncdp][nt] of anis2 functions 
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float *velint; array[nt] of vel for a particular trace 
    float *oa1t;	 array[nt] of anis1 for a particular trace 
    float *oa2t;	 array[nt] of anis2 for a particular trace 

The calling function requires this:

  float **oa1;
  float **oa2;
  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;
  float *oa1t;	
  float *oa2t;	

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv,oa1,oa2);

  interpovv(nt,ncdp,cdp,ovv,oa1,oa2,cdpgather,velint,oa1t,oa2t);

  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(velint);

  */

  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  int i, it;

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////


  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
    if (countparname("anis1")!=ncdp &&
	countparname("anis1")!=0)
      err("an anis1 array must be specified for each cdp, "
	  "or omitted at all");
    if (countparname("anis2")!=ncdp &&
	countparname("anis2")!=0)
      err("an anis2 array must be specified for each cdp, "
	  "or omitted at all");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
    if (countparname("anis1")>1)
      err("only one (or no) anis1 array must be specified");
    if (countparname("anis2")>1)
      err("only one (or no) anis2 array must be specified");    
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    nanis1 = countnparval(icdp+1,"anis1");
    nanis2 = countnparval(icdp+1,"anis2");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nanis1!=nvnmo && nanis1 != 0)
      err("number of vnmo and anis1 values must be equal");
    if (nanis2!=nvnmo && nanis2 != 0)
      err("number of vnmo and anis2 values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    if (nanis1==0) nanis1 = nvnmo;
    if (nanis2==0) nanis2 = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    anis1 = ealloc1float(nvnmo);
    anis2 = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    if (!getnparfloat(icdp+1,"anis1",anis1)) 
      for (i=0; i<nvnmo; i++) anis1[i] = 0.0;
    if (!getnparfloat(icdp+1,"anis2",anis2))
      for (i=0; i<nvnmo; i++) anis2[i] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis1,anis1[0],anis1[nanis1-1],1,&tn,&oa1[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis2,anis2[0],anis2[nanis2-1],1,&tn,&oa2[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
    free1float(anis1);
    free1float(anis2);
  }

  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    aoa1 = oa1[jcdp];
    aoa2 = oa2[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
      oa1[icdp+1] = oa1[icdp];
      oa2[icdp+1] = oa2[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
    oa1[icdp+1] = aoa1;
    oa2[icdp+1] = aoa2;
  } 

}








#include "su.h"
#include "radonhybrid.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, inv_par inv, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps, int symmetric)
{
  int it, ih;
  float *dtemp;
  float *Wd;
  int testadj=0;
  float dq1=0;
  float dq2=0;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2; 
  TRACE;
  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);
  
  fprintf(stderr,"fmax1=%f,fmax2=%f\n",fmax1,fmax2);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  }

  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  
  

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);
  
  dataweigths(h,nh,Wd,TRUE);
  TRACE;
  radon_param_2op(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2, nq1, nq2, depth1, depth2,
		   rtmethod1, rtmethod2, factor1, factor2,&dq1,&dq2, symmetric);
  TRACE;

  radon_wtcgls_2op(data,ph1,ph2,nh,t,nt,dt,model,q,nq,q1,nq1,q2,nq2,inv,Wd,testadj,fmax1,fmax2);
  TRACE;
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);  

  radoninv_2op(data,ph1,ph2,nh,t,nt,model,q,nq,q1,nq1,q2,nq2,fmax1,fmax2,filtering,npoly,ffilter,amps);
  
  /////////////////////////////

  plotgather(data[0],nt,nh,"data_before_inmo");

  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }    
  }

  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(Wd);
  free1float(dtemp);

  return;

}
























#include "su.h"
#include "radonhybrid.h"


/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/


void radon_wtcgls_2op(float **data, float *ph1, float *ph2, int nh, float *t, int nt, float dt, float **model, float *q, int nq, float *q1, int nq1, float *q2, int nq2, inv_par inv, float *Wd, int testadj, float fmax1, float fmax2)
{
  int iq,iter;
  complex czero;  czero.r=czero.i=0;
  float w=0,df;
  int freq, nf,nfft;
  int maxfreq1, maxfreq2;
  complex **m2;
  complex **d2;
  complex **L;
  complex **L1;
  complex **L2;
  float *Wm;
  float  *Jtot;
  float fmin=0; // Minimum freq to compute in Hz;
  float sigmam;
  float sigmad;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;

  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)
  Jtot=ealloc1float(20);
  Wm=ealloc1float(nq);

  //for (iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);

  zero_vector(Jtot,inv.iter_end);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
  
  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);
  

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
    
     //if (testadj) test=testadj_rad_f(L,LH);
    
    for (iter=1;iter<=inv.iter_end;iter++){
      weights_inv(m2[freq],nq,inv.norm,sigmam,Wm,iter);
      if (iter==1 && freq < 10 ){
	for (iq=0;iq<nq1;iq++) Wm[iq]*=2;
	for (iq=nq1;iq<nq;iq++) Wm[iq]*=0.5;
      }
      wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,0,inv.step,inv.itercg);
    }
    
    fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
  }
  //freqweighting(m2,nf,nq,df,fmin,fmax1);

  if (nq2>0){
    fprintf(stderr,"From now on one operator only\n");
    //////////////////////////////////////////////////////////////
    // For Ground roll,
    // Between maxfreq1 and freqmax2   we use only one operator
    complex **m2w=window(m2,0,nf-1,nq1,nq-1);
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      
      for (iter=1;iter<=inv.iter_end;iter++){
	if (iter==2) 
	  deviations(m2w[freq],nq2,d2[freq],nh,inv.norm,quantil1,quantil2,&sigmam,&sigmad);
	weights_inv(m2w[freq],nq2,inv.norm,sigmam,Wm,iter);   
	wtcgls(d2[freq],L2,m2w[freq],Wm,Wd,nh,nq2,0,inv.step,inv.itercg);
      }
      fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
    }
    freqweighting(m2w,nf,nq2,df,fmin,fmax2);
   }

   fprintf(stderr,"w=%f\n",w);
   
   
   for (iq=nq1;iq<nq;iq++) m2[freq][iq]/=nq;      
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  



   free1float(Wm);
   free1float(Jtot);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}








































#include "su.h"
#include "radonhybrid.h"

void radoninv_2op(float **data,  float *ph1, float *ph2, int nh, float *t, int nt, 
		  float **model, float *q,  int nq, float *q1, int nq1, float *q2, int nq2, 
		  float fmax1, float fmax2, int filter, int npoly, float *f, float *amps)
{
  int freq;
  complex **m2;
  complex **d2;
  complex czero; czero.r=czero.i=0;
  complex **L;
  complex **L1;
  complex **L2;
  float w,df;
  float dt=t[1]-t[0];
  int nfft;
  int nf;
  int maxfreq1, maxfreq2;
  float fmin=0;
  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh);  
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);

  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  if (filter==3) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);
  //if (filter==4) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
 
    Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
    fprintf(stderr,":");
    
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
  }
  if (nq2>0){ // Between maxfreq1 and freqmax2   we use only one operator
    fprintf(stderr,"From now on one operator only\n");
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      Atimesx(d2[freq],L2,&m2[freq][nq1],nh,nq2,FALSE);
      fprintf(stderr,".");      
    }
   }

  fprintf(stderr,"freq=%f\n",freq*df);      

  freqweighting(d2,nf,nh,df,fmin,fmax2);
  
  fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf); 

  TRACE;
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  TRACE;
  return;

}






























#include "su.h"
#include "stddef.h"

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod)
/* Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    Daniel Trad- UBC- 16-2-99
*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   if (rtmethod==2) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1) { //LRT
	  dq= 1/(fmax*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod,
float factor)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   
   if (rtmethod==2 || rtmethod==5) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1 || rtmethod==3) { //LRT
	  dq= 1/(fmax*fabs(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}
	















#include "su.h"
#include "stddef.h"



void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);
void interval(float *pos, int nh, float *dx_max, float *dx_av);


void radon_param_2op(float fmax1, float fmax2, float *h, int nh, float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1,
		      float factor2,float *pdq1, float *pdq2, int symmetricq1)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{	
  float qmaxt;
  float qmax;
  float dx_max;
  float dx_av;
  float dq1;
  float dq2;
  int iq;

  fprintf(stderr,"qmin1=%e,qmin2=%e\n",qmin1,qmin2);

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
  radon_param(fmax1,h,nh,dx_av,qmin1,&qmaxt,&qmax,&dq1,nq1,rtmethod1,factor1);
    
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax1,dq1);

  if (symmetricq1){
    for (iq=0;iq<nq1/2;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
    //fprintf(stderr,"q[%d]=%f\n",nq1/2-1,q1[nq1/2-1]);
    q1[nq1/2]=-q1[nq1/2-1];
    fprintf(stderr,"q[%d]=%f\n",nq1/2,q1[nq1/2]);
    for (iq=nq1/2+1;iq<nq1;iq++){
      q1[iq]=q1[iq-1]+dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
  }
  else{
    //qmin2=q1[nq1-1]+dq1;
    for (iq=0;iq<nq1;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q1[%d]=%f\n",iq,q1[iq]);
    }
  }
  
  radon_param(fmax2,h,nh,dx_av,qmin2,&qmaxt,&qmax,&dq2,nq1,rtmethod2,factor2);

  for (iq=0;iq<nq2;iq++) q2[iq]=qmin2+iq*dq2;

  for (iq=0;iq<nq1;iq++) q[iq]=q1[iq];
  for (iq=0;iq<nq2;iq++){
    q[iq+nq1]=q2[iq];
    fprintf(stderr,"q2[%d]=%e\n",iq,q2[iq]);
  }
  *pdq1=dq1;
  *pdq2=dq2;

  return;
}










#include "su.h"

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t)
{
	static int indx=0;
	int it;
	float a1,a2;

	/* if before first cdp, constant extrapolate */
	if (cdpt<=cdp[0]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[0][it];
			oa1t[it] = oa1[0][it];
			oa2t[it] = oa2[0][it];
		      };
	
	/* else if beyond last cdp, constant extrapolate */
	} else if (cdpt>=cdp[ncdp-1]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[ncdp-1][it];
			oa1t[it] = oa1[ncdp-1][it];
			oa2t[it] = oa2[ncdp-1][it];
		      };
	
	/* else, linearly interpolate */
	} else {
		xindex(ncdp,cdp,cdpt,&indx);
		a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
		a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
		for (it=0; it<nt; ++it) {
			ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
			oa1t[it] = a1*oa1[indx][it]+a2*oa1[indx+1][it];
			oa2t[it] = a1*oa2[indx][it]+a2*oa2[indx+1][it];
		      };
	}
}

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt)
{
  static int indx=0;
  int it;
  float a1,a2;
  
  /* if before first cdp, constant extrapolate */
  if (cdpt<=cdp[0]) for (it=0; it<nt; ++it) ovvt[it] = ovv[0][it];
  /* else if beyond last cdp, constant extrapolate */
  else if (cdpt>=cdp[ncdp-1]) for (it=0; it<nt; ++it) ovvt[it] = ovv[ncdp-1][it];
  /* else, linearly interpolate */
  else {
    xindex(ncdp,cdp,cdpt,&indx);
    a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
    a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
    for (it=0; it<nt; ++it) ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
  }
}











/*
  It computes the term Cm for the model weight 
  for the Hessian (L' Cd^{-1} L+ Cm^{-1})
  This term corresponds to the probability model, so that
  the distribution parameters sigma and norm are passed.
  Wm is a vector of dimension nx, but in fact is the diagonal 
  of the Wm matrix of size nx x nx.
  Inqut 
          m: model
          nx: number of model traces
          norm: implemented 1 Huber, 0 Cauchy, else L2
  Output
          Wm 
          eps1: standard deviation of the model

  Daniel Trad- 14 March 2000. UBC- Canada
  Based in Sacchi, 1996. phD thesis. UBC. Canada
*/
#include "su.h"
#include "math.h"
#include "Complex.h"

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=MAX(fabs(m[i]),sigmam);
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}
void weights_cgfft(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=sqrt(MAX(fabs(m[i]),sigmam));
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=sqrt(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter)
{ 
      int i;
      
      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
      
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=1./sqrt(abs(m[i])*sigmam);
      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=1./sqrt(sigmam*sigmam+abs(m[i]*m[i]));
      else if(norm==2){ // Mask
	for (i=0;i<nx;i++)
	  Wm[i]=1+200./(1+exp(1*(abs(m[i])-sigmam)+0.5));
      }
      return;
}

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter)
{ 
  /* 
     It computes the model weights using a window in the model space m(f,q).
     For example, Wm[iq] is a function of m(f,iq), m(f-1,iq), ..., m(f-buffer,iq)
     Particularly useful for the dealiased RT (Hermman et al.)

  */     
      int  iq, iw;
      float maveg;

      if (iter==1){ 
	for (iq=0;iq<nq;iq++) Wm[iq]=1;
	return;
      }
      
      if (norm==1){ 
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(fabs(maveg)*sigmam);
	}
      }
      else if(norm==0){
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(sigmam*sigmam+fabs(maveg*maveg));
	}
      }
      else if(norm==2){ // Mask
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1+200./(1+exp(1*(fabs(maveg)-sigmam)+0.5));
	}
      }
      return;
}

void modelweight(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=fabs(m[isamax(nx,m,1)]);
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=fabs(m[i]);
      }
      else if(norm==0){
	if (maxm>1e-4) 
	  for (i=0;i<nx;i++){
	    // Solved !!!!!!!!!!!!!!!!!
	    // The right Wm from Cauchy is 
	    // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // But if M^-1 ATA x = M^-1 AT b is solved instead
	    // of the satndard form M=WmT *Wm 
	    Wm[i]=(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // Actually it works even better with (I don't know why)
	    //Wm[i]=Wm[i]*Wm[i];
	    //if (Wm[i]>2) Wm[i]=2; 
	  }
	else for (i=0;i<nx;i++) Wm[i]=1e-3;
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}

void modelweight_inv(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=m[isamax(nx,m,1)];
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=1./MIN(fabs(m[i]),eps1);
      }
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
 	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=1./(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}





void modelweight(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
            
      if (norm==1)
	for (i=0;i<nx;i++) Wm[i]=abs(m[i]);

      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=eps1+pow(abs(m[i]),2.0);

      else for (i=0;i<nx;i++) Wm[i]=1.;     


      return;
}

void modelweight_inv(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
  float *mabs;
  int ix;
  float sigma;

  mabs=ealloc1float(nx);
  if (norm==0 || norm ==1){ 
    for (ix=0;ix<nx;ix++) mabs[ix]=abs(m[ix]);
    sigma=quest(eps1,nx,mabs);
    //sigma=MAX(sigma,1e-3);
    //sigma=eps1;
    if (norm==1) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(MAX(mabs[ix],sigma)))+1e-7;
    else if(norm==0) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(sigma*sigma+pow(mabs[ix],2.0)))+1e-7;
  }
  else for (ix=0;ix<nx;ix++) Wm[ix]=1.;     

  free1float(mabs);
  return;
}



void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=1;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=fabs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=fabs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);
      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);
      free1float(mabs);
      free1float(dabs);

      return;
}

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=0;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=abs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=abs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);

      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);

      free1float(mabs);
      free1float(dabs);

      return;
}


void dataweigths(float *pos, int nh, float *Wd, int add)
{
  int ih;
  float *dh;
  float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
  int verbose=0;

  dh=ealloc1float(nh);

  for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   

  dh[0]=pos[1]-pos[0];
  dh[nh-1]=pos[nh-1]-pos[nh-2];
  
  for (ih=0;ih<nh;ih++){
    if (add) Wd[ih]*=MIN(dh[ih],dhav);
    else  Wd[ih]=MIN(dh[ih],dhav);
    if (verbose) fprintf(stderr," Wd[%d]=%f\n",ih,Wd[ih]);
  }

  free1float(dh);
  return;

}









































//#include <math.h>
#include "su.h"

void nmo(float *d,float *m,float *t,float h,float *vel,int invert,int nt,float dt,float smute) 
{
  register int it;
  float moveout;
  float t2;
  float vel2;
  float dt2=dt*dt;
  float *ttn;
  float *tnt;
  //  float smute=1.5; /* zero samples with NMO stretch exceeding smute */
  float osmute;	   /* 1/smute */
  int lmute=25;    /* length in samples of linear ramp for mute */
  int itmute=0;	   /* zero samples with indices less than itmute */
  float *atn;	   /* amplitude a(tn) for NMO */
  float *at;	   /* amplitude a(t) for inverse NMO */
  int sscale=0;	   /* if non-zero, apply NMO stretch scaling */
  ttn=ealloc1float(nt);
  at = ealloc1float(nt);
  atn = ealloc1float(nt);
  tnt=ealloc1float(nt);

  if (invert) for(it=0;it<nt;it++) d[it]=0;
  else for (it=0;it<nt;it++) m[it]=0;

  for (it=0;it<nt;it++){
    vel2=vel[it]*vel[it];
    t2=t[it]*t[it];
    moveout=h*h/vel2;
    ttn[it]=sqrt(t2/dt2+moveout/dt2);
  }
  /* compute inverse of stretch factor a(tn) */
  atn[0] = ttn[1]-ttn[0];
  for (it=1; it<nt; ++it)
    atn[it] = ttn[it]-ttn[it-1];
  //fprintf(stderr,"smute=%f\n",smute);
  /* determine index of first sample to survive mute */
  osmute = 1.0/smute;
  for (it=0,itmute=0; it<nt && atn[it]<osmute; ++it)
    itmute++;
  
  if (invert){ 
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    /* adjust mute time */
    itmute = (int) (1.0+ttn[itmute]);
    itmute = MIN(nt-2,itmute);
    
    /* compute a(t) */
    if (sscale) {
      for (it=itmute+1; it<nt; ++it)
	at[it] = tnt[it]-tnt[it-1];
      at[itmute] = at[itmute+1];
    }

    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      d[it] = 0.0;
			
    /* if specified, undo NMO stretch factor scaling */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	d[it] *= at[it];


  }
  else{
    ints8r(nt,1.0,0,d,0.0,0.0,nt,ttn,m);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      m[it] = 0.0;
    
    /* apply linear ramp */
    for (it=itmute; it<itmute+lmute && it<nt; ++it)
      m[it] *= (float)(it-itmute+1)/(float)lmute;
    
    /* if specified, scale by the NMO stretch factor */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	m[it] *= atn[it];
  }

  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}










#include "su.h"

void radon_matrix(complex *R, complex **l,complex **lh,float *g,float *q,int nh,int nq,float w,float *dh)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];

        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;

        for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      lh[iq][ih]=(1./Aperture)*dh[ih]*co;
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
	for (iq=0;iq<nq;iq++){
	  R[iq].r=0;
	  R[iq].i=0;
	  for (ih=0;ih<nh;ih++)
	    R[iq]+=lh[0][ih]*l[ih][iq]; //Top row of LL=LH*L
	  //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
	}

        return;
}

void radon_matrix(complex **l, float *g,float *q,int nh,int nq,float w)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];


	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
        return;
}

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth)
{
  int ih;
  
  if (rtmethod==1) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih];
  else if (rtmethod==2) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih]*h[ih];
  else if (rtmethod==3) 
    for (ih=0;ih<nh;ih++) g[ih]=sqrt(h[ih]*h[ih]+depth*depth)-depth;   

  return;

}


void radon_matrix_irrq(complex **L, float *h, float *q,int nh,int nq,float w)
{
        
        register int ih;
	register int iq;  
        complex  arg;
	arg.r=0;
        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.i=w*h[ih]*q[iq];
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}


void radon_matrix_2op(complex **L, complex **L1, complex **L2, int nh, int nq1, int nq2)
{
  int ih, iq;

  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq1;iq++)
      L[ih][iq]=L1[ih][iq];
    for (iq=0;iq<nq2;iq++) 
      L[ih][iq+nq1]=L2[ih][iq];
  }	


  return;
}





#include "su.h"
void interval(float *x,int lx,float *pmx, float *pax)
{
   int i;
   float dx, mx, ax;		
   mx=0;
   ax=0;	
     for (i=1;i<lx;i++){
	 dx=fabs(x[i]-x[i-1]);
	 if (dx>mx) 
	        mx=dx;
	 ax=ax+dx;
     }
     ax=ax/(lx-1);
     *pmx=mx;
     *pax=ax;
     return;
}
#include "su.h"
int taper(float **data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;
  if (ntaper<=0) return EXIT_SUCCESS;
  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    //taper[k] = s*s;
    taper[k]= pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 

  if(flag==4){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  if(flag==5){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih][it]*=scalemute;
      }
    } 
  
  
  return EXIT_SUCCESS;
}

int taper(float *data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;

  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    taper[k] = pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih*nt+it]*=scalemute;
      }
    } 

  if(flag==5)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=0;   
  
  return EXIT_SUCCESS;
}


int taper(float *data, int nh, int ntaper)
{
  float *taper;
  int k, ih;
  float s;
  
  taper=ealloc1float(ntaper);

  for (k=0;k<ntaper;++k){ 
    s=sin(k*PI/(2*ntaper));
    taper[k]=s;
  }  
  taper[0]=1e-5; 
  /* Taper at the left end of the data set */
  for (ih = 0; ih < ntaper; ++ih) data[ih] /= taper[ih];

  /* Taper at the right end of the data set */
  for (ih=nh - ntaper ; ih < nh; ++ih)  data[ih] /= taper[nh - ih - 1];

  return EXIT_SUCCESS;

}










#include "su.h"
#include "clibrary.h"
#include <math.h>

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	    float *Wd,int nh,int nq, float tol, float step, int itercg)
{
  /* This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
     Notice that LH=FH WdT and L= Wd F
     so that Wd only needs to be used the first time 
  */
 
  float normb,dq,dq2,nit,beta, alpha, alphanum, alphaden;
  int j,in,num;
  complex czero;

  register int i;
  complex *q;
  complex *q1;
  complex *s;
  complex *x1;
  complex *z;
  complex *z1;
  complex *r;
  complex *Az;
  complex *xold;
  float *eta;
  float *rho;
  float *gcv;
  float J;
  complex *r2;
  complex *xtemp;

  q=ealloc1complex(nq); 
  q1=ealloc1complex(nq);
  s=ealloc1complex(nq);
  x1=ealloc1complex(nq);
  z=ealloc1complex(nq);
  z1=ealloc1complex(nq);
  r=ealloc1complex(nh);
  Az=ealloc1complex(nh);
  eta=ealloc1float(nq);
  rho=ealloc1float(nq);
  gcv=ealloc1float(nq);
  r2=ealloc1complex(nh);
  xtemp=ealloc1complex(nq);
  xold=ealloc1complex(nq);

  czero.r=czero.i=0;
  for (i=0;i<nq;i++) x[i]=czero;
  normb=sqrt(rcdot(nh,b,b));
  //xequaly(r,b,nh);
  for (i=0;i<nh;i++) r[i]=Wd[i]*b[i];
  for (i=0;i<nh;i++) r2[i]=r[i]*Wd[i];  
  Atimesx(r2,L,s,nh,nq,1);
  
  nit=MIN(itercg,nq);
  for(i=0;i<nq;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  xequaly(z,q,nq);
  dq=rcdot(nq,s,q);
  xequaly(z1,q1,nq);
  for(i=0;i<nq;i++) x1[i]=czero;       
  for (j=0;j<nit;j++){
    Atimesx(Az,L,z,nh,nq,0);            
    for (i=0;i<nh;i++) Az[i]*=Wd[i];  
    alphanum=dq;
    alphaden=rcdot(nh,Az,Az);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1.e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e\n",
	      alphanum,alphaden);
      break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
         
    //// Update model u and residuals

    for(i=0;i<nh;i++) r[i]-=alpha*Az[i];
    for(i=0;i<nq;i++) xold[i]=x[i];    
    for(i=0;i<nq;i++) x[i]+=alpha*z[i];    

    for(i=0;i<nh;i++) r2[i]=r[i]*Wd[i];
    Atimesx(r2,L,s,nh,nq,1);

    for(i=0;i<nq;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=rcdot(nq,s,q);
    beta=dq2/dq;
    dq=dq2;
    for (i=0;i<nq;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(rcdot(nq,s,s))/normb;
    //fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);
    for (i=0;i<nq;i++) xtemp[i]=x[i]/Wm[i];
    J=rcdot(nh,r2,r2)+rcdot(nq,xtemp,xtemp);
    for (i=0;i<nq;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(rcdot(nq,x1,x1));
    if ((tol==0) && (j>1)){ // GCV criteria
       in = j;
       for (i=0;i<in;i++){
       num=(nq-i)*(nq-i); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if ((gcv[j-2]<gcv[j-1])&&(gcv[j-1]<gcv[j])) { 
         if (1) fprintf(stderr,"GCV Criteria, iteration %d\n",j-1);
	 for(i=0;i<nq;i++) x[i]=xold[i];    
         nit = j-1;
         break;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        if (1) fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        break;
       }
    }          
  }
  if (0) fprintf(stderr,"j=%d\n",j);

  free1complex(xold);
  free1complex(xtemp);
  free1complex(r2);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1complex(Az);
  free1complex(r);
  free1complex(z1);
  free1complex(z);
  free1complex(x1);
  free1complex(s);
  free1complex(q1);
  free1complex(q);

 

  return(J);
}















#include "su.h"

 
void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv)
{
  /* This function reads the par file in the command line and computes the vector
  ovv. This vecotr can then be used with interpovv to get the interpolated velocities at the desired cdp
 input :
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 

The calling function requires this:

  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv);

  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  free1float(cdp);
  free2float(ovv);
  free1float(velint);

  */
  int it;
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////

  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  //fprintf(stderr,"countparname=%d\n",countparname("vnmo"));
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
  }
  
  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
  } 
}





void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1,
		   float **oa2)
{
  /* This function reads the par file in the command line and computes the vectors
  ovv, oa1, and oa2. These 3 vecotrs can then be sued with interpovv to get the interpolated velocities at the desired cdp
 input :
    float **oa1;	 array[ncdp][nt] of anis1 functions 
    float **oa2;	 array[ncdp][nt] of anis2 functions 
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float *velint; array[nt] of vel for a particular trace 
    float *oa1t;	 array[nt] of anis1 for a particular trace 
    float *oa2t;	 array[nt] of anis2 for a particular trace 

The calling function requires this:

  float **oa1;
  float **oa2;
  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;
  float *oa1t;	
  float *oa2t;	

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv,oa1,oa2);

  interpovv(nt,ncdp,cdp,ovv,oa1,oa2,cdpgather,velint,oa1t,oa2t);

  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(velint);

  */

  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  int i, it;

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////


  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
    if (countparname("anis1")!=ncdp &&
	countparname("anis1")!=0)
      err("an anis1 array must be specified for each cdp, "
	  "or omitted at all");
    if (countparname("anis2")!=ncdp &&
	countparname("anis2")!=0)
      err("an anis2 array must be specified for each cdp, "
	  "or omitted at all");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
    if (countparname("anis1")>1)
      err("only one (or no) anis1 array must be specified");
    if (countparname("anis2")>1)
      err("only one (or no) anis2 array must be specified");    
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    nanis1 = countnparval(icdp+1,"anis1");
    nanis2 = countnparval(icdp+1,"anis2");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nanis1!=nvnmo && nanis1 != 0)
      err("number of vnmo and anis1 values must be equal");
    if (nanis2!=nvnmo && nanis2 != 0)
      err("number of vnmo and anis2 values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    if (nanis1==0) nanis1 = nvnmo;
    if (nanis2==0) nanis2 = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    anis1 = ealloc1float(nvnmo);
    anis2 = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    if (!getnparfloat(icdp+1,"anis1",anis1)) 
      for (i=0; i<nvnmo; i++) anis1[i] = 0.0;
    if (!getnparfloat(icdp+1,"anis2",anis2))
      for (i=0; i<nvnmo; i++) anis2[i] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis1,anis1[0],anis1[nanis1-1],1,&tn,&oa1[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis2,anis2[0],anis2[nanis2-1],1,&tn,&oa2[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
    free1float(anis1);
    free1float(anis2);
  }

  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    aoa1 = oa1[jcdp];
    aoa2 = oa2[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
      oa1[icdp+1] = oa1[icdp];
      oa2[icdp+1] = oa2[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
    oa1[icdp+1] = aoa1;
    oa2[icdp+1] = aoa2;
  } 

}








