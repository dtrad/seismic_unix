/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADTD:  $Date: June 1999  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrarytd.h"
#include <time.h>
#include "dan.h"
#include "radonl1.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " suradonl1   Forward High Resolution Hyperbolic Radon transform      ", 
  "             using linear programming with interior point method 	",
  " 									",
  " suradonl1 < stdin > stdout modelfile= [optional parameters] 	",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=1            =0 test for the adjoint                     	",
  "                     =1 use slowness                                 ",
  "                     =1 use slowness                                 ",
  "                     =2 use Velocities                               ",
  " itercg = 200	number of internal iterations  		        ",
  " iter_end = 4	number of external iterations          	        ",
  " nq=nh		Number of Radon Traces                         	",
  " qmin=0              minimum velocity or slownes                     ",
  " qmax=0              maximum velocity or slowness                    ",
  " modelfile           name for output model file                      ", 
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : offset time domain and model file                          ",
  " Input : sudata file  (offset time domain)              		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradonl1 < data.su > data.rec modelfile=data.surad                 ", 
  "                                                                     ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/

segy tr,tri; 
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for header storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	
  FILE *myfilep; 
  cwp_String modelfile=""; /* output sufile for the model */  
  time_t start,finish;
  double elapsed_time;
  int iq,ih,ntr;
  int it;
  float **d;
  float **m;
  float *q;
  float *t;
  float *h;
  float eps;
  float qmin;
  float qmax;
  float t0;
  int nt;
  int nh;
  int nq;
  int itercg;
  int iter_end;
  int method;
  int nx;
  int ny;
  float dt,dq;
  float eps1;
  float eps2;
  float step;
  int testadj;
  
  
  //////////////////////////////////////////////
  fprintf(stderr,"**************************\n");

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 1;
  if (!getparfloat("eps1", &eps1))  eps1 = 1;
  if (!getparfloat("eps2", &eps2))  eps2 = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("step", &step))  step = .9;
  if (!getparfloat("qmin", &qmin))  qmin = 0;
  if (!getparfloat("qmax", &qmax))  qmax = 0;
  if (!getparint("nq", &nq))  nq = 0;
  if (!getparint("itercg", &itercg))  itercg = 200;
  if (!getparint("testadj", &testadj))  testadj =0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparfloat("t0",&t0)) t0=0;

  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  if (!nq) nq=nh;

  //  fprintf(stderr,"nh=%d,nt=%d,nq=%d,dt=%f,eps2=%e\n",nh,nt,nq,dt,eps2); 

  // Allocate memory for data and model
  d=ealloc2float(nt,nh);
  m=ealloc2float(nt,nq);
  q=ealloc1float(nq); 
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  for (ih=0;ih<nh;ih++) h[ih]=0;

  headerfp = etmpfile();
  tracefp = etmpfile();
  if (verbose) warn("using tmpfile() call");

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;
  dq=(qmax-qmin)/(nq-1);
  for(iq=0;iq<nq;iq++) q[iq]=qmin+iq*dq;
  
  ///////////////////////////////////////////////////////////////////////

  memset( (void *) d[0], (int) '\0', nh * nt *FSIZE);
  ntr = 0;
  do {
    ntr++;
    //    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);
  nh=ntr;


  for (ih=0;ih<nh;ih++){
    efread(tr.data, FSIZE, nt, tracefp);
    memcpy((void *) d[ih],(const void *) tr.data,nt*sizeof(float));
    efread(&tr,HDRBYTES,1,headerfp);
    h[ih]=tr.offset;
  }

  ny=nt*nh;
  nx=nt*nq;

  memset( (void *) m[0], (int) '\0', nq * nt *FSIZE);
  radon_PD_lp_bm(d,m,t,q,h,nt,nq,nh,iter_end,itercg,method);

  //  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  
  if ((myfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);


  for (iq=0;iq<nq;iq++){
    tr.f2=q[iq];  // copy velocity perturbation
    tr.ntr=nq;
    tr.tracl=iq;
    tr.tracr=iq;
    //    fprintf(stderr,"iq=%d\n",iq);
    memcpy((void *) tr.data,(const void *) m[iq],nt*sizeof(float));   
    fputtr(myfilep,&tr);

  }
  

  erewind(headerfp);  
  for(ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) d[ih],nt*sizeof(float));
    puttr(&tr);
  }

  if (verbose) fprintf(stderr,"nt=%d, nh=%d \n ",nt,nh);

  //  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free1float(t);
  free1float(h);
  free1float(q);
  free2float(m);
  free2float(d);
  efclose(headerfp);
  efclose(tracefp);
  efclose(myfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}












































