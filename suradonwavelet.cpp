/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADONWAVELET  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonwavelet.h"
#include <signal.h>
#include <math.h>
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONWAVELET - Wavelet domain Hyperbolic Radon Transform          ",
  " 	   								",
  " suradonwavelet < stdin > stdout [optional parameters]      		",
  "                                                                     ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  " quantil=0.90      defines filter.                                   ", 
  " scalefilt=6     smallest scale to be filtered                       ",
  "                                                                     ",
  " nq=80           Size of Radon EOM                                   ",
  " verbose=0       =1  Extra Information                               ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, ntr, dt, offset            */
/**************** end self doc ***********************************/


int main(int argc, char **argv)
{
  segy tr;
  cwp_String modelfile=""; /* output sufile for the model */ 
  FILE *modelfilep;       /* fp for model file output  		*/
  ///////////////////////////////////////////////////////////////
  time_t start,finish;

  double elapsed_time;

  int  ih, iq; // General counters 
  register int it;

  float **data;   /* Common Scattering gather */
  float **model;            // Temporal array for the model
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float *q;

  int nh;   // number of offset traces
  int nt=0;   // number of time samples after zero padding
  int nt0;  // Original number of samples
  int nq;  

  float dt; // Time sampling
  float dq;
  float t0;      // First useful time

  float qmin;
  float qmax;
  int verbose;		/* flag for echoing info		*/

  int wavn;
  float quantil;        /* Defines filter */
  int scalefilt;        /* Scales with smaller resolution than this are filtered */

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
  // Get info from first trace 
  if (!gettr(&tr)) err("cannot read first trace"); 
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  
  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt0 = (int) tr.ns;
  if (!(nh=tr.ntr)) err("***ntr must be set\n");
  
  fprintf(stderr,"nt=%d,nh=%d,dt=%f\n",nt,nh,dt);  

  nt=findpower2(nt0);
  if ((nt<nt0)||(!nt)) err("nt is wrong. Check tr.ns or findpower2\n");
  

  fprintf(stderr,"Original nt0=%d, new nt=%d \n",nt0,nt);

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparfloat("qmax", &qmax))  qmax = 5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparint("wavn",&wavn)) wavn=4;
  if (!getparint("scalefilt",&scalefilt)) scalefilt=6;
  if (!getparfloat("quantil",&quantil)) quantil=0.90;

  modelfilep=efopen(modelfile,"w");
  

  // Allocate memory for data and model
  data=ealloc2float(nt,nh);
  model=alloc2float(nt,nq);
  t=ealloc1float(nt); 
  q=ealloc1float(nq);
  h=ealloc1float(nh);

  // Velocity axis for Radon space

  dq = (qmax-qmin)/(nq-1);
  for (iq=0;iq<nq;iq++) q[iq] = qmin+iq*dq;

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ih = 0;
  do {
    memcpy((void *) data[ih],(const void *) tr.data,nt0*sizeof(float));
    h[ih]=tr.offset;
    ih++;
  } while (gettr(&tr));
  erewind(stdin);

  if (nh != ih) err("ntr was wrong\n");

  fprintf(stderr,"nh=%d\n",nh);

  radonwavelet0(data,t,nt,h,nh,model,q,nq,wavn,quantil,scalefilt);

  // Output trace
  for (ih=0;ih<nh;ih++){ 
    gettr(&tr);
    memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
    puttr(&tr);
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


  free1float(q);
  free1float(t);
  free1float(h);
  free2float(data);
  free2float(model);

  efclose(modelfilep);

  finish=time(0);

  elapsed_time=difftime(finish,start);

  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}
































