/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON00FORW  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "dan.h"
#include <signal.h>
#include <time.h>
#include <math.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUFFT2 - 2D FFT                                                     ",
  " 	   								",
  " sufft2  < stdin > stdout                             		",
  "                                                                     ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " offsetfile      ascii file that contains the desired output offset  ",
  "                 IMPORTANT: nxe has to be larger than the number of  ",
  "                 new offsets.                                        ",
  "                                                                     ", 
  " nk=80           Size of Radon EOM                                   ",
  " testadj=0       =1 Test adjoint with random numbers                 ",
  " smooth=0        =1 Pass triangular filter to smooth  output         ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " verbose=0       =1  Extra Information                               ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1.9   nmofactor * equiv offset is used for nmo		",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset            */
/**************** end self doc ***********************************/
/* Globals (so can trap signal) defining temporary disk files */
FILE *modelfilep;       /* fp for model file output  		*/

///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr;
  cwp_String modelfile=""; /* output sufile for the model */ 
  //cwp_String offsetfile=""; /* output sufile for the model */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nx;   // number of input offset traces
  int nt;   // number of time samples 
  float dt; // Time sampling
  int  ih, ik, iw, it; // General counters 
  float **dataxt;   /*  Input gather */
  complex **datakf;   /* Output gather */
  float  *t;     // time axis for input and output
  float  *h=0;      // input offset
  float dx; 
  //float t0;      // First useful time

  int verbose;		/* flag for echoing info		*/
  float *k;
  float dk;
  int nk;  
  //float fk;

  float *w;
  float dw;
  //float fw;
  int nw;  
  int ntr;
  //int ntfft;
  //int nxfft;

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

  nx=ntr;

  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);
  
  fksize(nt,nx,&nw,&nk);

  fprintf(stderr,"nw=%d,nk=%d\n",nw,nk);
  // Allocate memory for data and model
  dataxt=ealloc2float(nt,nx);
  datakf=alloc2complex(nw,nk);
  t=ealloc1float(nt);
  w=ealloc1float(nw);
  k=ealloc1float(nk);
  h=ealloc1float(nx);

  for(it=0;it<nt;it++) t[it]=0+it*dt;
  memset( (void *) dataxt[0], (int) '\0', nx * nt *FSIZE);

  ih = 0;
  do {
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    memcpy((void *) dataxt[ih],(const void *) tr.data,nt*sizeof(float));
    h[ih]=tr.offset;
    ih++;
  } while (gettr(&tr));
  dx=fabs(h[1]-h[0]);
  erewind(stdin);
  nx=ih;
  if (ntr!=nx) fprintf(stderr,"ntr=%d is different from nx=%d \n",ntr,nx);
  if (nx>ntr) err("Change ntr in header\n");

  fft2_parameters(nt,dt,nx,dx,&nw,&dw,&nk,&dk,w,k);
  fft2_xt2kf(dataxt, datakf,nt,nx,nw,nk);
  /* Apply any desired filter here*************/
  // example
  //zero_vector(datakf[(int) (nk/2)],nw);
  /********************************************/
  fft2_kf2xt(dataxt, datakf,nt,nx,nw,nk);

  for (ih=0;ih<nx;ih++){ 
    gettr(&tr);
    memcpy((void *) tr.data,(const void *) dataxt[ih],nt*sizeof(float));
    fputtr(modelfilep,&tr);    
  }

  for (ik=0;ik<nk;ik++){
    //memcpy((void *) tr.data,(const void *) model[ik],nt*sizeof(float));
    for (iw=0;iw<nw;iw++) tr.data[iw]=abs(datakf[ik][iw]);
    /* set header values */
    tr.tracl = ik + 1;
    tr.ns = nw;
    tr.ntr= nk;
    tr.dt = 0;  /* d1 is now the relevant step size */
    tr.trid = KOMEGA;
    tr.d1 = dw;
    tr.f1 = w[0];
    tr.d2 = dk;
    tr.f2 = k[0];		
    puttr(&tr); 
  }

  free1float(w);
  free1float(k);
  free1float(t);
  free1float(h);
  free2float(dataxt);
  free2complex(datakf);
  efclose(modelfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}






























