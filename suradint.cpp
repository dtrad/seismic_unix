/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADTD:  $Date: June 1999  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrarytd.h"
#include <time.h>

/* Maximum value an `unsigned short int' can hold.  (Minimum is 0).  */
#ifndef USHRT_MAX
#define USHRT_MAX 65535
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADINT Forward  High Resolution Hyperbolic Radon transform        ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradint< stdin > stdout [optional parameters]          		",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  "             0  Variable velocity at every q trace                   ",
  "                and variable delta q - Output is offset-space        ",
  " rtmethod=3  1-LRT 2-PRT 3-HRT  (only for method 2 and 5)            ",
  "                                                                     ",
  " eps1 =1		data variance (for Weight function in WTCGLS)   ",
  " eps2 =1             model variance (for Weight function in WTCGLS)  ",
  " eps=1e-7            small number for conjugate gradient             ",
  " step=0.9            the step for CG is shortened by step            ",
  " itercg = 10	        number of internal iterations  		        ",
  " iter_end =1		number of external iterations          	        ",
  " pervmin = 10        minimum perturbation on velocity                ",
  " dperv=0.1           defines the rate of incresing velocity spacing  ",
  " nq=20		Number of Radon Traces                         	",
  " factor=0.8		multiply dq critical by factor                	",
  " precond=0           Uses precondtioning                             ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : Iterpolated traces  (offset time domain)                   ",
  " Input : sudata file  (offset time domain)              		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradint pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudataint  ", 
  "                                                                     ",
  " If nh output is > nh input only the header words of the first nh    ",
  "  traces are preserved                                               ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/

segy tr; 
float factor;
int nt,nh,nq,nx,ny,rtmethod,itercg,iter_end,method,norm, reorth;
int nh2;
float dt,dh,dq,eps1,eps2, step, thres,theta;
int testadj;
int taper;
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	
  FILE *myfilep; 
  cwp_String modelfile=""; /* output sufile for the model */  
  time_t start,finish;
  double elapsed_time;
  int j,i,iq,ih;
  int it;
  float *d, *m;
  float *q, *t, *h, eps, qmin, qmax, fmax;
  float t0=0;
  int smooth;
  int model; 
  extern int nt, nh, nq,rtmethod, norm, itercg, iter_end, method, nx ,ny;
  extern int reorth;
  extern float dt,dh,dq, eps1,eps2, step, thres, theta;
  extern float factor; // multiply dq by factor
  extern int testadj;
  extern int taper;
  int precond;
  // Interpolated data
  float *dint;
  extern int nh2;
  float dh2;
  float h2min;
  float h2max;
  float hmin;
  float hmax;
  float *h2;
  /// Velocity Trend
  float *tvel;
  float *vel;
  float *velint;
  int ntvel;
  int nvel;
  int itv;
  float pervmin;
  float dperv;
  int rhofilter;
  float alum,blum,clum; // Lumley preconditioner parameters
  // smoothing
  int nl=3;  //  npoints left hand side
  int nr=3;  //  npoints left hand side
  int flag=2;  // 1 rectangular, 2 triangular
  //////////////////////////////////////////////
  fprintf(stderr,"**************************\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("eps1", &eps1))  eps1 = 1;
  if (!getparfloat("eps2", &eps2))  eps2 = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("step", &step))  step = .9;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("itercg", &itercg))  itercg = 10;
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparfloat("fmax", &fmax))  fmax =70; // max freq used for dq
  if (!getparint("norm", &norm))  norm =1; 
  if (!getparfloat("theta", &theta))  theta =1;
  if (!getparint("smooth", &smooth))  smooth =0;
  if (!getparint("testadj", &testadj))  testadj =0;
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparint("taper", &taper))  taper =0;
  if (!getparint("model",&model)) model = 1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("precond",&precond)) precond=0;
  if (!getparint("rhofilter",&rhofilter)) rhofilter=0;
  if (!getparfloat("alum",&alum)) alum=0.5;
  if (!getparfloat("blum",&blum)) blum=1.;
  if (!getparfloat("clum",&clum)) clum=1.;

  /* Introduce velocity trend to apply Hyp Vel Filtering */
  ntvel = countparval("tvel");
  if (ntvel==0) ntvel = 1;
  tvel = ealloc1float(ntvel);
  if (!getparfloat("tvel",tvel)) tvel[0] = 0.0;
  nvel = countparval("vel");
  if (nvel==0) nvel = 1;
  if (nvel!=ntvel) err("number of tmig and vmig must be equal");
  vel = ealloc1float(nvel);
  if (!getparfloat("vel",vel)) vel[0] = 2000.0;
  for (itv=1; itv<ntvel; ++itv)
    if (tvel[itv]<=tvel[itv-1])
      err("tvel must increase monotonically");
  ////////////////////////////////////////////////////////   
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  fprintf(stderr,"nh=%d,nt=%d,dt=%f,eps2=%e\n",nh,nt,dt,eps2); 

  //if ((nt*nq) > USHRT_MAX) err("nt*nq > maximum unsigned short\n");
  //if ((nt*nh) > USHRT_MAX) err("nt*nh > maximum unsigned short\n");

  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  
  if ((m=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
  
  if ((q=ealloc1float(nq))==NULL)
    fprintf(stderr,"***Sorry, space for q could not be allocated\n");
  
  if ((h=ealloc1float(nh))==NULL)
    fprintf(stderr,"***Sorry, space for h could not be allocated\n");

  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");    

  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"*Sorry, space for velint could not be allocated\n"); 

      	
  for (ih=0;ih<nh;ih++) h[ih]=0;

  headerfp = etmpfile();
  if (verbose) warn("using tmpfile() call");  

  // Loop over traces 
  ih=0;
  hmin=1e10;
  hmax=1e-10;
  do {
    register int i;
    efwrite(&tr,HDRBYTES,1,headerfp);    
    h[ih]=(float) tr.offset;
    if (h[ih]<hmin) hmin=h[ih];
    if (h[ih]>hmax) hmax=h[ih]; 
    for (i=0;i<nt;i++){
      d[ih*nt+i]=(float) tr.data[i];    //if sort by time-offset 
    }
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(headerfp);
  nh=ih;

  // Define new offset axis
  // The default is to keep same number of traces 
  if (1){
    if (!getparfloat("h2min",&h2min)) h2min=hmin;
    if (!getparfloat("h2max",&h2max)) h2max=hmax; 
    if (!getparfloat("dh2",&dh2)) dh2=(h2max-h2min)/(nh-1);
    nh2=(int) ((h2max-h2min)/dh2 + 1 );
    if ((h2=alloc1float(nh2))==NULL) err("Cannot allocate h2\n");
    for (h2[0]=h2min,ih=1;ih<nh2;ih++) h2[ih]=h2[ih-1]+dh2;
  }
  ////////////////////////////////
  else{
  // Keep axis for residuals study
    h2min=hmin;
    h2max=hmax; 
    nh2=nh;
    if ((h2=alloc1float(nh2))==NULL) err("Cannot allocate h2\n");
    for (ih=0;ih<nh2;ih++) h2[ih]=h[ih];
  }

  fprintf(stderr,"nh2=%d,nt=%d,dh2=%f\n",nh2,nt,dh2);
  if ((dint=alloc1float(nt*nh2))==NULL) err("Cannot allocate dint\n");
  /* Time axis */
  for (i=0;i<nt;i++) t[i]=t0+i*dt;

  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);


  if (verbose) fprintf(stderr,"Before radtd nq=%d, nt=%d, nh=%d, eps=%f \n",nq,nt,nh,eps);  

  ny=nt*nh;
  nx=nt*nq;

  radint(t,q,h,h2,m,d,dint,eps,qmin,qmax,fmax,velint,dperv,pervmin,
	 alum,blum,clum,norm);

  if (rhofilter) rhofilt(dint,t,h2,nt,nh2,dt);
      
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  if (model){
    if (smooth) smoothing(m,nt,nq,nl,nr,flag);
    if ((myfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);
    iq=0;
    do{
      tr.offset=(int) q[iq];  // copy velocity perturbation
      tr.ntr=nq;    
      for (i=0;i<nt;i++)
	tr.data[i]=m[iq*nt+i];    
      fputtr(myfilep,&tr);
      iq++;
    } while(iq<nq);
  }
  
  if (smooth) smoothing(dint,nt,nh2,nl,nr,flag);
  ih=0;
  do{
    if (ih<nh)  efread(&tr,HDRBYTES,1,headerfp);
    tr.offset=(int) h2[ih];
    tr.ntr=nh2;
    tr.f2=0;
    for (i=0;i<nt;i++)
      tr.data[i]=dint[i+ih*nt];
    puttr(&tr);
    ih++;
  } while(ih<nh2);
  
  if (verbose) fprintf(stderr,"nt%d, nh2 %d \n ",nt,nh2);

  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free1float(dint);
  free1float(h2);
  free1float(velint);  
  free1float(t);
  free1float(h);
  free1float(q);
  free1float(m);
  free1float(d);
  free1float(vel);
  free1float(tvel);
  efclose(headerfp);
  efclose(myfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}






























