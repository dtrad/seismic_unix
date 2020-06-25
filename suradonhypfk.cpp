/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADONHYPFK:  $Date: June 1999 - Last version October 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>
#include "dan.h"
#include "inversion_par.h"
//#include "radonhypfk.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONHYPFK Forward  Hyperbolic Radon transform  in the f-k domain ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradonhypfk < stdin > stdout [optional parameters]          	",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Input : sudata file  (offset time domain)              		",
  " Output : adjoint model                                              ",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradon1 pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudata     ", 
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
//inv_par inv;




int verbose;
int radonhypfk(float *t, float *q, float *h, float **model, float **data,
	       int nt, int nh, int nq, float dt, float t0, inv_par inv, int plot, 
	       float dx, int adj);

int main(int argc, char **argv)
{
  segy tr; 
  inv_par inv;
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep; 
  time_t start,finish;
  double elapsed_time;
  int it,iq,ih;
  float vmin,dv;
  float t0=0;
  float **data;
  float **model;
  float *q, *t, *h;
  int nt, nh, nq,rtmethod; 
  int method;
  int plot;
  float dt,dq,dx;
  float factor; // multiply dq by factor
  int adj;

  //////////////////////////////////////////////
  fprintf(stderr,"*******SURADONHYPFK*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("vmin", &vmin))  vmin = 1500;
  if (!getparfloat("dv", &dv))  dv = 100;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("dq",&dq)) dq =2e-8;
  if (!getparfloat("step", &inv.step))  inv.step =0.9;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparint("norm", &inv.norm))  inv.norm =1; 
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparfloat("dx",&dx)) dx=10; 
  if (!getparint("adj",&adj)) adj=1; 

  // Get info from first trace 

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;

  // Allocate memory for data and model

  data=ealloc2float(nt,nh);
  model=ealloc2float(nt,nq);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);
  
  // Loop over traces 
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  
  for (it=0;it<nt;it++) t[it]=0+it*dt;
  for (iq=0;iq<nq;iq++) q[iq]=vmin+iq*dv;
  radonhypfk(t,q,h,model,data,nt,nh,nq,dt,t0,inv,plot,dx,adj);
  
  //if (smooth) smoothing(model[0],nt,nq,nl,nr,flag);
  if ((modelfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);

  for (iq=0;iq<nq;iq++){
    memcpy((void *) tr.data,(const void *) model[iq],nt*sizeof(float));
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nq;
    tr.ns=nt;
    tr.tracl=iq+1;
    tr.tracr=iq+1;
    tr.f2=q[iq];
    fputtr(stdout,&tr);    
  }

  free1float(t);
  free1float(h);
  free1float(q);
  free2float(model);
  free2float(data);

  efclose(modelfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}













