/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADONDTD):  $Date: June 1999 - Last version October 2000  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include "interptx_win2.h"
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADINT Forward  High Resolution Hyperbolic Radon transform        ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradon1< stdin > stdout [optional parameters]          		",
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
  " dq=2e-8             defines the rate of constant slowness spacing   ",
  " nq=20		Number of Radon Traces                         	",
  " factor=0.8		multiply dq critical by factor                	",
  " precond=0           Uses precondtioning                             ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " parnmute=1      if q > parmute the radon trace is tapered           ",
  " resample=0      =0 keeps original offset                            ",
  "                 =1 defines a new offset axis with similar aperture  ",
  "                 that original one but regular and resamples the data",
  "                 =2 search for gaps and fills them keeping original  ",
  "                  traces.                                            ",
  "                 =3 resamples to the offset given by offsetfile      ",
  "                 =0 keeps original offset                            ", 
  " keepdata=0      when =1 if there are traces with equal offset than  ",
  "                 the interpolated ones, keep the original traces.    ",
  " flagvel=1       to use velocities for vgrid                         ",
  "        =0       to use square slowness like conventional Hyp RT     ",
  " centralq=0      the q trace where the velocity law lays             ",
  " dataprec =0      event with energy higher than quantil 99 are        ",
  "                 downweighted by Q50/Q99                             ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : Iterpolated traces  (offset time domain)                   ",
  " Input : sudata file  (offset time domain)              		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradon1 pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudataint  ", 
  "                                                                     ",
  " If nh output is > nh input only the header words of the first nh    ",
  "  traces are preserved                                               ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
//inv_par inv;
inv_par inv;
segy tr,tri; 

int verbose;

int main(int argc, char **argv)
{
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep; 
  time_t start,finish;
  double elapsed_time;
  int it,iq,ih;
  float **data;
  float **model;
  float **datap;
  float *q, *t, *h;
  float t0;
  int smooth;
  int nt, nh, nq,rtmethod; 
  int method;
  float dt,dq;
  float factor; // multiply dq by factor
  /* Velocity */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  float *velint;/* array[nt] of vel for a particular trace */
  float pervmin;
  float dperv;
  int centralq; /* q trace where we put the velocity law */
  int dataprec;  /*  dataprec =0      event with energy higher than quantil 99 are  
                    downweighted by Q50/Q99  */                         
  float cdpgather;
  int taperflag;
  // smoothing
  int nl=2;  //  npoints left hand side
  int nr=2;  //  npoints left hand side
  int flag=2;  // 1 rectangular, 2 triangular
  ////////////////////////
  int nw;
  float fpeak;
  int typewav;  // type of wavelet
  int LI;  // Linear interpolation 
  int nreg; // defines how many traces besides the velocity trend have reg sampling
  int plot;
  int outputmodel;
  float tm;
  int itm;
  float smute;
  int mute;
  float parmute;

  // data weights
  int lstaper=0;
  int lbtaper=0;  
  ////////////////

  //inv.restart=1;
  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  //////////////////////////////////////////////
  fprintf(stderr,"*******SURADONTD0*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("dq",&dq)) dq =2e-8;
  if (!getparfloat("step", &inv.step))  inv.step =0.9;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparint("norm", &inv.norm))  inv.norm =1; 
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparint("taperflag", &taperflag))  taperflag =0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("smooth",&smooth)) smooth=0;
  if (!getparint("centralq",&centralq)) centralq=0;
  if (!getparint("dataprec",&dataprec)) dataprec =0;
  if (!getparint("nw",&nw)) nw =0;
  if (!getparfloat("fpeak",&fpeak)) fpeak =25;
  if (!getparint("typewav",&typewav)) typewav = 1;
  if (!getparint("LI",&LI)) LI = 0;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("nreg",&nreg)) nreg = 5;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparfloat("tm",&tm)) tm=0;
  if (!getparfloat("smute",&smute)) smute=0;
  if (!getparint("outputmodel",&outputmodel)) outputmodel = 1;
  if (!getparfloat("parmute",&parmute)) parmute=1e-8;
  if (!getparint("mute",&mute)) mute=0;
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  if (!getparint("lstaper",&lstaper)) lstaper=0;
  if (!getparint("lbtaper",&lbtaper)) lbtaper=0;

  // Get info from first trace 

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");
  if (!tr.cdp) cdpgather=1;
  else cdpgather=tr.cdp;

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  nq=2*nh;
  // Time for the first multiple in index units
  itm=(int) (tm/dt+0.5);

  /* Introduce velocity trend to apply Hyp Vel Filtering */
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  fprintf(stderr,"ncdp=%d,dt=%f,nt=%d,cdpgather=%f\n",ncdp,dt,nt,cdpgather);

  getvelocities(dt,nt,ncdp,cdp,ovv);

  fprintf(stderr,"ncdp=%d\n",ncdp);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);
  //for (it=0;it<nt;it++) fprintf(stderr,"velint[%d]=%f\n",it,velint[it]);
  ////////////////////////////////////////////////////////   

  // Allocate memory for data and model

  data=ealloc2float(nt,nh);
  model=ealloc2float(nt,nq);
  datap=ealloc2float(nt,nq);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);
  
  // L
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  

  /* If offetfile name is given read it */
  memset( (void *) q, (int) '\0', nq * FSIZE);
  if (offsetfile) nq=read_ascii_file(offsetfile,q);  
  erealloc1float(q,nq);

  for (it=0;it<nt;it++) t[it]=0+it*dt;
 
  /* apply side and bottom tapers */
  for (ih=0; ih<nh; ++ih)  taper(lstaper,lbtaper,nh,ih,nt,data[ih],0);

  if (smute) smute_gather(data,nt,nh,t,h,velint,smute);
  interptx_sparse(t,q,h,model,data,datap,nt,nh,nq,dt,velint,dperv,pervmin,t0,inv,
		 centralq,dataprec,nw,fpeak,typewav,LI,nreg,parmute,mute,itm,plot);
  
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
    if (outputmodel) fputtr(modelfilep,&tr);    
  }

  //if (smooth) smoothing(data[0],nt,nh,nl,nr,flag);
  ih=0;

  if (0){
    for (ih=0;ih<nh;ih++){ 
      gettr(&tr);
      memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
      puttr(&tr);
    }
  }
  else{
    for (iq=0;iq<nq;iq++){ 
      gettr(&tr);
      memcpy((void *) tr.data,(const void *) datap[iq],nt*sizeof(float));
      tr.ntr=nq;
      tr.ns=nt;
      tr.tracl=iq+1;
      tr.tracr=iq+1;
      tr.offset=q[iq];
      puttr(&tr);
    }
  }

  free1float(cdp);
  free2float(ovv);
  free1float(velint);  
  free1float(t);
  free1float(h);
  free1float(q);
  free2float(datap);
  free2float(model);
  free2float(data);

  efclose(modelfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}













