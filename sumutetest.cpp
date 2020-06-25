/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADONWAVELET  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "./clibrary/radonclass.hpp"
#include "mutetest.h"
//#include <signal.h>
#include <math.h>
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUMUTETEST -  Program to test automatic mute                        ",
  " 	   								",
  " sumutetest   < stdin > stdout [optional parameters]      		",
  "                                                                     ",
  " verbose=0       =1  Extra Information                               ",
  NULL};
/* Credits:
 *	Daniel Trad.
 Example:
 sumutetest < hyperbolas3.su modelfile=model.su > pp qmin=-1e-5 
 qmax=1e-5 mute=1  nq=200 ihmin_m=10 ihmax_m=190 thres_m=5 slope_m=0 
 tmin_m=0.5 tmax_m=1
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
  float h0=0;
  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nq;  
  float dt; // Time sampling
  float dh;
  float dq;
  float t0;      // First useful time

  float qmin;
  float qmax;
  int verbose;		/* flag for echoing info		*/

  float *trace; // auxiliar trace

 



  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
  // Get info from first trace 
  if (!gettr(&tr)) err("cannot read first trace"); 
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (!(nh=tr.ntr)) err("***ntr must be set\n");
  fprintf(stderr,"nt=%d,nh=%d,dt=%f\n",nt,nh,dt);  


  /* for mute in the migrated space we need to set the geometry of the mask */
  int mute; 
  mutemask_par par;


  if (!getparint("mute",&mute)) mute = 0;
  if (!getparfloat("tmin_m",&par.tmin)) par.tmin = 0.; 
  if (!getparfloat("tmax_m",&par.tmax)) par.tmax = nt*dt;
  if (!getparint("ihmin_m",&par.ihmin)) par.ihmin = (int) (nh/2+2);
  if (!getparint("ihmax_m",&par.ihmax)) par.ihmin = (int) (nh);
  if (!getparfloat("slope_m",&par.slope)) par.slope = 3;     
  if (!getparfloat("thres_m",&par.threshold)) par.threshold = 0.2;    
 
  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparfloat("qmax", &qmax))  qmax = 5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";

  dq = (qmax-qmin)/(nq-1);
  dh=1;

  modelfilep=efopen(modelfile,"w");
  trace=ealloc1float(nt);
  TGather data(nt,nh,dt,dh,t0,h0);
  TGather model(nt,nq,dt,dq,t0,qmin);
   
  ih = 0;
  do {
    memcpy((void *) trace,(const void *) tr.data,nt*sizeof(float));
    data.settrace(trace,ih);
    data.setx((float) tr.offset,ih);
    ih++;
  } while (gettr(&tr));
  erewind(stdin);

  if (nh != ih) err("ntr was wrong\n");
  fprintf(stderr,"nh=%d\n",nh);
  
  radonhyper(data,model);
  data.plot("gather.su","offset");
  //model.plot("model","f2");
  if (mute) mutemask(model,par);
  model.plot("model","f2");

  // Output trace
  for (ih=0;ih<nh;ih++){ 
    gettr(&tr);
    data.gettrace(trace,ih);
    memcpy((void *) tr.data,(const void *) trace,nt*sizeof(float));
    puttr(&tr);
  }
  
  for (iq=0;iq<nq;iq++){
    model.gettrace(trace,iq);
    memcpy((void *) tr.data,(const void *) trace,nt*sizeof(float));
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nq;
    tr.ns=nt;
    tr.tracl=iq+1;
    tr.tracr=iq+1;
    tr.f2=model.getx(iq);
    fputtr(modelfilep,&tr);    
  }
 
  free1float(trace);
  
  efclose(modelfilep);
  finish=time(0);
  elapsed_time=difftime(finish,start);

  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}



void radonhyper(TGather &data, TGather &model)
{
  register int it;
  int ih,iq;
  int nt, nh, nq;
  int adj=1;
  float ftime;
  float time,hxh,pxhxh;
  int iqxnt,ihxnt;
  int itime;

  model.getsize(nt,nq);
  data.getsize(nt,nh);

  float *h=ealloc1float(nh);
  float *q=ealloc1float(nq);
  float *t=ealloc1float(nt);

  float **m=ealloc2float(nt,nq);
  float **d=ealloc2float(nt,nh);

  model.gettaxis(t);
  model.getxaxis(q);
  data.getxaxis(h);

  float dt=t[1]-t[0];

  
  if (adj) data.getdata(d);
  else model.getdata(m); 


  fprintf(stderr,"Crude approach\n");
  if (adj)  for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) m[iq][it]=0;
  else for(ih=0;ih<nh;ih++) for (it=0;it<nt;it++) d[ih][it]=0;   

  for (it=0;it<nt;it++){
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){    
	pxhxh=hxh*q[iq];
	iqxnt=iq*nt;
        time=sqrt(t[it]*t[it]+pxhxh);
	ftime=time/dt;
	itime=(int) floor(ftime+0.5);
	if (itime<nt){
	  if (adj) m[iq][it]+=d[ih][itime];
          else d[ih][itime]+=m[iq][it];
	}
      }            
    }
  } 


  if (adj) model.setdata(m);
  else data.setdata(d); 

  free2float(m);
  free2float(d);

  return;
}



//#define tmin 0.7
//#define threshold 0.5

void mutemask(TGather &model, mutemask_par par)
     /* 
	Function to compute a muting mask for elimination of 
	diffractions (kencmpwin.su example).
	structure mutemask contains:
	float tmin, float tmax, int ihmin, int ihmax, 
	float threshold, float slope 
     */

{
  int ih,it;
  int itmin;
  int itmax;
  int shift;
  int nt;
  int nh;
  //int ihmin;
  //int ihmax;

  model.getsize(nt,nh);

  float *h=ealloc1float(nh);
  float *t=ealloc1float(nt);
  float **d=ealloc2float(nt,nh);
  float **M=ealloc2float(nt,nh);

  model.gettaxis(t);
  model.getxaxis(h);
  model.getdata(d); 
  float dt=t[1]-t[0];


  itmin=(int) (par.tmin/dt);
  itmax=(int) (par.tmax/dt);
  //ihmin=(int) (nh/2+1);
  //ihmax=(int) (nh/2+50);
  fprintf(stderr,"ihmin=%d,ihmax=%d,itmin=%d,itmax=%d,slope=%f,threshold=%f,nt=%d,nh=%d\n",par.ihmin,par.ihmax,itmin,itmax,par.slope,par.threshold,nt,nh);  


  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++) M[ih][it]=0;


  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++)
      if ((it>itmin)&&(it<itmax)){
	shift=(int) ((it-itmin)*dt*par.slope);
	if ((ih>par.ihmin+shift)&&(ih<par.ihmax+shift))
	  if (fabs(d[ih][it])>par.threshold)
	    M[ih][it]=1;
      }
  
      
  // test
  //model.plot("model1","f2");
  AtimesB(d,M,nh,nt);
  model.setdata(d); 
  //model.plot("model2","f2");
  free2float(M);
  free2float(d);
  free1float(t);
  free1float(h);
  return;
}

void AtimesB(float **A, float **B, int nh, int nt)
{
  int ih,it;
  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++)
      A[ih][it]*=(B[ih][it]);

  return;
}
 














