#include "su.h"
#include "segy.h"
#include "radonline_stack.h"
#define plotdata 0 /* set to 1 to plot data after NMO */


/*
radonline_stack0 

Daniel Trad - January- 2002
*/


void radonline_stack0(float **data, complex **oper,float *trace, float *h, int nh,float *t, int nt, float dt, float *vel, float smute, float nmofactor,  float fmax, int nfft, float df)
{
  int it, ih;
  float *dtemp;
  complex **d=0;
  int maxfreq=(int) (fmax/df);
  complex czero; czero.r=czero.i=0;
  int nf;
  complex *sf;
  float *st;
  int freq;
  float w;
  float wa;
  int verbose=0;
  nf=(nfft/2+1);
  if (verbose) fprintf(stderr,"nf=%d,nfft=%d,nh=%d,nt=%d,maxfreq=%d\n",nf,nfft,nh,nt,maxfreq);

 
  d=ealloc2complex(nh,nt);  // data in the frequency domain
  sf=ealloc1complex(nt);  // stack in the freq domain
  st=ealloc1float(nt);    // stack in the time domain

  dtemp=ealloc1float(nt);
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  free1float(dtemp);

  if ((plotdata)&&(0)) plotgather_pipe(data,nh,nt,"After NMO");

  /* transsform to frequency and apply Stack+IRT+Mute+RT operator */
  fftgo_xt2fx(-1,data,d,nh,nt,dt,nfft,nf);

  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax); 
    sf[freq]=czero;
    for (ih=0;ih<nh;ih++) sf[freq]+=d[freq][ih]*oper[freq][ih];
    if ((wa<1)&&(wa>0)) sf[freq]*=wa;
  }
  sf[0]=czero;
  for (freq=maxfreq;freq<nt;freq++) sf[freq]=czero;

  fftback_fx2xt(1,trace,sf,nt,dt,nfft,nf);

  free1float(st);
  free1complex(sf);
  free2complex(d);

  return;

}


void line_stack0(float **data, complex **oper,float *trace, float *h, int nh,float *t, int nt, float dt, float *vel, float smute, float nmofactor,  float fmax, int nfft, float df)
{
  int it, ih;
  float *dtemp;
  complex **d=0;
  int maxfreq=(int) (fmax/df);
  complex czero; czero.r=czero.i=0;
  int nf;
  complex *sf;
  float *st;
  int verbose=0;
  float *tr;
  nf=(nfft/2+1);
  if (verbose) fprintf(stderr,"nf=%d,nfft=%d,nh=%d,nt=%d,maxfreq=%d\n",nf,nfft,nh,nt,maxfreq);

 
  d=ealloc2complex(nh,nt);  // data in the frequency domain
  sf=ealloc1complex(nt);  // stack in the freq domain
  st=ealloc1float(nt);    // stack in the time domain
  tr=ealloc1float(nt);
  if ((plotdata)&&(0)) plotgather_pipe(data,nh,nt,"Before NMO");
  dtemp=ealloc1float(nt);
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  free1float(dtemp);

  if ((plotdata)&&(0)) plotgather_pipe(data,nh,nt,"After NMO");

  for (it=0;it<nt;it++){
    st[it]=0;
    for (ih=0;ih<nh;ih++) st[it]+=(data[ih][it]/nh);
  }

  if (0){
    fftgo_xt2fx(-1,st,sf,nt,dt,nfft,nf);
    fftback_fx2xt(1,st,sf,nt,dt,nfft,nf);
  }

  for (it=0;it<nt;it++) trace[it]=st[it];
  if (0) plotgather_pipe(trace,1,nt,"trace");

  free1float(tr);
  free1float(st);
  free1complex(sf);
  free2complex(d);

  return;

}


void test2(segy tr, int nh, int nt)
{
  int it;
  TRACE;
  for (it=0;it<nt;it++){
    tr.data[it]=0;
    fprintf(stderr,"after tr.data[%d]=%f\n",it,tr.data[it]);
  }
  TRACE;
  return;
}


void radonline_LS(float **data, complex ***RT, complex ***L, float **model, float *trace,  float *h, int nh,float *t, int nt, float *q, int nq, float dt, float *vel, float smute, float nmofactor,  float fmax, int nfft, float df, float *mute, float *stackvector)
{
  int it, ih, iq;
  float *dtemp;
  complex **d=0;
  complex **m=0;
  complex *tracefreq;
  
  int maxfreq=(int) (fmax/df);
  complex czero; czero.r=czero.i=0;
  int nf;
  int freq;
  float w;
  float wa;
  int verbose=1;
  nf=(nfft/2+1);

  int save=1; /* =1 for saving data at different stages */
  int recoverdata=1; /* =1 to apply inverse nmo after RT */

  /* max index in frequency as a function of maximum frequency */
  if (fmax) maxfreq=(int) (fmax/df);
  else maxfreq=(int) (0.9*nf);  

  if (verbose) fprintf(stderr,"nf=%d,nfft=%d,nh=%d,nt=%d,maxfreq=%d\n",nf,nfft,nh,nt,maxfreq);

  if ((plotdata)&&(0)) plotgather_pipe(data,nh,nt,"Before NMO");
  if (save){
    save_gather(data,nh,h,nt,dt,"before");
    system("suxwigb key=offset < before title=\"NO NMO\" perc=99 & \n");
  }
  /* Apply NMO correction to the time domain */
  dtemp=ealloc1float(nt);
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }


 
  /* Arrays to work in the frequency domain */
  d=ealloc2complex(nh,nt);  // data in the frequency domain
  m=ealloc2complex(nq,nt);  // Radon domain in the frequency domain
  tracefreq=ealloc1complex(nt);
  if ((plotdata)&&(0)) plotgather_pipe(data,nh,nt,"After NMO");
  if (save){
    save_gather(data,nh,h,nt,dt,"after");
    system("suxwigb key=offset < after title=\"NMO\" perc=99 & \n");
  }

  /* transsform to frequency and apply Stack+IRT+Mute+RT operator */
  fftgo_xt2fx(-1,data,d,nh,nt,dt,nfft,nf);

  /* Loop in frequency */
  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax); 
    Atimesx(m[freq],RT[freq],d[freq],nq,nh);
    if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) m[freq][iq]*=wa;
    xtimesy(m[freq],m[freq],mute,nq);
  }
  /* Set to zero uncalculated frequencies */
  for (iq=0;iq<nq;iq++) m[0][iq]=czero;
  for (freq=maxfreq;freq<nt;freq++) for (iq=0;iq<nq;iq++) m[freq][iq]=czero;
  
  fftback_fx2xt(1,model,m,nq,nt,dt,nfft,nf);

  if (0) plotgather_pipe(model,nq,nt,"Radon Domain");
  if (save){
    save_gather(model,nq,q,nt,dt,"RT");
    system("suxwigb key=f2 < RT title=\"Radon\" perc=99 & \n");
    system("cat RT >> radondata  & \n");
  }
  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax); 
    Atimesx(d[freq],L[freq],m[freq],nh,nq);
    if ((wa<1)&&(wa>0)) for (ih=0;ih<nh;ih++) d[freq][ih]*=wa;
    tracefreq[freq]=czero;
    for (ih=0;ih<nh;ih++) tracefreq[freq]+=(stackvector[ih]*d[freq][ih]);
  }
  /*  Set to zero uncalculated frequencies */
  for (ih=0;ih<nh;ih++) d[0][ih]=czero;
  for (freq=maxfreq;freq<nt;freq++) for (ih=0;ih<nh;ih++) d[freq][ih]=czero;
  tracefreq[0]=czero;
  for (freq=maxfreq;freq<nt;freq++) tracefreq[freq]=czero;
  /*****************************************/
  fftback_fx2xt(1,data,d,nh,nt,dt,nfft,nf);
  fftback_fx2xt(1,trace,tracefreq,nt,dt,nfft,nf);

 
 if (0) plotgather_pipe(data,nh,nt,"Recovered time Domain");
 if (recoverdata){
   for (ih=0;ih<nh;ih++){
     for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
     nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
   }    
 }
 if ((save)&&(recoverdata)){
    save_gather(data,nh,h,nt,dt,"demultiple");
    system("suxwigb key=offset < after title=\"Demultiple\" perc=99 & \n");
 }
 free1float(dtemp);
 free1complex(tracefreq);
 free2complex(m);
 free2complex(d);

  return;

}



















