#include "su.h"

int read_ascii_file(const char *name,float *x);

void plotvector(float *d, int nt, const char *s);

float dot(int n, float *a, float *b);

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

int get_wavelet(float *wavelet,const char *name,int nw, int type, float dt, float fpeak)
{
  
  char buf[80];
  int iw;
  int plot=1;

  if (type==1){
    nw=read_ascii_file(name,wavelet);
    if (plot){ 
      sprintf(buf,"xgraph < wavelet n=%d pairs=2 d1=1 style=normal title=\"BP_wavelet\" &",nw);
      plotvector(wavelet,nw,buf);
    }
  }   


  if (type==2){
    fprintf(stderr,"nw=%d,fpeak=%f\n",nw,fpeak);
    ricker1_wavelet (nw,dt,fpeak,wavelet);
    iw=5;while ((fabs(wavelet[iw]) > 1e-5)&&(iw<nw)) iw++;
    nw=iw;
    if(plot){
      sprintf(buf,"xgraph < wavelet n=%d pairs=2 d1=1 style=normal title=\"wavelet\" & ",nw);
      plotvector(wavelet,nw,buf);
    }
  }


  if (type==3){
    for (iw=0;iw<nw;iw++) wavelet[iw]=1;
    /*
      wavelet[0]=-0.25;
      wavelet[1]=1;
      wavelet[2]=-0.25;
      nw=3;
    */
  }

  if (type==4){
    for (iw=0;iw<nw;iw++) wavelet[iw]=0;
    wavelet[0]=1;
    //wavelet[1]=1;
    //wavelet[2]=0;
    nw=1;
  }

  return(nw);
}


