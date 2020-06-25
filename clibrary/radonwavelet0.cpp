#include "su.h"
#include "radonwavelet.h"

#define GO 1
#define BACK -1
/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonwavelet0(float **data, float *t, int nt, float *h, int nh, float **model, float *q, int nq, int wavn, float quantil, int scalefilt)
{
  int iw;
  float **modelw=NULL; // window for the model
  float **dataw=NULL;  // window for the data
  int *itd=NULL;   // axis with dyad sizes
  float *td=NULL;  // time axis for every dyad
  float dtd;       // time sampling for every dyad
  int ntw;         // Dimensions (nt,nh,nq) for a window; 
  int nhw;
  int nqw;
  int nw;          // Number of dyadics
  int fw;
  float dt=t[1]-t[0];
  float t0=t[0];
  float **datatemp;
  int scale;
  char buf[80];

  nhw=nh;
  nqw=nq;
  itd=ealloc1int(20);
  nw=dyad(nt,itd); // Estimate dyadic sizes
  for (iw=0;iw<nw;iw++) fprintf(stderr,"itd[%d]=%d\n",iw,itd[iw]);
  pwtset(wavn);  // set wavelet coefficients
   
  plotgather(data[0],nt,nh,"data");

  if (0){ // test for window_cpy
    dataw=ealloc2float(nt/2,nh);
    window_cpy(data,0,nh,nt/2-1,nt,dataw,TRUE);
    plotgather(dataw[0],nt/2,nh,"data_test_before");
    zero_vector(dataw[20],nt/2);
    window_cpy(data,0,nh,nt/2-1,nt,dataw,FALSE);
    plotgather(data[0],nt,nh,"data_test_after");
  }

  dwt(data,nh,nt,GO); // Wavelet transform
  plotgather(data[0],nt,nh,"data_go");

  fw=4;
  // Loop in dyadics for time, offset domain
  if (1){
    datatemp=ealloc2float(nt,nh);
    char buf2[7];
    for (iw=fw;iw<nw;iw++){
      scale=iw;
      mrgathers(data,nt,nh,scale,itd,datatemp);
      sprintf(buf2,"scale%d",scale);
      save_gather(datatemp,nh,h,nt,dt,buf2);
      //plotgather_pipe(datatemp,nh,nt,buf);
    }

    free2float(datatemp);
  }

  // Loop in dyadics for wavelet coefficients. 
  if (0) 
    for (iw=fw;iw<nw;iw++){
      ntw=itd[iw];intprint(ntw);
      dataw=ealloc2float(ntw,nh);
      window_cpy(data,0,nh,ntw-1,2*ntw,dataw,TRUE);      
      mrcoefficients(dataw,ntw,nh,iw,quantil,scalefilt);
      window_cpy(data,0,nh,ntw-1,2*ntw,dataw,FALSE);
      free2float(dataw);
  }

 // Loop in dyadics for wavelet coefficients. 
  if (0) 
    for (iw=fw;iw<nw;iw++){
      TRACE;
      ntw=itd[iw];
      fprintf(stderr,"ntw=%d\n",ntw);
      dataw=ealloc2float(ntw,nh);
      modelw=ealloc2float(ntw,nq);
      td=ealloc1float(ntw);
      dtd=dyad_timeaxis(nt,ntw,td,dt,t0);
      window_cpy(data,0,nh,ntw-1,2*ntw,dataw,TRUE);
      //plotgather_pipe(dataw,nh,ntw,"data_scale");
      radonwindow(dataw,td,ntw,h,nh,modelw,q,nq,iw,quantil,scalefilt);
      window_cpy(data,0,nh,ntw-1,2*ntw,dataw,FALSE);
      free2float(dataw);
      free2float(modelw);
      free1float(td);
  }

  dwt(data,nh,nt,BACK); // Inverse Wavelet transform

  plotgather(data[0],nt,nh,"data_back");


  free1int(itd);

  return;

}



void dwt(float **data, int nh, int nt, int sign)
{
  
  int ih;
  
  for(ih=0;ih<nh;ih++) wt1(data[ih]-1,nt,sign,pwt);
  return;
} 


int dyad(int nt, int *id)
{
  int nn, j;
  for (j=0,nn=1;nn<=nt;nn<<=1,j++) ;
  intprint(j);
  //id=ealloc1int(j);
  for (j=0,nn=1;nn<nt;nn<<=1,j++) id[j]=nn;

  return(j);
}


float dyad_timeaxis(int nt, int ntd, float *td, float dt, float t0)
{
  int it;
  float dtd=dt*ntd/nt;
  
  for (it=0;it<ntd;it++) td[it]=t0+it*dtd;

  
  return(dtd);
}

void mrcoefficients(float **data, int nt, int nh, int scale, float quantil, int scalefilt)
{
  
  char buf[20];

  sprintf(buf,"scale=%d\n",scale);
  plotgather_pipe(data,nh,nt,buf);
  if (scale<=scalefilt) filtoutliers(data,nh,nt,quantil);
  return;

}



void mrgathers(float **data, int nt, int nh, int scale, int *itd, float **datatemp)
{
  float **dataw;
  int ntw=itd[scale];

  dataw=ealloc2float(ntw,nh);
  zero_array(datatemp,nh,nt);
  window_cpy(data,0,nh,ntw-1,2*ntw,dataw,TRUE);
  window_cpy(datatemp,0,nh,ntw-1,2*ntw,dataw,FALSE);
  dwt(datatemp,nh,nt,BACK);
  free2float(dataw);

  return;

}



void radonwindow(float **data, float *t, int nt, float *h, int nh, float **model, float *q, int nq, int scale, float quantil, int scalefilt)
{
  
  char buf[20];

  //radontd_sparse(t,q,h,m,d,nt,nh,nq,dt,vel,dperv,pervmin,t0,inv_par inv,centralq,dataprec,nw,
  //		 fpeak, typewav, LI, nreg);  
  sprintf(buf,"scale=%d\n",scale);
  plotgather_pipe(data,nh,nt,buf);
  if (scale<=scalefilt) filtoutliers(data,nh,nt,quantil);
  

  return;

}






