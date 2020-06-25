#include "dan.h"
#include "su.h"
#define SLOPE 2

int taper(float **data, int nt, int nh, int ntaper, int flag);

/* Function to compute and subtract the multiples generated after muting the primaries 
   in the RT space.
   Because the forward operator can take different shapes I need to write different 
   routines, but the only difference is the forward operator used. 

   Some experimental variables.
   plot: =1 plot gathers before and after a
   itm:  First multiple is below this time

   Daniel Trad UBC - June 2001

*/

void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int, float *, int), 
		    float **model, float **data, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float *wavelet, int nw, float parmute, 
		    float *q, float *h, int itm, int plot)
{
  int side=1;  // side 2 --> right mute ; side 1 left mute
  float **datarec=0; // multiples
  float dpd;
  float dpdp;
  float scale;
  float **modeltemp=0;
  int iq,ih,it;
  int nmute=0;
  int slope=SLOPE; // Defines a slope in the muting

  modeltemp=ealloc2float(nt,nq);
  datarec=ealloc2float(nt,nh);
  

  side=1;
  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;
  fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);

  for (iq=0;iq<nq;iq++)
    memcpy(modeltemp[iq],model[iq],nt*FSIZE);
  
  taper(modeltemp,nt,nq,nmute,side);
  // Temporal change to leave primaries untouched at the shallow times
  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(itm-(iq-nmute)*slope);it++) 
      modeltemp[iq][it]=0;   

  if (plot){
    save_gather(modeltemp,nq,q,nt,0.004,"mutedmodel.su");
    system("suxwigb < mutedmodel.su perc=95 key=f2 title=mutedmodel & ");
  }

  oper(modeltemp[0],datarec[0],index,0,nt,nh,nq,nsparse,wavelet,nw);
  
  if (plot){
    save_gather(datarec,nh,h,nt,0.004,"multiples.su");
    system("suxwigb < multiples.su perc=95 key=offset  title=multiples & ");
  }

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datarec[0]);
  dpdp=dot(nh*nt,datarec[0],datarec[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale 2===>%f\n",scale);

  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datarec[ih][it];
  
  
  free2float(modeltemp);
  free2float(datarec);
  return;
}




void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int), 
		    float **model, float **data, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float parmute, float *q, float *h, int itm, int plot)
{
  int side=1;  // side 2 --> right mute ; side 1 left mute
  float **datarec=0; // multiples
  float dpd;
  float dpdp;
  float scale;
  float **modeltemp=0;
  int iq,ih,it;
  int nmute=0;
  int slope=SLOPE; // Defines a slope in the muting
  
  modeltemp=ealloc2float(nt,nq);
  datarec=ealloc2float(nt,nh);
  

  side=1;
  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;
  fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);

  for (iq=0;iq<nq;iq++)
    memcpy(modeltemp[iq],model[iq],nt*FSIZE);
  
  taper(modeltemp,nt,nq,nmute,side);
  // Temporal change to leave primaries untouched at the shallow times
  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(itm-(iq-nmute)*slope);it++)  // 0 should be it0 
      modeltemp[iq][it]=0; 

  if ((plot==2)||(plot==3)){
    //system("rm mutedmodel.su");
    save_gather(modeltemp,nq,q,nt,0.004,"mutedmodel.su");
    system("suxwigb < mutedmodel.su perc=97 key=f2 title=mutedmodel ");
  }

  oper(modeltemp[0],datarec[0],index,0,nt,nh,nq,nsparse);

  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(datarec,nh,h,nt,0.004,"multiples.su");
    system("suxwigb < multiples.su perc=97 key=offset  title=multiples ");
  }

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datarec[0]);
  dpdp=dot(nh*nt,datarec[0],datarec[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale 2===>%f\n",scale);

  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datarec[ih][it];
  
  
  free2float(modeltemp);
  free2float(datarec);
  return;
}
