#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "su.h"

void radonopi(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	time=sqrt(t[itau]*t[itau]+pxhxh);
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
        dtemp[itau]=d[ihxnt+itau];
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
 	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            
      }
    }
  }
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}

void radonopi(float *m, float *t, float *h, float *q, float *d, float *ww) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");  
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	//if (fabs(ww[k])>1e-2){
	time=sqrt(t[itau]*t[itau]+pxhxh);
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
	dtemp[itau]=d[ihxnt+itau];
	//}
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            	  
      }
    }
  } 
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}

void radonopi_old(float *m, float *t, float *h, float *q, float *d, float *ww) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	if (fabs(ww[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+pxhxh);
	  /*time=(t[itau]*t[itau])+pxhxh;*/
	  it=time/dt;
          itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
          //a1=1; a2=0;
	  if ((it!=nt)&&(j<(ny-1))&&(k<nx)) 
	    m[k]=m[k]+a1*d[j]+a2*d[j+1];/**sincin(it-floor(it));*/
	  
	}
      }
    } 
  }
  return;
}


void radonopi_old(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	time=sqrt(t[itau]*t[itau]+pxhxh);
	/*time=(t[itau]*t[itau])+pxhxh;*/
	it=time/dt;
	itint=(int) floor(it);
	j=ihxnt+itint;
	a1=1-(it-itint);
	a2=it-itint;
      
	if ((it!=nt)&&(j<ny-1)&&(k<nx)) 
	  m[k]=m[k]+a1*d[j]+a2*d[j+1];/**sincin(it-floor(it));*/	  
      }      
    }
  }
  return;
}


void radonopi( float *m, float *t, float *h, float *q, float *d, float *ww, float theta) 
{
  int i,k,ih,iq,itau,itint;;
  unsigned int j;
  float time,it,hxh,qxhxh,qxhxhxs,sintheta,a1,a2;
  int iqxnt,ihxnt;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  theta=theta/180.*acos(-1.);
  sintheta=sin(theta)*sin(theta);    
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      qxhxh=hxh*q[iq];
      qxhxhxs=qxhxh*sintheta;
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	if (fabs(ww[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+qxhxh-qxhxhxs);
	  /*time=(t[itau]*t[itau])+pxhxh;*/
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
	  a1=1-(it-itint);
	  a2=it-itint;
	  if ((it!=nt)&&(j<ny-1)&&(k<nx))
	    m[k]=m[k]+a1*d[j]+a2*d[j+1];/**sincin(it-floor(it));*/
	  //m[k]=m[k]+d[j];/**sincin(it-floor(it));*/	  
	}      
      }
    }
  }
  return;
}

void radonopi( float *m, float *t, float *h, float *q, float *d, float *ww, char flag) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh,pxhxh;
  int iqxnt,ihxnt;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
 	pxhxh=hxh*q[iq];
	iqxnt=iq*nt;
	for (itau=0;itau<nt;itau++){
	  k=iqxnt+itau;
	  if (fabs(ww[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+pxhxh);
	  /*time=(t[itau]*t[itau])+pxhxh;*/
	  it=time/dt;
	  j=ihxnt+(int) floor(it);
	  if ((it!=nt)&&(j<ny)&&(k<nx)) 
	    m[k]=m[k]+d[j];/**sincin(it-floor(it));*/
	  
	}
      }
      
    }
  }
  return;
}


// Radon operator with input driven estimation
// Idea based in Wisecup(1998) Geophysics.
// The idea here is to store all the d at their exact location
// (no interpolation) and the time at which they will contribute 
// for a given slope. Finally this irregularly sampled output is 
// manipulated in some way (moving average,  perhaps), to 
// produce a free alias model.
// The problem is how to combine the list of times and models
// (is like stack of irregular traces)
  
void radonopi_id(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau,ktime,itime,*iktime;
  unsigned int j;
  float *tau,*mint,*dtemp,time,hxh,pxhxh,tt,*tau2,*dtemp2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  float kk,kkk;

  if ((tau=alloc1float(nh*nt))==NULL)
    err("cannot allocate memory for tau \n");
  if ((tau2=alloc1float(nh*nt))==NULL)
    err("cannot allocate memory for tau2 \n");
  if ((mint=alloc1float(nt))==NULL)
    err("cannot allocate memory for mint \n");
  if ((dtemp=alloc1float(nt*nh))==NULL)
    err("cannot allocate memory for dtemp \n");
  if ((dtemp2=alloc1float(nt*nh))==NULL)
    err("cannot allocate memory for dtemp \n");
  if ((iktime=alloc1int(nt*nh))==NULL)
    err("cannot allocate memory for iktime \n");
  //rwa_smoothing_filter (2, int nl, int nr, float *filter);
  for (i=0;i<(nq*nt);i++) m[i]=0;
  for (iq=0;iq<nq;iq++){
    pxhxh=hxh*q[iq];
    iqxnt=iq*nt;
    for (ktime=0,ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;    
      for (itime=0;itime<nt;itime++)
	if ((tt=t[itime]*t[itime]-pxhxh) >= 0.0){
          ktime+=1;
          //iktime[ktime]=ktime;
	  tau[ktime]=sqrt(tt);
	  dtemp[ktime]=d[ihxnt+ktime];          
	}
    }

    /* Comments
    Here I have a list of tau and d that contributes to m[tau,iq]
    This is a irregular sampling for m[tau,iq]
    How to stack them???????????
    here I need a function that performs interpolation from
    irregular to regular.
    Example the sinc8 interpolation 
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    goes from regular grid in
    dtemp, to a irregular grid dint[itime]. 
    The opposite is to go from dtemp[tau]
    */


    //qkisort (ktime,tau,iktime);
    //for (i=0;i<ktime;i++){
    //  tau2[i]=tau[iktime[i]];
    //  dtemp2[i]=dtemp[iktime[i]];
      //fprintf(stderr,"tau2[%d]=%e,dtemp2[%d]=%e\n",i,tau2[i],i,dtemp2[i]);
    //}

    for (itau=0;itau<nt;itau++){
      k=iqxnt+itau;kk=itau*dt;kkk=kk+dt;
      for (i=0;i<ktime;i++) 
	if ((tau[i] > kk) && (tau[i] < kkk))  
	  m[k]+=dtemp[i];
    }
    //intlin(ktime,tau2,dtemp2,dtemp2[0],dtemp2[ktime-1],nt,t,mint);     
    //for (itau=0;itau<nt;itau++){
    //  k=iqxnt+itau;        
    //  m[k]=mint[itau];            
    //}    
  }
  free1float(dtemp);
  free1float(tau);
  free1float(dtemp2);
  free1float(tau2);
  free1float(mint);
  free1int(iktime);
  return;
}

// Routine used to plot the LH operator in time domain.
// Present form is not useful, just for test on operator LH

void radonopi_id2(float *m, float *t, float *h, float *q, float *d) 
{
  int i,j,k,ih,iq,itau;
  float time,hxh,qxhxh,a1,a2;
  float *LH, *mtemp, *dtemp;
  complex *LHC, *dc, *mc;
  int nfft;
  int iqxnt,ihxnt,itint;
  int itime;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  nfft=npfar(2*ny); 
  if ((LH=ealloc1float(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for LH\n");
  //for (i=0;i<(nx);i++) for (j=0;j<(ny);j++) LH[i][j]=0;

    
  if ((LHC=ealloc1complex(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for LHC\n");
  if ((dc=ealloc1complex(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for dc\n");   
   if ((mc=ealloc1complex(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for mc\n"); 
   if ((mtemp=ealloc1float(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for mc\n"); 
    if ((dtemp=ealloc1float(nfft))==NULL)  //==> LH(nq*nt x nh*nt)
    err("cannot allocate memory for mc\n");   
  for (i=0;i<(nq*nt);i++) m[i]=0;
  //for (iq=0;iq<nq;iq++){

  for (iq=0;iq<1;iq++){
    iqxnt=iq*nt;
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      qxhxh=hxh*(q[2]-q[1]);
      for (itau=0;itau<nt;itau++){
	time=t[itau]+qxhxh;
	itime=int (time/dt);
	//if (itime<(nt-1)) LH[itau*nq+iq][itime*nh+ih]=1;
	//if (itime<(nt-1)) LH[iqxnt+itau][ihxnt+itime]=1;
        LH[ihxnt+itime]=1;
      }
    } 
  }
 
  /*
  FILE *myfilep;
  if((myfilep=fopen("LH.bin","w"))==NULL)
    err("cannot open file=%s\n","LH.bin");
  //for (i=0;i<(nx);i++)
  //  fwrite(LH[i],sizeof(float),ny,myfilep);
  //fclose(myfilep);

  /* 
  for (i=0;i<(nx);i++){
    for (j=0;j<(ny);j++){
      fprintf(myfilep,"%2.0f",LH[i][j]);
    }
    fprintf(myfilep,"\n");
  }
		 
  fclose(myfilep);
		 */

  
  for (i=nfft;i>0;i--) LH[i]=LH[nfft-i+1];
  
  for (i=ny;i<nfft;i++) dtemp[i]=0;
  for (i=0;i<ny;i++) dtemp[i]=d[i];
  
  pfarc(1,nfft,LH,LHC);
  pfarc(1,nfft,dtemp,dc);  
  for (i=0;i<nfft;i++) mc[i]=LHC[i]*dc[i];
  pfacr(-1,nfft,mc,mtemp);
  for (i=0;i<nx;i++) m[i]=mtemp[i]/nfft;
  
  /*  
  for (iq=0;iq<nq;iq++){
    iqxnt=iq*nt;
    for (itau=0;itau<nt;itau++){
      k=iqxnt+itau;
      for (i=0;i<nh*nt;i++)        
	m[k]=m[k]+LH[k][i]*d[i];            
    }     
  }
  */

 	 
  free1float(mtemp);
  free1float(LH);
  free1complex(dc);
  free1complex(mc);
  free1complex(LHC);
  free1float(dtemp);
  fprintf(stderr,"Here we go +++++++++++++=%d,%d,%d\n",i,nx,ny);
  return;
}



void radonopi_lin(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxh=h[ih]*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	time=t[itau]+pxh;
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
        dtemp[itau]=d[ihxnt+itau];
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
 	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            
      }
    }
  }
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}




void radonopi_par(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	time=t[itau]+pxhxh;
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
        dtemp[itau]=d[ihxnt+itau];
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
 	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            
      }
    }
  }
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}









