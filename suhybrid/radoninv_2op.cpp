#include "su.h"
#include "radonhybrid.h"

void radoninv_2op(float **data,  float *ph1, float *ph2, int nh, float *t, int nt, 
		  float **model, float *q,  int nq, float *q1, int nq1, float *q2, int nq2, 
		  float fmax1, float fmax2, int filter, int npoly, float *f, float *amps)
{
  int freq;
  complex **m2;
  complex **d2;
  complex czero; czero.r=czero.i=0;
  complex **L;
  complex **L1;
  complex **L2;
  float w,df;
  float dt=t[1]-t[0];
  int nfft;
  int nf;
  int maxfreq1, maxfreq2;
  float fmin=0;
  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh);  
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);

  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  if (filter==3) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);
  //if (filter==4) filtering(m2,nf,nq1,nfft,dt,f,amps,npoly);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2); 
    Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
    fprintf(stderr,":");
    
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
  }
  if (nq2>0){ // Between maxfreq1 and freqmax2   we use only one operator
    fprintf(stderr,"From now on one operator only\n");
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      Atimesx(d2[freq],L2,&m2[freq][nq1],nh,nq2,FALSE);
      fprintf(stderr,".");      
    }
   }

  fprintf(stderr,"freq=%f\n",freq*df);      
  freqweighting(d2,nf,nh,df,fmin,fmax2);
  fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf); 
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  
  return;

}






























