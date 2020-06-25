#include "eomig0freq.h"

void stack_rteom_operator(complex **F, int maxfreq, float df, float *h, int nh, 
		       float *q, int nq, float *t, int nt, float eps, float parmute1, 
		       float parmute2)
{
  /* Variables with capital letters corresponds to matrices */

  int ih,freq;
  float w=0;
  float *g=0;  
  int rtmethod=2;
  float depth=4000; // Not used
  complex **L;    // RT oprator
  complex **RT=0; // LS RT operator
  float *mute=0; // mute operator 
  complex **IRT=0; // LS RT operator after mute 
  float *stackop=0;
  complex *finalop;

  g=ealloc1float(nh);
  L=ealloc2complex(nq,nh); 
  RT=ealloc2complex(nh,nq);
  IRT=ealloc2complex(nq,nh);
  mute=ealloc1float(nq);
  stackop=ealloc1float(nh);
  finalop=ealloc1complex(nh);

  radon_moveout(h,g,nh,rtmethod,depth); 
  mute_vector(mute,q,nq,parmute1,parmute2);
  fprintf(stderr,"generate stack+mute+RT operator\n");
  stack_vector(stackop,nh);
  fprintf(stderr,"eps=%f\n",eps);
  for (freq=1;freq<maxfreq;freq++){
    fprintf(stderr,"freq=%d\n",freq);
    w=2*PI*freq*df;
    radon_matrix(L,g,q,nh,nq,w);
    LSmatrix(L,RT,nh,nq,eps); 
    diagxtimesA(RT,mute,RT,nq,nh);
    AtimesBm(IRT,L,RT,nh,nq,nh);
    xtimesA(finalop,stackop,IRT,nh,nh);
    for (ih=0;ih<nh;ih++) F[freq][ih]=finalop[ih];
  }
  
  free1complex(finalop);
  free1float(stackop);
  free1float(mute);
  free2complex(IRT);
  free2complex(RT);
  free2complex(L);
  free1float(g);

  return;
}

void stack_rteom_LS_operators(complex ***LS, complex ***LT, int maxfreq, float df, 
			   float *h, int nh, float *q, int nq, float *t, int nt, float eps, 
			   float parmute1, float parmute2, float *mute, float *stackop)
{
  /* Variables with capital letters corresponds to matrices */
  /* This routine returns not only the mute+stack+operaor but also the 
     Radon LS operator (required only for testing purposes */

  int freq;
  float w=0;
  float *g=0;  
  int rtmethod=3;
  float depth=4000; // Not used
  complex **L, **LH;    // RT oprator
  complex **RT=0; // LS RT operator

  g=ealloc1float(nh);
  L=ealloc2complex(nq,nh); 
  RT=ealloc2complex(nh,nq);
  LH=ealloc2complex(nh,nq); 

  radon_moveout(h,g,nh,rtmethod,depth); 
  mute_vector(mute,q,nq,parmute1,parmute2);

  stack_vector(stackop,nh);

  fprintf(stderr,"eps=%f\n",eps);
  for (freq=1;freq<maxfreq;freq++){
    fprintf(stderr,"freq=%d\n",freq);
    w=2*PI*freq*df;
    radon_matrix(L,g,q,nh,nq,w);
    LSmatrix(L,RT,nh,nq,eps); 
    transpose(L,LH,nh,nq); 
    AequalB(LS[freq],RT,nh,nq);
    AequalB(LT[freq],L,nq,nh);
  }

  free2complex(RT);
  free2complex(LH);
  free2complex(L);
  free1float(g);

  return;
}


void testop(complex **F, int maxfreq, int nh)
{
  /* Create a unit test operator */

  int ih,freq;
  for (freq=1;freq<maxfreq;freq++){
    for (ih=0;ih<nh;ih++){
      F[freq][ih].r=1;
      F[freq][ih].i=0;
    }
  }
  return;
}


void LSmatrix(complex **F, complex **FLS, int nh, int nk, float eps)
{
  int ik;
  complex **FF=0;
  complex **FT=0;

  FT=ealloc2complex(nh,nk);
  FF=ealloc2complex(nk,nk);


  AHtimesB(FF,F,F,nh,nk,nh,nk);
  // regularization
  //  for (ik=0;ik<nk;ik++) FF[ik][ik]+=(eps*abs(FF[ik][ik]));
  for (ik=0;ik<29;ik++) FF[ik][ik]+=(eps*abs(FF[ik][ik]));
  for (ik=29;ik<36;ik++) FF[ik][ik]+=(eps*abs(FF[ik][ik]));
  for (ik=36;ik<nk;ik++) FF[ik][ik]+=(eps*abs(FF[ik][ik]));
  transpose(F,FT,nh,nk);
  
  inverse_matrix_multiply(nk,FF,nh,nk,FT,FLS);  

  free2complex(FT);
  free2complex(FF);


  return;
}

void mute_vector(float *mute,float *q, int nq, float parmute1, float parmute2)
{
  int iq;
  for (iq=0;iq<nq;iq++)
    if ((q[iq]>=parmute1)&&(q[iq]<=parmute2)) mute[iq]=1;  
    else mute[iq]=0;
  
  return;
}


void stack_vector(float *stackop, int nh)
{
  int ih;
  for (ih=0;ih<nh;ih++) stackop[ih]=1./nh;
  return;
} 










