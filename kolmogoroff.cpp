#include "su.h"
#include "Complex.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void kolmogoroff(int nt, float *x) 
{
  register int it;
  int nfft;
  complex czero;
  complex *cx;

  czero.r=czero.i=0;
  nfft=npfa(nt);
  

  if ((cx=alloc1complex(nfft))==NULL)
    err("cannot allocate memory for cxz");

  for (it=0;it<nt;it++) {cx[it].r=x[it];cx[it].i=0;}
  for (it=nt;it<nfft;it++) cx[it]=czero;

  pfacc(1,nfft,cx);  
  for (it=0;it<nfft;it++) cx[it]*=sqrt(1./nfft);
   
  for (it=0;it<nfft;it++) cx[it]=cx[it]*conjg(cx[it]); 

  for (it=0;it<nfft;it++) cx[it]=log(cx[it]);

  pfacc(-1,nfft,cx);
  
  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);

  cx[0]/=2.;

  cx[nfft/2]/=2.;

  for (it=nfft/2+1;it<nfft;it++) cx[it]=czero;

  pfacc(1,nfft,cx);  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);

  for (it=0;it<nfft;it++)  cx[it]=exp(cx[it]);
   
  pfacc(-1,nfft,cx);  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);

  for (it=0;it<nt;it++) x[it]=cx[it].r; 

  free1complex(cx);

  return;
}







