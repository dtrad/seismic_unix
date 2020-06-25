#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "cwp.h"
#include "radoncgfft.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void radon_matrix_cgfft(float **L, float *RC, int nh, int nq, int nf2)

{
  int ih,iq;
  register int i;
  complex czero;czero.r=czero.i=0;
  complex *R;

  R=ealloc1complex(nq);
  for (iq=0;iq<nq;iq++){
    R[iq]=czero;
    for (ih=0;ih<nh;ih++) R[iq]+=lh[ih][0]*l[ih][iq]; //Top row of LL=LH*L
    //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
  }

  // Construct the circulant matrix//////////////////////////////////////
  for (i=0; i<nq;i++) RC[i]=conjg(R[i]);       
  for (i=nq; i<nf2;i++) RC[i]=czero; // pad with zeros the autocorrelation
  for (i=1;i<=nq-1;i++) RC[nf2-i]=conjg(RC[i]);      //use the DFT

  pfacc(-1,nf2,RC);
  ///////////////////////////////////////////////////////////////////////
  free1complex(R);

  return;
}



















