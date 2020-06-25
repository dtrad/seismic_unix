#include "dft2.h"

/*
Daniel Trad - December - 2000
*/

void dft_matrix(complex *R, complex **F,complex **FH,float *h,float *k,int nh,int nk, float *dh)
{

  /******************************************************************** 
     FFT Transformation matrices. F, FH and R= top row of LH*L
     Input parameters:
     Out parameter:
     Notes:
     Daniel Trad- 22-02-99
  *********************************************************************/
   
     int ih, ik;
     complex  co;
     complex  dco;
     complex phase, dphase;
     float dk=k[1]-k[0];
     float Aperture=0.0;

     for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

     for (ih=0;ih<nh;ih++){
       phase.r=dphase.r=0;
       phase.i=(h[ih]*(k[0]-dk));
       dphase.i=(h[ih]*dk);
       co=exp(phase);
       dco=exp(dphase);
       for (ik=0;ik<nk;ik++){
	 co*=dco;
	 F[ih][ik]=conjg(co);
	 FH[ik][ih]=(1./Aperture)*dh[ih]*co;
	 
       }
     }
  	      
     for (ik=0;ik<nk;ik++){
       R[ik].r=0;
       R[ik].i=0;
       for (ih=0;ih<nh;ih++)
	 R[ik]+=FH[0][ih]*F[ih][ik]; //Top row of LL=LH*L
       //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
     }
     
     return;
}

void dft_matrix(complex **F,float *h,float *k,int nh,int nk)
{

  /******************************************************************** 
     FFT Transformation matrices. F, FH and R= top row of LH*L
     Input parameters:
     Out parameter:
     Notes:
     Daniel Trad- 22-02-99
  *********************************************************************/
   
     int ih, ik;
     complex  co;
     complex  dco;
     complex phase, dphase;
     float dk=k[1]-k[0];

     for (ih=0;ih<nh;ih++){
       phase.r=dphase.r=0;
       phase.i=(h[ih]*(k[0]-dk));
       dphase.i=(h[ih]*dk);
       co=exp(phase);
       dco=exp(dphase);
       for (ik=0;ik<nk;ik++){
	 co*=dco;
	 F[ih][ik]=conjg(co);
       }
     }
  	      
     return;
}
