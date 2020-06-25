#include "su.h"
#include "segy.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#include "Complex.h"

void fftgo(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0)
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   int nfft;		/* transform length			*/
   int nf;		/* number of frequencies		*/
   float d1;		/* output sample interval in Hz		*/   
   //   int *nf0; is the (address of) number of frequencies
   //              used to go back in fftback

   nfft = npfaro(nt, LOOKFAC * nt);
   if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
   nf = nfft/2 + 1;
   d1 = 1.0/(nfft*dt);

   if ((rt = ealloc1float(nfft))==NULL)
     err("cannot allocate memory for rt\n");
   if ((ct = ealloc1complex(nf))==NULL)
     err("cannot allocate memory for ct\n"); 

   fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){

       for (it=0;it<MIN(nt,nfft);it++)
	 rt[it]=d[it][ih];
       
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++)
           m[it][ih]=ct[it];
       /*
       if (nf<nt){
          for (it=nf;it<nt;it++){
	    //fprintf(stderr,"In fftgo  it=%d \n",it);
            m[it][ih].r=0;
	    m[it][ih].i=0;
          }
       }
       */
    }
   *nf0=nfft+2;
   fprintf(stderr,"nf0=%d\n",*nf0);   
   
   free1complex(ct);
   free1float(rt);
   
   return;

}


