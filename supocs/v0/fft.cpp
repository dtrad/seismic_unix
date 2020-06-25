#include "su.h"
#include "segy.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#include "Complex.h"

void fft_parameters(int nt, float dt, int *pnfft, int *pnf, float *pdf)
{
  /********** fft_parameters *************
    Inputs: Time parameters
       nt
       dt
     Outputs: Freq parameters
     *nf  (address of number of freq) 
     *df  (address of df)
  Daniel Trad - Dec 2000   
  ***************************************/
  int nfft;
  nfft = npfaro(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt
  if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);

  *pnfft=nfft;
  *pdf=1.0/(nfft*dt);
  *pnf = nfft/2 + 1;

  return;
}

void fftgo_xt2fx(int sign,float **d,complex  **m, int nh, int nt, float dt, int nfft, int nf)
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   // input int nfft;	      transform length			
   // int nf;		      number of frequencies		
   complex czero; czero.r=czero.i=0;
   int verbose=0;

   if (nf > nt) err("nf must <= nt\n");
   
   rt = ealloc1float(nfft);
   ct = ealloc1complex(nf);

   if (verbose) fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
       for (it=0;it<nt;it++) rt[it]=d[ih][it];
       for (it=nt;it<nfft;it++) rt[it]=0;
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++) m[it][ih]=ct[it];
       for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }
    
   free1complex(ct);
   free1float(rt);
   
   return;

}

int fftback_fx2xt(int sign,float **d,complex  **m, int nh, int nt, float dt, int nfft, 
		  int nf)
{
  register complex *ct;	/* complex input trace			*/
  register float *rt;	/* real output trace			*/
  float onfft;		/* 1.0/nfft				*/
  int ih;                 /*Index for traces */
  int verbose=0;

  onfft = 1.0/nfft;

  if (verbose) fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nh=%d nt=%d\n",nfft,
		       nf,onfft,sign,nh,nt);

  /* Allocate fft arrays */
  ct   = ealloc1complex(nf);
  rt   = ealloc1float(nfft);

  /* Main loop over traces */
  for (ih=0;ih<nh;++ih) {
    register int it;
    /* Load traces into ct (pfa fills in negative freqs) */
    for (it = 0; it < nf; ++it) ct[it] = m[it][ih];
    /* Inverse FFT */
    pfacr(sign, nfft, ct, rt);
    /* Load back and scale for inverse fft */
    for (it = 0; it < nt; it++) d[ih][it] = rt[it] * onfft;
    
  }
  free1complex(ct);
  free1float(rt);
  return EXIT_SUCCESS;
}



void fftgo_xt2fx(int sign,float *d,complex  **m, int nh, int nt, float dt, int nfft, int nf)
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   // input int nfft;	      transform length			
   // int nf;		      number of frequencies		
   complex czero; czero.r=czero.i=0;
   int verbose=0;

   if (nf > nt) err("nf must <= nt\n");
   
   rt = ealloc1float(nfft);
   ct = ealloc1complex(nf);

   if (verbose) fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
       for (it=0;it<nt;it++) rt[it]=d[ih*nt+it];
       for (it=nt;it<nfft;it++) rt[it]=0;
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++) m[it][ih]=ct[it];
       for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }
    
   free1complex(ct);
   free1float(rt);
   
   return;

}

int fftback_fx2xt(int sign,float *d,complex  **m, int nh, int nt, float dt, int nfft, 
		  int nf)
{
  register complex *ct;	/* complex input trace			*/
  register float *rt;	/* real output trace			*/
  float onfft;		/* 1.0/nfft				*/
  int ih;                 /*Index for traces */
  int verbose=0;

  onfft = 1.0/nfft;

  if (verbose) fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nh=%d nt=%d\n",nfft,
	  nf,onfft,sign,nh,nt);

  /* Allocate fft arrays */
  ct   = ealloc1complex(nf);
  rt   = ealloc1float(nfft);

  /* Main loop over traces */
  for (ih=0;ih<nh;++ih) {
    register int it;
    /* Load traces into ct (pfa fills in negative freqs) */
    for (it = 0; it < nf; ++it) ct[it] = m[it][ih];
    /* Inverse FFT */
    pfacr(sign, nfft, ct, rt);
    /* Load back and scale for inverse fft */
    for (it = 0; it < nt; it++) d[ih*nt+it] = rt[it] * onfft;
    
  }
  free1complex(ct);
  free1float(rt);
  return EXIT_SUCCESS;
}

/************* Fourier transform for a single trace ***********/
/**** Forward ******/
void fftgo_xt2fx(int sign,float *d,complex  *m, int nt, float dt, int nfft, int nf)
{
   int it;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   // input int nfft;	      transform length			
   // int nf;		      number of frequencies		
   complex czero; czero.r=czero.i=0;
   int verbose=0;

   if (nf > nt) err("nf must <= nt\n");
   
   rt = ealloc1float(nfft);
   ct = ealloc1complex(nf);

   if (verbose) 
     fprintf(stderr,"In fftgo nfft=%d, nf=%d nt=%d\n",nfft,nf,nt);   

   
   for (it=0;it<nt;it++) rt[it]=d[it];
   for (it=nt;it<nfft;it++) rt[it]=0;
   pfarc(sign, nfft, rt, ct);
   //for (it=0;it<nf;it++) ct[it]/=nfft;
   for (it=0;it<nf;it++) m[it]=ct[it];
   for (it=nf;it<nt;it++) m[it]=czero;       

    
   free1complex(ct);
   free1float(rt);
   
   return;

}
/**** Inverse ******/
int fftback_fx2xt(int sign,float *d,complex  *m, int nt, float dt, int nfft, int nf)
{
  register complex *ct;	/* complex input trace			*/
  register float *rt;	/* real output trace			*/
  float onfft;		/* 1.0/nfft				*/
  int verbose=0;

  onfft = 1.0/nfft;

  if (verbose) 
    fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nt=%d\n",
	    nfft,nf,onfft,sign,nt);

  /* Allocate fft arrays */
  ct   = ealloc1complex(nf);
  rt   = ealloc1float(nfft);

  /* Main loop over traces */
  register int it;
  /* Load traces into ct (pfa fills in negative freqs) */
  for (it = 0; it < nf; ++it) ct[it] = m[it];
  /* Inverse FFT */
  pfacr(sign, nfft, ct, rt);
  /* Load back and scale for inverse fft */
  for (it = 0; it < nt; it++) d[it] = rt[it] * onfft;
    
  
  free1complex(ct);
  free1float(rt);
  return EXIT_SUCCESS;
}


























