#include "su.h"
#include "segy.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#include "Complex.h"

void fftgo(int sign,float **d,complex  **m, int nh, int nt, float dt, 
	    int *nf0)
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   int nfft;		/* transform length			*/
   int nf;		/* number of frequencies		*/
   float d1;		/* output sample interval in Hz		*/   
   complex czero; czero.r=czero.i=0;

   nfft = npfaro(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt

   if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
   nf = nfft/2 + 1;
   d1 = 1.0/(nfft*dt);
   if (nf > nt) err("nf must <= nt\n");

   if ((rt = ealloc1float(nfft))==NULL)
     err("cannot allocate memory for rt\n");
   if ((ct = ealloc1complex(nf))==NULL)
     err("cannot allocate memory for ct\n"); 

   fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
       for (it=0;it<nt;it++) rt[it]=d[it][ih];
       for (it=nt;it<nfft;it++) rt[it]=0;
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++) m[it][ih]=ct[it];
       for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }
   *nf0=(2*nf); // number of frequencies 
   fprintf(stderr,"nf0=%d\n",*nf0);   
   
   free1complex(ct);
   free1float(rt);
   
   return;

}

/* 
Now the input is an array of traces, 
the output is an array of frequency slices
*/

void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, 
	    int *nf0)
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   int nfft;		/* transform length			*/
   int nf;		/* number of frequencies		*/
   float d1;		/* output sample interval in Hz		*/   
   complex czero; czero.r=czero.i=0;

   nfft = npfaro(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt

   if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
   nf = nfft/2 + 1;
   d1 = 1.0/(nfft*dt);
   if (nf > nt) err("nf must <= nt\n");

   if ((rt = ealloc1float(nfft))==NULL)
     err("cannot allocate memory for rt\n");
   if ((ct = ealloc1complex(nf))==NULL)
     err("cannot allocate memory for ct\n"); 

   fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
       for (it=0;it<nt;it++) rt[it]=d[ih][it];
       for (it=nt;it<nfft;it++) rt[it]=0;
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++) m[it][ih]=ct[it];
       for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }
   *nf0=(2*nf); // number of frequencies 
   fprintf(stderr,"nf0=%d\n",*nf0);   
   
   free1complex(ct);
   free1float(rt);
   
   return;

}


/* 
Now the input is an array of traces, arranged as a column vector 
the output is an array of frequency slices
*/


void fftgo1(int sign,float *d,complex  **m, int nh, int nt, float dt, 
	    int nf0) // <== 1st difference
{
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   int nfft;		/* transform length			*/
   int nf;		/* number of frequencies		*/
   float d1;		/* output sample interval in Hz		*/   
   complex czero; czero.r=czero.i=0;

   nfft= nf0 - 2; // Number of frequencies (pos and neg)
   nf= nf0/2;     // Number of positive freq required for computations

   /*
   nfft = npfaro(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt
   if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
   nf = nfft/2 + 1;
   */
   
   d1 = 1.0/(nfft*dt);
   if (nf > nt) err("nf must <= nt\n");

   if ((rt = ealloc1float(nfft))==NULL)
     err("cannot allocate memory for rt\n");
   if ((ct = ealloc1complex(nf))==NULL)
     err("cannot allocate memory for ct\n"); 

   //   fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
     for (it=0;it<nt;it++) rt[it]=d[ih*nt+it];  // <=== 2nd difference
     for (it=nt;it<nfft;it++) rt[it]=0;
     pfarc(sign, nfft, rt, ct);
     //for (it=0;it<nf;it++) ct[it]/=nfft;
     for (it=0;it<nf;it++) m[it][ih]=ct[it];
     for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }

   /*
   *nf0=(2*nf); // number of frequencies 
   fprintf(stderr,"nf0=%d\n",*nf0);   
   */

   free1complex(ct);
   free1float(rt);
   
   return;

}

