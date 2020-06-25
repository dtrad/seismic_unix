#include "su.h"
#include "segy.h"
#include "Complex.h"
#define PFA_MAX	720720		/* Largest allowed fft	*/

// Based on suiift.c from Colorado Scholl of mines
// April 4, 1999 --- Daniel Trad. UBC
int fftback(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0)
{
	register complex *ct;	/* complex input trace			*/
	register float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	float onfft;		/* 1.0/nfft				*/
        int ih;                 /*Index for traces */
        
	nfft = nf0 - 2; /* nf0=2*number of frequencies */
	nf = nf0/2;

	onfft = 1.0/nfft;
	if (nfft < nt) err("nfft must >= nt\n");
 	if (nf > nt) err("nf must <= nt\n");

        fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nh=%d nt=%d\n",nfft,nf,onfft,sign,nh,nt);

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
		for (it = 0; it < nt; it++) d[it][ih] = rt[it] * onfft;

        }
	free1complex(ct);
	free1float(rt);
	return EXIT_SUCCESS;
}

/* The Following variant does the same but the output is a array of traces
instead of an array of time slices 
*/

int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0)
{
	register complex *ct;	/* complex input trace			*/
	register float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	float onfft;		/* 1.0/nfft				*/
        int ih;                 /*Index for traces */
        
	nfft = nf0 - 2; /* nf0=2*number of frequencies */
	nf = nf0/2;

	onfft = 1.0/nfft;
	if (nfft < nt) err("nfft must >= nt\n");
 	if (nf > nt) err("nf must <= nt\n");

        fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nh=%d nt=%d\n",nfft,nf,onfft,sign,nh,nt);

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
		/*Note the change here */
		for (it = 0; it < nt; it++) d[ih][it] = rt[it] * onfft;

        }
	free1complex(ct);
	free1float(rt);
	return EXIT_SUCCESS;
}


/* The Following variant does the same but the output is a array of traces
instead of an array of time slices 
*/

int fftback1(int sign,float *d,complex  **m, int nh, int nt, float dt, int nf0)
{
	register complex *ct;	/* complex input trace			*/
	register float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	float onfft;		/* 1.0/nfft				*/
        int ih;                 /*Index for traces */
        
	nfft = nf0 - 2; /* nf0=2*number of frequencies */
	nf = nf0/2;

	onfft = 1.0/nfft;
	if (nfft < nt) err("nfft must >= nt\n");
 	if (nf > nt) err("nf must <= nt\n");

	//        fprintf(stderr,"in fftback nfft=%d, nf=%d; onfft=%f, sign=%d nh=%d nt=%d\n",nfft,nf,onfft,sign,nh,nt);

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







