#include "su.h"
#include "segy.h"
#include "clibrarytd.h"

/******************************************************************************

        Subroutine to compute the rho filter in frequenccy
         domain for the time domain inverse slant stack

******************************************************************************/
void rho_filter(int npoints, int nt, float dt, float *rho)
/******************************************************************************
Input Parameters:
npoints         number of point for the rho filter
nt              number of time samples
dt              time sampling interval

Output Parameters:
rho             1-D array of filter points

Credits:
written by Gabriel Alvarez, CWP
corrected by Bjoern E. Rommel, IKU, Petroleumsforskning
******************************************************************************/
{
	const int maxnpfa = 720720;
				/* maximum number that npfa can handle */
                                /* (see function npfa) */
        int it,iff; 		/* loop counters */
        int nfh,nfhp;  		/* half number of frequency coefficients */
	int ntfft;		/* time samples to pad */
	int nph=npoints/2;
        float f, df;            /* frequency and frequency sampling interval */
        complex *cx;            /* array of frequencies */

	/* oddness of npoints */
	if (2 * nph == npoints)   err ("npoints must be odd!\n");

	/* compare filter length with number of time samples */
	if (nph > nt)   
		err ("filter length larger than number of time samples!");

	/* compute padding factor */
	if (nt > maxnpfa)   
		err ("number of time samples too large for npfa!");
	ntfft = npfa(nt);

	/* allocate working space */
	cx = alloc1complex(ntfft);

	/* define constants */
	nfh = ntfft/2;
	nfhp = nfh+1;
	df = 1.0/(2*dt*nfh);

	/* compute filter coefficients */
	cx[0] = cmplx (0.0, 0.0);
	for (iff = 1, f = df; iff < nfhp; f = ++iff * df)
		cx[iff] = cx[ntfft-iff] = cmplx (f, 0.0);
  
	/* inverse Fourier transform from f to t */
	pfacc(-1,ntfft,cx);

	/* normalized output filter coefficients */
	rho[nph] = 1.0;
	for (it = 0; it < nph; it++)
		rho[nph-it-1] = rho[nph+it+1] = cx[it+1].r / cx[0].r;
  
	/* clean up */
	free1complex(cx);
}
