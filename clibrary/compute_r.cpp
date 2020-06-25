#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
#include "Dcomplex.h"
void compute_r( float w, int nx, float *g, int np, float dp, complex *r)
/*******************************************************************
Compute the top row of the Hermitian Toeplitz Matrix
			+
		  R = B B

		  i w p g(x)
where B = (1/np) e	    for equal increments in p as

     +           -i w p g(x)
and B = (1/nx) e 

as used for the Discrete Radon Transform computation for
linear or parabolic tau-p.


		 nx-1	i w j dp g(x )
r[j] = 1/(nx*np) Sum	e	    k
		 k=0
						  2
g(x ) is initialized to  x  for linear tau-p or x   for the parabolic transform
   k		          k		         k
prior to calling this routine.  The use of g is intended to emphasize that the
spatial locations do not have to be equally spaced for either method.
In general, this routine can be called for g specified as any function
of spatial position only.  For a more general function of x, dp will
not correspond to an increment in slowness or slowness squared but
rather to a more general parameter.

******************************************************************
Function parameters:

float w	input as angular frequency component of interest
int   nx      number of spatial positions stored in g
float g[]     spatial function for this Radon Transform
int   np      number of slowness (or slowness squared) components
float dp      increment in slownes (or slowness squared)
float r[]     output vector of { real r0, imaginary r0, real r1, 
	      imaginary r1, ...}
******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
******************************************************************/
{
	int j,k;
	float rsum, isum, fac;

	fac= 1./(nx*np);
	for(j=0;j<np;j++) {
		rsum=0.;
		isum=0.;
		for(k=0;k<nx;k++) {
			rsum = rsum+cos( w*j*dp*g[k] );
			isum = isum+sin( w*j*dp*g[k] );	
		}
		r[j].r    = fac*rsum;
		r[j].i    = fac*isum;
	}
}










