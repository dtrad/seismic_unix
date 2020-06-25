#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
#include "Dcomplex.h"
void compute_rhs( float w, int nx, float *g, complex *data, int np, 
		float pmin, float dp, complex *rhs)
/*********************************************************************
				     +
Compute the right-hand-side vector  B  data(x)

	+	    -i w p g(x)
where B   = (1/nx) e	        for equal increments in p as
used for the Discrete Radon Transform computation for
linear or parabolic tau-p.

Function parameters:

float w	input angular frequency of interest
int   nx	number of spatial positions ( defines length of g and data )
float g[]      spatial function corresponding to spatial locations of data
complex data[] data as a function of spatial position for a single
		angular frequency w as complex values 
int   np	number of output slownesses p (may be slowness squared
		or a more general function)
float pmin     starting value of output p
float dp	increment in output p
complex rhs[]  np complex values for the result 
*********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*********************************************************************/

{
	int ip, ix;
	float p, rsum, isum, dr, di, tr, ti, fac;
	fac=1./nx;
        //floatprint(w);
        //expfprint(dp);
        //intprint(np);
        //expfprint(pmin);
	for(ip=0;ip<np;ip++) {
		p = pmin + ip*dp;
		rsum = isum = 0.;
		for(ix=0;ix<nx;ix++) {
			tr = cos(w*p*g[ix]);
			ti = sin(w*p*g[ix]);
			dr = data[ix].r;
			di = data[ix].i;
			rsum += tr*dr - ti*di;
			isum += tr*di + ti*dr;
		}
		rhs[ip].r   = fac*rsum;
		rhs[ip].i   = fac*isum;
	}
}



