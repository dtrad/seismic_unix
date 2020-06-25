#include "radonsolver.h"
void radon_matrix_toep(complex **L, complex *R, int nh, int nq);
int ctoeplitz( int n, complex *r, complex *a, complex *b, complex *f, complex *g );

float radon_toeplitz(complex *d, complex **L, complex *m, float sigmad, int nh, int nq) 
{
  complex *rtoep, *ftoep, *gtoep, *madj, *dp, *res;
  int iter, i;
  float resid;
  

  madj=ealloc1complex(nq);
  rtoep=ealloc1complex(nq);
  ftoep=ealloc1complex(nq);
  gtoep=ealloc1complex(nq); 
  dp=ealloc1complex(nh); 
  res=ealloc1complex(nh);  

  radon_matrix_toep(L,rtoep,nh,nq);
  Atimesx(d,L,madj,nh,nq,TRUE);
  //fprintf(stderr,"rtoep[0].r=%f, rtoep[0].i=%f \n",rtoep[0].r,rtoep[0].i);
  rtoep[0].r*=(1.+sigmad);
  iter=ctoeplitz(nq,rtoep,m,madj,ftoep,gtoep);
  if (iter!=nq) warn("titer != nq ");

  Atimesx(dp,L,m,nh,nq,FALSE);
  
  for (i=0;i<nh;i++) res[i]=d[i]-dp[i];
  resid=rcdot(nh,res,res);

  free1complex(res);
  free1complex(dp);
  free1complex(gtoep);
  free1complex(ftoep);
  free1complex(rtoep);
  free1complex(madj); 

  
  return(resid);

}

void radon_matrix_toep(complex **L, complex *R, int nh, int nq)
{
  register int ih;
  int iq;
  complex czero; czero.r=czero.i=0;
  for (iq=0;iq<nq;iq++){
    R[iq]=czero;
    for (ih=0;ih<nh;ih++) R[iq]+=conjg(L[ih][0])*L[ih][iq]/nh; //Top row of LL=LH*L
    //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
  }

  return;
}


int ctoeplitz( int n, complex *r, complex *a, complex *b,
		 complex *f, complex *g )
/***********************************************************************	
Complex Hermitian Toeplitz Solver for

N-1
Sum  R	     A  = B      for i=0,1,2,...,N-1
j=0   (i-j)   j    i

where R is Hermitian Toeplitz and A and B are complex.  For
an example 4 x 4 system,  A returns as the solution of


   R0  R1  R2  R3	A0	     B0

     *
   R1  R0  R1  R2	A1	     B1
				=    
     *   *
   R2  R1  R0  R1	A2	     B2

     *   *   *
   R3  R2  R1  R0	A3	     B3

and


   R0  R1  R2  R3	F0	     1

     *
   R1  R0  R1  R2	F1	     0
				=    
     *   *
   R2  R1  R0  R1	F2	     0

     *   *   *
   R3  R2  R1  R0	F3	     0


***********************************************************************	
where the function parameters are defined by

n     dimension of system
*r    provides the top row of the Hermitian Toeplitz matrix R 
*a    returns the complex solution vector A
*b    input as complex vector B (not changed during call)
*f    returns the complex spiking filter F
      (may be needed later for Simpson's sideways recursion
      if do search for optimum filter lag)
*g    work space of length n complex values

The function value returns as the number of successfully
computed complex filter coefficients (up to n) if successful or
0 if no coefficients could be computed.
***********************************************************************	
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
***********************************************************************/

{
	float er, ei, vr, vi, cr, ci, vsq;
	int j;	  	/*  for the jth iteration, j=0,n-1 	*/
	int k;			/*  for the kth component, k=0,j-1 	*/
	int jmk;		/*  j-k 				*/
	if (r[0].r==0.) return 0;

	f[0].r = 1.0/r[0].r;
	f[0].i = 0.;
	a[0].r = b[0].r/r[0].r;
	a[0].i = b[0].i/r[0].r;
	vr=1.;
	vi=0.;
	vsq=1.;

	for(j=1;j<n;j++) {     	/* iteration loop for iteration j	*/
	/*  	Compute spiking filter that outputs {v,0,0,...} 
		for this iteration step j			*/
		f[j].r=0.;
		f[j].i=0.;
		er=ei=0.;
		for(k=0;k<j;k++) {
			jmk=j-k;
			er+=r[jmk].r*f[k].r+r[jmk].i*f[k].i;
			ei+=r[jmk].r*f[k].i-r[jmk].i*f[k].r;
		}
		cr  = (er*vr - ei*vi)/vsq;
		ci  = (er*vi + ei*vr)/vsq;
		vr  = vr - (cr*er+ci*ei);
		vi  = vi + (cr*ei-ci*er);
		vsq =  vr*vr + vi*vi;
		if (vsq <= 0.) break;
		for(k=0;k<=j;k++) {
			jmk=j-k;
			g[k].r = f[k].r - cr*f[jmk].r - ci*f[jmk].i;
			g[k].i = f[k].i + cr*f[jmk].i - ci*f[jmk].r;		
		}
		for(k=0;k<=j;k++) {
			f[k]=g[k];
		}

		/*  Compute shaping filter for this iteration */
		a[j].r=0.;
		a[j].i=0.;
		er=ei=0.;
		for(k=0;k<j;k++) {
			jmk=j-k;
			er+=r[jmk].r*a[k].r+r[jmk].i*a[k].i;
			ei+=r[jmk].r*a[k].i-r[jmk].i*a[k].r;
		}
		er  = er-b[j].r;
		ei  = ei-b[j].i;
		cr  = (er*vr - ei*vi)/vsq;
		ci  = (er*vi + ei*vr)/vsq;
		for(k=0;k<=j;k++) {
			jmk=j-k;
			a[k].r += - cr*f[jmk].r - ci*f[jmk].i;
			a[k].i += + cr*f[jmk].i - ci*f[jmk].r;		
		}	
	}

	/* Properly normalize the spiking filter so that R F = {1,0,0,...} */
	/* instead of {v,0,0,...}.  To be accurate, recompute vr,vi,vsq */ 
	vr=vi=0.;
	for(k=0;k<j;k++) {
		vr+=r[k].r*f[k].r-r[k].i*f[k].i;
		vi+=r[k].r*f[k].i+r[k].i*f[k].r;
	}

	vsq = vr*vr + vi*vi;

	/*  Compute (er,ei) = 1./(vr,vi)   */
	er = vr/vsq;
	ei = -vi/vsq;
	for(k=0;k<j;k++) {
		f[k].r = er*f[k].r - ei*f[k].i;
		f[k].i = er*f[k].i + ei*f[k].r;	
	}
	return (j);
}
