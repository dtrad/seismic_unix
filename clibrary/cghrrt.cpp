#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"

int cghrrt( int niter, int n, complex **a, complex *x, complex *y, 
	complex *s, complex *ss, complex *g, complex *rr)

/*********************************************************************

Hestenes and Stiefel conjugate gradient algorithm 
specialized for solving Hermitian Toeplitz
system.  a[] is input as a vector defining the only the
top row of A.  x[] is the solution vector returned.
y[] is input.  niter is the maximum number of conjugate 
gradient steps to compute.  The function returns as
the number of steps actually computed.  The other 
vectors provide workspace.

Complex Hermitian Toeplitz Solver for

N-1
Sum  A	     x  = y      for i=0,1,2,...,N-1
j=0   (i-j)   j    i

where A is Hermitian Toeplitz and x and y are complex.  For
an example 4 x 4 system,  x returns as the solution of


   A0  A1  A2  A3	x0	     y0

     *
   A1  A0  A1  A2	x1	     y1
				=    
     *   *
   A2  A1  A0  A1	x2	     y2

     *   *   *
   A3  A2  A1  A0	x3	     y3

********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*********************************************************************/
/* Added Notes;
g: is the gradient= A'*Ax-A'*y)=A'(Ax-y)=A'rr
s: is the conjugate gradient s=g+beta*s
rr: is the updated residuals (Ax_i+1-y)= r-alfa*A*s
ss: is the auxiliar vector As
rsq: norm  of residuals
*/
{
	int j, iter;
	complex czero;
	float alpha, beta, gamma, gammam, rsq, rp, test;
	float eps=1.0e-6;
	float rcdot(int n, complex *a, complex *b);

	rp   = rcdot(n,y,y);
	test = n*eps*eps*rp;
	czero.r=czero.i=0.;

	for(j=0;j<n;j++) {
		x[j]=czero;
		rr[j]=y[j];
	}

	Atimesx(g,a,rr,n,n);	   /*  adjoint matrix multiply */

	for(j=0;j<n;j++) s[j]=g[j];
	gammam=rcdot(n,g,g);

	for(iter=0;iter<niter;iter++) { /* forward matrix multiply  */
		  
        	Atimesx(ss,a,s,n,n);	   /*  adjoint matrix multiply */
		alpha  = gammam/rcdot(n,ss,ss);
		for(j=0;j<n;j++) {
			x[j] +=alpha*s[j];
			rr[j]-=alpha*ss[j];
		}
		rsq = rcdot(n,rr,rr);
		if ( iter>0 && ( rsq==rp || rsq<test ) ) return(iter-1);
		rp = rsq;
        	Atimesx(g,a,rr,n,n);	   /*  adjoint matrix multiply */
		gamma  = rcdot(n,g,g);
		//if (gamma<eps) break;
		beta   = gamma/gammam;
		gammam = gamma;			

		for(j=0;j<n;j++) {
			s[j] =g[j]+beta*s[j];
		}	
	}
return(iter);
}









