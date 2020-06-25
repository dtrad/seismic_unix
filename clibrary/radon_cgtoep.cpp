#include "radonsolver.h"
void radon_matrix_toep(complex **L, complex *R, int nh, int nq);
void htmul(int n, complex *a, complex *x, complex *y);
int ctoephcg( int niter, int n, complex *a, complex *x, complex *y, 
	      complex *s, complex *ss, complex *g, complex *rr);

float radon_cgtoep(complex *d, complex **L, complex *m, float sigmad, int nh, int nq) 
{
  complex *rtoep, *s, *ss, *madj, *g, *rr;
  int iter;
  float resid;
  

  madj=ealloc1complex(nq);
  rtoep=ealloc1complex(nq);
  s=ealloc1complex(nq);
  ss=ealloc1complex(nq); 
  g=ealloc1complex(nq); 
  rr=ealloc1complex(nq);  

  radon_matrix_toep(L,rtoep,nh,nq);
  Atimesx(d,L,madj,nh,nq,TRUE);
  //fprintf(stderr,"rtoep[0].r=%f, rtoep[0].i=%f \n",rtoep[0].r,rtoep[0].i);
  rtoep[0].r*=(1.+sigmad);

  iter=ctoephcg(nq/7,nq,rtoep,m,madj,s,ss,g,rr);

  resid=rcdot(nq,rr,rr);

  free1complex(rr);
  free1complex(g);
  free1complex(ss);
  free1complex(s);
  free1complex(rtoep);
  free1complex(madj); 

  
  return(resid);

}

int ctoephcg( int niter, int n, complex *a, complex *x, complex *y, 
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
{
	int j, iter;
	complex czero;
	float alpha, beta, gamma, gammam, rsq, rp, test;
	float eps=1.0e-12;
	float rcdot(int n, complex *a, complex *b);

	rp   = rcdot(n,y,y);
	test = n*eps*eps*rp;
	czero.r=czero.i=0.;

	for(j=0;j<n;j++) {
		x[j]=czero;
		rr[j]=y[j];
	}

	htmul(n,a,rr,g);	   /*  adjoint matrix multiply */

	for(j=0;j<n;j++) s[j]=g[j];
	gammam=rcdot(n,g,g);

	for(iter=0;iter<niter;iter++) { /* forward matrix multiply  */
		htmul(n,a,s,ss);  
		alpha  = gammam/rcdot(n,ss,ss);
		for(j=0;j<n;j++) {
			x[j] +=alpha*s[j];
			rr[j]-=alpha*ss[j];
		}
		rsq = rcdot(n,rr,rr);
		if ( iter>0 && ( rsq==rp || rsq<test ) ) return(iter-1);
		rp = rsq;

		htmul(n,a,rr,g);   /*  adjoint matrix multiply  */
		gamma  = rcdot(n,g,g);
		if (gamma<eps) break;
		beta   = gamma/gammam;
		gammam = gamma;			

		for(j=0;j<n;j++) {
			s[j] =g[j]+beta*s[j];
		}	
	}
	return(iter);
}

void htmul(int n, complex *a, complex *x, complex *y)

/*******************************************************************
   Hermitian Toeplitz matrix multiply

     solve for y = A x   where A is Hermitian Toeplitz

     and defined by the vector a giving the top row of A.
     x is input.  y is output. 
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	int j,irow;
	complex czero;
	czero.r=czero.i=0.;

	for(irow=0;irow<n;irow++) {
		y[irow]=czero;
		for(j=0;j<irow;j++)
			y[irow] += conjg(a[irow-j])*x[j];
		for(j=irow;j<n;j++)
			y[irow] += a[j-irow]*x[j];
	}
}

