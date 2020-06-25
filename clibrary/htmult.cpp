#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
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





