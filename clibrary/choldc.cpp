//#include "math.h"
#include "su.h"
#include "segy.h"
#include "Complex.h"
// Note that always j>i>k. Hence we access upper part of a and return
// L in the lower part of a. ( a=L*LT )
// See Numerical recipes

void choldc(complex **a, int n, double p[])
{
	
	int i,j,k;
	complex sum;
	//float sum;

	for (i=0;i<n;i++) {
	    for (j=i;j<n;j++) {
	        for (sum=conjg(a[i][j]),k=i-1;k>=0;k--) 
		          sum -=conjg(a[i][k])*a[j][k]; // actually L*LT
		       if (i == j) {
			  if (sum.r <= 0.0)
			     fprintf(stderr,"choldc failed");
			  p[i]=sqrt(abs(sum)); // diagonal term should be 
			        // real because a(i,i)=conjg(a(i,i))
		       } else a[j][i]=sum/p[i]; // This is in fact L
	     }
	}
	return;
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

#include <math.h>

void choldc(float **a, int n, float p[])
{
	void nrerror(char error_text[]);
	int i,j,k;
	float sum;

	for (i=1;i<=n;i++) {
		for (j=i;j<=n;j++) {
			for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					nrerror("choldc failed");
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */








