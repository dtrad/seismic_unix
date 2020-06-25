#include <math.h>
#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-20;

void ludcmp_c(complex **a, int n, int *indx, complex *d)
{
	int i,imax,j,k;
	complex dumc,sum;
	float *vv,big,dumr,temp;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=abs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dumr=vv[i]*abs(sum)) >= big) {
				big=dumr;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dumc=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dumc;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dumc=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dumc;
		}
	}
	free_vector(vv,1,n);
}
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
