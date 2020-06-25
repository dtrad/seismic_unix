#include "su.h"
#include "segy.h"
#include "Complex.h"
void cholsl(complex **a, int n,double p[],complex b[],complex x[])
{
	int i,k;
	complex sum;

	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= conjg(a[k][i])*x[k];
		x[i]=sum/p[i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
void cholsl(float **a, int n, float p[], float b[], float x[])
{
	int i,k;
	float sum;

	for (i=1;i<=n;i++) {
		for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n;i>=1;i--) {
		for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
		x[i]=sum/p[i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
