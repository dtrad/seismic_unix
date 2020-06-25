/* Driver for routine wt1 */
#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"


int wavfilt(int *data, int nmax)
{
	unsigned long i,nused;
	int itest,k;
	float *u,*v,*w,frac,thresh,tmp;

	u=vector(1,nmax);
	v=vector(1,nmax);
	w=vector(1,nmax);

	for (i=1;i<=nmax;i++) w[i]=v[i]=data[i-1];
	
	if (0){
	  printf("Enter k (4, -4, 12, or 20) and frac (0.0 to 1.0):\n");
	  //if (scanf("%d %f",&k,&frac) == EOF) break;
	  frac=FMIN(1.0,FMAX(0.0,frac));
	  itest=(k == -4 ? 1 : 0);
	  if (k < 0) k = -k;
	  //if (k != 4 && k != 12 && k != 20) continue;
	}

	k=4;
	itest=1;
	frac=0.001;

	if (!itest) pwtset(k);
	wt1(v,nmax,1,itest ? daub4 : pwt);
	for (i=1;i<=nmax;i++) u[i]=fabs(v[i]);
	thresh=select((int)((1.0-frac)*nmax),nmax,u);
	nused=0;
	for (i=1;i<=nmax;i++) {
	  if (fabs(v[i]) <= thresh)
	    v[i]=0.0;
	  else
	    nused++;
	}
	wt1(v,nmax,-1,itest ? daub4 : pwt);
	for (thresh=0.0,i=1;i<=nmax;i++)
	  if ((tmp=fabs(v[i]-w[i])) > thresh) thresh=tmp;
	printf("k,nmax,nused= %d %d %d\n",k,nmax,nused);
	printf("discrepancy= %12.6f\n",thresh);
	for (i=1;i<=nmax;i++) data[i-1]=v[i];
	
	free_vector(w,1,nmax);
	free_vector(v,1,nmax);
	free_vector(u,1,nmax);
	return 0;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */





