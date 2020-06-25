/* Driver for routine lubksb */
#include "c:\daniel\radon\flibrary\complex.h"
#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "c:\daniel\cpress\util\nr.h"
#include "c:\daniel\cpress\util\nrutil.h"
#include "c:\daniel\cpress\util\nrutil.c"
#include "c:\daniel\radon\flibrary\lubksb_c.c"
#include "c:\daniel\radon\flibrary\ludcmp_c.c"
#define NP 20
#define MAXSTR 80

int main(void)
{
	int j,k,l,m,n,*indx;
	complex p,*x,**a,**b,**c;
	char dummy[MAXSTR];
	FILE *fp;

	indx=ivector(1,NP);
//	complex x[NP+1],a[NP+1][NP+1],b[NP+1][NP+1],c[NP+1][NP+1];
	x=vectorc(1,NP);
	a=matrixc(1,NP,1,NP);
	b=matrixc(1,NP,1,NP);
	c=matrixc(1,NP,1,NP);
	if ((fp = fopen("c:/daniel/cpress/util/matrx1c.dat","r")) == NULL)
		nrerror("Data file matrx1.dat not found\n");
	while (!feof(fp)) {
		fgets(dummy,MAXSTR,fp);
		fgets(dummy,MAXSTR,fp);
		fscanf(fp,"%d %d ",&n,&m);
		fgets(dummy,MAXSTR,fp);
		for (k=1;k<=n;k++)
			for (l=1;l<=n;l++) {
				fscanf(fp,"%f %f",&a[k][l].r,&a[k][l].i);
				//printf("a[%d][%d]=%f%f\n",k,l,a[k][l].r,a[k][l].i);
			}
		fgets(dummy,MAXSTR,fp);
		for (l=1;l<=m;l++)
			for (k=1;k<=n;k++) {
				fscanf(fp,"%f ",&b[k][l].r);
				b[k][l].i=0;
				//printf("b[%d][%d]=%f%f\n",k,l,b[k][l].r,b[k][l].i);
			}
		b[1][1].i=1;
		/* Save matrix a for later testing */
		for (l=1;l<=n;l++)
			for (k=1;k<=n;k++) c[k][l]=a[k][l];
		/* Do LU decomposition */
		ludcmp_c(c,n,indx,&p);
		/* Solve equations for each right-hand vector */
		for (k=1;k<=m;k++) {
			for (l=1;l<=n;l++) x[l]=b[l][k];
			lubksb_c(c,n,indx,x);
			/* Test results with original matrix */
			printf("right-hand side vector:\n");
			for (l=1;l<=n;l++){
				printf("%12.6f",(b[l][k].r));
				printf("%12.6f",(b[l][k].i));
				printf("\n");
			}
			printf("\n%s%s\n","result of matrix applied",
				" to sol'n vector");
			for (l=1;l<=n;l++) {
				b[l][k]=0.0;
				for (j=1;j<=n;j++)
					b[l][k] += (a[l][j]*x[j]);
			}
			for (l=1;l<=n;l++){
				printf("%12.6f",b[l][k].r);
				printf("%12.6f",b[l][k].i);
				printf("\n");
			}
			printf("\n*********************************\n");
		}
		printf("press RETURN for next problem:\n");
		(void) getchar();
	}
	fclose(fp);
	free_matrixc(c,1,NP,1,NP);
	free_matrixc(b,1,NP,1,NP);
	free_matrixc(a,1,NP,1,NP);
	free_vectorc(x,1,NP);
	free_ivector(indx,1,NP);
	return 0;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
