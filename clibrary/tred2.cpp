#include <math.h>
//#include <clibrary.h>
//#include <Complex.h>
#include "su.h"
//#include "/usr/local/cwp/include/segy.h"
#include "Complex.h"
#include "clibrary.h"
//#include "Dcomplex.h"
void tred2(complex **b, int n, complex d[],complex e[])
{
	int l,k,j,i;
        float scale;
	complex hh,h,g,f;
        complex **a, czero;
        czero.r=czero.i=0;
	// Conversion from 0 index to 1  
        if ((a=(complex**) malloc((unsigned) n*sizeof(complex*)))==NULL)
             err("cannot allocate memory for a\n");
        
        for (i=1;i<=n;i++) *a[i]=*b[i-1];
        d=d-1;
	e=e-1;
        ////////
	for (i=n;i>=2;i--) {
		l=i;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += abs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(abs(f) >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=czero;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=czero;
	e[1]=czero;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i].r!=0) {
			for (j=1;j<=l;j++) {
				g=czero;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i].r=1.0;
                a[i][i].i=0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=czero;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */











