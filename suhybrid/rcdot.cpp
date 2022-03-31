#include "su.h"
#include "segy.h"
#include "Complex.h"
float rcdot(int n, complex *a, complex *b)
/********************************************************************  
return the real part of a complex dot product where
    the first vector is the one complex conjugated
*********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*********************************************************************/
{
	int j;
	float sum=0.;
	for(j=0;j<n;j++) sum += a[j].r * b[j].r + a[j].i * b[j].i;
	return(sum);
}
