#include "su.h"

//#define	MAX(x,y) ((x) > (y) ? (x) : (y))

float normalize(int n, float *a)
/********************************************************************  
return the  normalized vecotr
*********************************************************************/
{
	int j;
	float max=0.;
	for(j=0;j<n;j++) max = MAX(max,fabs(a[j]));
	if (max) for(j=0;j<n;j++) a[j]/=max;
	return(max);
}

float normalize(int n, complex *a)
/********************************************************************  
return the  normalized vecotr
*********************************************************************/
{
	int j;
	float max=0.;
	for(j=0;j<n;j++) max = MAX(max,abs(a[j]));
	//fprintf(stderr,"max=%f\n",max);
	if (max) for(j=0;j<n;j++) a[j]/=max;
	return(max);
}

float max(int n, complex *a)
/********************************************************************  
return the  max
*********************************************************************/
{
	int j;
	float max=0.;
	for(j=0;j<n;j++) max = MAX(max,abs(a[j]));
	fprintf(stderr,"max=%f\n",max);

	return(max);
}

float max(int n, float *a)
/********************************************************************  
return the  max
*********************************************************************/
{
	int j;
	float max=0.;
	for(j=0;j<n;j++) max = MAX(max,fabs(a[j]));
	fprintf(stderr,"max=%f\n",max);

	return(max);
}

float min(int n, complex *a)
/********************************************************************  
return the  min
*********************************************************************/
{
	int j;
	float min=FLT_MAX;
	for(j=0;j<n;j++) min = MIN(min,abs(a[j]));
	fprintf(stderr,"min=%f\n",min);

	return(min);
}

float min(int n, float *a)
/********************************************************************  
return the  min
*********************************************************************/
{
	int j;
	float min=FLT_MAX;
	for(j=0;j<n;j++) min = MIN(min,fabs(a[j]));
	fprintf(stderr,"min=%f\n",min);

	return(min);
}

float min_with_sign(int n, float *a)
/********************************************************************  
return the  min
*********************************************************************/
{
	int j;
	float min=FLT_MAX;
	for(j=0;j<n;j++) min = MIN(min,a[j]);
	fprintf(stderr,"min=%f\n",min);

	return(min);
}


void maxmin(float *x,int lx, float *pmin, float *pmax)
{
/********************************************************************  
return the  min and max of increments
*********************************************************************/

  /* 
     Given a vector x of length lx, with irregular increments 
     find the max and min of the increments.
  */
   int i;
   float dx, max, min;		
   max=0;
   min=1e10;	
   for (i=1;i<lx;i++){
     dx=x[i]-x[i-1];
     if (dx>max) max=dx;
     if (dx<min) min=dx;
   }
   *pmax=max;
   *pmin=min;
   return;
}






