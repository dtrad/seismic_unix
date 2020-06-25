double dot(int n, double *a, double *b)
/********************************************************************  
return the  dot product
*********************************************************************/
{
	int j;
	double  sum=0.;
	for(j=0;j<n;j++) sum += a[j]*b[j];
	return(sum);
}
