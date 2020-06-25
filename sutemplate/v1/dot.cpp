float dot(int n, float *a, float *b)
/********************************************************************  
return the  dot product
*********************************************************************/
{
	int j;
	float sum=0.;
	for(j=0;j<n;j++) sum += (a[j]*b[j]);
	return(sum);
}
