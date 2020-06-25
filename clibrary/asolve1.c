void asolve1(double **A,double b[],double x[],int n)
{
	int i;
        
	for(i=1;i<=n;i++) x[i]=(A[i][i] != 0.0 ? b[i]/A[i][i] : b[i]);
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
