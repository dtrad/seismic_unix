#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)
void radtd(double *t, double *q, double *h, double *x,double *y,double eps, 
	   double qmin, double qmax);
void rstack(double *t,double *q,double *h,double *m,double *d,int conj);
double dot(int n, double *a, double *b);
void cgtdhr(double *t, double *q,double *h, double *x, double *y,double eps);
void cgtd(double *t, double *q,double *h, double *x, double *y, double eps);
void wtcgls(double *t, double *qaxis, double *h,double *x,double *b,double tol);
void wtcgls3(double *t, double *qaxis, double *h,double *x,double *b,double tol);
void lsqr(double *t,double *qaxis,double *h,double *x,double *d,double  tol,
	  int reorth);
