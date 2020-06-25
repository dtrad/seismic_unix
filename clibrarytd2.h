#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)
void radtd(double *t, double *q, double *h, double *x,double *y,double eps, 
double qmin, double qmax);
void rstack(double *t,double *q,double *h,double *m,double *d,int conj);
double dot(int n, double *a, double *b);
void cgtd(double *t, double *q,double *h, double *x, double *y,
	  double *r,double *g, double *s, double *ss, double eps);
void cgtdhr(double *t, double *q,double *h, double *x, double *y,
	  double *r,double *g, double *s, double *ss, double *D, double eps);

