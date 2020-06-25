/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrary.h"
#include <time.h>

void hrrtf(float *pos, float **d ,float dt,float eps1,float eps2, float qmin, 
	   float **m, float *q, float dq, float freq, float eps,
	   float *bad, int taper);

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUHRRT2 -Forward  High Resolution Parabolic Radon transform           ", 
"	   Program in development                                 	",
" 	   								",
" suhrrt2 < stdin > stdout [optional parameters]          		",
" 									",
" 									",
" Optional parameters:							",
" method=1                                                      	",
" 		1 Normal RT (Gauss- Gauss) Levinson                     ",
"               2 HRRT      (L1-Gauss or Cauchy-Gauss)			",
"                   WTCGLS  solution  (Fast HRRTF)                      ",
"               3 HRRT      (L1-Gauss or Cauchy-Gauss)			",
"                   Overdetermined solution with Cholesky               ",
"               4 HRRT      (L1-Gauss or Cauchy-Gauss)			",
"                   Underdetermined solution with Cholesky              ",
"               5 Normal RT (Gauss- Gauss) with CGLS     		",
"                                                                       ",
"               14 HRRT      (L1-Gauss or Cauchy-Gauss)			",
"                   WTCGLS  solution  (Fast HRRTF)                      ",
"               10 CG computed with circulant matrix algorihtm          ",
"                                                                       ",
"		All other methods are exeprimental versions that        ",
"               do not work very well yet.                              ",
"                                                                       ",
" eps1 =3		noise level in percent for data Covariance      ",
" eps2 =1e-3            variance to use in model Covariance matrix.     ",
"                       Not necesary in most cases                      ",
" eps=1e-7              small number to avoid division by zero          ",
"                       if = 0 with method 2 uses GCV as stopping      ",
"                       criteria  (approximated GCV)      		",
" iter = 3		number of external iterations         	        ",
" itercg= 10            number of internal iterations (for CG methods)  ", 
" qmin =0               minimum Radon parameter in sec/offset^2         ",
" 		                           				",
" nq=150		Number of Radon Traces                         	",
" rtmethod = 2          1-LRT 2-PRT                                     ",
" freq                  defines the freq max value used to compute      ",
"                       dq according to the formula                     ",
"                  PRT:    dq=0.8/(fmax*(xmax-xmin)^2)                  ",
"                  LRT:    dq=0.8*(1/fmax*(xmax-xmin))                  ",
"                       qmax is estimated as qmin+(nq-1)*dq             ",
" 									",
" step=0.5  (between 0-1)  step size to use in method 14 (CG).          ",
" costflag = 0            = 1 Performs extra computations to computing  ",
"                         the cost functions (slower).                  ", 
"                                                                       ",
" bad=none  Number of bad traces to attenuate with data Covariance      ",
" factor=.8  times dq critical for sampling Radon space                 ",
" Required parameters:		[None]					",
"                                                                	",
" Output : traces with header in the Radon domain.                     	",
" Input : sudata file                                    		",
"		                                        		",
" Example: 		                                                ",
" #Forward  Fast High Resoution Parabolic  Radon transform              ",
" suhrrtf method=2 rtmethod=2 eps1=1e-3 step=0.5 iter=3 itercg=10      ",
"         qmin=1.e-8  nq=100 freq=70   < sudata > sudatarad             ", 
"                                                                       ",
" key=f2 contains the radon parameter                                   ",
" If nh < np ( as usual) header words of first nh traces are preserved  ",
" and radon parameter is kept in f2 header word. Hence, when suhrrti2   ",
" is used to recover the data the header is preserved.                  ",
" To increase qmax increase nq or  qmin because dq is computed from     ",
" fmax, Xmax, Xmin.   	                                                ",
" If you are new to this method try first with method =3 (Levinson)     ",
" and then try to improve with High resolution methods 7 (slow Cholesky ",
" or 14 (WTCGLS) faster.                                                ",
" **************Important for method 2 ******************************** ", 
" Note: method 2  (wtcgls) usually very well but step parameter is      ",
"  critical. In the future it will be computed automatically.           ",
" If the results are bad compared with method 1                         ",
" (Cholesky) try reduce the step size. You will see this problem when   ",
" because rho increases instead of decreasing. The right step size is   ",
" when eps=0 GCV gives about 5-20 iterations before stopping.           ",
" If it only gives 2 iterations the step is too large.                  ", 
"                                                                       ",
NULL};
/* Credits:
 * 	Daniel Trad
 *
 * Trace header fields accessed: ns, dt, key=offset
 * Trace header fields modified f2
 
     method=1 Least Square Radon with Levinson
     method=2 Weighted Conjugate Gradient Least Squares
     method=3 Cholesky overdetermined
     method=4 Cholesky undertermined
     method=5 CGLS for Toeplitz (SU)
     method=6 SVD overdetermined  
     method=7 CG non linear with linear search
     method=8 LSQR 
     method=9 CG (Mauricio's method)
     method=10 Very fast CG with circulant matrices
     method=11 BCG from NR 
     method=12 Simplest CGLS   */

/**************** end self doc ***********************************/

segy tr;
complex czero=(0.0);
int nt, nh, nq, method, iter_end, rtmethod, norm, itercg, freqflag, costflag;
int reorth;
int nbad;
float step;
float factor;
int testadj;
float depth;
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	
	FILE *myfilep;
	time_t start,finish;
	double elapsed_time;
	int ih,iq,i;
	float **d, **m;
        float *q, eps1,eps2, eps, qmin, dq; 
        float *pos, dt, freq;
	extern int nt, nh, nq, method, iter_end, rtmethod, norm, 
	               itercg, freqflag, costflag, reorth;
        extern float step;
	float *bad;  /* pointer to list of traces to reject  */
	extern int nbad;   /* number of rejected traces		*/
        extern float factor;
	extern int testadj;
        int taper;
	// Foster and Mosher offset function
	extern float depth;
	float *gFM; 


	// Initialize 
	initargs(argc, argv);
	requestdoc(1);

	start=time(0);

        if((myfilep=fopen("radonpar","w"))==NULL)
                        err("cannot open file=%s\n","radonpar");

	// Get parameters 
	if (!getparint("costflag", &costflag))  costflag = 0; 
	if (!getparint("freqflag", &freqflag))  freqflag = 0; //=1 goes to Nyq
	if (!getparint("method", &method))  method = 1;
	if (!getparfloat("eps1", &eps1))  eps1 = 1e-2;
	if (!getparfloat("eps2", &eps2))  eps2 = 1e-2;
	if (!getparfloat("eps", &eps))  eps = 1e-7;
	if (!getparint("iter_end", &iter_end))  iter_end = 3;
	if (!getparfloat("qmin", &qmin))  qmin = -1e-7;
	if (!getparint("nq", &nq))  nq = 150;
	if (!getparint("itercg", &itercg))  itercg = 10;
	if (!getparfloat("freq", &freq))  freq =70; // freq=0 use Fnyquist
	if (!getparint("rtmethod", &rtmethod))  rtmethod =2; // PRT default
	if (!getparint("norm", &norm))  norm =1; // PRT default
	if (!getparfloat("step", &step))  step = .95;	
	if (!getparint("reorth", &reorth))  reorth =0;
	if (!getparfloat("factor", &factor))  factor = .8;
	if (!getparint("taper", &taper))  taper = 0;
	if (!getparint("testadj", &testadj))  testadj = 0;
	if (!getparfloat("depth", &depth))  depth = 1000;
	if (!getparint("verbose", &verbose))  verbose =0;
	/* Getpar the reject vector for covariance matrix*/
       	nbad = countnparval(1,"bad");
	if ((bad=ealloc1float(nbad))==NULL)
	  fprintf(stderr,"***Space for bad could not be allocated\n");
	if (!getnparfloat(1,"bad",bad)) bad[0]=999;        

        for (i=0;i<nbad;i++) fprintf(stderr,"bad[%d]=%f\n",i,bad[i]);

        if ((method!=2)&&(eps==0)) eps=1e-7;
        //fprintf(stderr,"eps=%e\n",eps);							
	// Get info from first trace 

	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ns) err("ns header field must be set");
	if (!tr.ntr) err("ntr header field must be set");
	dt   = ((float) tr.dt)/1000000.0;
	nt = (int) tr.ns;
        nh= (int) tr.ntr;
        if (verbose) fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt); 

	// Allocate memory for data and model

	if ((d=ealloc2float(nh,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
 
	if ((m=ealloc2float(nq,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");

	if ((q=ealloc1float(nq))==NULL)
	  fprintf(stderr,"***Sorry, space for q could not be allocated\n");

	if ((pos=ealloc1float(nh))==NULL)
	  fprintf(stderr,"***Sorry, space for pos could not be allocated\n");
       	
	for (ih=0;ih<nh-1;++ih) pos[ih]=0;

	headerfp = etmpfile();
	if (verbose) warn("using tmpfile() call");  

	ih=0;
	// Loop over traces 
	do {
   		register int i;
		efwrite(&tr,HDRBYTES,1,headerfp);    
      		pos[ih]=(float) tr.offset;
		for (i=0;i<nt;i++){
		        d[i][ih]=(float) tr.data[i];
		}
      		ih++;
	} while (gettr(&tr));
	erewind(headerfp); 

	nh=ih;	
	hrrtf(pos,d,dt,eps1,eps2,qmin,m,q,dq,freq,eps,bad,taper);
	if (verbose) fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",
			     nq,nt,nh);
	
	iq=0;
	do{
	  if (iq<nh)  efread(&tr,HDRBYTES,1,headerfp);
      	  tr.f2=q[iq];  // copy radon parameter in tr.f2
       	  tr.ntr=nq; 
	  for (i=0;i<nt;i++)
	    tr.data[i]=m[i][iq];
	  puttr(&tr);
	  iq++;
	} while(iq<nq);
 
	// Binary file with radon parameter for other non su programs
	
       	fwrite(q,sizeof(float),nq,myfilep);
        if (verbose) fprintf(stderr,"nq=%d,nt=%d,nh=%d,iq=%d\n",nq,nt,nh,iq);
	fclose(myfilep);

	free1float(bad);
	free1float(pos);
	free1float(q);
	free2float(m);
	free2float(d);
	efclose(headerfp);       

	finish=time(0);
	elapsed_time=difftime(finish,start);
	fprintf(stderr,"Total time required=====>: %f \n", elapsed_time);

	return EXIT_SUCCESS;
}




















