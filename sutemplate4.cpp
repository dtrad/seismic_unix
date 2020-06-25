
#include "Complex.h"
#include "clibrarytd.h"
#include <time.h>

/* Maximum value an `unsigned short int' can hold.  (Minimum is 0).  */
#ifndef USHRT_MAX
#define USHRT_MAX 65535
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADINT Forward  High Resolution Hyperbolic Radon transform        ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradint< stdin > stdout [optional parameters]          		",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  "             0  Variable velocity at every q trace                   ",
  "                and variable delta q - Output is offset-space        ",
  " rtmethod=3  1-LRT 2-PRT 3-HRT  (only for method 2 and 5)            ",
  "                                                                     ",
  " eps1 =1		data variance (for Weight function in WTCGLS)   ",
  " eps2 =1             model variance (for Weight function in WTCGLS)  ",
  " eps=1e-7            small number for conjugate gradient             ",
  " step=0.9            the step for CG is shortened by step            ",
  " itercg = 10	        number of internal iterations  		        ",
  " iter_end =1		number of external iterations          	        ",
  " pervmin = 10        minimum perturbation on velocity                ",
  " dperv=0.1           defines the rate of incresing velocity spacing  ",
  " nq=20		Number of Radon Traces                         	",
  " factor=0.8		multiply dq critical by factor                	",
  " precond=0           Uses precondtioning                             ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " parnmute=1      if q > parmute the radon trace is tapered           ",
  " keepdata=1      if there are traces with equal offset than the      ",
  "                 interpolated ones, keep the original traces.        ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : Iterpolated traces  (offset time domain)                   ",
  " Input : sudata file  (offset time domain)              		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradint pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudataint  ", 
  "                                                                     ",
  " If nh output is > nh input only the header words of the first nh    ",
  "  traces are preserved                                               ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/

segy tr,tri; 
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for header storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	
  FILE *myfilep; 
  cwp_String modelfile=""; /* output sufile for the model */  
  time_t start,finish;
  double elapsed_time;
  int j,i,iq,ih;
  int it;
  float **d;
  float **m;
  float *q;
  float *t;
  float *h;
  float eps;
  float qmin;
  float qmax;
  float t0;
  int model; 
  int nt;
  int nh;
  int nq;
  int rtmethod;
  int norm;
  int itercg;
  int iter_end;
  int method;
  int nx;
  int ny;
  float dt,dh,dq, eps1,eps2, step, thres, theta;
  float factor; // multiply dq by factor
  int testadj;
  
  
  //////////////////////////////////////////////
  fprintf(stderr,"**************************\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("eps1", &eps1))  eps1 = 1;
  if (!getparfloat("eps2", &eps2))  eps2 = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("step", &step))  step = .9;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("itercg", &itercg))  itercg = 10;
  if (!getparint("testadj", &testadj))  testadj =0;
  if (!getparfloat("factor", &factor))  factor =0.8;

  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;

  if (!getparfloat("t0",&t0)) t0=0;

  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  fprintf(stderr,"nh=%d,nt=%d,dt=%f,eps2=%e\n",nh,nt,dt,eps2); 

  // Allocate memory for data and model
  d=ealloc2float(nt,nh);
  m=ealloc2float(nt,nq);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  for (ih=0;ih<nh;ih++) h[ih]=0;

  headerfp = etmpfile();
  tracefp = etmpfile();
  if (verbose) warn("using tmpfile() call");

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////////////////////////
  
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ntr = 0;
  do {
    ntr++;
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);
  ih=0;
  nh=ntr;

  for (ih=0;ih<nh;ih++){
    efread(tr.data, FSIZE, nt, tracefp);
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    efread(&tr,HDRBYTES,1,headerfp);
    h[ih]=tr.offset;
  }
  int cdpgather=tr.cdp;

  ny=nt*nh;
  nx=nt*nq;

  radint(t,q,h,h2,m[0],d[0],dint[0],eps,qmin,qmax,fmax,velint,dperv,pervmin,
	 alum,blum,clum,norm,t0,mute,parmute);

  if (rhofilter) rhofilt(dint[0],t,h2,nt,nh2,dt);
      
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);

  
  if ((myfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);
  iq=0;

  do{
    tr.f2=q[iq];  // copy velocity perturbation
    tr.ntr=nq; 
    memcpy((void *) tr.data,(const void *) model[iq],nt*sizeof(float));   
    fputtr(myfilep,&tr);
    iq++;
  } while(iq<nq);
  
  
  for(ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
    puttr(&tr);
  }
  
  if (verbose) fprintf(stderr,"nt=%d, nh=%d \n ",nt,nh);

  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free1float(t);
  free1float(h);
  free1float(q);
  free2float(m);
  free2float(d);
  efclose(headerfp);
  efclose(tracefp);
  efclose(myfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}












































