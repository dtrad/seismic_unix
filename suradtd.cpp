/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADTD:  $Date: June 1999  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrarytd.h"
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADTD -Forward  High Resolution Hyperbolic Radon transform        ", 
  "	   Program in development                                 	",
  " 	   								",
  " suhrrt2 < stdin > stdout [optional parameters]          		",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  " 		0 Inverse Radon Stack	                                ",
  "		1  WTCGLS with Semblance Mask                           ",
  "             2  WTCGLS with Cm function to get Sparse Transform      ",
  "             3  LSQR  (Does not work well yet                        ",
  "             4  Tests                                                ",
  "             5  rho filter with WTCGLS for iterations                ",
  "                (use iter_end=0) for only rho filter                 ",
  "             6  Hilbert Transform for Cm instead of abs(m)           ",
  " rtmethod=3  1-LRT 2-PRT 3-HRT  (only for method 2 and 5)            ",
  "                                                                     ",
  " eps1 =1		data variance (for Weight function in WTCGLS)   ",
  " eps2 =1             model variance (for Weight function in WTCGLS)  ",
  " eps=1e-7            small number for conjugate gradient             ",
  " step=0.9            the step for CG is shortened by step            ",
  " itercg = 10	        number of internal iterations  		        ",
  " iter_end =1		number of external iterations          	        ",
  " qmin =0             minimum Radon parameter in (sec/offset)^2       ",
  " qmax =1e-6          maximum Radon parameter in (sec/offset)^2       ",
  " thres=0.1           threshold for mask function with Semblance      ",
  "                     Semblance is computed and only those values     ",
  "                     greater than thres take part in the computation ",
  " nq=100		Number of Radon Traces                         	",
  " factor=0.8		multiply dq critical by factor                	",
  " crude=0             crude=1 test simple midpoint Radon (no sinc)    ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : traces with header in the Radon domain.                  	",
  " Input : sudata file                                    		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradtd method=1  qmin=0.e-8 qmax=1e-8 nq=100 < sudata > sudatarad  ", 
  "                                                                     ",
  " key=f2 contains the radon parameter                                 ",
  " If nh < np ( as usual) header words of first nh traces are preserved",
  " and radon parameter is kept in f2 header word. Hence, when suradtdi ",
  " To see what is the mask for method 1 type:                          ",
  " ximage < semblance n1=$NT                                           ",
  " where semblance is an output file created when method 2 is selected ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 * Trace header fields modified f2
 */
/**************** end self doc ***********************************/

segy tr; 
float factor;
int nt,nh,nq,nx,ny,rtmethod,itercg,iter_end,method,norm, reorth;
float dt,dh,dq,eps1,eps2, step, thres,theta;
int testadj;
int taper;
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;
int crude;

int main(int argc, char **argv)
{
	
  //FILE *myfilep;
  int j,i,iq,ih;
  int it;
  float *d, *m;
  float *q, *t, *h, eps, qmin, qmax, fmax;
  float t0=0;
  int smooth; 
  extern int nt, nh, nq,rtmethod, norm, itercg, iter_end, method, nx ,ny;
  extern int reorth;
  extern float dt,dh,dq, eps1,eps2, step, thres, theta;
  extern float factor; // multiply dq by factor
  extern int testadj;
  extern int taper;
  extern int crude;
  /// Velocity Trend
  float *tvel;
  float *vel;
  float *velint;
  int ntvel;
  int nvel;
  int itv;
  float pervmin;
  float dperv;
                    	
  // smoothing
  int nl=3;  //  npoints left hand side
  int nr=3;  //  npoints left hand side
  int flag=2;  // 1 rectangular, 2 triangular
  //////////////////////////////////////////////
  fprintf(stderr,"**************************\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
    
  // Get parameters 
  if (!getparint("method", &method))  method = 1;
  if (!getparfloat("eps1", &eps1))  eps1 = 1;
  if (!getparfloat("eps2", &eps2))  eps2 = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("qmin", &qmin))  qmin = 0;
  if (!getparfloat("qmax", &qmax))  qmax = 1e-6;
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("step", &step))  step = .9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("itercg", &itercg))  itercg = 10;
  if (!getparfloat("thres", &thres))  thres = 0.1;  
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparfloat("fmax", &fmax))  fmax =0; // max freq used for dq
  if (!getparint("norm", &norm))  norm =10; 
  if (!getparint("reorth", &reorth))  reorth =1; 
  if (!getparfloat("theta", &theta))  theta =1;
  if (!getparint("smooth", &smooth))  smooth =0;
  if (!getparint("testadj", &testadj))  testadj =0;
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparint("taper", &taper))  taper =0;    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("crude",&crude)) crude=0;

  /* Introduce velocity trend to apply Hyp Vel Filtering */
  ntvel = countparval("tvel");
  if (ntvel==0) ntvel = 1;
  tvel = ealloc1float(ntvel);
  if (!getparfloat("tvel",tvel)) tvel[0] = 0.0;
  nvel = countparval("vel");
  if (nvel==0) nvel = 1;
  if (nvel!=ntvel) err("number of tmig and vmig must be equal");
  vel = ealloc1float(nvel);
  if (!getparfloat("vel",vel)) vel[0] = 2000.0;
  for (itv=1; itv<ntvel; ++itv)
    if (tvel[itv]<=tvel[itv-1])
      err("tvel must increase monotonically");
  ////////////////////////////////////////////////////////   
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
  
  if ((d=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  
  if ((m=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
  
  if ((q=ealloc1float(nq))==NULL)
    fprintf(stderr,"***Sorry, space for q could not be allocated\n");
  
  if ((h=ealloc1float(nh))==NULL)
    fprintf(stderr,"***Sorry, space for h could not be allocated\n");

  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");    

  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"*Sorry, space for velint could not be allocated\n"); 
      	
  for (ih=0;ih<nh-1;++ih) h[ih]=0;

  headerfp = etmpfile();
  if (verbose) warn("using tmpfile() call");  

  ih=0;
  // Loop over traces 
  do {
    register int i;
    efwrite(&tr,HDRBYTES,1,headerfp);    
    h[ih]=(float) tr.offset;
    for (i=0;i<nt;i++){
      d[ih*nt+i]=(float) tr.data[i];    //if sort by time-offset 
    }
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));

  erewind(headerfp);
  nh=ih;

  /* Time axis */
  for (i=0;i<nt;i++) t[i]=t0+i*dt;


  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);

  //if (verbose) for (it=0;it<nt;it++) fprintf(stderr,"velint[%d]=%f,t[%d]=%f\n",it,velint[it],it,t[it]);

  if (verbose) fprintf(stderr,"Before radtd nq=%d, nt=%d, nh=%d, eps=%f \n",nq,nt,nh,eps);  
  ny=nt*nh;
  nx=nt*nq;
  radtd2(t,q,h,m,d,eps,qmin,qmax,fmax,velint,dperv,pervmin);

  
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  
  if (smooth) smoothing(m,nt,nq,nl,nr,flag);
   
  iq=0;
  do{
    if (iq<nh)  efread(&tr,HDRBYTES,1,headerfp);
    tr.f2=q[iq];  // copy radon parameter in tr.f2
    //if (method==8) tr.f2=vgrid[iq][0];
    tr.ntr=nq;
    
    for (i=0;i<nt;i++)
      tr.data[i]=m[iq*nt+i];    
    puttr(&tr);
    iq++;
  } while(iq<nq);
  
  // Binary file with radon parameter for other non su programs

  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free1float(velint);  
  free1float(t);
  free1float(h);
  free1float(q);
  free1float(m);
  free1float(d);
  free1float(vel);
  free1float(tvel);
  efclose(headerfp);
  
  return EXIT_SUCCESS;
}




















