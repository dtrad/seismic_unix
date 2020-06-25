/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUDMORAD -Radon DMO                                                 ", 
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

segy tr,*trr; // for malloc use *trr

int nt; //Number of time samples
int nh;
int nq;
int nx;
int ny;
int rtmethod,itercg,iter_end,method,norm, reorth;
float dt,dh,dq,eps1,eps2, step, thres,theta;
int main(int argc, char **argv)
{
	
  //FILE *myfilep;
  int j,i, maxtr;
  float *d, *m;
  float *q, *t, *h, eps, qmin, qmax; 
  extern int nt, nh, nq,rtmethod, norm, itercg, iter_end, method, nx ,ny;
  extern int reorth;
  extern float dt,dh,dq, eps1,eps2, step, thres, theta;	
  fprintf(stderr,"**************************\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  
  //if((myfilep=fopen("radonpar","w"))==NULL)
  //  err("cannot open file=%s\n","radonpar");
  
  // Get parameters 
  if (!getparint("method", &method))  method = 1;
  if (!getparfloat("eps1", &eps1))  eps1 = 1;
  if (!getparfloat("eps2", &eps2))  eps2 = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("qmin", &qmin))  qmin = 0;
  if (!getparfloat("qmax", &qmax))  qmax = 1e-6;
  if (!getparfloat("step", &step))  step = .9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("itercg", &itercg))  itercg = 10;
  if (!getparfloat("thres", &thres))  thres = 0.1;  
  // freq=0 use Fnyquist form data
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparint("norm", &norm))  norm =10; 
  if (!getparint("reorth", &reorth))  reorth =1; 
  if (!getparfloat("theta", &theta))  theta =1;    
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
  
  //fprintf(stderr,"&d=%p,dd=%g\n",&d,d[15][20]);
  if ((m=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
  
  if ((q=ealloc1float(nq))==NULL)
    fprintf(stderr,"***Sorry, space for q could not be allocated\n");
  
  if ((h=ealloc1float(nh))==NULL)
    fprintf(stderr,"***Sorry, space for h could not be allocated\n");

  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");    
  // Because we want to use same struct array for data and model 
  // the maximun number of traces between them will be taken.
  
  maxtr= (nq>nh) ? nq : nh; 
  fprintf(stderr,"maxtr=%d\n",maxtr); 
  
  if ((trr=(segy*) malloc(maxtr*sizeof(segy)))==NULL)
    fprintf(stderr,"**Sorry, space for traces could not be allocated\n");
  
  //  if ((trr=(segy*) ealloc1(maxtr,sizeof(segy)))==NULL)
  //  fprintf(stderr,"**Sorry, space for traces could not be allocated\n");	
  for (j=0;j<nh-1;++j) h[j]=0;
  
  j=0;
  // Loop over traces 
  do {
    register int i;
    
    trr[j]=tr;
    h[j]=(float) tr.offset;
    //fprintf(stderr,"h[j]=%f\n",h[j]);
    for (i=0;i<nt;i++){
      ///d[j+nh*i]=(float) tr.data[i]; //if sort by offset-time
      d[i+nt*j]=(float) tr.data[i];    //if sort by time-offset 
      //	fprintf(stderr,"i=%d\n",i);
    }
    j++;
  } while (gettr(&tr));
  fprintf(stderr,"Before radtd nq=%d, nt=%d, nh=%d\n",nq,nt,nh);  
  nh=j;
  ny=nt*nh;
  nx=nt*nq;
  radtd(t,q,h,m,d,eps,qmin,qmax);
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  
  /*fwrite(model,sizeof(float),ntf*nq,stdout);*/
  
  j=0;
  do{
    
    if (j>=nh)  trr[j]=trr[0];  // otherwise replace it with first one
    trr[j].f2=q[j];  // copy radon parameter in tr.f2
    trr[j].ntr=nq;
    
    for (i=0;i<nt;i++)
      //trr[j].data[i]=m[j+nq*i];
      trr[j].data[i]=m[i+nt*j];
    
    puttr(&trr[j]);
    j++;
  } while(j<nq);
  
  // Binary file with radon parameter for other non su programs
  // fwrite(q,sizeof(float),nq,myfilep);
  fprintf(stderr," nq=%d, nt=%d, nh=%d, j=%d\n",nq,nt,nh,j);
  //fclose(myfilep);
  
  free(trr);
  free1float(t);
  free1float(h);
  free1float(q);
  free1float(m);
  free1float(d);
  
  return EXIT_SUCCESS;
}




















