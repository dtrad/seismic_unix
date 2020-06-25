/* Copyright (c) University of Alberta, 2012.*/
/* All rights reserved.                       */
/* supocs5d  :  $Date: July     2012- Last version July    2012  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>


#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif



void save_gather(float **d, int nh, int nt, float dt, const char* name);
bool process(float** datain, float** dataout, int nt, int  nh, float dt, int x1_num, int x2_num, int x3_num, int x4_num, int dx1, int dx2, int dx3, int dx4);
bool process2(float** datain, float** dataout, int nt, int  nh, float dt);
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUPOCS5D  Reconstruction of 5D data using Projection Onto           ",
  "           Convex Sets. The input should be regularly sampled and    ",
  "           sorted by cdp, offset, gelev, selev.                      ",
  "           Make sure ntr is set.                                     ",
 NULL};
/* Credits:
 *	Aaron Stanton.
 * Trace header fields accessed: ns, dt, ntr, cdp, offset, gelev, selev
 * Last changes: July : 2012 
 */
/**************** end self doc ***********************************/



int main(int argc, char **argv)
{
  float *h;   
  int *x1;	    
  int *x2;
  int *x3;	    
  int *x4;   
  int verbose;
  segy tr; 
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep;

  time_t start,finish;
  double elapsed_time;
  int it,ih;
  float t0=0;
  float **datain    = 0;
  float **dataout   = 0;
  float **Wd        = 0;
  float *t;

  int nt, nh; 
  int method;
  int plot; 
  float dt;
  int dx1,dx2,dx3,dx4;
  float kmin;
  float fmax;
  char buf[80];
  ////////////////
    
  fprintf(stderr,"*******SUPOCS*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =1;
  if (!getparint("plot",&plot)) plot = 0;

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");
  //if (!tr.cdp) err("cdp header field must be set");
  //if (!tr.offset) err("offset header field must be set");
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh = (int) tr.ntr;
  if (!getparfloat("fmax",&fmax)) fmax = 0.5/dt;
  fmax = MIN(fmax,0.5/dt);

  // Allocate memory for data and model
  datain   = ealloc2float(nt,nh);
  dataout  = ealloc2float(nt,nh);
  Wd=ealloc2float(nt,nh);
  h=ealloc1float(nh);
  x1=ealloc1int(nh);
  x2=ealloc1int(nh);
  x3=ealloc1int(nh);
  x4=ealloc1int(nh);
  t=ealloc1float(nt);
  memset( (void *) h, (int) '\0', nh * FSIZE);
  // Loop over traces
  // Dead traces are marked with trid = 2. They are weighted to zero
  ih=0;

  do {
    if (tr.trid==2) for(it=0;it<nt;it++) Wd[ih][it]=0;
    else  for(it=0;it<nt;it++) Wd[ih][it]=1;
    x1[ih]=(int) tr.cdp;
    x2[ih]=(int) tr.offset;
    x3[ih]=(int) tr.gelev;
    x4[ih]=(int) tr.selev;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  
  if (verbose) fprintf(stderr,"processing %d traces \n", nh);
  for (it=0;it<nt;it++) t[it]=t0+it*dt;  /* Not implemented for t0 != 0  */
  
  // Figure out min/max values of x1, x2, x3, and x4 and the dimensions.
  int x1_min = x1[0];
  int x1_max = x1[0];
  int x2_min = x2[0]; 
  int x2_max = x2[0];
  int x3_min = x3[0]; 
  int x3_max = x3[0];
  int x4_min = x4[0]; 
  int x4_max = x4[0];
  int x1_num = 1;
  int x2_num = 1;
  int x3_num = 1;
  int x4_num = 1;
  int x1_prev = x1[0];
  int x2_prev = x2[0];
  int x3_prev = x3[0];
  int x4_prev = x4[0];
    
  for (ih=0;ih<nh;ih++){
    //    fprintf(stderr,"x1_prev=%d, x1[%d]=%d, x2_prev=%d, x2[%d]=%d, x3_prev=%d, x3[%d]=%d, x4_prev=%d, x4[%d]=%d\n",x1_prev,ih,x1[ih],x2_prev,ih,x2[ih],x3_prev,ih,x3[ih],x4_prev,ih,x4[ih]);
    if (x1[ih]<x1_min) x1_min = x1[ih];
    if (x1[ih]>x1_max) x1_max = x1[ih];
    if (x2[ih]<x2_min) x2_min = x2[ih];
    if (x2[ih]>x2_max) x2_max = x2[ih];
    if (x3[ih]<x3_min) x3_min = x3[ih];
    if (x3[ih]>x3_max) x3_max = x3[ih];
    if (x4[ih]<x4_min) x4_min = x4[ih];
    if (x4[ih]>x4_max) x4_max = x4[ih];
    if (x1[ih]!=x1_prev){ 
      x1_num = x1_num + 1;
      x2_num = 1;
      x3_num = 1;
      x4_num = 1;
    }
    if(x2[ih]!=x2_prev){ 
      //      fprintf(stderr,"=====================> x2_prev=%d, x2[%d]=%d\n",x2_prev,ih,x2[ih]);
      x2_num = x2_num + 1;
      x3_num = 1;
      x4_num = 1;
    }
    if(x3[ih]!=x3_prev){ 
      //      fprintf(stderr,"=====================> x3_prev=%d, x3[%d]=%d\n",x3_prev,ih,x3[ih]);
      x3_num = x3_num + 1;
      x4_num = 1;
    }
    if(x4[ih]!=x4_prev){ 
      //      fprintf(stderr,"=====================> x4_prev=%d, x4[%d]=%d\n",x4_prev,ih,x4[ih]);
      x4_num = x4_num + 1;
    }
    x1_prev = x1[ih]; 
    x2_prev = x2[ih];
    x3_prev = x3[ih];
    x4_prev = x4[ih];
  }

  MARK;
  fprintf(stderr,"x1_min=%d, x1_max=%d, x1_num=%d\n",x1_min,x1_max,x1_num);
  fprintf(stderr,"x2_min=%d, x2_max=%d, x2_num=%d\n",x2_min,x2_max,x2_num);
  fprintf(stderr,"x3_min=%d, x3_max=%d, x3_num=%d\n",x3_min,x3_max,x3_num);
  fprintf(stderr,"x4_min=%d, x4_max=%d, x4_num=%d\n",x4_min,x4_max,x4_num);

  fprintf(stderr,"x1_num=%d, x2_num=%d, x3_num=%d, x4_num=%d\n",x1_num,x2_num,x3_num,x4_num);

  dx1 = 1;
  dx2 = 1;
  dx3 = 1;
  dx4 = 1;
  if (x1_num > 1) dx1=(x1_max-x1_min)/(x1_num-1);
  if (x2_num > 1) dx2=(x2_max-x2_min)/(x2_num-1);
  if (x3_num > 1) dx3=(x3_max-x3_min)/(x3_num-1);
  if (x4_num > 1) dx4=(x4_max-x4_min)/(x4_num-1);
  MARK;
  
  if (plot){ /* additional plots only for debugging  */  
    save_gather(datain,nh,nt,dt,"datain.su");
    system("suximage < datain.su perc=98 key=offset title=datain xbox=%d curve=curve1 npair=5 &");  
  }

 fprintf(stderr,"nt=%d, nh=%d, dt=%f, x1_num=%d, x2_num=%d, x3_num=%d, x4_num=%d\n",nt,nh,dt,x1_num,x2_num,x3_num,x4_num);
 process(datain,dataout,nt,nh,dt,x1_num,x2_num,x3_num,x4_num,dx1,dx2,dx3,dx4);
  
  rewind(stdin);
  for (ih=0;ih<nh;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    tr.offset=(int) h[ih];
    tr.ntr=nh;
    if (Wd[ih][0]==0) tr.trid=2; // dead trace
    fputtr(stdout,&tr);
  }
  /******** End of interpolated output **********/


  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}










