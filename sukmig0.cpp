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
  " SUKMIG0 - First attempt for Kirchhoff migration                     ", 
  "	                                                          	",
  " 	   								",
  " sukmig0 < stdin > stdout [optional parameters]          		",
  " 									",
  " 									",
  "                                      				",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/

segy tr; // for malloc use *trr
void kmig0(float **d, float **m, float *t, float *h, float *x, float vel, float cdp) ;

int nt,nh,nx; // Size of volume x,h,t
float dt,dh,dx; // Sampling
int main(int argc, char **argv)
{
  int j,i, k, nhmax;
  register int it;
  float **d, **m;
  float *x, *t, *h, cdp, vel, cdpmin, cdpmax; 
  extern int nt, nh, nx;
  extern float dt,dh,dq;	
  
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
    
  // Get parameters 
  if (!getparfloat("vel", &vel))  vel = 2000;
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparint("nx", &nx))  nx = 100;
  if (!getparint("nhmax", &nhmax))  nhmax = 100;  
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
 
  fprintf(stderr,"nhmax=%d,nt=%d,dt=%f,nx=%d,vel=%f\n",nhmax,nt,dt,nx,vel); 

  // Allocate memory for data and model
  
  if ((d=ealloc2float(nhmax,nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
 
  if ((m=ealloc2float(nx,nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");
 
  if ((h=ealloc1float(nhmax))==NULL)
    fprintf(stderr,"***Sorry, space for h could not be allocated\n");
 
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 

   
  for (k=0;k<nx;k++) for (i=0;i<nt;i++) m[i][k]=0;


  dx=(cdpmax-cdpmin)/(nx-1);      
  for(i=0;i<nx;i++) x[i]=cdpmin+dx*i;
  for(it=0;it<nt;it++) t[it]=it*dt;
  fprintf(stderr,"**************************\n");  
  k=0;  
  do {   // Loop over midpoint  
    j=0;
    fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    do {  // Loop over traces  
      register int i;    
      h[j]=(float) tr.offset;h[j]/=2;  // halfoffset
      cdp=(float) tr.cdp;
      //fprintf(stderr,"h[%d]=%f,cdp=%f,k=%d\n",j,h[j],cdp,k);
      for (i=0;i<nt;i++)
	d[i][j]=(float) tr.data[i];
	      
      j++;
      gettr(&tr);
      
    } while (((float) tr.cdp==cdp)&&((nh=(j+1))<nhmax));
    kmig0(d,m,t,h,x,vel,cdp);
    //for (int ik=0;ik<nx;ik++) for (i=0;i<nt;i++) 
    //if (m[i][ik] > 1e-8 ) fprintf(stderr,"m[%d][%d]=%f\n",i,ik,m[i][ik]);
    k++;    
  } while (tr.cdp<cdpmax);
  for (j=0,k=0;k<nx;k++){
    j++;
    for (i=0;i<nt;i++)
      tr.data[i]=m[i][k];
    tr.cdp=(int) x[k];
    tr.dt=(int) (dt*1e6);
    tr.ntr=nx;
    tr.ns=nt;
    tr.tracl=j;
    tr.tracr=j;
    tr.offset=0;
    tr.sx=(int) x[k];
    tr.gx=(int) x[k];
    puttr(&tr);
  }   
  free1float(t);
  free1float(h);
  free1float(x);
  free2float(m);
  free2float(d);
  
  return EXIT_SUCCESS;
}




















