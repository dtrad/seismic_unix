/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUSYNTH:  $Date: December 1999  */

#include "su.h"
#include "segy.h"
#include "math.h"

void kolmogoroff(int nt, float *x); 
void ricker1_wavelet (int nt, float dt, float fpeak, float *wavelet);
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUMINPHASE  > output                                                ",
  " Required parameters  ... at least one                               ",
  " Optional parameters                                                 ",
  "  nh                                                                 ", 
  "  nx                                                                 ",
  "  dh                                                                 ",
  "  dx                                                                 ",
  "  nt                                                                 ",
  "  dt                                                                 ",
  "  wave                                                               ",
  "  vel                                                                ",
  "  tau                                                                ",
  "  coeff                                                              ",
  "  avo                                 				",
  "  avof                                                               ", 
  "  avop                                                               ",
  "  dh                                                                 ",
  "  dx                                                                 ",
  "  nt                                                                 ",
  "  dt                                                                 ",
  "  wave=0      0 Ricker 1 MP Ricker                                   ",
  "  nw                                                                 ",
  "  vel                                                                ",
  "  tau                                                                ",
  "  coeff                                                              ",
  "  avo                                 				",
  "  shape=3     1 Line 2 Parabola 3 Hyperbola 				",
  "  EXAMPLE                                                            ",
  "     susynth nh=30 nx=3 dh=100 dx=250 fpeak=20 nt=1024               ",
  "     midpoint=0 vel=1500,2000,3000  tau=0.5,1,2 coef=1,-0.5,1        ",
  "     wave=0 avof=1.1 avop=1 avo=2,1,0 | suxwigb  &                   ",
  "                                                                     ",     
  NULL};
/* Credits:
 *	Daniel Trad.
/**************** end self doc ***********************************/


int main(int argc, char **argv)
{
  int ih, it, iv, ix;
  int itint;
  float itfloat,a1,a2;
  int nt;
  int ntr;
  float *wavelet;
  float fpeak;
  float dt;
  float *data;
  float offset;
  float dh;
  float dx;
  float tau0;
  int itime;
  float hnear;
  float midpoint;
  float time;
  float veloc;
  float hoff;
  int nvel;
  float *vel;
  int ntau;
  float *tau;
  int ncoef;
  float *coef;
  int wave;
  int nx;
  int nh;
  float avof; // frequency of AVO 
  float avop; // AVO phase
  int *avo;
  int navo;
  int iavo;
  float amplitude;
  float avoscale;
  int nw;
  int shape; // 1 Line 2 Parabola 3  Hyperbola
  
  segy tr;

  // Initialize 

  initargs(argc, argv);
  if (argc==1) requestdoc(1);
  fprintf(stderr,"arc=%d\n",argc);

  if (!getparfloat("fpeak", &fpeak))  fpeak = 25;
  if (!getparint("nt", &nt))  nt = 1024;
  if (!getparint("nh", &nh))  nh = 20;
  if (!getparint("nx", &nx))  nx = 1;
  if (!getparfloat("dt", &dt))  dt = 0.004;  
  if (!getparfloat("dh", &dh))  dh = 25;
  if (!getparfloat("dx", &dx))  dx = 250;
 
  if (!getparfloat("hnear", &hnear))  hnear = (nh*dh/2.);  
  if (!getparfloat("midpoint", &midpoint))  midpoint = 0; 
  if (!getparint("wave", &wave))  wave = 0;
  if (!getparint("nw", &nw))  nw = 100;    
  if (!getparfloat("avof", &avof))  avof = 0;
  if (!getparfloat("avop", &avop))  avop = 0;
  if (!getparint("shape", &shape))  shape = 3;    
  ntr=nx*nh;
  //////////////////////////////////////////////////////////////////
  nvel = countnparval(1,"vel");if (nvel==0) nvel=1;
  if ((vel=ealloc1float(nvel))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");
  if (!getnparfloat(1,"vel",vel)) vel[0]=1500;
  /////////////////////////////////////////////////////////////////        
  ntau = countnparval(1,"tau");if (ntau==0) ntau=1;
  if ((tau=ealloc1float(ntau))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");
  if (!getnparfloat(1,"tau",tau)) tau[0]=0.5;        
  if (ntau!=nvel) err("nvel and ntau must be equal\n");
  //////////////////////////////////////////////////////////////////
  ncoef = countnparval(1,"coef");
  if (ncoef!=nvel) warn("***, ncoef differs from nvel. used coeff=1 instead");
  if ((coef=ealloc1float(nvel))==NULL)
    fprintf(stderr,"***Space for coef could not be allocated\n");
  if (!getnparfloat(1,"coef",coef)|| ncoef!=nvel) 
    for (iv=0;iv<nvel;iv++) coef[iv]=1; 
  //////////////////////////////////////////////////////////////////
  navo = countnparval(1,"avo");
  if ((avo=ealloc1int(navo))==NULL)
    fprintf(stderr,"***Space for avo could not be allocated\n");
  if (!getnparint(1,"avo",avo)) 
    for (iv=0;iv<navo;iv++) avo[iv]=0; 
  //////////////////////////////////////////////////////////////////
     
  for (it=0;it<nvel;it++) fprintf(stderr,"vel[%d]=%f,tau[%d]=%f\n",it,vel[it],it,tau[it]);

  for (iavo=0;iavo<navo;iavo++) fprintf(stderr,"avo[%d]=%d\n",iavo,avo[iavo]);

  if ((wavelet=alloc1float(nw))==NULL)
    err("cannot allocate memory for wavelet");

  if ((data=alloc1float(nt))==NULL)
    err("cannot allocate memory for data");

  if (wave==0 || wave==1) ricker1_wavelet (nw,dt,fpeak,wavelet);
  if (wave==1) kolmogoroff(nw,wavelet);

  fprintf(stderr,"nt=%d,dt=%f;nw=%d,ntr=%d,midpoint=%f\n",nt,dt,nw,ntr,midpoint);
  // Allocate memory for data and model
  fprintf(stderr,"**************************\n");    
  
  for (ix=0;ix<nx;ix++){    // loop over midpoints
  
    for (ih=0;ih<nh;ih++) {   // Loop over halfoffsets
      offset=ih*dh-hnear;hoff=offset/2.;
      
      for (it=0;it<nt;it++) data[it]=0;
      
      for (iv=0;iv<nvel;iv++){

	veloc=vel[iv];
	tau0=tau[iv];
	if (shape==3)
	  time=sqrt(tau0*tau0/4+((midpoint+hoff)*(midpoint+hoff))
		    /(veloc*veloc))+sqrt(tau0*tau0/4+
		    ((midpoint-hoff)*(midpoint-hoff))/(veloc*veloc));
	else if (shape==2)
	  time=(tau0/2+((midpoint+hoff)*(midpoint+hoff))
		    /(veloc*veloc))+(tau0/2+
		    ((midpoint-hoff)*(midpoint-hoff))/(veloc*veloc));	
	else if (shape==1)
	  time=tau0+(fabs(2*hoff)/veloc);
	
      
	itfloat=time/dt;
	itint=(int) floor(itfloat);
	a1=1-(itfloat-itint);
	a2=itfloat-itint;
	amplitude=coef[iv];
        for (iavo=0;iavo<navo;iavo++) 
	  if (avo[iavo]==iv) {
	    avoscale=cos((avof/1000.)*offset-(avof/10.));
	    amplitude=coef[iv]*avoscale;
            //fprintf(stderr,"avoscale=%f\n",avoscale);
	  }
	if (itint<nt) data[itint]=data[itint]+a1*amplitude;
	if ((itint+1)<nt) data[itint+1]=data[itint+1]+a2*amplitude;
      }

      conv(nt,0,data,nw,0,wavelet,nt,0,tr.data);
      
      tr.cdp=(int) midpoint;    
      tr.offset=(int) offset;
      tr.ns=(unsigned short) nt;
      tr.dt=(unsigned short) (dt*1e6);
      tr.ntr=(int) ntr;
      puttr(&tr);
    }
    midpoint+=dx;
  }
  free1float(data); 
  free1float(wavelet);
  free1int(avo);
  free1float(coef);
  free1float(tau);
  free1float(vel);  
  
  return EXIT_SUCCESS;
}




















