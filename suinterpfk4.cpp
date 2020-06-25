/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADONHYPFK:  $Date: June 1999 - Last version October 2000  */

#include "interpfk.h"
#include "segy.h"
#include "header.h"
#include <time.h>


//#include "radonhypfk.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONHYPFK Forward  Hyperbolic Radon transform  in the f-k domain ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradonhypfk < stdin > stdout [optional parameters]          	",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Input : sudata file  (offset time domain)              		",
  " Output : adjoint model                                              ",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradon1 pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudata     ", 
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
//inv_par inv;

int verbose;

int main(int argc, char **argv)
{
  segy tr; 
  inv_par inv;
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep;
  //  FILE *offsetfile; 
  time_t start,finish;
  double elapsed_time;
  int it,ih, ik;
  float vel,dv;
  float t0=0;
  float **datain=0;
  float **dataout=0;
  float **datains=0;
  float **dataouts=0;
  float **dataout0;
  float *t, *h, *h2, *k;
  complex **F, **FLS, **F2, **FLS2;
  int nt, nh, nh2, nk; 
  int method;
  int plot;
  float dt,dk,dh2;
  float kmin;
  int testadj;
  int dft;
  float epsfft;
  float fmax;
  float ascale;
  int option;

  // Velocity law and stretching
  float *tmig;
  float *vmig;
  int ntmig;
  int nvmig;
  int itmig;
  float smig;
  float vscale,vstolt,vmin,vmax,ft,du,*v,*ut,*tu;
  int nu;
  
  // interpolation by zero padding or zero traces */
  float **dataout1;
  float interpfact;
  float dh;

  // data weights
  float **Wd=0;
  float **Wds=0;
  
  ////////////////
    
  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  //////////////////////////////////////////////
  fprintf(stderr,"*******SURADONHYPFK*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("vel", &vel))  vel = 1500;
  if (!getparfloat("dv", &dv))  dv = 100;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparfloat("epsfft", &epsfft))  epsfft = 1e-3;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("step", &inv.step))  inv.step =1;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 3;
  if (!getparint("norm", &inv.norm))  inv.norm =0; 
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparint("testadj",&testadj)) testadj=0; 
  if (!getparint("dft",&dft)) dft=1; 
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  if (!getparfloat("kmin",&kmin)) kmin=0;
  if (!getparfloat("interpfact",&interpfact)) interpfact=2;
  if (!getparint("option",&option)) option=3;
 
  // Velocity law
  if (!getparfloat("vscale",&vscale)) vscale = 1.0;
  if (!getparfloat("ascale",&ascale)) ascale = 1.0;
  ntmig=countparval("tmig");
  if (ntmig==0) ntmig=1;
  tmig=ealloc1float(ntmig);
  if (!getparfloat("tmig",tmig)) tmig[0]=0.0;

  nvmig=countparval("vmig");
  if (nvmig==0) nvmig=1;
  vmig=ealloc1float(nvmig);
  if (!getparfloat("vmig",vmig)) vmig[0]=2000.0;
  
  if (ntmig!=nvmig) err("number of tmig and vmig must be equal");
  for (itmig=1;itmig<ntmig;++itmig)
    if (tmig[itmig]<=tmig[itmig-1])
      err("tmig must increase monotonically");
  if (!getparfloat("smig",&smig)) smig=1.0;
  // Wavelet for the RT operator
  int nw;        // number of point for the wavelet
  float fpeak;   // peak frequency for the wavelet
  int typewav;   // type of wavelet
  float *wavelet=0;
  void (*convop)  (int, int, int, float *,int, float *,float *);

  if (!getparint("nw",&nw)) nw =0;
  if (!getparfloat("fpeak",&fpeak)) fpeak =25;
  if (!getparint("typewav",&typewav)) typewav = 1;
  convop=contruc_2;

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  nh2=10*nh;
  ft = tr.delrt/1000.0;
  if (ft!=0.0) err("cannot handle non-zero time of first sample");

  if (!getparfloat("fmax",&fmax)) fmax = 0.5/dt;
  fmax = MIN(fmax,0.5/dt);
  
  /* make uniformly sampled rms velocity function of time */
  makev(ntmig,tmig,vmig,vscale,nt,dt,ft,&v,&vmin,&vmax);
  
  /* Stolt migration velocity is the minimum velocity */
  vel = vstolt = vmin;

  /* make u(t) and t(u) for Stolt stretch */
  makeut(vstolt,fmax,v,nt,dt,&ut,&nu,&du,&tu);
  free1float(v);
  

  /******************************************************/
  // Generate a wavelet for the forward and adjoint operator and plot it
  if (nw){
    /*
    typewav=1==>   Read a wavelet generated by Matlab (stored in an ascii file)
    typewav=2==>   Generate a Ricker with SU
    typewav=3==>   Design a simple wavelet by hand
    */
    //if ( typewav==1 || typewav == 2 ) nw=50;
    wavelet=ealloc1float(nw);
    nw=get_wavelet(wavelet,"BP_wavelet",nw,typewav,du,fpeak);
    wavelet=erealloc1float(wavelet,nw);
    fprintf(stderr,"Test dot product for contran\n");
    //test=testadjop_conv(convop,nw,wavelet,nu,nu+nw-1);
  }      
  // Get info from first trace 
  /******************************************************/
  //fprintf(stderr,"nu=%d,du=%f,tu[1]-tu[0]=%f,ut[1]-ut[0]=%f\n",nu,du,tu[1]-tu[0],ut[1]-ut[0]);  
  fprintf(stderr,"nu=%d,du=%f,option=%d\n",nu,du,option);  
  //fprintf(stderr,"tu[1]=%f\n",tu[1]);  
  // Allocate memory for data and model

  datain=ealloc2float(nt,nh);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);
  /* If offetfile name is given read it */
  h2=ealloc1float(nh2); // allocate more to play safe
  memset( (void *) h2, (int) '\0', nh2 * FSIZE);


  if (option==3) if (offsetfile) nh2=read_ascii_file(offsetfile,h2); 
  else nh2=nh;

  dataout0=ealloc2float(nt,nh);
  // Loop over traces 
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
    fprintf(stderr,"ih=%d\n",ih);
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;

  /* Time axis */
  for (it=0;it<nt;it++) t[it]=t0+it*dt;

  if ((option==1)||(option==2)){
    dataout1=ealloc2float(nt,200);
    Wd=ealloc2float(nt,nh*(int) interpfact);
  }
  else if (option==3){
    dataout1=ealloc2float(nt,nh2);
    Wd=ealloc2float(nt,nh2);
  }
  else{
    dataout1=ealloc2float(nt,nh);
    Wd=ealloc2float(nt,nh);
    
  } 

  /* initial interpolation by zero padding in f-k */
  if (1){ 
    if (0){
      save_gather(datain,nh,nt,0.004,"datain.su");
      system("suxwigb < datain.su perc=100 title=datain &");  
    }
    if (option==1) nh2=fft2_zeropad(datain,dataout1,nt,nh,interpfact);
    else if (option==2){  nh2=add_zerotraces(datain,dataout1,Wd,nt,nh,interpfact);}
    else if (option==3){ nh2=add_zerotraces(datain,dataout1,Wd,nt,nh,nh2,h,h2);}
    if (1){
      save_gather(dataout1,nh2,nt,0.004,"dataout1.su");
      system("suxwigb < dataout1.su perc=100 title=dataout2 &");    
    }
    else for (ih=0;ih<nh2;ih++) h2[ih]=h[ih];
  }
  if ((option==1)||(option==2)){
    /* Create new offset axis */
    dh=(h[nh-1]-h[0])/(nh-1);
    dh2=dh/interpfact;
    for (ih=0;ih<nh2;ih++) h2[ih]=h[0]+ih*dh2;
  }
  else if (option==3){
    dh2=(h2[nh2-1]-h2[0])/(nh2-1);
  }
  else dh2=(h2[nh2-1]-h2[0])/(nh2-1);

  
  /* Change all offset axis (to keep unchanged rest of the code*/
  h=realloc1float(h,nh2);
  for (ih=0;ih<nh2;ih++) h[ih]=h2[ih];
  nh=nh2;
  
  /* Change original data */
  free2float(datain);
  datain=ealloc2float(nt,nh2);
  for (ih=0;ih<nh2;ih++) for (it=0;it<nt;it++) datain[ih][it]=dataout1[ih][it]; 

  /* Temporal array dataout1 is no longer needed */
  free2float(dataout1);

  datains=ealloc2float(nu,nh);
  dataouts=ealloc2float(nu,nh2);
  Wds=ealloc2float(nu,nh2);
  
  stretch(datain,datains,nt,nu,nh,t,tu,ut,dt,du,1);
  stretch(Wd,Wds,nt,nu,nh,t,tu,ut,dt,du,1);

  free2float(Wd);
  
  if (nu==nt){
    tu=ealloc1float(nt);
    ut=ealloc1float(nt);
    for (it=0;it<nt;it++) tu[it]=ut[it]=t[it];
  }
 
  kaxis(nh2,vel,du,nu,dh2,&nk,&dk);
  fprintf(stderr,"nk=%d,dk=%f,dh2=%f,nx=%d\n",nk,dk,dh2,nh2);

  
  k=ealloc1float(nk);
  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;

  /* We need the operator F that performs d=Fm
     and the least squares operator FLS=(F^T F)^-1 F^T
     that finds m=FLS d 
     For regular sampling they are equal F=FLS */
  
  F=ealloc2complex(nk,nh);
  F2=ealloc2complex(nk,nh2);
  FLS=ealloc2complex(nk,nh);
  FLS2=ealloc2complex(nk,nh2);
  
  FTmatrix(F,FLS,h,k,nh,nk,epsfft);
  FTmatrix(F2,FLS2,h2,k,nh2,nk,epsfft);
  
  if (0){
    save_gather(datain,nh,nt,0.004,"datain.su");
    system("suxwigb < datain.su perc=100 title=datain &");
  }	

  if (testadj) adjteststoltz(nu,nh,nh2,tu,h,h2,vel,F,F2,wavelet,nw);
  
  if (0){
    stoltzop2(datains,dataouts,nu,nh,nh2,tu,h,h2,vel,F,F2,wavelet,nw,1);
    save_gather(dataouts,nh2,nu,0.004,"dataout.su");
    system("suxwigb < dataout.su perc=100 title=dataout &");
  }	
   
  if (0){
    stoltzop2(datains,dataouts,nu,nh,nh2,tu,h,h2,vel,F,F2,wavelet,nw,0);
    save_gather(datains,nh,nu,0.004,"datain.su");
    system("suxwigb < datain.su perc=100 title=datain &");
  }	 

  if (0){
    save_gather(Wds,nh,nu,0.004,"Wd");
    system("suxwigb < Wd perc=100 title=Wd &");
  }

  fprintf(stderr,"nt=%d,nu=%d\n",nt,nu);
  stoltz_wtcgls(datains,dataouts,Wds,h,nh,tu,nu,h2,nh2,vel,inv,F,F2,wavelet,nw);
  // data weights no need any longer 
  free2float(Wds);

  if (0){
    save_gather(dataouts,nh2,nu,0.004,"migrated.su");
    system("suxwigb < migrated.su perc=100 title=migrated &");
  }	 


  /* Create output array */
  dataout=ealloc2float(nt,nh2);  
  
  stretch(dataout,dataouts,nt,nu,nh2,t,tu,ut,dt,du,-1);
  
  if ((modelfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);

  rewind(stdin);
  for (ih=0;ih<nh2;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    /* There is an error in the scale
       I apply this correction for now until finding the problem */
    for (it=0;it<nt;it++) tr.data[it]*=ascale;
    tr.offset=(int) h2[ih];
    tr.ntr=nh2;
    fputtr(modelfilep,&tr);
  }
  efclose(modelfilep);
  
  stoltzopinv2(dataouts,dataouts,nu,nh2,tu,h2,vel,F2,wavelet,nw);
  
  if (0){
    save_gather(dataouts,nh2,nu,0.004,"demigrated.su");
    system("suxwigb < demigrated.su perc=100 title=demigrated &");
  }	

  stretch(dataout,dataouts,nt,nu,nh2,t,tu,ut,dt,du,-1);
  
  rewind(stdin);
  for (ih=0;ih<nh2;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    /* There is an error in the scale
       I apply this correction for now until finding the problem */
    for (it=0;it<nt;it++) tr.data[it]*=ascale;
    tr.offset=(int) h2[ih];
    tr.ntr=nh2;
    fputtr(stdout,&tr);
  }
  
  if (nw) free1float(wavelet); 
  free2complex(FLS2);
  free2complex(FLS);
  free2complex(F2);
  free2complex(F);
  free1float(k);
  free1float(h2);
  free1float(t);
  free1float(h);
  free2float(datain);
  free2float(dataout);
  free2float(datains);
  free2float(dataouts);
  free2float(dataout0);

  free1float(tmig);
  free1float(vmig);
  if (nt==nu){
    free1float(tu);
    free1float(ut);
  }

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}


static void stretch(float **data, float **datas, int nt, int nu, int nh, float *t, float *tu, float *ut,
	     float dt, float du, int str)
{
  int ih, it, iu;

  if (str==1){/* stretch */
    if (nu!=nt){// && tu!=NULL) {
      float *temp=ealloc1float(nu);
      for (ih=0; ih<nh; ++ih) {
	ints8r(nt,dt,0.0,data[ih],0.0,0.0,nu,tu,temp);
	for (iu=0; iu<nu; ++iu) datas[ih][iu] = temp[iu];
      }
      free1float(temp);
    }
    else{
      for (ih=0; ih<nh; ++ih) 
	for (it=0; it<nt; ++it) datas[ih][it] = data[ih][it];
    }
  }
  else if (str==-1){   /* unstretch */
    if (nu!=nt){// && ut!=NULL) {
      float *temp=ealloc1float(nt);
      for (ih=0; ih<nh; ++ih) {
	ints8r(nu,du,0.0,datas[ih],0.0,0.0,nt,ut,temp);
	for (it=0; it<nt; ++it) data[ih][it] = temp[it];
      }
      free1float(temp);
    }
    else{
      for (ih=0; ih<nh; ++ih) 
	for (it=0; it<nt; ++it) data[ih][it] = datas[ih][it];
      
    }
  }
  return;
}


static void makev (int nmig, float *tmig, float *vmig, float vscale,
	    int nt, float dt, float ft, float **v, float *vmin, float *vmax)
/*****************************************************************************
make uniformly sampled rms velocity function v(t) for migration
******************************************************************************
Input:
nmig		number of tmig,vmig pairs
tmig		array[nmig] of times
vmig		array[nmig] of rms velocities
vscale		velocity scale factor
nt		number of time samples
dt		time sampling interval
ft		first time sample

Output:
v		array[nt] of rms velocities
vmin		minimum velocity
vmax		maximum velocity
******************************************************************************
Author:	 Dave Hale, Colorado School of Mines, 10/22/91
*****************************************************************************/
{
	int it;
	float t,*vel,velmin,velmax,(*vmigd)[4];
	
	vmigd = (float(*)[4])ealloc1float(nmig*4);
	cmonot(nmig,tmig,vmig,vmigd);
	vel = ealloc1float(nt);
	for (it=0,t=ft; it<nt; ++it,t+=dt)
		intcub(0,nmig,tmig,vmigd,1,&t,&vel[it]);
	for (it=0; it<nt; ++it)
		vel[it] *= vscale;
	for (it=1,velmin=velmax=vel[0]; it<nt; ++it) {
		velmin = MIN(velmin,vel[it]);
		velmax = MAX(velmax,vel[it]);
	}
	free1float((float*)vmigd);
	*v = vel;
	*vmin = velmin;
	*vmax = velmax;
}

static void makeut (float vstolt, float fmax, float *vrms,
	     int nt, float dt, float **ut, int *nu, float *du, float **tu)
/*****************************************************************************
Compute u(t) and t(u) for Stolt stretch
******************************************************************************
Input:
vstolt		Stolt migration velocity
fmax		maximum frequency
vrms		array[nt] of RMS velocities
nt		number of t samples
dt		t sampling interval (first t assumed to be zero)

Output
ut		array[nt] of u(t); NULL if constant velocity
nu		number of u samples
du		u sampling interval (first u assumed to be zero)
tu		array[nu] of t(u); NULL if constant velocity
*****************************************************************************/
{
	int it;
	
	/* check for constant velocity */
	for (it=1; it<nt; ++it)
		if (vrms[it]!=vrms[0]) break;
		
	/* if constant velocity */
	if (it==nt) {
		*ut = NULL;
		*tu = NULL;
		*nu = nt;
		*du = dt;

	/* else if velocity not constant */
	} else {
		int it,nuu;
		float duu,delu,umax,*u,*t;
		
		/* u(t) */
		u = alloc1float(nt);
		makeu(vstolt,vrms,nt,dt,u);

		/* smallest du and maximum u */
		duu = FLT_MAX;
		for (it=1; it<nt; ++it) {
			delu =	u[it]-u[it-1];
			if (delu<duu) duu = delu;
		}
		umax = u[nt-1];

		/* u sampling */
		duu = duu/(2.0*fmax*dt);
		nuu = 1+NINT(umax/duu);

		/* t(u) */
		t = alloc1float(nuu);
		yxtoxy(nt,dt,0.0,u,nuu,duu,0.0,0.0,(nt-1)*dt,t);

		/* set output parameters before returning */
		*ut = u;
		*tu = t;
		*nu = nuu;
		*du = duu;
	}
}

static void makeu (float vstolt, float *v, int nt, float dt, float *u)
/*****************************************************************************
Compute				     t
	u(t) = sqrt( 2/vstolt^2 * Integral ds*s*(v(s)^2) )
				     0
via the trapezoidal rule.
******************************************************************************
Input:
vstolt		Stolt migration velocity
v		array[nt] of RMS velocities
nt		number of t samples
dt		t sampling interval

Output
u		array[nt] of u(t)
*****************************************************************************/
{
	int it;
	float t,scale,sum;

	scale = 2.0/(vstolt*vstolt);
	u[0] = sum = 0.0;
	for (it=1,t=dt; it<nt; ++it,t+=dt) {
		sum += 0.5*dt*(t*v[it]*v[it]+(t-dt)*v[it-1]*v[it-1]);
		u[it] = sqrt(scale*sum);
	}
}

















