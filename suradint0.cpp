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

int main(int argc, char **argv)
{
	
  //FILE *myfilep;
  int j,i,iq,ih;
  int it;
  float *d, *m;
  float *q, *t, *h, eps, qmin, qmax, fmax;
  float t0=0;
  int smooth;
  int model; 
  extern int nt, nh, nq,rtmethod, norm, itercg, iter_end, method, nx ,ny;
  extern int reorth;
  extern float dt,dh,dq, eps1,eps2, step, thres, theta;
  extern float factor; // multiply dq by factor
  extern int testadj;
  extern int taper;
  // Interpolated data
  float *dint;
  int nh2;
  float dh2;
  float h2min;
  float h2max;
  float *h2;
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
  if (!getparint("model",&model)) model = 1;    
  if (!getparint("verbose", &verbose))  verbose =0;
 
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
  
  if ((dint=ealloc1float(nt*nh))==NULL)
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
    if (fabs(h[ih])<fabs(hmin)) hmin=h[ih];
    if (fabs(h[ih])>fabs(hmax)) hmax=h[ih]; 
    for (i=0;i<nt;i++){
      d[ih*nt+i]=(float) tr.data[i];    //if sort by time-offset 
    }
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(headerfp);
  nh=ih;

  // Define new offset axis
  // The default is to keep same number of traces with
  if (!getparint("dh2",&dh2)) dh2=(h2max-h2min)/(nh-1);
  if (!getparfloat("h2min",&h2min)) h2min=hmin;
  if (!getparfloat("h2max",&h2max)) h2max=hmax; 
  nh2=(h2max-h2min)/dh2 + 1;
  if ((h2=alloc1float(nh2))==NULL) err("Cannot allocate h2\n");
  for (h2[0]=h2min,ih=1;ih<nh2;ih++) h2[ih]=h2[ih-1]+dh2;
  ////////////////////////////////
  if ((dint=alloc1float(nt*nh2))==NULL) err("Cannot allocate dint\n");
  /* Time axis */
  for (i=0;i<nt;i++) t[i]=t0+i*dt;

  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);

  //if (verbose) for (it=0;it<nt;it++) fprintf(stderr,"velint[%d]=%f,t[%d]=%f\n",it,velint[it],it,t[it]);

  if (verbose) fprintf(stderr,"Before radtd nq=%d, nt=%d, nh=%d, eps=%f \n",nq,nt,nh,eps);  
  ny=nt*nh;
  nx=nt*nq;

  radint(t,q,h,h2,m,d,dint,eps,qmin,qmax,fmax,velint,dperv,pervmin);
  
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  if (model){
    if (smooth) smoothing(m,nt,nq,nl,nr,flag);
    FILE *myfilep;
    if ((myfilep=fopen("model.su","w"))==NULL)
      err("cannot open file=%s\n","model.su");
    iq=0;
    do{
      tr.f2=q[iq];  // copy radon parameter in tr.f2
      //if (method==8) tr.f2=vgrid[iq][0];
      tr.ntr=nq;
    
      for (i=0;i<nt;i++)
	tr.data[i]=m[iq*nt+i];    
      fputtr(&tr);
      iq++;
    } while(iq<nq);
    fclose(myfilep);
  }
  
  if (smooth) smoothing(d,nt,nh,nl,nr,flags);
  ih=0;
  do{
    if (ih<nh)  efread(&tr,HDRBYTES,1,headerfp);
    tr.offset=(int) h2[ih];
    tr.ntr=nh2;
    tr.f2=0;
    for (i=0;i<nt;i++)
      tr.data[i]=dint[i+ih*nt];
    puttr(&tr);
    ih++;
  } while(ih<nh2);
  
  if (verbose) fprintf(stderr,"nt%d, nh2 %d \n ",nt,nh2);

  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free1dint(dint);
  free1float(h2);
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


void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin)

  /*
    subroutine: radtd.cpp
    RADON TRANSFORM IN TIME DOMAIN
    Three possible curves, linear, parabolic and hyperbolic
    are implemented in the file radon.cpp
    Inside this file the functions are called radonlin, radonpar, radonhyp
    At the beginning a pointer to function is declared and 
    associated with the chosen function. 
    This pointer is used as a reference to the operator 
    in WTCGLS, TESTADJ or any other function. 

    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
  */
{
  extern int nt,nh,nq, nx, ny, method, iter_end, reorth, rtmethod;
  extern float dt,dh,dq,eps1,eps2,thres,theta;
  register int i;
  int  j, iter, ih, iq;
  float *Wm;
  float *Wd;
  float dx_max;
  float dx_av;
  float qmaxt;
  float test;
  extern float factor; // Used to compute dq, i.e., dq=factor*dq;
  extern float step;
  extern int itercg;
  extern int testadj;
  extern int taper;
  //////////////////////////////////////////////////////////////////////////
  void (*radon) (float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon2) (float *,float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon3) (float *,float *,float *,float *,float *,float **,int ,int ,int,int);
  
  if (rtmethod==3) {
    radon=radonhyp;
    radon2=radonhyp;
    radon3=radonhyp;
  }
  else if (rtmethod==2) radon=radonpar;
  else if (rtmethod==1) radon=radonlin;
  
  /////////////////////////////////////////////////////////////////////////
 
  if ((Wm=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n"); 
  if ((Wd=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n"); 
  ///////////////////////////////////////////////////////////////////
  // Compute inv(CD)=WdT Wd       
  for (i=0;i<ny;i++) Wd[i]=1.;
  if (taper==1)
    for (ih=0;ih<5;ih++) 
      for (i=0;i<nt;i++){
	Wd[ih*nt+i]=1-exp(-(ih+.3));
	Wd[(nh-ih-1)*nt+i]=Wd[ih*nt+i];
      }
  
  //   Velocity axis //////////////////////////////////////////////////
  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f, eps=%f\n", dx_max, dx_av, eps);

  if (fmax==0) fmax=0.6/(2*dt); // sinc8 is valid only until 0.6 Nyquist 
  if ((rtmethod==1)||(rtmethod==2)) 
    radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  
  else dq = (qmax-qmin)/(nq-1);

  for (i=0;i<nq;i++) q[i] = qmin+i*dq;
  float **vgrid;
  float *perv;

  if (method==8){ // Irregular velocity grid
    // Define irregular spacing along the horizontal
    // Here q  is perturbartion around the central velocity law
    // dperq is a parameter used to generate the increasing space
    // perqmin defines the minimum distance to the perturbation
    int nqh=nq/2;
   
    
    unsigned int it;
    if ((perv=alloc1float(nq))==NULL) err(" Cannot allocate perv");
    if ((vgrid=alloc2float(nt,nq))==NULL) err("Cannot allocate vgrid");
    perv[nqh]=0;
    for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
    for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
    for (iq=0;iq<nq;iq++) fprintf(stderr,"perv[%d]=%f,pervmin=%f,dperv=%f\n",iq,perv[iq],pervmin,dperv);
    // Define irregualr qgrid as a sumation of 
    for (iq=0;iq<nq;iq++){
      for (it=0;it<nt;it++){
	vgrid[iq][it]=vel[it]*vel[it]+perv[iq]*perv[iq]+2*vel[it]*perv[iq];
      }
    } 
  }

  ///////////////////////////////////////////////////////////////////////  

  fprintf(stderr,"Inside radtd nq=%d, nt=%d, nh=%d dq=%f\n",nq,nt,nh,dq);
  
  for (i=0;i<nx;i++) m[i]=0.;
 
  if (method==0){ // Adjoint and tests
    radon(m,t,h,q,d,1,nt,nh,nq); 
    if (testadj) for (i=0;i<iter_end;i++) test=testadjop(radon,t,h,q,nt,nh,nq);
  }  

  /* Marfurt's method with Semblance analysis */
  else if (method==1) {
    //Semblance computation
    FILE *myfilep;
    if((myfilep=fopen("semblance","w"))==NULL)
      err("cannot open file=%s\n","semblance");    
    semblance(Wm,t,h,q,d,nt,nh,nq,dt);
    fwrite(Wm,sizeof(float),nx,myfilep);
    fclose(myfilep);
    
    // Mask computation
    if ((thres>=0.7)||(thres<=0.01)) thres=0.3;
    FILE *myfilep2;
    if((myfilep2=fopen("mask","w"))==NULL)
      err("cannot open file=%s\n","mask");
    /*Mask is given by losigm function*/
    for (i=0;i<nx;i++) Wm[i]=1+1./(1+exp(100*(Wm[i]-thres)+0.5));
    fwrite(Wm,sizeof(float),nx,myfilep2);
    fclose(myfilep2);

    for (iter=1;iter<=iter_end;iter++){
      fprintf(stderr,"iter_end=%d;iter=%d;\n",iter_end,iter);
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);
      if (iter<iter_end) {
        //thres=0.04; 
	for (i=0;i<nx;i++) 
	  Wm[i]=1+1./(1+exp(100*(m[i]-thres)+0.5));
      }
    }
  }
  /* Sparse approach (Thompson, 1984) 
     First it is computed the adjoint model, and  a Covariance matrix
     is computed from it. WTCGLS solves the system 
     (LH*L + Cm) m = LH d
  */
  else if (method==2) {
    radon(m,t,h,q,d,1,nt,nh,nq);
    if (testadj) test=testadjop(radon,t,h,q,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);  
      }
  }
  /* LSQR does not work properly */
  else if (method==3) {
    lsqr(t,q,h,m,d,eps,reorth);
  }
  /* Method=4 is used for test */ 
  else if (method==4)  radonopi_id(m,t,h,q,d); 
  else if (method==5) {

    // Method 5 is similar to Thorson but the rho filter is applied as
    // a preconditioner to the adjoint model before start WTCGLS
    int lfilter=51;
    int np2=lfilter/2;	
    float *ptrace, *ptrace2;
    float *rho;
    if ((ptrace=alloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for ptrace could not be allocated\n");
    if ((ptrace2=alloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for ptrace could not be allocated\n");   
    if ((rho=ealloc1float(lfilter))==NULL)
      fprintf(stderr,"***Sorry, space for rho could not be allocated\n"); 
    rho_filter(lfilter,nt,dt,rho);
    //for (i=0;i<lfilter;i++) fprintf(stderr,"rho[%d]=%f\n",i,rho[i]);
    fprintf(stderr,"RHO filter......................\n");

    radon(m,t,h,q,d,1,nt,nh,nq);
    
    for (int iq=0;iq<nq;iq++){ 
      for (i=0;i<nt;i++) ptrace[i]=m[iq*nt+i];
      conv(nt,-np2,ptrace,lfilter,0,rho,nt,0,ptrace2);
      for (i=0;i<nt;i++) m[iq*nt+i]=ptrace2[i];
    }    
  }
  else if (method==6) {
    radon(m,t,h,q,d,1,nt,nh,nq);
    hilbert(nt*nq,m,Wm);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+Wm[i])+eps1;
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);  
    }
  }
  else if (method==7) {
    radon2(m,t,h,q,d,vel,1,nt,nh,nq);
    if (testadj) test=testadjop(radon2,t,h,q,vel,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon2,nt,nh,nq,t,h,q,m,d,Wm,Wd,vel,eps,step,eps1,eps2,itercg);  
    }
  }
  else if (method==8) {
    radon3(m,t,h,q,d,vgrid,1,nt,nh,nq);
    // if (testadj) test=testadjop(radon2,t,h,q,vel,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon3,nt,nh,nq,t,h,q,m,d,Wm,Wd,vgrid,eps,step,eps1,eps2,itercg); 
    }
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq);
  }
  if (method==8){
    free2float(vgrid);  
    free1float(perv);
  }
  free1float(Wm);
  free1float(Wd);  
  return;
  
}


















