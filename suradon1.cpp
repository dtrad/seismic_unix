/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADTD:  $Date: June 1999  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrarytd.h"
#include <time.h>
#include "inversion_par.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADINT Forward  High Resolution Hyperbolic Radon transform        ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradon1< stdin > stdout [optional parameters]          		",
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
  " dq=2e-8             defines the rate of constant slowness spacing   ",
  " nq=20		Number of Radon Traces                         	",
  " factor=0.8		multiply dq critical by factor                	",
  " precond=0           Uses precondtioning                             ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " parnmute=1      if q > parmute the radon trace is tapered           ",
  " resample=0      =0 keeps original offset                            ",
  "                 =1 defines a new offset axis with similar aperture  ",
  "                 that original one but regular and resamples the data",
  "                 =2 search for gaps and fills them keeping original  ",
  "                  traces.                                            ",
  "                 =3 resamples to the offset given by offsetfile      ",
  "                 =0 keeps original offset                            ", 
  " keepdata=0      when =1 if there are traces with equal offset than  ",
  "                 the interpolated ones, keep the original traces.    ",
  " flagvel=1       to use velocities for vgrid                         ",
  "        =0       to use square slowness like conventional Hyp RT     ",
  " centralq=0      the q trace where the velocity law lays             ",
  " filtout =0      event with energy higher than quantil 99 are        ",
  "                 downweighted by Q50/Q99                             ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Output : Iterpolated traces  (offset time domain)                   ",
  " Input : sudata file  (offset time domain)              		",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradon1 pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudataint  ", 
  "                                                                     ",
  " If nh output is > nh input only the header words of the first nh    ",
  "  traces are preserved                                               ",
  NULL};

/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
//inv_par inv;
inv_par inv;
segy tr,tri; 
options_par opt;

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
  int j,i,iq,ih, ih2;
  int it;
  float **d, **m;
  float *q, *t, *h, eps, qmin, qmax, fmax;
  float t0;
  int smooth;
  int model; 
  int nt, nh, nq,rtmethod; 
  int method, nx ,ny;
  float dt,dh,dq,thres;
  float factor; // multiply dq by factor
  // Interpolated data
  float **dint;
  int nh2;
  float dh2;
  float dhmin;
  float h2min;
  float h2max;
  float hmin;
  float hmax;
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
  int resample; // =1 produces data with a new offet axis given by h2max h2min
  int keepdata;
  int flagvel; /* Set on this flag to use velocities for vgrid 
		    or =0 to use square slowness like conventional Hyperb RT */
  int centralq; /* q trace where we put the velocity law */
  int filtout;  /*  filtout =0      event with energy higher than quantil 99 are  
                    downweighted by Q50/Q99  */                         
  // smoothing
  int nl=3;  //  npoints left hand side
  int nr=3;  //  npoints left hand side
  int flag=2;  // 1 rectangular, 2 triangular
  ////////////////////////
  inv.restart=1;

  //////////////////////////////////////////////
  fprintf(stderr,"**************************\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("dq",&dq)) dq =2e-8;
  if (!getparfloat("step", &inv.step))  inv.step =0.9;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparint("norm", &inv.norm))  inv.norm =1; 
  if (!getparint("testadj", &opt.testadj))  opt.testadj =0;
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparint("taperflag", &opt.taperflag))  opt.taperflag =0;
  if (!getparint("model",&model)) model = 1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("resample",&resample)) resample=0;
  if (!getparfloat("dhmin",&dhmin)) dhmin=20;
  if (!getparfloat("parmute",&opt.parmute)) opt.parmute=1;
  if (!getparint("mute",&opt.mute)) opt.mute=0;
  if (!getparint("keepdata",&keepdata)) keepdata=1;
  if (!getparint("smooth",&smooth)) smooth=0;
  if (!getparint("flagvel",&flagvel)) flagvel=0;
  if (!getparint("centralq",&centralq)) centralq=0;
   if (!getparint("filtout",&filtout)) filtout =0;
  if (!flagvel) dperv=dq;

  /* Introduce velocity trend to apply Hyp Vel Filtering */
  ntvel = countparval("tvel");
  if (ntvel==0) ntvel = 1;
  tvel = ealloc1float(ntvel);
  if (!getparfloat("tvel",tvel)) tvel[0] = 0.0;
  nvel = countparval("vel");
  if (nvel==0) nvel = 1;
  if (nvel!=ntvel) err("number of tvel and vel must be equal");
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
  fprintf(stderr,"nh=%d,nt=%d,dt=%f,eps2=%e\n",nh,nt,dt,inv.eps2); 

  //if ((nt*nq) > USHRT_MAX) err("nt*nq > maximum unsigned short\n");
  //if ((nt*nh) > USHRT_MAX) err("nt*nh > maximum unsigned short\n");

  // Allocate memory for data and model
  
  d=ealloc2float(nt,nh);
  m=ealloc2float(nt,nq);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  t=ealloc1float(nt);
  velint=ealloc1float(nt);
      	
  memset( (void *) h, (int) '\0', nh * FSIZE);

  headerfp = etmpfile();
  tracefp = etmpfile();
  if (verbose) warn("using tmpfile() call");


  // Loop over traces 
  ih=0;
  hmin=1e10;
  hmax=1e-10;
  do {
    register int i;
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);      
    h[ih]=(float) tr.offset;
    if (h[ih]<hmin) hmin=h[ih];
    if (h[ih]>hmax) hmax=h[ih]; 
    for (i=0;i<nt;i++){
      d[ih][i]=(float) tr.data[i];    //if sort by time-offset 
    }
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(headerfp);
  erewind(tracefp);
  nh=ih;
  
  // Define new offset axis
  // The default is to keep same number of traces 
  
  if (resample==1){
    if (!getparfloat("h2min",&h2min)) h2min=hmin;
    if (!getparfloat("h2max",&h2max)) h2max=hmax; 
    if (!getparfloat("dh2",&dh2)) dh2=fabs(h2max-h2min)/(nh-1);
    nh2=(int) ((h2max-h2min)/dh2 + 1 );
    if ((h2=alloc1float(nh2))==NULL) err("Cannot allocate h2\n");
    for (h2[0]=h2min,ih=1;ih<nh2;ih++) h2[ih]=h2[ih-1]+dh2;
  }
  ////////////////////////////////
  else if(resample==2){
    h2=ealloc1float(5*nh);
    nh2=findgaps(h,h2,nh,dhmin,hmin,hmax);
    for (ih=0;ih<nh;ih++) fprintf(stderr,"h2[%d]=%f,h[%d]=%f\n",ih,h2[ih],ih,h[ih]);
    fprintf(stderr,"+nh2=%d\n",nh2);
  }
  else if(resample==3){
    FILE *offsetfilep;
    cwp_String offsetfile=""; /* file containing positions */
    int nn;
    if (getparstring("offsetfile",&offsetfile)){
      h2=ealloc1float(nh*2);
      fprintf(stderr,"New offset given by %s\n",offsetfile);  
      offsetfilep=efopen(offsetfile,"r");
      ih=0;
      do{
	nn=fscanf(offsetfilep,"%f",&h2[ih]); 
	ih++;
      }while(nn==1);
      nh2=ih-1;
      fprintf(stderr,"nh2= %d\n",nh2);
      efclose(offsetfilep);
    }
    else err("No offset file\n");
  }
  else{
  // Keep axis for residuals study
    h2min=hmin;
    h2max=hmax; 
    nh2=nh;
    if ((h2=alloc1float(nh2))==NULL) err("Cannot allocate h2\n");
    for (ih=0;ih<nh2;ih++) h2[ih]=h[ih];
  }

  fprintf(stderr,"nh2=%d,nt=%d\n",nh2,nt);
  if ((dint=alloc2float(nt,nh2))==NULL) err("Cannot allocate dint\n");
  /* Time axis */
  for (i=0;i<nt;i++) t[i]=0+i*dt;

  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);


  if (verbose) fprintf(stderr,"Before radtd nq=%d, nt=%d, nh=%d, eps=%f \n",nq,nt,nh,eps);  

  ny=nt*nh;
  nx=nt*nq;

  radint_sparse4(t,q,h,h2,m,d,dint,nt,nh,nq,nh2,dt,velint,dperv,pervmin,t0,opt,inv,
		 flagvel,centralq,filtout);
 
  fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  if (model){
    if (smooth) smoothing(m[0],nt,nq,nl,nr,flag);
    if ((myfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);
    iq=0;
    do{
      //tr.offset=(int) q[iq];  // copy velocity perturbation
      tr.f2=q[iq];
      tr.ntr=nq;    
      for (i=0;i<nt;i++)
	tr.data[i]=m[iq][i];    
      fputtr(myfilep,&tr);
      iq++;
    } while(iq<nq);
  }
  
  if (smooth) smoothing(dint[0],nt,nh2,nl,nr,flag);
  ih=0;


  int *hflag=ealloc1int(nh2);
  memset( (void *) hflag, (int) '\0', nh2 * FSIZE);

  if ((!opt.mute)&&(keepdata)) sethflag(hflag,nh2,h,nh,h2,nh2);

  ih=0;
  for(ih2=0;ih2<nh2;ih2++){
    if ((hflag[ih2]==1)&&(ih<nh)){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      ih++;
    }
    else{
      tr.offset=(int) h2[ih2];
      tr.f2=0;
      for (i=0;i<nt;i++)  tr.data[i]=dint[ih2][i];
    }
    //tr.offset=(int) h2[ih2];
    tr.ntr=nh2;
    puttr(&tr);
  }
  free1int(hflag);
  
  if (verbose) fprintf(stderr,"nt%d, nh2 %d \n ",nt,nh2);

  fprintf(stderr," nq=%d, nt=%d, nh=%d, iq=%d\n",nq,nt,nh,iq);

  free2float(dint);
  free1float(h2);
  free1float(velint);  
  free1float(t);
  free1float(h);
  free1float(q);
  free2float(m);
  free2float(d);
  free1float(vel);
  free1float(tvel);
  efclose(headerfp);
  efclose(tracefp);
  efclose(myfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}
