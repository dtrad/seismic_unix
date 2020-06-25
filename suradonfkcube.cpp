/* Copyright (c) University of British Columbia, 2002.*/
/* All rights reserved.                       */
/* suradonfk0  :  $Date: April    2002- Last version July    2002  */

#include "radonfkcube.h"
#include "segy.h"
#include "header.h"
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONFK0CUBE  RT for  apex shifted hyperbolas  with LS square     ",
  "                 time migration Stolt fk operator and multiple vel   ",
  " 	   								",
  " suradonfk0  < stdin > stdout [optional parameters]          	",
  " 									",
  " Prediction  of apex shifted hyoerbolas by using least squares       ",
  " Stolt migration.                                                    ",
  " The input file is a cdp file. It is migrated such that hyperbolas   ",
  " collapse to small areas. By using a pass mute window, undesired     ",
  " hyperbolic coherent noise can be predicted and subtracted from the  ",
  " the data.                                                           ",
  "                                                                     ",
  " Since the code uses 2d FFT then regular sampling is required.       ",
  " If the original gather has gaps then the program sufill             ",
  " can be used to add zero traces into these gaps.                     ",
  " If the data are irregularly sampled then you need to use the        ",
  " program suinterpfk0, that used DFT instead.                         ",
  "                                                                     ",
  " Note that this code is very similar to suinterpfk0 but it is not    ",
  " intended for interpolation. Hence all the complexities of adding    ",
  " zero traces or using DFT are left out. Instead some facilities for  ",
  " the mute filter are added, including plotting of the pass mute      ",
  " window and plots iteration for iteration.                           ",
  "                                                                     ",
  " A sparseness constraint can be enforced (in the time domain)        ",
  " by using more than one iteration (iter_end > 1).                    ",
  " This improves resolution but prediction is worse than without       ",
  " sparseness.                                                         ",
  "                                                                	",
  " A very powerful way of controlling sparseness is to set a           ",
  " upper limit (with Wmthreshold). If Wmthreshold is low then the      ",
  " model weights are not that different between different regions.     ",
  "                                                                	",
  " In this version the model weights are calculated using windows      ",
  " such that sparseness is enforced in regions rather than points      ",
  " The size of this window is hard code for now, but they can be       ",
  " changed in the code radonfk_wtcgls.cpp, function weight_test.       ",
  " the size for now is 5x3 (5 units in time x 3 in offset              ",
  "                                                                     ",  
  " This code has the capability of plotting iteration per iteration    ",
  " such that convergence and artifact creation can be studied carefully",
  "                                                                	",
  " To do this you need to change one line in radonfk_wtcgls.cpp        ",
  " In line two change #define plot 0 to #define plot 1                 ",
  "                                                                     ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Standard input : data file  (offset time domain)          		",
  " Standard output : Interpolated data file                            ",
  "                                                                     ",
  " nv=5            number of velocity laws in the model space          ",
  " dv=200          increment between velocity laws                     ",
  " modelfile=      migrated domain (sufile)                            ",
  " vmig=2000       Velocities                                          ",
  " tmig=0          Times for velocities                                ",
  " eps1=0.8        Quantil to use in the numerator hyperparameter      ",
  "                 for Wm (not used if iter_end=1)                     ",
  " eps1=0.8        Quantil to use in the denomirator hyperparameter    ",
  "                 for Wm (not used if iter_end=1)                     ",
  " itercg=3        Internal iterations for CG (2-5 are enough)         ",
  " iter_end=1      External iterations for CG (>1 for sparseness)      ",
  " step=1          step scale for CG (use 0.8-0.9 for more stable      ",
  " norm=0          0 Cauchy 1 L1                                       ",
  " testadj=0       Calculate adjoint test (Claerbout,1992)             ",
  " verbose=1       =1  verbose output                                  ",
  " plot=0          =1  additional plots (debug and testing)            ",
  " plot2=0         =1  additional plots for testing migration operator ",
  "                 without LS                                          ",
  " tmin_m=0        min time for pass mute window                       ",
  " tmax_m=nt*dt    max time for pass mute window                       ",
  " ihmin_m=nh/2+2  first trace for pass mute window                    ",
  " ihmax_m=nh      last trace number for pass mute window              ",
  " thres_m=0.2    values less than threshold are removed inside window ",
  " slope_m=3      slope in pass mute window (ih top / ih bottom)       ",
  " lstaper=0	    length of side tapers (# of traces)		        ",
  " lbtaper=0	    length of bottom taper (# of samples)  		",
  " Wmthreshold=1000 upper limit for the model weights                  ",
  "  		                                        		",
  " Example:                                                            ",
  " Step 1:   Fill gaps with zero traces                                ",
  " sufill  < filein > tempfile;                                        ",
  " ( Zero traces are identified with tr.trid=2, because inside         ",
  "   suradonfk0 Wd needs to be set to 0 for these added traces         ",
  "   and after suradonfk0 we need to remove dead traces and set ntr to ",
  "   the original )                                                    ",
  "                                                                     ",
  " Step 2: filter and remove zero traces                               ",
  "                                                                     ",
  "   suradonfkcube < tempfile vmig=3600,3000 tmig=0,4 plot=1           ",
  "   itercg=5 iter_end=1 step=0.98 testadj=1 eps1=9.8e-1 eps2=9.8e-1   ",
  "   modelfile=muted_model ascale=1                                    ",
  "   tmin_m=0.8 tmax_m=2 ihmin_m=128 ihmax_m=158 thres_m=0 slope_m=27  ",
  "   mute=1 norm=0 | suremovedead hmax=800 | susetntr > fileout        ",
  "   cp migrated.su non_muted_model                                    ",
  "                                                                     ",
  " Step 3: the predicted traces in fileout can be subtracted           ",
  "                                                                     ",
  "   suadapdiff filein fileout > fileclean                             ",
  "                                                                     ",
 NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 * Last changes: July : 2002 
 */
/**************** end self doc ***********************************/



int main(int argc, char **argv)
{
  int verbose;
  segy tr; 
  inv_par inv;
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep;

  time_t start,finish;
  double elapsed_time;
  int it,ih,iv;
  float *velarray=0;
  float vel; // after time stretching vel represents constant velocity
  float t0=0;

  float **datain=0;
  float **datains=0;
  float **dataout=0;
  float **dataouts=0;

  float ***model=0;
  float ***models=0;

  float *t, *h;

  int nt, nh, nv; 
  int method;
  int plot; 
  int plot2;
  float dt;
  float kmin;
  int testadj;
  float fmax;
  float ascale;

  // Velocity law and stretching
  float *tmig;
  float *vmig;
  int ntmig;
  int nvmig;
  int itmig;
  float smig;
  float vscale,vstolt,vmin,vmax,ft=0,du,*v,*ut,*tu;
  int nu;
  float dv=0; // velocity increment
  
  float dh=0;

  // data weights
  float **Wd=0;
  float **Wds=0;
 
  float Wmthreshold; // upper limit for the model weights
  float threshold;   // Factor to reduce Wmthreshold and use it as mask
  int   smoothMute;  // =1 apply smoothing into pass-mute zone

  // Number of elements to taper on side and bottom
  int lstaper=0;
  int lbtaper=0;
  ////////////////
    
  fprintf(stderr,"*******SURADONFK*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 0.8;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 0.8;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("step", &inv.step))  inv.step =1;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 3;
  if (!getparint("norm", &inv.norm))  inv.norm =0; 
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =1;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparint("plot2",&plot2)) plot2 = 0;
  if (!getparint("testadj",&testadj)) testadj=0; 
  if (!getparfloat("kmin",&kmin)) kmin=0;
  if (!getparint("nv",&nv)) nv=5;
  if (!getparfloat("dv",&dv)) dv=200;
  if (!getparint("lstaper",&lstaper)) lstaper=0;
  if (!getparint("lbtaper",&lbtaper)) lbtaper=0;
  if (!getparfloat("Wmthreshold",&Wmthreshold)) Wmthreshold=100;
  if (!getparfloat("threshold",&threshold)) threshold=0.8;
  if (!getparint("smoothMute",&smoothMute)) smoothMute=0;

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

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;

  /* for mute in the migrated space we need to set the geometry of the mask */
  int mute; 
  mutemask_par par;

  if (!getparint("mute",&mute)) mute = 0;
  if (!getparfloat("tmin_m",&par.tmin)) par.tmin = 0.; 
  if (!getparfloat("tmax_m",&par.tmax)) par.tmax = nt*dt;
  if (!getparint("ihmin_m",&par.ihmin)) par.ihmin = (int) (nh/2+2);
  if (!getparint("ihmax_m",&par.ihmax)) par.ihmax = (int) (nh);
  if (!getparfloat("slope_m",&par.slope)) par.slope = 3;     
  if (!getparfloat("thres_m",&par.threshold)) par.threshold = 0.2;     
  

  /**************************************************************************/

  if (!getparfloat("fmax",&fmax)) fmax = 0.5/dt;
  fmax = MIN(fmax,0.5/dt);

  /* Correct for a factor 2 the given velocities */
  for (itmig=0;itmig<ntmig;++itmig) vmig[itmig]*=2;

  /* make uniformly sampled rms velocity function of time */
  makev(ntmig,tmig,vmig,vscale,nt,dt,ft,&v,&vmin,&vmax);
  
  /* Stolt migration velocity is the minimum velocity */
  vel = vstolt = vmin;

  /* make u(t) and t(u) for Stolt stretch */
  makeut(vstolt,fmax,v,nt,dt,&ut,&nu,&du,&tu);
  free1float(v);

  if (verbose)
    fprintf(stderr,"New time axis after stretching: nt=%d ,dt=%f \n",nu,du);  

  // Allocate memory for data and model before stretching

  datain=ealloc2float(nt,nh);
  dataout=ealloc2float(nt,nh);
  model=ealloc3float(nt,nh,nv);  // model before stretching
  Wd=ealloc2float(nt,nh);
  velarray=ealloc1float(nv);      
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  /* Create the velocity axis */
  dv*=2; // take care of the 2 factor in Stolt equation 
  velarray[0]=vel;  for (iv=1;iv<nv;iv++) velarray[iv]=velarray[iv-1]+dv;

  memset( (void *) h, (int) '\0', nh * FSIZE);

  // Loop over traces
  // Dead traces are marked with trid = 2. They are weighted to zero
  ih=0;
  do { 
    if (tr.trid==2) for(it=0;it<nt;it++) Wd[ih][it]=0; 
    else  for(it=0;it<nt;it++) Wd[ih][it]=1;
    h[ih]=(float) tr.offset;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  
  if (verbose) fprintf(stderr,"processing %d traces \n", nh);
  /* Time axis */
  /* If the first sample is not zero the code assumes it is zero. */
  /* The only consequence is that the required velocity is a parameter 
     without physical meaning anyway */

  for (it=0;it<nt;it++) t[it]=t0+it*dt;  /* Not implemented for t0 != 0  */

  /***************************************************************** 
     From this point the program has to create two new data spaces
     first: original sampling in offset to new offset sampling (interpolation) 
     second: original time sampling to new time sampling (stretching)
  *******************************************************************/

  dh=(h[nh-1]-h[0])/(nh-1);

 if (plot){ /* additional plots only for degugging  */  
    save_gather(datain,nh,h,nt,dt,"datain.su");
    system("suximage < datain.su perc=98 key=offset title=datain xbox=%d curve=curve1 npair=5 hbox=900 wbox=700 &");  
  }
 /* apply side and bottom tapers */
 if ( (lstaper) || (lbtaper) )
      for (ih=0; ih<nh; ++ih)
	   taper(lstaper,lbtaper,nh,ih,nt,datain[ih],0);


 /**************** Time stretching  ******************/
  /*** The new data space has different nt (due to the time stretching) */
  if (verbose) fprintf(stderr,"Time stretching: old nt=%d, New nt=%d\n",nt,nu);

  /* define arrays for data and model after stretching */
  datains=ealloc2float(nu,nh);
  dataouts=ealloc2float(nu,nh);
  models=ealloc3float(nu,nh,nv);  // model after stretching
  Wds=ealloc2float(nu,nh);
  
  stretch(datain,datains,nt,nu,nh,t,tu,ut,dt,du,1);
  stretch(Wd,Wds,nt,nu,nh,t,tu,ut,dt,du,1);
  // test to see the effect of the data weights
  //for (ih=0;ih<nh;ih++) for (it=0;it<nu;it++) Wds[ih][it]=1;
  
  /* If the velocity is constant the time axis was not created before */
  if (nu==nt){
    tu=ealloc1float(nt);
    ut=ealloc1float(nt);
    for (it=0;it<nt;it++) tu[it]=ut[it]=t[it];
  }

  /* Muting mask ***********************************************/
  /* Apply correction to the mute time by the stretching factor*/
  float stfact=(1.0*nu)/(1.0*nt);
  fprintf(stderr,"stfact=%f\n",stfact);
  par.tmin*=stfact;
  par.tmax*=stfact;
  /* Write the new mute parameters to a file for plots */
  write_curve_mute(par);
  
  float **  M  = ealloc2float(MAX(nu,nt),nh);  // Mask (pass mute)
  float *** Wm = ealloc3float(MAX(nu,nt),nh,nv);  // Model Weigths for CG
  
  /************* Test for operators ***********/  
  if (testadj) adjteststoltzcube(nu,nh,nv,tu,h,velarray);

  if (plot2){
    stoltzop2(datains,models,nu,nh,nv,tu,h,velarray,1);
    save_gather(models,nv,nh,nu,0.004,"modeladj.su");
    system("suximage < modeladj.su perc=100 title=modeladj &");
  }	
   
  if (plot2){
    stoltzop2(datains,models,nu,nh,nv,tu,h,velarray,0);
    save_gather(datains,nh,nu,0.004,"datain.su");
    system("suximage < datain.su perc=100 title=datain &");
  }	 

  if (plot2){
    save_gather(Wds,nh,nu,0.004,"Wd");
    system("suximage < Wd perc=100 title=Wd &");
  }
  /********************************************************/
  TRACE;
  stoltz_wtcgls_cube(datains,models,Wds,h,nh,tu,nu,velarray,nv,inv,Wmthreshold,Wm);
  TRACE;
  // data weights no need any longer 
  free2float(Wds);
  
  /******* Migrated output **************/
  /* array for output and plots */
  if (plot){
       for (iv=0;iv<nv;iv++)    
	    stretch(model[iv],models[iv],nt,nu,nh,t,tu,ut,dt,du,-1);   
       save_gather(model,nv,velarray,nh,h,nt,t,"migrated.su");
       system("suximage < migrated.su perc=98 title=migrated xbox=%d curve=curve1 npair=5 hbox=900 wbox=700 &");
  }	 
  TRACE;
  if (mute){ 
       for (iv=0;iv<nv;iv++) mutemask(M,models[iv],nh,nu,du,par);

       if (0){
	    save_gather(M,nh,h,nu,du,"mask.su");
	    system("suximage < mask.su perc=100 title=mask &");
       }
       TRACE;
       if (1){
	    save_gather(Wm,nv,velarray,nh,h,nu,du,"Wmthres0.su");
	    system("suximage < Wmthres0.su perc=100 legend=1 title=Wmthresh0 &");
       }
       TRACE;
       for (iv=0;iv<nv;iv++) makeMask(Wm[iv],nh,nu,Wmthreshold*threshold);
       
       if (1){
	    save_gather(Wm,nv,velarray,nh,h,nu,du,"Wmthres.su");
	    system("suximage < Wmthres.su perc=100 legend=1 title=Wmthresh &");
       }
       TRACE;
       if (smoothMute){
	    // smoothing
	    int nl=11;  //  npoints left hand side
	    int nr=11;  //  npoints left hand side
	    int flag=2;  // 1 rectangular, 2 triangular
	    for (iv=0;iv<nv;iv++)  smoothing(Wm[iv],nu,nh,nl,nr,flag);
       }
       TRACE;
       if (1){
	    save_gather(Wm,nv,velarray,nh,h,nu,du,"WmthresSmooth.su");
	    system("suximage < WmthresSmooth.su perc=100 legend=1 title=WmthreshSmooth &");
       }
       
       for (iv=0;iv<nv;iv++) AtimesB(Wm[iv],M,nh,nu);       
       for (iv=0;iv<nv;iv++) AtimesB(models[iv],Wm[iv],nh,nu);

  }
  TRACE;
  free3float(Wm);       
  free2float(M);
  /**************************************************************/

  /* unstretched  migrated model and plot */
  for (iv=0;iv<nv;iv++){
    stretch(model[iv],models[iv],nt,nu,nh,t,tu,ut,dt,du,-1);
    if (0)  xplotgather(model[iv],nh,nu,du,"models.su",0,"perc=97 xbox=600 legend=1&");
  }
  TRACE;
  modelfilep=efopen(modelfile,"w");
  TRACE;
  for (iv=0;iv<nv;iv++){
    erewind(stdin);
    for (ih=0;ih<nh;ih++){ 
      fgettr(stdin,&tr);
      memcpy((void *) tr.data,(const void *) model[iv][ih],nt*sizeof(float));
      /* If the scale is wrong (sometimes for dft=1) 
	 apply a scale correction ascale*/
      //for (it=0;it<nt;it++) tr.data[it]=ascale*model[iv][ih][it];
      tr.offset=(int) h[ih];
      tr.f2=velarray[iv];
      tr.ntr=nh*nv;
      fvputtr(modelfilep,&tr);
    }
  }
  TRACE;
  efclose(modelfilep);

  if (0){
  /** plot */ 
    char buf[80];
    sprintf(buf,"suxwigb < %s key=offset perc=99  title=%s & \n",modelfile,modelfile);
    system(buf);  
  }
  /******** End of migrated output **********/

  /******** Interpolated output *************/
  TRACE;
  stoltzop2(dataouts,models,nu,nh,nv,tu,h,velarray,0);
  stretch(dataout,dataouts,nt,nu,nh,t,tu,ut,dt,du,-1);
  TRACE;
 /* apply side and bottom tapers */
 if (( (lstaper) || (lbtaper) )&&(0))
   for (ih=0; ih<nh; ++ih)
     taper(lstaper,lbtaper,nh,ih,nt,dataout[ih],1);

 
  rewind(stdin);
  for (ih=0;ih<nh;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    for (it=0;it<nt;it++) tr.data[it]*=ascale;
    tr.offset=(int) h[ih];
    tr.ntr=nh;
    if (Wd[ih][0]==0) tr.trid=2; // dead trace
    fputtr(stdout,&tr);
  }

  /******** End of interpolated output **********/

  free1float(velarray);
  free2float(Wd);
  free1float(t);
  free1float(h);
  free2float(datain);
  free2float(dataout);
  free3float(model);

  free2float(datains); 
  free2float(dataouts);
  free3float(models);

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


















