/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUEOMIG - First attempt for Equivalent Offset migration             ",
  "                                                                     ", 
  "	                                                          	",
  " 	   								",
  " sueomig0 < stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " vel=2000        Velocity of migration (in meters) 			",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

segy tr,trf; // for malloc use *trr
void eomig1(float **csp, float *he, float *m, float *t, float *vel) ;
void rjwfilter(float **d,int nt,int nh, float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

int nt,nh,nx,nhmax; // Size of volume x,h,t
float dt,dh,dx,dhe; // Sampling
int main(int argc, char **argv)
{
  int j, i, k, ih, ix, itr; // General counters 
  register int it;
  float *d;      /* single trace */
  float **csp;   /* Common Scattering gather */
  float  *m;    /* Final migrated trace   */
  float *x;      /* axis for CDPs */
  float *he;     /* axis for equivalent offset */
  float  *t;     // time axis for input and output
  float  h;      // halfoffset
  float  cdp;    // midpoint coordinate

  /* Velocity */
  float *vel;    // Array for velocity function
  float *time;   // Array for time in velocity function 
  int   nvel;    // Number of velocities
  int   ntime;   // Number of time points for Velocities
  float *velint; // Array for velocity function after interpolation
  float vel1;    // Temporal variable for velocity

  int ntr;     // Total number of input traces
  extern float dt,dh,dx,dhe;	
  const double  pi=acos(-1.); 
  float xm;    // midpoint distance from trace to scp = cdp_tr - x_scp
  extern int  nhmax;      //Maximum number of traces for the CSP
  extern int  nh;
  extern int  nt;
  extern int  nx;         // number of output cdps
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  hemin;          //limit for the equivalent offset
  float  hemax;          //limit for the equivalent offset     
  float  *hec;           /* equivalent offset for every trace */ 
  int    ihec;           /* equivalent offset index */
  float  hec2;           
  float  *tc;            /* time when equivalent offset bin increases */
  int    itc;            /* corresponding index for tc[i] */   
  int    nhec;           /* Number of offset bins at every trace */
  int verbose;		/* flag for echoing info		*/
  char *tmpdir;		/* directory path for tmp files		*/
  int output;           /* =0 output are CSP gathers,           */ 
			/* =1 output are  migrated traces       */
                        /* =2 output are original traces, used for test */
  float xm2;             // xm*xm
  float h2;              // h*h
  float crossterm;       // 4*xm^2*h*2/vel
  float *fold;           /* Number of traces added at each offset of a CSP */
  float fmax;
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  
  // Get parameters 
  //if (!getparfloat("vel", &vel))  vel = 2000;
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparint("nx", &nx))  nx = 100;
  if (!getparfloat("hemin", &hemin))  hemin = 0;
  if (!getparfloat("hemax", &hemax))  hemax = 1000;
  if (!getparint("nhmax", &nhmax))  nhmax = 50;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparint("output", &output)) output = 1;

  nvel = countnparval(1,"vel");
  ntime = countnparval(1,"time");
  intprint(nvel);
  intprint(ntime);
       
  if (nvel!=ntime)
    err("number of vel and time values must be equal");  


  
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  fmax=0.8/(2*dt);
  
  fprintf(stderr,"nt=%d,dt=%f,nx=%d\n",nt,dt,nx); 
  
  // Store traces in the disk for later use
  tracefp = etmpfile();
  headerfp = etmpfile();
  if (verbose) warn("using tmpfile() call");
  ntr = 0;
  do {
    ntr++;
    //filt(tr.data,nt,dt,fmax,0,50,trf.data);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data, FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);

  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");

  if ((csp=ealloc2float(nhmax,nt))==NULL)
    fprintf(stderr,"***Sorry, space for csp could not be allocated\n");
  
  if ((m=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");

  if ((he=ealloc1float(nhmax))==NULL)
    fprintf(stderr,"***Sorry, space for he could not be allocated\n"); 
   
  if ((fold=ealloc1float(nhmax))==NULL)
    fprintf(stderr,"***Sorry, space for fold could not be allocated\n");   

  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");

  if ((tc=ealloc1float(nhmax))==NULL)
    fprintf(stderr,"***Sorry, space for tc could not be allocated\n");

  if ((hec=ealloc1float(nhmax))==NULL)
    fprintf(stderr,"***Sorry, space for hec could not be allocated\n");

  // arrays for nmo model
  
  if ((time=ealloc1float(ntime))==NULL)
    fprintf(stderr,"***Space for time could not be allocated\n");
  
  if ((vel=ealloc1float(nvel))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");

  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  if (!getnparfloat(1,"vel",vel)) vel[0] = 2000.0;
  if (!getnparfloat(1,"time",time)) time[0] = 0.0;        
  for (it=0; it<ntime; ++it){
    fprintf(stderr,"vel[%d]=%e\n",it,vel[it]);
    fprintf(stderr,"time[%d]=%e\n",it,time[it]);
  }
  for (it=1; it<ntime; ++it)
    if (time[it]<=time[it-1])
      err("tnmo values must increase monotonically");




  /* Create axis for output cpds, equivalent offset and time */
  
  dx=(cdpmax-cdpmin)/(nx);
  dhe=(hemax-hemin)/(nhmax);      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;
  for(ih=0;ih<nhmax;ih++) he[ih]=hemin+ih*dhe;  
  
  /* Create axis for velocities */
  intlin(ntime,time,vel,vel[0],vel[nvel-1],nt,t,velint);


  fprintf(stderr,"**************************\n");  
  /* Loop to compute one CSP gather */
  for (ix=0;ix<nx;ix++){
    
    /* Initialize to zero the ouput trace and the csp gather corresp to ix */

    for (ih=0;ih<nhmax;ih++) for (it=0;it<nt;it++) csp[it][ih]=0; // CSP trace
    for (ih=0;ih<nhmax;ih++) fold[ih]=1; 

    if (output==1) for (it=0;it<nt;it++) m[it]=0; // Output trace

    /* Loop for mapping one trace to the CSP gather */
    for (itr=0;itr<ntr;itr++){     
      efread(tr.data, FSIZE, nt, tracefp);
      efread(&tr,HDRBYTES,1,headerfp);
      if (output==2){
	fprintf(stderr,"tr.tracr=%d\n",tr.tracr);
	puttr(&tr);      
      }
      
      h=(float) tr.offset;
      h/=2;  // halfoffset
      cdp=(float) tr.cdp;
      xm=fabs(cdp-x[ix]);
     
      xm2=xm*xm;
      h2=h*h;
      //crossterm=4*xm2*h2/(vel*vel);
     
      /* Precompute times for changing offset bin 
      hec[i] is the equivalent offset at bin i
      ihec   is the corresponding index
      tc[i]  is the time at which euivalent offset bin changes
      itc is the corresponding offset
      hec2 is hec[i]*hec[i]
      xm2  is xm*xm
      h2   is h*h
      crossterm is 4*xm2*h2/(vel*vel);
      */
      float hecmax;
      if ((xm!=0)&&(h!=0)){
	i=0;
	tc[i]=2*xm/velint[0];  // First time to use
	itc=int (tc[i]/dt);if (itc<nt) vel1=velint[itc];else vel1=velint[nt];
        hec[i]=xm;
        hecmax=sqrt(xm2+h2);    
	do{
	  i++;
	  //hec2=xm2+h2-(crossterm/(tc[i-1]*tc[i-1]));
          hec[i]=hec[i-1]+dhe;
      	  hec2=hec[i]*hec[i];
          /* First estimate of time uses previous velocity */
	  tc[i]=2*xm*h/(vel1*sqrt(xm2+h2-hec2));
	  itc=int (tc[i]/dt);
          if (itc<nt) vel1=velint[itc];
          else vel1=velint[nt];
	}while((tc[i] <= t[nt-2]) && (i<nhmax) && (hec[i]<=hecmax));
	nhec=i++;

      }
      else{
        tc[0]=2*(xm+h)/velint[0];
        hec[0]=sqrt(xm2+h2);
	nhec=1;
      }
      tc[nhec]=t[nt-2];
      hec[nhec]=hecmax;
      //for (i=0;i<=nhec;i++) fprintf(stderr,"tc[%d]=%f,hec[%d]=%f,nhec=%d,hecmax=%f\n",i,tc[i],i,hec[i],nhec,hecmax);

      it=int (tc[0]/dt);
      // it=0; 
      for (i=0;i<nhec;i++){
        itc=(int) (tc[i+1]/dt);
	ihec=(int) (hec[i]/dhe);
	if (ihec < nhmax){
	  fold[ihec]+=1;    
	  do{
	    it++;
	    csp[it][ihec]+=tr.data[it];	    
	  }while(it<=itc);
	}
        else it=itc;
      }
    }  
    /* the loop for the input traces finishes here */
  
    //for (ih=0;ih<nhmax;ih++) for (it=0;it<nt;it++) csp[it][ih]/=fold[ih];
    
    if (output==1) eomig1(csp,he,m,t,velint);
    
    if (output==0){
      for (itr=0;itr<nhmax;itr++){
	for (it=0;it<nt;it++) tr.data[it]=csp[it][itr];
	tr.cdp=(int) x[ix]; // front of the CSP
	tr.dt=(int) (dt*1e6);       
	tr.ntr=nhmax;
	tr.ns=nt;
	tr.tracl=itr+1;
	tr.tracr=itr+1;
	tr.offset=(int) (2*he[itr]);
	tr.sx=(int) x[ix];
	tr.gx=(int) x[ix];
	puttr(&tr);
	//fprintf(stderr,"he[%d]=%e,tr.offset=%d\n",itr,he[itr],tr.offset);    
      }
    }         
    erewind(tracefp);
    erewind(headerfp);    
    
    
    //for (ix=0;ix<nx;ix++) 
    //  for (it=0;it<nt;it++)
    //    m[it][ix]/=(vel*sqrt(pi*MAX(t[it],1e-3)));
    //rjwfilter(m,nt,nx,dt);
    if (output==1){
      for (it=0;it<nt;it++)
	m[it]/=(vel[0]*sqrt(pi*MAX(t[it],1e-2)));
      for (it=0;it<nt;it++)
	tr.data[it]=m[πœØå——ÁÎÎø˜ÏÛ+ÆhÍ˚¸Œ|≥ØkÒD=øH=}±Gá'	¸ÿÜmz}ˇ ?∂>Cw‘ÌçÀ≈˝€V%£svë∫”≤·∑.X˜`aå’≤G±2}Î+{|05y\Àïq˜kª#±.}ˇ0.(ÚÎûÌ40;ˆ—˙µ(ˆ≈>HπËS,F⁄Œæü¡¢–r‰=˚`?∫u±}ªÛ1ı>ø6"πµı´Ù6$∑“v‹ÉÅ…>]™Ô7íGu6ÙÅ·8qá„JÙÅëÄ9Ø¶¬[WVôbéﬁ˜Ú=Ωﬂ¨AñÀ,Ú˜Pπ}Â~Jπ/•Ì˚©%>∑éØ^œ~]Ï˚/{D¸∫˜˚ÒW)˝˜P˙≈=ÈW≠—“g¨(ßÿ^{FÂÆÌníE∏k¶j√ﬁ˙,€ÆÌµÔÎvYTª∂è~Æ[é¨RÆÖ–0ãg?^ﬁ‰oÛÏg˛∑yˆW±z‘¶∏gü¢övqØñSÑ|ÿf∑\∫™œØÇg˜Ég˜h±Jû=WZ√‡Ÿ„lé<Û£œy±Î±}ˆº.n=W~∫uçÆÌÜá(±Ë—$∫∂gSÎèòbù≤±–`Ñ[œôbó¸ã´gmÊ⁄Fd\0øÛD°ßÇæz¨iD†kªÁh1Sπ˛µw´&b*∂˚¸w˛ºU˙;»ÛÁ´/õ<?JUıté˝™¢nR5Á≥é©Ôm);
  free2float(csp);
  free1float(d);
  efclose(headerfp);
  efclose(tracefp);  
  return EXIT_SUCCESS;
}




















