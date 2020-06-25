#include "dan.h"
#include "su.h"

/* SImple mute routine. 
   The parmute parameter defines the boundary between muted and non muted spaces
   A taper is applied for smoothness but it may be too smooth in many cases,
   By default is off.
   the minimum time to mute follows a slope give by it0-(iq-nmute)*slope
   This slope is a function of nt/nq. In general values between 2-5 must work fine */

#define VERBOSE 1
#define SLOPE 2
#define TAPER 0


typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;

int taper(float **data, int nt, int nh, int ntaper,int flag);

void mutemask(float **M, float **d, int nh, int nt, float dt, mutemask_par par)
     /* 
	Function to compute a muting mask for elimination of 
	diffractions (kencmpwin.su example).
	structure mutemask contains:
	float tmin, float tmax, int ihmin, int ihmax, 
	float threshold, float slope
 
	to use threshold, the mask is data dependent
     */

{
  int ih,it;
  int itmin;
  int itmax;
  int shift;
  //int ihmin;
  //int ihmax;

  itmin=(int) (par.tmin/dt);
  itmax=(int) (par.tmax/dt);
  //ihmin=(int) (nh/2+1);
  //ihmax=(int) (nh/2+50);
  fprintf(stderr,"ihmin=%d,ihmax=%d,itmin=%d,itmax=%d,slope=%f,threshold=%f\n",
	  par.ihmin,par.ihmax,itmin,itmax,par.slope,par.threshold);  

  for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)
	if ((it>itmin)&&(it<itmax)){
	  shift=(int) ((it-itmin)*dt*par.slope);
	  if ((ih>par.ihmin+shift)&&(ih<par.ihmax+shift))
	    if (fabs(d[ih][it])>par.threshold)
	      M[ih][it]=1;
	      //M[ih][it]=M[nh-ih][it]=1;
	}
  else
    M[ih][it]=0;
  return;

}

void mutemask2(float **M, int nh, int nt, float dt, mutemask_par par)
     /* 
	Function to compute a muting mask for elimination of 
	diffractions (kencmpwin.su example).
	structure mutemask contains:
	float tmin, float tmax, int ihmin, int ihmax, 
	float threshold, float slope 
	
	this version is not data dependent (it does not use threshold)
     */

{
  int ih,it;
  int itmin;
  int itmax;
  int shift;
  //int ihmin;
  //int ihmax;

  itmin=(int) (par.tmin/dt);
  itmax=(int) (par.tmax/dt);
  //ihmin=(int) (nh/2+1);
  //ihmax=(int) (nh/2+50);
  fprintf(stderr,"ihmin=%d,ihmax=%d,itmin=%d,itmax=%d,slope=%f,threshold=%f\n",
	  par.ihmin,par.ihmax,itmin,itmax,par.slope,par.threshold);  

  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++)
      if ((it>itmin)&&(it<itmax)){
	shift=(int) ((it-itmin)*dt*par.slope);
	if ((ih>par.ihmin+shift)&&(ih<par.ihmax+shift))
	  M[ih][it]=1;
      }
      else
	M[ih][it]=0;
  return;
  
}

void makeMask(float **A, int nh, int nt, float threshold)
{
  int ih,it;
  float expo;
  float tolx  =  3;
  float ntolx = -3;
  float debug =  0;
  
  fprintf(stderr,"start makeMask\n");
  for (ih=0;ih<nh;ih++)
       for (it=0;it<nt;it++){
	    if (1){
		 /* use factor = 0.1 for best results */
		 expo= 0.01 * ( fabs(A[ih][it])- threshold  );
		 if      ( expo > tolx) A[ih][it] = 1.;
		 else if ( expo < ntolx) A[ih][it] = 0.;
		 else{
		      A[ih][it] = 1 - exp(-exp(expo));
		      if (debug)
			   fprintf(stderr,"expo = %f , exp(expo) = %f, exp( - exp(expo) =%f \n",
				   expo, exp(expo), exp(-exp(expo)));
		 }
	    }
	    else{
		 if  ( fabs(A[ih][it]) > threshold ) A[ih][it]=1;
		 else A[ih][it] = 0;
	    }
       }
  fprintf(stderr,"End makeMask\n");
  return;
}
 


void AtimesB(float **A, float **B, int nh, int nt)
{
  int ih,it;
  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++)
      A[ih][it]*=(B[ih][it]);
  fprintf(stderr,"AtimesB\n");
  return;
}
 
void write_mutemask(mutemask_par par, char* fileName)
{
  
  FILE *fp;

  int ihmax2;
  int ihmin2;
  int shift;

  /* Calculate the corners of the mute window at the bottom */

  shift=(int) ((par.tmax-par.tmin)*par.slope);
  ihmax2=par.ihmax+shift;
  ihmin2=par.ihmin+shift;

  fp=efopen(fileName,"w");
  
  fprintf(fp,"%f %d\n",par.tmin, par.ihmin);
  fprintf(fp,"%f %d\n",par.tmin, par.ihmax);
  fprintf(fp,"%f %d\n",par.tmax, ihmax2);
  fprintf(fp,"%f %d\n",par.tmax, ihmin2);
  fprintf(fp,"%f %d\n",par.tmin, par.ihmin);

  efclose(fp);

  fprintf(stderr,"%f %d\n",par.tmin, par.ihmin);
  fprintf(stderr,"%f %d\n",par.tmin, par.ihmax);
  fprintf(stderr,"%f %d\n",par.tmax, ihmax2);
  fprintf(stderr,"%f %d\n",par.tmax, ihmin2);
  fprintf(stderr,"%f %d\n",par.tmin, par.ihmin);


  
  return;
}

void write_curve_mute(mutemask_par par)
{
  
  FILE *fp;
  fp=efopen("curve1","w");
  
  fprintf(fp,"%f %i\n",par.tmin, par.ihmin);
  fprintf(fp,"%f %i\n",par.tmin, par.ihmax);
  fprintf(fp,"%f %i\n",par.tmax, par.ihmax);
  fprintf(fp,"%f %i\n",par.tmax, par.ihmin);
  fprintf(fp,"%f %i\n",par.tmin, par.ihmin);

  efclose(fp);
  
  return;
}

void muting(float **model, int nq, int nt, float *q, float *t, float parmute, float t0, int side)
{

  // side 2 --> right mute ; side 1 left mute
  int iq, it;
  int nmute=0;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt); // define a minimum time for mute2
  int slope=SLOPE; // Defines a slope in the muting

  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;

  if (VERBOSE) fprintf(stderr,"MUTING at nmute=%d \n",nmute);
  
  if (TAPER) taper(model,nt,nq,nmute,side);
  else{
    if (side==1)  
      for (iq=0;iq<MIN(nq,nmute);iq++) for (it=0;it<nt;it++) model[iq][it]=0; 
    else   
      for (iq=MIN(nq,nmute);iq<nq;iq++) for (it=0;it<nt;it++) model[iq][it]=0; 
  }

  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(it0-(iq-nmute)*slope);it++) 
      model[iq][it]=0; 
  
  return;
}

int taper(float **data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;
  if (ntaper<=0) return EXIT_SUCCESS;
  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    //taper[k] = s*s;
    taper[k]= pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 

  if(flag==4){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  if(flag==5){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih][it]*=scalemute;
      }
    } 
  
  
  return EXIT_SUCCESS;
}

int taper(float *data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;

  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    taper[k] = pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih*nt+it]*=scalemute;
      }
    } 

  if(flag==5)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=0;   
  
  return EXIT_SUCCESS;
}


int taper(float *data, int nh, int ntaper)
{
  float *taper;
  int k, ih;
  float s;
  
  taper=ealloc1float(ntaper);

  for (k=0;k<ntaper;++k){ 
    s=sin(k*PI/(2*ntaper));
    taper[k]=s;
  }  
  taper[0]=1e-5; 
  /* Taper at the left end of the data set */
  for (ih = 0; ih < ntaper; ++ih) data[ih] /= taper[ih];

  /* Taper at the right end of the data set */
  for (ih=nh - ntaper ; ih < nh; ++ih)  data[ih] /= taper[nh - ih - 1];

  return EXIT_SUCCESS;

}




 
 


