#include "su.h"
//#include "interpfk.h"

//#define tmin 0.7
//#define threshold 0.5
typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;


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













