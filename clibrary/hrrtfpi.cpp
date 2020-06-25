#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
void conjgradrt3pi(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax);
void cholundpi(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
void hrrtfpi(float *pos, float **d ,float dt,float eps1,float eps2,float qmin,
	   float **m,float *q,float dq,float freq, float eps)

	/*
       subroutine: hrrt.cpp
       HIGH RESOLUTION PARABOLIC TRANSFORMS
       Mauricio D. Sacchi, Daniel Trad
       E-mail: dtrad@geop.ubc.ca
       *********************************************************
       GENERAL DESCRIPTION:
       Given the seimic data  computes a parabolic stack. The P.Stack 
       is then maped back to the data space. Undesired components
       can be masked in the velocity gather before the reconstruction 
       of the data. Doing so, we can isolate and keep only
       part of the information contain in the original data.
	*/
{

     int i, nqt;
     float  qmax, t0, dqt, qmaxt, fmax;
     float  dx_av, dx_max, dh;
     //float  *pdq, *pqmaxt, *pqmax;
     extern int nt, nh, nq, method, iter_end, rtmethod, norm;
     fprintf(stderr,"nq=%d, nt=%d, nh=%d\n, rtmethod=%d\n",nq,nt,nh,rtmethod);
     
     t0=0;
                
     interval(pos,nh,&dx_max,&dx_av);
     fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
     // maximum and minimum q parameter to invert

     if (freq==0)
           fmax=1/(2*dt);
     else
           fmax=freq;

     radon_param(fmax,pos[0],pos[nh-1],dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod);
         
     fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
     fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
     for (i=0;i<nq;i++)  q[i]=qmin+i*dq;

     if (rtmethod==1)
       fprintf(stderr,"LRT........\n");
     else if (rtmethod==2)
       fprintf(stderr,"PRT........\n");

     if (method==2) conjgradrt3pi(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==3) toepradon(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==4) cholover(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==5) toepradoncg(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==9) choleskyrt(pos,d,dt,m,q,dq,eps1,eps2,fmax);
     else if (method==7) cholundpi(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==10) conjgradrt4(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==11) svdover(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==12) conjgrad3r(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);   
     else if (method==13) cgls0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);   
     else if (method==14) wtcgls0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);  
     else if (method==15) cgnl0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     //else p_stack(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else fprintf(stderr,"method no implemented/n");
     return;
}
 











