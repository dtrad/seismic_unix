#include "su.h"
#include "clibrary.h"

void hrrtf(float *pos, float **d ,float dt,float eps1,float eps2,float qmin,
	   float **m,float *q,float dq,float freq, float eps, float *bad, 
	   int taper)


  /*       subroutine: hrrt.cpp
       HIGH RESOLUTION PARABOLIC TRANSFORMS
     method=1 Least Square Radon with Levinson
     method=2 Weighted Conjugate Gradient Least Squares
     method=3 Cholesky overdetermined
     method=4 Cholesky undertermined
     method=5 CGLS for Toeplitz (SU)
     method=6 SVD overdetermined  
     method=7 CG non linear with linear search
     method=8 LSQR 
     method=9 CG (Mauricio's method)
     method=10 Very fast CG with circulant matrices
     method=11 BCG from NR 
     method=12 Simplest CGLS   
    
     Daniel Trad
     E-mail: dtrad@geop.ubc.ca
  */
{

     int i, nqt, ih;
     float  qmax, t0, dqt, qmaxt, fmax;
     float  dx_av, dx_max, dh;
     float *Wd; 
     //float  *pdq, *pqmaxt, *pqmax;
     extern int nt, nh, nq, method, iter_end, rtmethod, norm, nbad;
     extern float factor;

     if ((Wd=alloc1float(nh))==NULL)
       err("cannot allocate memory for Wd\n");


     fprintf(stderr,"nq=%d, nt=%d, nh=%d\n, rtmethod=%d\n",nq,nt,nh,rtmethod);
     
     t0=0;
                
     interval(pos,nh,&dx_max,&dx_av);
     fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
     // maximum and minimum q parameter to invert

     if (freq==0)
           fmax=1/(2*dt);
     else
           fmax=freq;
     
     radon_param(fmax,pos,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
     dq=1.778932e-06;  //dq=1.190476e-05; // Test remove later...........................
     fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
     fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
     ///////////////////////////////////////////////////////////////////
     // Compute inv(CD)=WdT Wd       
     for (ih=0;ih<nh;ih++) Wd[ih]=1.;
     if (taper==1)
       for (ih=0;ih<5;ih++) {
	 Wd[ih]=1-exp(-(ih+.3));
	 Wd[nh-ih-1]=Wd[ih];
       }
   
     for (ih=0;ih<nh;ih++)
       for (i=0;i<nbad;i++) 
	 if ((ih+1) == int(bad[i])) Wd[ih]=0;

     // for (ih=0;ih<nh;ih++) fprintf(stderr,"Wd[%d]=%f\n",ih,Wd[ih]);
     ///////////////////////////////////////////////////////////////////

     fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
     fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
     for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
     if (rtmethod==1)
       fprintf(stderr,"LRT........\n");
     else if (rtmethod==2)
       fprintf(stderr,"PRT........\n");


     if (method==1) toepradon(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==2) wtcgls0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax,Wd);
     else if (method==3) cholover(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==4) cholund(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==5) toepradoncg(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==6) svdover(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax); 
     else if (method==7) cgnl0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==8) lsqr0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==9) conjgradrt3(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);
     else if (method==10) conjgradrt4(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax,Wd);
     else if (method==11) conjgrad3r(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);   
     else if (method==12) cgls0(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax);   
     else if (method==13) rad_wpcgnr(pos,d,dt,m,q,dq,eps1,eps2,eps,fmax,Wd);
     else fprintf(stderr,"method no implemented/n");

     free1float(Wd);
     return;
}
 











