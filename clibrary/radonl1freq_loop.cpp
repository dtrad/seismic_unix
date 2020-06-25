#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "radonl1freq.h"

/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void radonl1freq_loop(float **data, float **model, float *t, float *q, float *h, 
		      int nt, int nq, int nh, float dt, int itercg, int iter_end, 
		      op_param op2)
{
  int ih,iq;
  int freq, nf, nfft;
  complex **m2=0;
  complex **d2=0;
  complex czero;
  complex **L=0;
  complex **L1=0;
  complex **L2=0;
  float w,wa,df;
  float dq1;
  float dq2;
  float *q1=0;
  float *q2=0;
  float *ph1=0;
  float *ph2=0;
  int nq1=op2.nq1;
  int nq2=op2.nq2;
  //const double  pi=acos(-1.);
  float *d2r=0;
  float *m2r=0;
  float **LR=0;
  int maxfreq1;
  int maxfreq2;
 
  czero.r=czero.i=0;
  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh); 
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);

  ////////////////////////////////////////////////////

  radon_param_beam(op2.fmax1,op2.fmax2,h,nh,q, nq, op2.qmin1,op2.qmin2,q1,q2,ph1,ph2,
		   nq1,nq2, op2.depth1,op2.depth2,op2.rtmethod1,op2.rtmethod2,op2.factor1,
		   op2.factor2,&dq1,&dq2);

  ///////////////////////////////////////////////////
  fft_parameters(nt,dt,&nfft,&nf,&df);
  fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);

  floatprint(df);

  maxfreq1=(int) (op2.fmax1/df);
  maxfreq2=(int) (op2.fmax2/df);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d, , dt=%f, df=%f\n",maxfreq1,
	  maxfreq2,dt,df);
  w=0;
  for (freq=1;freq<maxfreq1;freq++){
     w=2*PI*freq*df;
     wa=freqweight(freq,df,op2.fmax1-10,op2.fmax1);
     
     if (op2.rtmethod1==3) matrix(L1,ph1,q1,nh,nq1,w,dq1,op2.rtmethod1);
     else  matrix(L1,h,q1,nh,nq1,w,dq1,op2.rtmethod1);
     if (op2.rtmethod2==3) matrix(L2,ph2,q2,nh,nq2,w,dq2,op2.rtmethod2);
     else matrix(L2,h,q2,nh,nq2,w,dq2,op2.rtmethod2);

     for (ih=0;ih<nh;ih++){
       for (iq=0;iq<nq1;iq++) L[ih][iq]=L1[ih][iq];
       for (iq=0;iq<nq2;iq++) L[ih][iq+nq1]=L2[ih][iq];
     }	 
     
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d===>2 oper\n",freq);

     // Transform the complex system of equation in a real one 
     d2r=ealloc1float(2*nh);
     m2r=ealloc1float(2*nq);
     LR=ealloc2float(2*nq,2*nh);

     complex2real(d2[freq],L,d2r,LR,nh,nq);
  
     // test
     if (0) matmul0 (d2r,LR,m2r,1,2*nh,2*nq);
     
     radonl1freq_karmarkar(d2r,LR,m2r,2*nq,2*nh,iter_end,itercg);  
     matmul0 (d2r,LR,m2r,0,2*nh,2*nq);
     real2complex(d2[freq],d2r,m2[freq],m2r,nh,nq);     

     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq1;iq++)  m2[freq][iq]*=wa;
  }


  for (freq=maxfreq1;freq<nf;freq++){
    for (iq=0;iq<nq1;iq++) m2[freq][iq]=czero;
  }     

  for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  

  for (freq=maxfreq2;freq<nf;freq++) for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
       
  fprintf(stderr,"w=%f\n",w);

  fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf); 

  char buf[128];
  //plotgather(model,nq,nt,0.004,"suxwigb perc=99 title=\"model\"");
  sprintf(buf,"xwigb perc=99 n1=%d n2=%d title=\"model\"",nt,nq);
  plotgather(model,nq,nt,0.004,buf);

  if (0){
    save_gather(model,nq,nt,0.004,"model");
    system("sufft  < model | suamp | suxwigb title=\"d\" perc=99 & ");
  }

  free1float(d2r);
  free1float(m2r);
  free2float(LR);
  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  return;
}




void matrix(complex **L, float *h, float *q,int nh,int nq,float w, float dq, int rtmethod)
{

        register int ih;
	register int iq;  
        complex  arg;

        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.r=0;
              if ((rtmethod==1)||(rtmethod==3)) arg.i=w*h[ih]*q[iq];
              else if (rtmethod==2) arg.i=w*h[ih]*h[ih]*q[iq];
	      
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}













































