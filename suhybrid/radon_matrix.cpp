#include "su.h"

void radon_matrix(complex *R, complex **l,complex **lh,float *g,float *q,int nh,int nq,float w,float *dh)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];

        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;

        for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      lh[iq][ih]=(1./Aperture)*dh[ih]*co;
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
	for (iq=0;iq<nq;iq++){
	  R[iq].r=0;
	  R[iq].i=0;
	  for (ih=0;ih<nh;ih++)
	    R[iq]+=lh[0][ih]*l[ih][iq]; //Top row of LL=LH*L
	  //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
	}

        return;
}

void radon_matrix(complex **l, float *g,float *q,int nh,int nq,float w)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        int ih, iq;
        complex  co;
	complex  dco;
	complex phase, dphase;
	float dq=q[1]-q[0];


	for (ih=0;ih<nh;ih++){
	  phase.r=dphase.r=0;

	  phase.i=(w*g[ih]*(q[0]-dq));
	  dphase.i=(w*g[ih]*dq);

	  co=exp(phase);
	  dco=exp(dphase);

	  for (iq=0;iq<nq;iq++){
	      co*=dco;
    	      l[ih][iq]=conjg(co);
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",iq,ih,lh[iq][ih].r,lh[iq][ih].i);
	  }
	}
  	      
        return;
}

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth)
{
  int ih;
  
  if (rtmethod==1) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih];
  else if (rtmethod==2) 
    for (ih=0;ih<nh;ih++) g[ih]=h[ih]*h[ih];
  else if (rtmethod==3) 
    for (ih=0;ih<nh;ih++) g[ih]=sqrt(h[ih]*h[ih]+depth*depth)-depth;   

  return;

}


void radon_matrix_irrq(complex **L, float *h, float *q,int nh,int nq,float w)
{
        
        register int ih;
	register int iq;  
        complex  arg;
	arg.r=0;
        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.i=w*h[ih]*q[iq];
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}


void radon_matrix_2op(complex **L, complex **L1, complex **L2, int nh, int nq1, int nq2)
{
  int ih, iq;

  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq1;iq++)
      L[ih][iq]=L1[ih][iq];
    for (iq=0;iq<nq2;iq++) 
      L[ih][iq+nq1]=L2[ih][iq];
  }	


  return;
}





