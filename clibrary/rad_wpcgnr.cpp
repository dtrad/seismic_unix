#include "su.h"
#include "clibrary.h"
#include <math.h>
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

   IMPORTANT NOTE ABOUT CG SCALING
   Every multiplication by L means summation over nq
   and multiplication by LH means summation over nh
   Hence mult by L requires ==> /nq
   and   LH==> /nh
   Because the CG requires symmetry between L and LH then L/sqrt(nq*nh) 
   such that LH*L==> /(nh*nq)
	 	   

*/
void matrixcg(complex **l,float *pos,float *q,int nh,int nq,float w,float *dh, float *Wd, float dq, int rtmethod, float *g);


#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void rad_wpcgnr(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax, float *Wd)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2, **d2, czero;
  complex **L, *R;
  float w,wa,*dh,*dh2,df,*Jtot, Jdata, Jmod, normu;
  complex *d,*u;
  extern int nt,nh,nq,method,iter_end,rtmethod,norm,itercg,freqflag,costflag;
  complex *dtemp, *dc;
  float *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float extern step;  
  float power, Jtotlast, Jtotprev, bb, tempfloat; 
  float *uaux;
  extern int testadj;
  // Foster and Mosher offset function
  extern float depth;
  float *gFM;  

  czero.r=czero.i=0;
  eps1/=100;
   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   if ((d2=alloc2complex(nh,nt+100))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,nt+100))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((L=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n");
   if ((dh2=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh2\n");
   if ((u=alloc1complex(nq))==NULL)
     err("cannot allocate memory for u\n");
   if ((d=alloc1complex(nh))==NULL)
     err("cannot allocate memory for d\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((Wm=alloc1float(nq))==NULL)
     err("cannot allocate memory for Wm\n");
   if ((dc=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dc\n");
   if ((R=alloc1complex(nq))==NULL)
     err("cannot allocate memory for R\n");
   if ((uaux=alloc1float(nq))==NULL)
     err("cannot allocate memory for uaux\n");


   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];

   if (rtmethod==3){
     if ((gFM=alloc1float(nh))==NULL)
       err("cannot allocate memory for g\n");
     for (i=0;i<nh;i++) gFM[i]=sqrt(pos[i]*pos[i]+depth*depth)-depth;
   }
   //for (ih=0;ih<nh;ih++) dh2[ih]=sqrt(fabs(dh[ih]));

   for (i=0;i<=iter_end;i++) Jtot[i]=0.;
   fftgo(-1,data,d2,nh,nt,dt,&nf0);
   //fftback(1,data,d2,nh,nt,dt,nf0);
   //savesudata(data,nh,nt,dt,"data.su",1);
   //fftgo(-1,data,d2,nh,nt,dt,&nf0);

   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f\n",maxfreq,dt,df);
   const double  pi=acos(-1.);
   float test;
   int cgiter;
   void (*oper) (complex *,complex *, complex  **, int, int, int);
   oper=radon_freq;
   for (iq=0;iq<nq;iq++) Wm[iq]=1;
   for (iq=0;iq<nh;iq++) Wd[iq]=1;
   for (freq=1;freq<maxfreq;freq++){
     w=2*pi*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     //for (iq=0;iq<nq;iq++)  Wm[iq]=eps1;       
     for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
     matrixcg(L,pos,q,nh,nq,w,dh,Wd,dq,rtmethod,gFM);
     //power=rcdot(nh,d2[freq],d2[freq]);
     radon_freq(m2[freq],d2[freq],L,1,nh,nq);
     for (iter=1;iter<=iter_end;iter++){
       modelweight(m2[freq],nq,norm,eps1,Wm);
       cgiter=wpcgnr(oper,nh,nq,m2[freq],d2[freq],Wd,Wm,L,eps,step,itercg,1);
       fprintf(stderr,"freq=%d,it=%d\n",freq,cgiter);
     }
     /////////////////////////////////////////////////
     for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;      
     if ((wa<1)&&(wa>0))
	for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;

   }
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");

if (rtmethod==3) free1float(gFM);
free1float(uaux);
free1complex(R);
free1complex(dc);
free1float(Wm);
free1float(Jtot);
free1complex(d);
free1complex(u);
free1float(dh2);
free1float(dh);
free2complex(L);
free2complex(m2);
free2complex(d2);
return;
}


void matrixcg(complex **l,float *pos,float *q,int nh,int nq,float w,float *dh, 
	      float *Wd, float dq, int rtmethod, float *g)
{
        register int i;
	register int j;  
        complex  arg;
        //float Aperture=0.0;
        //for (i=0;i<nh;i++) Aperture+=dh[i]; 
        //float scale=1/sqrt(fabs(Aperture)*nq);
	/* Every multiplication by L means summation over nq
	   and multiplication by LH means summation over nh
	   Hence L ==> /nq
	   and   LH==> /nh
	   Because the adjoint requires symmetry L/sqrt(nq*nh) 
	   such that LH*L==> /(nh*nq)
	 */	   
        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
	      else if (rtmethod==3)         
		arg.i=w*g[i]*q[j];
	      
    	      //l[i][j]=scale*dh2[i]*exp(-arg);
	      //lh[j][i]=scale*dh2[i]*exp(arg);
	      l[i][j]=Wd[i]*exp(-arg)/sqrt(nq*nh);

	  }
	}
	
        return;
}

















































