#include "radonwtcgls.h"

void wtcgls(complex *b,complex **L, complex *x,complex *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol, float step, int itercg);

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



void wtcgls0(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,iter;
  complex **m2, **d2, czero;
  complex **L;
  float w=0,wa,*dh,*dh2,df,*Jtot;
  complex *dc; 
  float  *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float power; 
  float *uaux;
  float **J;
  // Foster and Mosher offset function
  float depth=10; // Not used
  float *g=0;  
  float quantil1=eps1;
  float quantil2=eps2;
  float sigmam;
  float sigmad;   

  int freqflag=0;  

  czero.r=czero.i=0;
  
   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh); 
   dh=ealloc1float(nh);
   dh2=ealloc1float(nh);
   Jtot=ealloc1float(20);
   qcg=ealloc1complex(nq); 
   q1cg=ealloc1complex(nq);
   scg=ealloc1complex(nq);
   x1cg=ealloc1complex(nq);
   zcg=ealloc1complex(nq);
   z1cg=ealloc1complex(nq);
   rcg=ealloc1complex(nh);
   Azcg=ealloc1complex(nh);
   eta=ealloc1float(nq);
   rho=ealloc1float(nq);
   gcv=ealloc1float(nq);
   Wm=ealloc1float(nq);
   dc=ealloc1complex(nh);
   uaux=ealloc1float(nq);
   g=ealloc1float(nh);

   

   for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   
   dh[0]=pos[1]-pos[0];
   dh[nh-1]=pos[nh-1]-pos[nh-2];
   float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
   for (ih=0;ih<nh;ih++){
     Wd[ih]*=MIN(dh[ih],dhav);
     fprintf(stderr," Wd[%d]=%f\n",ih,Wd[ih]);
   }


   if (rtmethod==1) for (ih=0;ih<nh;ih++) g[ih]=pos[ih];
   else if (rtmethod==2) for (ih=0;ih<nh;ih++) g[ih]=pos[ih]*pos[ih];
   else if (rtmethod==3) for (ih=0;ih<nh;ih++) g[ih]=sqrt(pos[ih]*pos[ih]+depth*depth)-depth;   

     
   for (i=0;i<=iter_end;i++) Jtot[i]=0.;
   fftgo0(-1,data,d2,nh,nt,dt,&nf0);
   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f\n",maxfreq,dt,df);
   /***********PLot Cost***********/
   int ncurves=20;
   int countcurves=0;float J0;
   J=ealloc2float(iter_end,ncurves);
   /***********PLot Cost***********/
   for (freq=1;freq<maxfreq;freq++){
     w=2*PI*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     for (iq=0;iq<nq;iq++)  Wm[iq]=eps1;
     //for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
     radon_matrix(L,g,q,nh,nq,w);
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);
     Atimesx(d2[freq],L,m2[freq],nh,nq,1);     
     for (iter=1;iter<=iter_end;iter++){
     
       if (iter==2) 
	 deviations_td(m2[freq],nq,d2[freq],nh,norm,quantil1,quantil2,&sigmam,&sigmad);
       fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
       if (iter>1)  weights_inv(m2[freq],nq,norm,sigmam,Wm);    
      
       
       if (0)
	 if (iter>1) modelweight_inv(m2[freq],nq,norm,eps1,Wm);
	 else modelweight_inv(m2[freq],nq,3,eps1,Wm);

       
       /*
	 modelweight_inv(m,nq,inv.norm,inv.eps1,Wm);   
	 for (iq=0;iq<nq;iq++)  {
	 if ((iter!=1)){
	 float qmodel=quest(.7,nq,mtemp);
	 //qmodel=1;
	 normu=abs(m2[freq][iq])/qmodel+eps2;
	 Wm[iq].r=1/normu+eps1;
	 }
	 else if (0){	   	   
	 normu=abs(uaux[iq])+eps2;
	 Wm[iq].r=1/normu+eps1;
	 }
	 else{ 
	 Wm[iq].r=eps1;
	 Wm[iq].i=0.;
	 }
	 }
       */

       J0=wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,step,itercg);

       if (countcurves < ncurves ) J[countcurves][iter-1]=(J0); 
       

       //fprintf(stderr,"J=%f\n",J);
     }
     for (iq=0;iq<nq;iq++) uaux[iq]=abs(m2[freq][iq]);            
     /////////////////////////////////////////////////
     for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];          
     if (countcurves < ncurves ) normalize(iter_end,J[countcurves]);
     countcurves++;
     
   }

   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++)  for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");

   if (1) plotcurves(J[0],iter_end,ncurves,"Cost");    


   free2float(J);
   free1float(g);
   free1float(uaux);
   free1complex(dc);
   free1float(Wm);
   free1float(gcv);
   free1float(rho);
   free1float(eta);
   free1complex(Azcg);
   free1complex(rcg);
   free1complex(z1cg);
   free1complex(zcg);
   free1complex(x1cg);
   free1complex(scg);
   free1complex(q1cg);
   free1complex(qcg);
   free1float(Jtot);
   free1float(dh2);
   free1float(dh);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);

   return;
}
























