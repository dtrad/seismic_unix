#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
#include "Dcomplex.h"
dcomplex sin(dcomplex);
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void p_stack(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float fmax)
{
   int it, ih, iq, freq, maxfreq, nfreq, exit, nf0, i,j, tempint, iter;
   complex **m2, **d2,ctemp, czero;
   complex **l, **lh, **ll, *Qp, *ll0, **ll2;
   float w,*dh,df,*diagll,eps=1e-4,testcholres,*Jtot,Jmod,Jdata,powd, *diagll2;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep, *ftoep, *gtoep, noise;
   extern int nt, nh, nq, method, iter_end, rtmethod, norm, itercg, freqflag,
          costflag;             
   complex *ss, *ss2, *gg, *rr; // for ctoephcg
   dcomplex aa,pp;
   aa.r=PI;
   aa.i=1;
   pp=exp(aa);
   fprintf(stderr,"pp.r=%g,pp.i=%g\n",pp.r,pp.i);
 
   czero.r=czero.i=0;
   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   if ((d2=alloc2complex(nh,(nt+100)))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,(nt+100)))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((l=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((lh=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((ll=alloc2complex(nq,nq))==NULL)
     err("cannot allocate memory for ll\n");
   if ((ll0=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ll0\n"); 
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((maux=alloc1complex(nq))==NULL)
     err("cannot allocate memory for maux\n");
   if ((maux2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for maux2\n");
   if ((diagll=alloc1float(nq))==NULL)
     err("cannot allocate memory for diagll\n");
   if ((daux=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((daux2=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   if ((ftoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ftoep\n");
   if ((gtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtoep\n");
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");

   if ((method==5)||(method==6)){
   if ((ss=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ss\n");
   if ((ss2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ss2\n");
   if ((gg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gg\n");
   if ((rr=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rr\n");
   }

   if (method==7){
   if ((ll2=alloc2complex(nh,nh))==NULL)
     err("cannot allocate memory for ll2\n");
   //if ((ll02=alloc1complex(nh))==NULL)
   //err("cannot allocate memory for ll02\n");
   if ((diagll2=alloc1float(nh))==NULL)
     err("cannot allocate memory for diagll2\n");
   }


   for (ih=1;ih<nh;ih++)
       dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1]; 
  
   fprintf(stderr,"In p_stack1 nh=%d, nt=%d eps1=%f method=%d\n",
   nh,nt,eps1,method);
   fftgo(-1,d,d2,nh,nt,dt,&nf0);
   //fprintf(stderr,"In p_stack2 nh=%d, nt=%d nf0=%d\n",nh,nt,nf0);   
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

   for (freq=1;freq<maxfreq;freq++){
       w=2*PI*freq*df;
       for (ih=0;ih<nh;ih++)
          daux[ih]=d2[freq][ih];
 
       matrix_2(l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       
       if (method==1){
       fprintf(stderr,"Cholesky..............\n");
       AtimesBm(ll,lh,l,nq,nh,nq,'t');
       
       Atimesx(maux,lh,daux,nq,nh);
       for (i=0;i<nq;i++)
         ll[i][i]=ll[i][i]+(eps1/100)*ll[0][0];
       choldc(ll,nq,diagll);
       testcholres=testchol(ll,diagll,nq,nq);
       fprintf(stderr,"testcholres=%e\n",testcholres);
       cholsl(ll,nq,diagll,maux,maux2);
       }
       else if(method==2){
       fprintf(stderr,"Conjugate Gradients Before.........\n");
       conjgrad(l,lh,daux,maux2,eps,dh,dq,fmax,Jtot);
       fprintf(stderr,"Conjugate Gradients After........\n");
       }
       else if(method==3){
       Atimesx(maux,lh,daux,nq,nh);
       for (i=0;i<nq;i++){
          rtoep[i]=czero;
          for (j=0;j<nh;j++) 
              rtoep[i]+=lh[i][j]*l[j][0];
       }
       noise=((eps1/100)*rtoep[0].r);
       rtoep[0].i=0;
       fprintf(stderr,"rtoep[0]=%e, noise=%e\n",rtoep[0].r,noise.r);
       rtoep[0]+=noise; 
       tempint=ctoeplitz(nq,rtoep,maux2,maux,ftoep,gtoep);
       //fprintf(stderr,"Levinson, nstep=%d\n",tempint);
       }
       else if (method==4){
       fprintf(stderr,"Cholesky.HR.............\n");
       AtimesBm(ll,lh,l,nq,nh,nq,'t'); // 't' ==> only top triang is computed
       // add noise Covarianza
       for (powd=0,j=0;j<nh;j++) powd+=real(daux[j]*conjg(daux[j]));
       powd*=((eps1/100)/nh);
       powd=1;
       expfprint(powd);
       for (i=0;i<nq;i++) ll0[i]=ll[i][i]; 
       Atimesx(maux,lh,daux,nq,nh);
       //for (i=0;i<nq;i++) {maux[i].r=1e-5;maux[i].i=1e-5;} 
       iter=0;
       while (iter<iter_end){
         iter++;
	 if (iter==1) eps2=modgrad(maux,nq,norm,powd,Qp);
         if (iter!=1) eps2=modgrad(maux2,nq,norm,powd,Qp);
         
   	 for (i=0;i<nq;i++)
	   ll[i][i]=ll0[i]+Qp[i];

       
	 choldc(ll,nq,diagll);
	 cholsl(ll,nq,diagll,maux,maux2);
        
         if (costflag==1){
	   Atimesx(daux2,l,maux2,nh,nq);
	   Jmod=modnorm(maux2,eps2,1,nq);
	   Jdata=misfit(daux2,daux,powd,dh,nh);
           //displayA(daux2,nh);
           //displayA(daux,nh);
	   fprintf(stderr,"Iter %d, Jd=%e, Jm=%e, J=%e, eps2=%e\n",
		   iter,Jdata,Jmod,(Jdata+Jmod),eps2);
         }
       }
       } 
       else if(method==5){
       Atimesx(maux,lh,daux,nq,nh);
       for (i=0;i<nq;i++){
          rtoep[i]=czero;
          for (j=0;j<nh;j++) 
              rtoep[i]+=lh[i][j]*l[j][0];
       }
       noise=((eps1/100)*rtoep[0].r);
       rtoep[0].i=0;
       //fprintf(stderr,"rtoep[0]=%e, noise=%e\n",rtoep[0].r,noise.r);
       rtoep[0]+=noise;
       tempint=ctoephcg(itercg,nq,rtoep,maux2,maux,ss,ss2,gg,rr);
       fprintf(stderr,"Levinson, nstep=%d\n",tempint);
       }
       else if(method==6){
       fprintf(stderr,"Cholesky.HR.....Index=%d\n",freq);
       AtimesBm(ll,lh,l,nq,nh,nq);
       for (i=0;i<nq;i++) ll0[i]=ll[i][i];       
       Atimesx(maux,lh,daux,nq,nh);
       iter=0;
       while (iter<=iter_end){
         iter++;
	 if (iter==1) for (i=0;i<nq;i++)
	   Qp[i]=(eps1/100)*ll[0][0];         
         if (iter!=1) modgrad(maux2,nq,norm,1,Qp);
   	 for (i=0;i<nq;i++)
	   ll[i][i]=ll0[i]+Qp[i];
           tempint=cghrrt(itercg,nq,ll,maux2,maux,ss,ss2,gg,rr);
           fprintf(stderr,"CG, nstep=%d\n",tempint);
         }
       }
       else if (method==7){
       fprintf(stderr,"Cholesky.HR2.......Index =%d\n",freq);

       // add noise Covarianza
       //for (powd=0,j=0;j<nh;j++) powd+=real(daux[j]*conjg(daux[j]));
       //powd*=((eps1/100)/nh);
       powd=1;
       expfprint(powd);
       
       Atimesx(maux,lh,daux,nq,nh);
       //for (i=0;i<nq;i++) {maux[i].r=1e-5;maux[i].i=1e-5;} 
       iter=0;
       while (iter<iter_end){
         iter++;
	 if (iter==1) eps2=modgrad(maux,nq,norm,powd,Qp);
         if (iter!=1) eps2=modgrad(maux2,nq,norm,powd,Qp);
         AtimesBm(ll2,l,lh,nh,nq,nh,Qp); // 't' ==> only top triang is computed
   	 for (i=0;i<nh;i++)
	   ll2[i][i]=ll2[i][i]*(1+eps1/100)+1e-4;

       
	 choldc(ll2,nh,diagll2);
	 cholsl(ll2,nh,diagll2,daux,daux2);
 	 Atimesx(maux2,lh,daux2,nq,nh);
         xtimesy(maux2,maux2,Qp,nq);       
         if (costflag==1){
	   Atimesx(daux2,l,maux2,nh,nq);
	   Jmod=modnorm(maux2,eps2,1,nq);
	   Jdata=misfit(daux2,daux,powd,dh,nh);
           //displayA(daux2,nh);
           //displayA(daux,nh);
	   fprintf(stderr,"Iter %d, Jd=%e, Jm=%e, J=%e, eps2=%e\n",
		   iter,Jdata,Jmod,(Jdata+Jmod),eps2);
         }
       }
       } 
/////////////////////////////////////////////////
      for (iq=0;iq<nq;iq++)
          m2[freq][iq]=maux2[iq]; 
   }
   for (iq=0;iq<nq;iq++)     m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,m,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);
   
   if (method==7){
     free1float(diagll2);
     //free1complex(ll02);
     free2complex(ll2);
   }


   if ((method==5)||(method==6)){ 
     free1complex(rr);
     free1complex(gg);
     free1complex(ss2);
     free1complex(ss);
   }


free1complex(Qp);
free1complex(gtoep);
free1complex(ftoep);
free1complex(rtoep);
free1float(Jtot);
free1complex(daux2);
free1complex(daux);
free1float(diagll);
free1complex(maux2);
free1complex(maux);
free1float(dh);
fprintf(stderr,"After freemem  \n");
free1complex(ll0);
free2complex(ll);
free2complex(lh);
free2complex(l);
free2complex(m2);
free2complex(d2);

}







