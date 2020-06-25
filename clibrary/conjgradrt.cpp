#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
#include "Dcomplex.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void conjgradrt(float *pos, float **data, float dt,float **model, float *q, float dq,
	     float eps1, float fmax)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2, **d2, czero;
  complex **L, **LH;
  float w,wa,*dh,df,eps=1e-4,*Jtot;
  complex *d,*u;
  extern int nt,nh,nq,method,iter_end,rtmethod,norm,itercg,freqflag,costflag;             
  complex *uold, *dtemp, *dc;;
  complex *g1, *g1old, *g2,  *Qp;
  complex *gtemp, *gtemp2, *gtemp3, *gtemp4, *gtemp5;
  float alfanum, alfaden, alfa, betanum, betaden, beta;
  float power, resid, eps2, residold, Jtotlast, Jtotprev, bb;
  czero.r=czero.i=0;

   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   if ((d2=alloc2complex(nh,(nt+100)))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,(nt+100)))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((L=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((LH=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((u=alloc1complex(nq))==NULL)
     err("cannot allocate memory for u\n");
   if ((d=alloc1complex(nh))==NULL)
     err("cannot allocate memory for d\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((uold=alloc1complex(nq))==NULL) 
     err("cannot allocate memory for uold\n");
   if ((g1=alloc1complex(nq))==NULL)
     err("cannot allocate memory for g1\n");
   if ((g1old=alloc1complex(nq))==NULL)
     err("cannot allocate memory for g1old\n");
   if ((g2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for g2\n"); 
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((gtemp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtemp\n");
   if ((gtemp2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtemp2\n");
   if ((gtemp3=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtemp3\n");
   if ((gtemp4=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtemp4\n");
   if ((gtemp5=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtemp5\n");
   if ((dc=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dc\n");
   if ((dtemp=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dtemp\n");


   for (ih=1;ih<nh;ih++)
       dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1]; 
  
   fprintf(stderr,"In conjgradrt nh=%d, nt=%d eps1=%f\n",nh,nt,eps1);
   fftgo(-1,data,d2,nh,nt,dt,&nf0);
 
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
   const double  pi=acos(-1.);
    
   for (freq=1;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       
       for (ih=0;ih<nh;ih++)
          d[ih]=d2[freq][ih];
         
       matrix_2(L,LH,pos,q,nh,nq,w,dh,dq,rtmethod);
       fprintf(stderr,"freq=%d\n",freq);
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Unconstrained gradient step

//-----Apply weights LH
      // AtimesDiag(LH,LH,dh,nq,nh);
//       dtemp=LH*d
      Atimesx(gtemp2,LH,d,nq,nh);
      bb=rcdot(nq,gtemp2,gtemp2);
      xequaly(u,gtemp2,nq);  // Initial guess
       
//-   Minimum for u=(eps,eps)
      for(i=0;i<nq;i++) if(abs(u[i]) < sqrt(eps)) u[i].r=u[i].i=eps;
//---     gtemp=LH*L*u
      Atimesx(dtemp,L,u,nh,nq);
      Atimesx(gtemp,LH,dtemp,nq,nh);      
//--      g1=LH*d-LH*L*u   
      xminusy(g1,gtemp2,gtemp,nq);

      iter=1; Jtotlast=0; Jtotprev=10; flag=0; //This flag is set to 1 when J increases
    
      while ((iter<iter_end)&&(flag==0)){

	 eps2=modgrad(u,nq,norm,1,Qp);
         //fprintf(stderr,"iter=%d,eps2=%f\n",iter,eps2);
         xequaly(uold,u,nq);
         
//- Gradient steps
//- NQ loop

         resid=rcdot(nq,g1,g1);
         k=0; power=1; 

         while ((k!=nq)&&(sqrt(resid)>(eps*bb))&&(k<itercg)){
            k++;
	    ///// Compute Beta
            if (k==1) xequaly(g2,g1,nq);
            else{ 
	      betanum=rcdot(nq,g1,g1);
	      betaden=rcdot(nq,g1old,g1old);
              //fprintf(stderr,"betaden=%f\n",betaden);
              if (betaden<eps) betaden=eps;
              beta=betanum/betaden;
	      //fprintf(stderr,"beta=%f\n",beta);
              for(i=0;i<nq;i++) g2[i]=g1[i]+beta*g2[i];
            }
            ///// Compute alfa
            alfanum=rcdot(nq,g1,g1);
            Atimesx(dtemp,L,g2,nh,nq); 
            Atimesx(gtemp,LH,dtemp,nq,nh);
            xtimesy(gtemp4,Qp,g2,nq);
            xplusy(gtemp5,gtemp,gtemp4,nq);
            alfaden=rcdot(nq,g2,gtemp5);
            if (alfaden<eps) alfaden=eps; 
            alfa=alfanum/alfaden;
            //fprintf(stderr,"alfa=%f, k=%d\n",alfa,k);
            //fprintf(stderr,"alfanum=%f, k=%d\n",alfanum,k);
            //fprintf(stderr,"alfaden=%f, k=%d\n",alfaden,k);
            //// Update model u and residuals
            for(i=0;i<nq;i++){
                u[i]=u[i]+alfa*g2[i];
                g1old[i]=g1[i];  // Kip previous residual
                g1[i]=g1[i]-alfa*gtemp5[i];  
	    }
            residold=resid;
            resid=alfanum;
            power=rcdot(nq,u,u);            
         } 
         /*    //!loop for k
         Atimesx(dc,L,u,nh,nq);
         Jtot[iter]=costfunc(d,dc,u,nh,nq,eps2,norm);
         fprintf(stderr,"Jiter=%f,%d,%d,%f\n",Jtot[iter],iter,k,freq);
         Jtotlast=Jtot[iter];
         if (iter==1) Jtotprev=Jtot[iter]+10; else Jtotprev=Jtot[iter-1];
         for(i=iter;i<iter_end;i++){
	   Jtot[i]=Jtot[iter] ;//In case of stop iterating keep the last J

	   if ((Jtotlast>(1.0*Jtotprev))&&(iter!=1)){
             flag=1;
             for(i=0;i<nq;i++)
               u[i]=uold[i];
      	   }
	   }*/
         iter++;
       }  //loop for niter        
      /////////////////////////////////////////////////

      if ((wa<1)&&(wa>0))
	for (iq=0;iq<nq;iq++)  u[iq]*=wa;
      for (iq=0;iq<nq;iq++)
          m2[freq][iq]=u[iq];

        
   }
   for (iq=0;iq<nq;iq++)     m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);

free1complex(dtemp);
free1complex(dc);
free1complex(gtemp5);
free1complex(gtemp4);
free1complex(gtemp3);
free1complex(gtemp2);
free1complex(gtemp);
free1complex(Qp);
free1complex(g2);
free1complex(g1old);
free1complex(g1);
free1complex(uold);      
   
free1float(Jtot);
free1complex(d);
free1complex(u);
free1float(dh);
free2complex(LH);
free2complex(L);
free2complex(m2);
free2complex(d2);
return;
}























