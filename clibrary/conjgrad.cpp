#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
void conjgrad(complex **L, complex **LH, complex *d, complex *u,
float eps,float *dh, float dq, float freq, float *Jtot)
{   
  int i, k, iter, laststep, flag;
  complex *uold;
  complex *g1, *g1old, *g2,  *Qp;
  complex *gtemp, *gtemp2, *gtemp3, *gtemp4, *gtemp5;
  complex *dtemp, *dc;
  float alfanum, alfaden, alfa;
  float betanum, betaden, beta;
  extern int  nh, nq, iter_end, norm;
  float bb;
  float power, resid;
  float pmin, eps2, residold, Jtotlast, Jtotprev;
  float *power2;

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
   if ((power2=alloc1float(nq))==NULL)
     err("cannot allocate memory for dtemp\n");

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Unconstrained gradient step

//-----Apply weights LH
      // AtimesDiag(LH,LH,dh,nq,nh);
      
//       dtemp=LH*d
      Atimesx(gtemp2,LH,d,nq,nh);  
      xequaly(u,gtemp2,nq);
      
//-   Minimum for u=(eps,eps)
      bb=0;
      for(i=0;i<nq;i++){
            if(abs(u[i]) < eps){ 
	      u[i].r=u[i].i=eps;
	      power2[i]=real(conjg(u[i])*u[i]);
            }
            bb=bb+power2[i];
       }
      
//-       Sigma=min(u[i]**2)

      pmin=power2[1];
      for(i=0;i<nq;i++)
              if (power2[i]<pmin) 
                     pmin=power2[i];
     
      pmin=pmin*1e-3;
      eps2=(pmin > 1.e-6)?pmin:1e-6; 

//---     gtemp=LH*L*u
      Atimesx(dtemp,L,u,nh,nq);
      Atimesx(gtemp,LH,dtemp,nq,nh);      

//--      g1=LH*d-LH*L*u
      xminusy(g1,gtemp2,gtemp,nq);

//-Estimate Qp

      iter=1;
      Jtotlast=0;
      Jtotprev=10;
      flag=0; //This flag is set to 1 when J increases
      
      while ((iter<iter_end)&&(flag==0)){
         fprintf(stderr,"iter=%d\n",iter);
	 modgrad(u,nq,norm,1,Qp);
         //displayA(Qp,nq);
         //fprintf(stderr,"maxQp=%f\n",(maxmax(Qp,nq)));
         xequaly(uold,u,nq);
         
//- Gradient steps
//- NQ loop
         resid=rcdot(nq,g1,g1);
         k=0;
         power=1;
         laststep=10;
         if (iter==iter_end) laststep=30;

         while ((k!=nq)&&(sqrt(resid)>(eps*bb))&&(k<laststep)){
            k++;
            if (k==1){ 
               xequaly(g2,g1,nq);
               betanum=rcdot(nq,g1,g1);
               betaden=rcdot(nq,g1old,g1old);
               fprintf(stderr,"betaden=%f\n",betaden);
               if (betaden<eps) betaden=eps;
               beta=betanum/betaden;
	       fprintf(stderr,"beta=%f\n",beta);
               for(i=0;i<nq;i++)
                  g2[i]=g1[i]+beta*g2[i];
               
            }
            alfanum=rcdot(nq,g1,g1);

            Atimesx(dtemp,L,g2,nh,nq); 
            Atimesx(gtemp,LH,dtemp,nq,nh);
            xtimesy(gtemp4,Qp,g2,nq);
            xplusy(gtemp5,gtemp,gtemp4,nq);
            
            alfaden=rcdot(nq,g2,gtemp5);
            
            if (alfaden<eps) alfaden=eps; 
                               
            alfa=alfanum/alfaden;
            fprintf(stderr,"alfa=%f, k=%d\n",alfa,k);
            fprintf(stderr,"alfa=%f\n",alfa);
            fprintf(stderr,"alfanum=%f, k=%d\n",alfanum,k);
            fprintf(stderr,"alfaden=%f, k=%d\n",alfaden,k);

            for(i=0;i<nq;i++){
                u[i]=u[i]+alfa*g2[i];
                g1old[i]=g1[i];
                g1[i]=g1[i]-alfa*gtemp5[i];
	    }
            residold=resid;
            resid=sqrt(alfanum);
            power=rcdot(nq,u,u);            
         } 
             //!loop for k
         Atimesx(dc,L,u,nh,nq);
         Jtot[iter]=costfunc(d,dc,u,nh,nq,eps2,norm);
         fprintf(stderr,"Jiter=%f,%d,%d,%f\n",Jtot[iter],iter,k,freq);
         Jtotlast=Jtot[iter];
         if (iter==1) 
                 Jtotprev=Jtot[iter]+10;
         else
                 Jtotprev=Jtot[iter-1];
         

         for(i=iter;i<iter_end;i++){
	   Jtot[i]=Jtot[iter] ;//In case of stop iterating keep the last J
	   /*         
c         if ((Jtotlast.gt.(1.1*Jtotprev)).and.(iter.ne.1)) then   
c            for(i=0;i<nq;i++)
c               u[i]=0.3*u[i]+0.7*uold[i]
c            enddo   
c            fprintf(stderr,"  \n", ); 'cost function increases',freq,iter
c         endif
	   */
	   if ((Jtotlast>(1.0*Jtotprev))&&(iter!=1)){
             flag=1;
             for(i=0;i<nq;i++)
               u[i]=uold[i];
      	   }
         }
         iter++;
         fprintf(stderr,"newiter=%d\n",iter);
       }  //loop for niter 
      free1complex(dtemp);
      free1complex(dc);
      free1float(power2);
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
return;
}



















