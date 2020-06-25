#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void conjgrad3r(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag,nn;
  complex **m2, **d2, czero;
  complex **L, **LH, *R, **LL;
  float w,wa,*dh,df,*Jtot, Jdata, Jmod, *Cd;
  double **U,*RHS,*x;
  complex *d, *u, *adj, *dc, *Qp;
  extern int nt,nh,nq,method,iter_end,rtmethod,norm,itercg,freqflag,costflag;
  float Rorig;
  float power, resid, Jtotlast, Jtotprev, bb, tempfloat;
  double ress,tol;
  nn=2*nq;   	
  FILE *myfile3; 
  czero.r=czero.i=0;
  tol=(double) eps;
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
   if ((LL=alloc2complex(nq,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for LL\n");
   if ((R=alloc1complex(nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for R\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((u=alloc1complex(nq))==NULL)
     err("cannot allocate memory for u\n");
   if ((adj=alloc1complex(nq))==NULL)
     err("cannot allocate memory for adj\n");
   if ((d=alloc1complex(nh))==NULL)
     err("cannot allocate memory for d\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((dc=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dc\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Cd\n");
   U=dmatrix(1,nn,1,nn);
   RHS=dvector(1,nn);
   x=dvector(1,nn);

   for (ih=0;ih<nh;ih++) Cd[ih]=1.;
   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   //for (ih=1;ih<nh;ih++) dh[ih]=1.;
   dh[0]=dh[1]; 
   for (i=0;i<=iter_end;i++) Jtot[i]=0.;
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
       
       for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
       matrix_3(R,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod,Cd);

       //AtimesBm(ll,lh,l,nq,nh,nq);
       fprintf(stderr,"freq=%d\n",freq);
       power=rcdot(nh,d,d); 
       //power=1; 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Unconstrained gradient step

//-----Apply weights LH
//       dtemp=LH*d
       //for (i=0;i<nq;i++) R[i]=conjg(R[i]);
      Rorig=R[0].r;
      //R[0].r*=(1.+eps1/100.);
      Atimesx(adj,LH,d,nq,nh);
      bb=rcdot(nq,adj,adj);
      for (i=0;i<nq;i++) 
	 for (j=0;j<nq;j++)
             if (j>i) LL[i][j]=R[j-i];
	     else LL[i][j]=conjg(R[i-j]);
           
//-   Minimum for u=(eps,eps)
      //for(i=0;i<nq;i++) if(abs(u[i]) < sqrt(eps)) u[i].r=u[i].i=eps;
      iter=1; Jtotlast=0; Jtotprev=10; flag=0; //This flag is set to 1 when J increases
      while ((iter<=iter_end)&&(flag==0)){
         if (iter==1) eps2=modgrad(adj,nq,norm,1.,Qp);
         if (iter!=1) eps2=modgrad(u,nq,norm,1.,Qp);
         if (norm==1) {
         eps2=1.;
         for (i=0;i<nq;i++) Qp[i]=1./(abs(u[i])+1e-9)+eps;
         }
         //if (iter!=1) R[0].r=Rorig+1e-7;
         for (i=0;i<nq;i++)  LL[i][i]=R[0].r+Qp[i];
         for (i=0;i<nq;i++) 
	   for (j=0;j<nq;j++){ 
	     U[i+1][j+1]=LL[i][j].r;
	     U[nq+i+1][j+1]=LL[i][j].i;
	     U[i+1][nq+j+1]=-LL[i][j].i;
	     U[nq+i+1][nq+j+1]=LL[i][j].r;
           } 
         for (i=0;i<nq;i++){
	   RHS[i+1]=adj[i].r;
           RHS[nq+i+1]=adj[i].i;
         }
       
       	 linbcg1(U,x,RHS,nn,1,tol,itercg,&k,&ress);
         fprintf(stderr,"k=%d,ress=%e\n\n",k,ress); 
 
         for (i=0;i<nq;i++){
	   u[i].r=x[i+1];
	   u[i].i=x[nq+i+1];           
         } 

         if (costflag==1){
	   Atimesx(dc,L,u,nh,nq);
	   Jmod=modnorm(u,1.,dq,nq);
	   Jdata=misfit(dc,d,power,dh,nh);
           //fprintf(stderr,"Jdata=%f,Jmod=%f\n",Jdata,Jmod);
           Jtot[iter-1]+=log10(Jdata+Jmod);
         }
         iter++;
       }  //loop for niter        
      /////////////////////////////////////////////////
      
      if ((wa<1)&&(wa>0))
	for (iq=0;iq<nq;iq++)  u[iq]*=wa;
      for (iq=0;iq<nq;iq++)
          m2[freq][iq]=u[iq];        
   }
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);
   if((myfile3=fopen("Jtottemp","w"))==NULL)
        err("cannot open file=%s\n","Jtottemp");    
   fwrite(Jtot,sizeof(float),iter_end,myfile3);
   fclose(myfile3);
free_dvector(x,1,nn);
free_dvector(RHS,1,nn);
free_dmatrix(U,1,nn,1,nn);     
free1float(Cd);
free1complex(dc);
free1complex(Qp);
free1float(Jtot);
free1complex(d);
free1complex(adj);
free1complex(u);
free1float(dh);
free1complex(R);
free2complex(LL);
free2complex(LH);
free2complex(L);
free2complex(m2);
free2complex(d2);
return;
}























