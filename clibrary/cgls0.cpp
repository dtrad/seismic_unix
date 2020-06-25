#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void cgls0(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2, **d2, czero;
  complex **L, **LH, *R;
  float w,wa,*dh,df,*Jtot, Jdata, Jmod, *Cd;
  complex *d,*u;
  extern int nt,nh,nq,method,iter_end,rtmethod,norm,itercg,freqflag,costflag;
  complex *dtemp, *dc, *Qp;
  complex *pcg, *rcg, *rauxcg, *scg, *qcg;
  float power, Jtotlast, Jtotprev, bb, tempfloat;
   FILE *myfile3; 
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
   if ((pcg=alloc1complex(nq))==NULL) 
     err("cannot allocate memory for pcg\n");
   if ((rcg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rcg\n");
   if ((scg=alloc1complex(nh))==NULL)
     err("cannot allocate memory for scg\n");
   if ((qcg=alloc1complex(nh))==NULL)
     err("cannot allocate memory for qcg\n"); 
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((rauxcg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rauxcg\n");
   if ((dc=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dc\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Cd\n");
   if ((R=alloc1complex(nq))==NULL)
     err("cannot allocate memory for R\n");
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
       for (iq=0;iq<nq;iq++)  {Qp[iq].r=eps1;Qp[iq].i=0.0;}
       //Qp[5]=.1;Qp[16]=.1;
       //Qp[11]=eps1*1e-3;Qp[31]=eps1*1e-3;        
       for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
       matrix_3(R,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod);
       
       fprintf(stderr,"freq=%d,R[0].r=%e\n",freq,R[0].r);
       power=rcdot(nh,d,d); 
       for (iter=0;iter<iter_end;iter++){
	 if (iter==0) for (iq=0;iq<nq;iq++)  {Qp[iq].r=eps1;Qp[iq].i=0.0;}
	 if (iter!=0) for (iq=0;iq<nq;iq++){
	   Qp[i]=1./(abs(u[i])+1e-9)+eps1/100.;
	   Qp[iq].i=0.0;
	 }  
	 cgls(d,L,LH,u,Qp,nh,nq,qcg,rcg,rauxcg,scg,pcg,eps1,eps2,eps);
       }
       //displayA(u,nq);       
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
   if((myfile3=fopen("Jtottemp","w"))==NULL)
        err("cannot open file=%s\n","Jtottemp");    
   //fwrite(Jtot,sizeof(float),iter_end,myfile3);
   fclose(myfile3);
free1complex(R);
free1float(Cd);
free1complex(dc);
free1complex(rauxcg);
free1complex(Qp);
free1complex(qcg);
free1complex(scg);
free1complex(rcg);
free1complex(pcg);      
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























