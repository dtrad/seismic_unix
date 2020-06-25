#include "radonsolver.h"

/* 
if GCV == TRUE then use GCV function,
otherwise uses residuals
*/

#define GCV 1

void gcvfunction(float *gcv, float *rho, float normdata, int nh, int nq, int iter);

float radon_cgfft(complex *d2, complex **L, complex *RC, complex *m2, int nh, int nq, int nf2, float *Wm, float *Cd, float eps, int itercg, float step)
{
  int k,iter;
  register int i;
  complex czero;  czero.r=czero.i=0;
  complex *fc, *gc;
  complex *u;
  complex *g1;    // residuals r0 
  complex *g2;    // Conjugate gradient (search direction)
  complex *Qp;    // Regularization
  complex *gtemp, *gtemp2, *dtemp;  // Temporal vectors required to compute Ap
  complex *Ag;  // Matrix times conjugate gradient
  complex *rhs; // Right hand side LH*d
  complex *xold;
  float *gcv;   // gcv function
  float *rho;   // normalized residuals
  float alphanum, alphaden, alpha, betanum, betaden, beta;
  float power, normrhs;
  float resold=0;
  float resid;
  

  // option=0 ------> circulant matrix-vector multplication
  // option=1 ------> Conventional matrix vector multiplication
  // option=2 ------> Both methods: prints out the numbers to check.


  // Note that c allocates memory per column
  // Hence 2 D array with alloc2type requires reversing dimension
  // i.e., first columns, then rows.

  g1=ealloc1complex(nq);
  g2=ealloc1complex(nq);
  u=ealloc1complex(nq);
  Qp=ealloc1complex(nq);
  gtemp=ealloc1complex(nq);
  gtemp2=ealloc1complex(nq);
  dtemp=ealloc1complex(nh);
  rhs=ealloc1complex(nq);
  Ag=ealloc1complex(nq);
  fc=ealloc1complex(nf2);
  gc=ealloc1complex(nf2);
  gcv=ealloc1float(itercg);
  rho=ealloc1float(itercg);
  xold=ealloc1complex(nq);

  power=rcdot(nh,d2,d2);
  
  //power=1;

  
  xtimesy(dtemp,Cd,d2,nh);
  Atimesx(dtemp,L,rhs,nh,nq,TRUE);

  ///////////////////////////////////////////////////////////////////////
  // Compute the model preconditioner
  // Setting power to larger values increase the non linearity
  // using smalll values make more least squares results. It is equivalent
  // to the Lagrange multiplier.      

  for (i=0;i<nq;i++){
    // Qp[i].r=sqrt(power)*Wm[i];
    Qp[i].r=1./Wm[i];
    Qp[i].i=0;
  }

  normalize(nq,Qp);

  if (0){
    max(nq,Qp);
    min(nq,Qp);
  }
  ///////////////////////////////////////////////////////////////////////

  xtimesy(rhs,Qp,rhs,nq);

  normrhs=rcdot(nq,rhs,rhs);
  iter=1;    

  for (i=0;i<nq;i++) u[i]=czero;
    
  xequaly(g1,rhs,nq);  // resid (gradient) g1 = rhs  
  xequaly(g2,rhs,nq);  // search direction g2 = rhs
 
  resid=rcdot(nq,g1,g1);
  normrhs=resid;
  //plotcurves(Wm,nq,1,"Wm");
  k=0;resold=2*resid;
  while ( k < itercg ){
    k++;
    
    ///// Compute alpha
    alphanum=rcdot(nq,g1,g1);
    xtimesy(gtemp,Qp,g2,nq);  // Wm^{-1} m
    circ_mult(nq,RC,gtemp,gtemp2,nf2,fc,gc); // A^T A Wm^{-1} m
    xtimesy(Ag,Qp,gtemp2,nq); // Wm^{-T} A^T A Wm^{-1} m
    //xtimesy(gtemp2,Qp,g2,nq);
    //xplusy(Ag,gtemp,gtemp2,nq);
    
    alphaden=rcdot(nq,g2,Ag);

    if (0)
      for (i=0;i<nq;i++) 
	fprintf(stderr,"gtemp[%d]=%f,gtemp2[%d]=%f\n",i,abs(gtemp[i]),i,abs(gtemp2[i]));

    //if (alphaden<eps) break;
    alpha=step*alphanum/alphaden;
    //fprintf(stderr,"alpha=%f\n",alpha);

    for(i=0;i<nq;i++){
      xold[i]=u[i];
      u[i]=u[i]+alpha*g2[i];
      g1[i]=g1[i]-alpha*Ag[i];  
    }
    
    ///// Compute Beta
    betanum=rcdot(nq,g1,g1);
    betaden=alphanum;              
    //if (betaden<eps) break;
    beta=betanum/betaden;	      
    //fprintf(stderr,"beta=%f\n",beta);
    for(i=0;i<nq;i++) g2[i]=g1[i]+beta*g2[i];
	  
    resold=resid;
    resid=betanum;

    // Stop if residuals increases and GCV is not used
    if (!GCV) if (1.5*resold < resid) break;

    // GCV function
    rho[k-1] = sqrt(resid)/normrhs;
    if (GCV){
      gcvfunction(gcv,rho,normrhs,nh,nq,k);
      if (k>1){
	if ( gcv[k-2] < gcv[k-1] ) { 
	  if (1) fprintf(stderr,"GCV Criteria, iteration %d\n",k);
	  for(i=0;i<nq;i++) u[i]=xold[i]; 
	  break;
	}
      }
    }

  } 

  if (1) fprintf(stderr,"k=%d\n",k);

  for(i=0;i<nq;i++) m2[i]=Qp[i]*u[i];

  free1complex(xold);
  free1float(rho);
  free1float(gcv);
  free1complex(gc);
  free1complex(fc);
  free1complex(Ag);
  free1complex(gtemp2);
  free1complex(rhs);
  free1complex(dtemp);
  free1complex(gtemp);
  free1complex(Qp);
  free1complex(g2);
  free1complex(u);
  free1complex(g1);
  
  // circ_mult


  /////////
  
  return(resid);
}

void radon_matrix_cgfft(complex **L, complex *RC, int nh, int nq, int nf2, float *Cd)

{
  int ih,iq;
  register int i;
  complex czero;czero.r=czero.i=0;
  complex *R;

  R=ealloc1complex(nq);
  for (iq=0;iq<nq;iq++){
    R[iq]=czero;
    //for (ih=0;ih<nh;ih++) R[iq]+=conjg(L[ih][0])*L[ih][iq]/nh; //Top row of LL=LH*L
    for (ih=0;ih<nh;ih++) R[iq]+=conjg(L[ih][0])*Cd[ih]*L[ih][iq]/nh; //Top row of LL=LH*L
    //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
  }

  // Construct the circulant matrix//////////////////////////////////////
  for (i=0; i<nq;i++) RC[i]=conjg(R[i]);       
  for (i=nq; i<nf2;i++) RC[i]=czero; // pad with zeros the autocorrelation
  for (i=1;i<=nq-1;i++) RC[nf2-i]=conjg(RC[i]);      //use the DFT
  //fprintf(stderr,"rtoep[0].r=%f, rtoep[0].i=%f \n",RC[0].r,RC[0].i);
  pfacc(-1,nf2,RC);
  ///////////////////////////////////////////////////////////////////////
  free1complex(R);
  
  return;
}

void  circ_mult(int n,complex *RC,complex *f,complex *g,int nf,
		complex *fc, complex *gc)
{

  /*
     Perform the multiplication of a Toeplitz form with a vector.
     The procedure is done using a circulant matrix derived from
     the first row of the Toeplitz form. The Toeplitz form, in fact,
     is not needed. Everything you need is the firt row of the matrix.

     INPUT  RC(n): first row of the autocorrelation matrix (Toeplitz)
            f(n): a filter
              nf: number of freq. samples to perform the circular conv.
            fc(nf): a pointer to an already allocated vector
            gc(nf): a pointer to an already allocated vector

     OUTPUT g(n): the product of the toeplitz form with f

     NOTE   that T(1) must be real in order to  have an Hermitian 
           form T

     Mauricio Sachi 

     adapted for SU and translated to C: Daniel Trad
  */
      
  complex czero;
  int i;
  czero.r=czero.i=0.;

  for (i=0;i<n;i++) fc[i]=f[i];  
  for (i=n;i<nf;i++) fc[i]=czero;    // pad with zeros the autocorrelation
 
  pfacc(-1,nf,fc);
  
  for (i=0;i<nf;i++) gc[i]=RC[i]*fc[i];
  
  pfacc (1,nf,gc);
  
  for (i=0;i<n;i++) g[i]=gc[i];

  for (i=0;i<n;i++) g[i]/=nf;            

  return;
}

/*
      ///////////////////////////////////////////////////////////////////
      // The following options are only to test multiplication with FFT
      // In general option=0
      //////////////////////////////////////////////////////////////////
      if (option==2){
	circ_mult(nq,RC,u,gtemp,nf2,fc,gc);
	for (i=0;i<10;i++){
	  fprintf(stderr,"Circ_matrix gtemp[%d].r=%e\n",i,gtemp[i].r);
	  fprintf(stderr,"Circ_matrix gtemp[%d].i=%e\n",i,gtemp[i].i);
	}
      }

      if (option==1){
	AtimesBm(LL,LH,L,nq,nh,nq);
	Atimesx(gtemp,LL,u,nq,nq);
      }
     
      if (option==2){
	AtimesBm(LL,LH,L,nq,nh,nq);
	Atimesx(gtemp,LL,u,nq,nq);
	for (i=0;i<10;i++){
	  fprintf(stderr,"LL*U gtemp[%d].r=%e\n",i,gtemp[i].r);
	  fprintf(stderr,"LL*U gtemp[%d].i=%e\n",i,gtemp[i].i);
	} 
      }
      //////////////////////////////////////////////////////////////////////

*/




void gcvfunction(float *gcv, float *rho, float normdata, int nh, int nq, int iter)
{
  int num;
  int i=iter-1;
  
  num=(nq-i)*(nq-i);
  gcv[i]=(rho[i]*rho[i])/num;
  if (0) fprintf(stderr,"gcv[%d]=%e\n",i,gcv[i]);

  return;

}
 




























