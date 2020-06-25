#include "radonlinetfd.h"

void radoncgfft_tfd_conv(float *pos, int nh, float **data, float *t, int nt, float dt,
		    float **model, float *q, int nq, float dq,  float fmax, float **Wd, 
		    int testadj, float quantil, int rtmethod, float depth, inv_par inv, 
		    float *wavelet, int nw)
{
  int verbose=0;
  register int it;
  int  j, iter, ih, iq;
  float **Wm; // Model weights
  //float **Wd;// Data weights
  int ntaper=5; // Number of lateral traces to taper 
  int nx=nt*nq;
  int ny=nh*nt;
  int freq;
  float *g=0; /* function of offset for moveot */
  complex ***L;
  complex **RC;
  complex **d2;
  complex **m2;
  float *Cd;
  int adj;
  int nfft;
  int nf;
  float df;
  int maxfreq;
  //float fmax=1./(2*dt)*0.7;
  int filtout=0;
  int taperflag=0;
  double pi=PI;
  complex arg;
  float w;
  int mute=0;
  int parmute=0;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;
  float sigmam;
  float sigmad;
  float *J;
  int nf2=npfa(2*nq);  
  float **data2;
  //////////////////////////////////////////////////////////////////////////
  // Set pointer to function for Conjugate Gradients
  void (*radon) (float *model, float *data,  complex ***L,int nt, int nh, int nq, 
		 float dt, complex **m2, complex **d2, int nfft, float fmax, int adj);
  
  // Set Fourier transform parameters 
  fft_parameters(nt,dt,&nfft,&nf,&df);

  maxfreq=(int) (fmax/df+0.5);

  if (verbose){
    fprintf(stderr,"radon:nfft=%d, nf=%d df=%f, maxfreq=%d\n",nfft,nf,df,maxfreq);   
    fprintf(stderr,"fmax=%f,rtmethod=%d,norm=%d\n",fmax,rtmethod,inv.norm);   
  }

  radon=radonf;

  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc3complex(nq,nh,nf);
  RC=ealloc2complex(nf2,nf);
  g=ealloc1float(nh);
  Wm=ealloc2float(nt,nq);
  data2=ealloc2float(nt,nh);
  J=ealloc1float(inv.iter_end+1);
  Cd=ealloc1float(nh);
  zero_array(d2,nt,nh);
  zero_array(m2,nt,nq);

  //Wd=ealloc2float(nt,nh);
  
  //if (taperflag==1) taper(data,nt,nh,ntaper,0);

  for (iq=0;iq<nq;iq++) memset((void *)model[iq],(int)'\0',nt*FSIZE);
  fprintf(stderr,"++++++++++++++++++1\n");
  // This is a filter for outliers or dominant bad data
  if (filtout){
    float qup=quest(0.99,nh*nt,data[0]);
    float qmean=quest(0.50,nh*nt,data[0]);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++) 
	if (fabs(data[ih][it])>qup) Wd[ih][it]=qmean/qup;
	else Wd[ih][it]=1.0;
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1.0;
  for (ih=0;ih<nh;ih++) Cd[ih]=1;
  radon_moveout(pos,g,nh,rtmethod,depth);  

  // Construct the Radon operator in Frequency
  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    radon_matrix(L[freq],g,q,nh,nq,w);
    radon_matrix_cgfft(L[freq],RC[freq],nh,nq,nf2,Cd);
  }

  testadjop(radon,L,nt,nh,nq,dt,m2,d2,nfft,fmax,adj);

  // Compute the adjoint 
  conv_corr(data[0],data2[0],nt,nh,wavelet,nw,TRUE);
  radon(model[0],data2[0],L,nt,nh,nq,dt,m2,d2,nfft,fmax,1);
  if (0) plotgather(model,nq,nt,dt,"suxwigb perc=99");
  fprintf(stderr,"quantil1=%f, quantil2=%f\n",quantil1,quantil2);

  for (j=1;j<=inv.iter_end;j++){
    if (0)  if (j==inv.iter_end) inv.taperflag=2; // change stop criteria at the last iteration
    if (j==2) deviations(model[0],nq*nt,data[0],nh*nt,inv.norm,quantil1,quantil2,
			 &sigmam,&sigmad);
    fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
    weights_cgfft(model[0],nt*nq,inv.norm,sigmam,Wm[0],j);    
    J[j-1]=cgfft(data[0],L,model[0],d2,RC,m2,nt,nh,nq,nf2,dt,nfft,fmax,Wm[0],Wd[0],inv,wavelet,nw);
  }
  for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) model[iq][it]/=nh;
  
  if (0) plotgather(model,nq,nt,dt,"suxwigb perc=99");
  if (0) plotcurves(J,inv.iter_end,1,"Cost");    
  if (mute){
    int nmute;
    iq=0; while(fabs(q[iq])<parmute) iq++; 
    //nmute=nq-iq;
    fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
    //taper(m,nt,nq,nmute,2);
    taper(model,nt,nq,iq,1);
  }
  
  //radon(model[0],data[0],L,nt,nh,nq,dt,d2,m2,nfft,fmax,0);
  free1float(J);
  free1float(Cd);
  free2float(Wm);
  free2float(data2);
  free1float(g);
  free2complex(RC);
  free3complex(L);
  free2complex(d2);
  free2complex(m2);

  return;
  
}

void radonf(float *model, float *data, complex ***L,int nt, int nh, int nq, 
	    float dt,complex **m2, complex **d2, int nfft, float fmax, int adj)
{
  int freq, it, iq, ih; // indexes
  int maxfreq;
  float w;
  float wa;
  complex czero;  czero.r=czero.i=0;  
  int nf=nfft/2+1;
  float df=1./(nfft*dt);

  maxfreq=(int) (fmax/df+0.5);

  if (adj){
    fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);

    for (freq=1;freq<maxfreq;freq++){
      w=2*PI*freq*df;
      wa=freqweight(freq,df,fmax-10,fmax);
      Atimesx(d2[freq],L[freq],m2[freq],nh,nq,adj);
      if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) m2[freq][iq]*=wa;
    }

    for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nf;freq++){
      for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
    }

    fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);

  }
  else{
    fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);
    for (freq=1;freq<maxfreq;freq++){
      w=2*PI*freq*df;
      wa=freqweight(freq,df,fmax-10,fmax);
      Atimesx(d2[freq],L[freq],m2[freq],nh,nq,adj);
      if ((wa<1)&&(wa>0)) for (ih=0;ih<nh;ih++) m2[freq][ih]*=wa;
    }
    
    for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nf;freq++){
      for (ih=0;ih<nh;ih++) d2[freq][ih]=czero;
    }
    
    fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf);
  }

  return;
}




//////////////////////////////////////////////////////////////////////////////////////
/* 
     This function solves the system of equations 
     (WmT FH WdT M^{-1} Wd F Wm) m =  FH WdT Wd M^{-1} d
     oper is an operator implemented with the function
     where M is a  preconditioner acting on data space.
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);
     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )
     Wd allows to downweight bad data by using large values. 
     In general is a diagonal size(ndata) 
     M is the preconditioner. Diagonal size(nmodel) Large values
     penalize, small values focus the solution to a desired model value.
     M changes the null space, W does not. 
     Hence prior information about the model must be implemented by M. 
     
     Taken from Yousef Saad, pag 260: 
     Iterative methods for sparse linear systems
     W has bee added to the original algorihm 
     Daniel Trad - UBC March-2000
    
  */

/* 
if gcvbool == 1 then use GCV function,
otherwise uses residuals
*/

float cgfft(float *b,  complex ***L, float *x, complex **d2, complex **RC, complex **m2, 
	    int nt, int nh, int nq, int nf2, float dt,  int nfft, float fmax, 
	    float *Wm, float *Wd, inv_par inv, float *wavelet, int nw)
 
{
  float normb, beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function and square residuals

  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // s  Search direction
  // w  Conjugate search direction
  // Wm model weights
  // Wd data weights
  // rc Weighted residuals to compute cost function
  // xc Weighted solution to compute cost function

  float *r,*g,*s,*stemp,*stemp2,*stemp3,*w,*rho,*eta,*gcv,*rc, *xc;
  complex **m2aux;
  int nx=nt*nq;
  int ny=nt*nh;  
  int stopc=inv.taperflag;

  g=ealloc1float(nx);  // Gradient vector (residuals in model space)
  s=ealloc1float(nx);  // Search direction
  stemp=ealloc1float(nx); // Temporal vector for the search direction
  stemp2=ealloc1float(nx); // Temporal vector for the search direction
  stemp3=ealloc1float(nx);  // Temporal vector for the search direction
  r=ealloc1float(ny);  // Residuals in data space
  w=ealloc1float(nx);
  rho=ealloc1float(inv.itercg+1);
  eta=ealloc1float(inv.itercg+1);
  gcv=ealloc1float(inv.itercg+1);
  xc=ealloc1float(nx); 
  rc=ealloc1float(ny);
  m2aux=ealloc2complex(nq,nt);


  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);
  int restart=1;



  // If restart (initial model = 0) then compute r = Wd (b - A x) = Wd b 
  if (restart){
    for (i=0;i<nx;i++) x[i]=0.;
    for (i=0;i<ny;i++) r[i]=Wd[i]*b[i];
  }
  else{
    radonf(x,r,L,nt,nh,nq,dt,m2,d2,nfft,fmax,0);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }
  // Compute Wmi L Wd b (rhs of the LS system)   
  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  radonf(g,r,L,nt,nh,nq,dt,m2,d2,nfft,fmax,1);
  conv_corr(g,stemp,nt,nq,wavelet,nw,FALSE);
  for(i=0;i<nx;i++) g[i]=stemp[i]*Wm[i];


  normb=dot(nx,g,g);

  for(i=0;i<nx;i++) s[i]=g[i];
  k=0;
  while(k<inv.itercg){
    // Calculate alpha num
    alphanum=dot(nx,g,g);
    // Calculate alpha den = ( L^T L s, s )
    xtimesy(stemp,Wm,s,nx);
    if (0) plotgather_pipe(stemp,nq,nt,"stemp");
    conv_corr(stemp,stemp3,nt,nq,wavelet,nw,TRUE);    
    Atimesx_by_fft_loop(stemp2,RC,stemp3,m2,m2aux,nt,nq,dt,nfft,fmax,nf2);
    conv_corr(stemp2,stemp3,nt,nq,wavelet,nw,FALSE);    
    if (0) plotgather_pipe(stemp3,nq,nt,"stemp2");    
    xtimesy(w,Wm,stemp3,nx);
    alphaden=dot(nx,w,s);

    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=inv.step*alphanum/alphaden;

    for(i=0;i<nx;i++) xc[i]=x[i]; // save last solution
    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<nx;i++) g[i]-=alpha*w[i];  
    // Calculate beta
    // beta numerator = alpha numerator = (g^T,g) with g_old
    // beta denominator = (g^T,g) with g new
    betanum=dot(nx,g,g);
    beta=betanum/alphanum;
    // New search direction
    for (i=0;i<nx;i++) s[i]=g[i]+beta*s[i];
    // model norm 
    eta[k]=dot(nx,x,x);
    // residuals norm 
    rho[k]=betanum;
    // If stopping criteria is true, copy back the last result and stop
    if (stopcriteria(gcv,rho,normb,nh,nq,k,stopc)){
      for(i=0;i<nq;i++) x[i]=xc[i];
      break;
    }
    k++;
  }
  // Keep the calculated model to compute the cost function
  for (i=0;i<nx;i++) xc[i]=x[i];
  // The calculated model x = Wm x_true  
  for (i=0;i<nx;i++) x[i]*=Wm[i];

  if (1) save_vector(&rho[0],k,"rhofile");


  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2
  
  /*****Cost Function *******/
  // J= Jm + Jd = ( x^T , x) +  (Wd^T ( b - L Wmi x )^T) (Wd ( b - L Wmi x))
  // Calculate Jd
  conv_corr(r,x,nt,nq,wavelet,nw,TRUE);    
  radonf(x,r,L,nt,nh,nq,dt,m2,d2,nfft,fmax,0);
  for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  
  for (i=0;i<ny;i++) r[i]*=Wd[i];

  J=dot(nx,xc,xc)+dot(ny,r,r);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);

  free2complex(m2aux);
  free1float(xc);
  free1float(rc);
  free1float(gcv);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(stemp);
  free1float(stemp2);
  free1float(stemp3);
  free1float(s);
  free1float(g);
  return(J);
}


void Atimesx_by_fft_loop(float *y, complex **A, float *x, complex **y2, complex **x2,
			 int nt, int n, float dt, int nfft, float fmax, int nf2)
{
  // Matrix vector multiplication ( y = A x ) using ffts.
  // The matrix is made circular and convolution in the time domain is performed 
  int freq, it, iq, ih; // indexes
  complex *yc, *xc; // auxiliar vectors to compute the FFT inside the frequency loop
  int maxfreq;
  float w;
  float wa;
  complex czero;  czero.r=czero.i=0;  
  int nf=nfft/2+1; // Only half of the fft samples are required
  float df=1./(nfft*dt); // sampling interval
  int nq=n;

  maxfreq=(int) (fmax/df+0.5);
  if (0) fprintf(stderr,"n=%d,nf=%d,nf2=%d,maxfreq=%d\n",n,nf,nf2,maxfreq);


  yc=ealloc1complex(nf);
  xc=ealloc1complex(nf);

  fftgo_xt2fx(-1,x,x2,n,nt,dt,nfft,nf);
  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax);
    Atimesx_by_fft(y2[freq],A[freq],x2[freq],n,nf2,yc,xc); // A^T A Wm^{-1} m 
    if (0) for (iq=0;iq<nq;iq++) y2[freq][iq]=x2[freq][iq];
    if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) y2[freq][iq]*=wa;
  }

  for (iq=0;iq<nq;iq++) y2[0][iq]=czero;  //dc can not be recovered  
  for (freq=maxfreq;freq<nf;freq++){
    for (iq=0;iq<nq;iq++) y2[freq][iq]=czero;
  }

  fftback_fx2xt(1,y,y2,n,nt,dt,nfft,nf);

  free1complex(yc);
  free1complex(xc);

  return;
}


void  Atimesx_by_fft(complex *y,complex *A,complex *x, int n, int nf, complex *yc, complex *xc)
{

  /*
     Perform the multiplication of a Toeplitz form with a vector.
     The procedure is done using a circulant matrix derived from
     the first row of the Toeplitz form. The Toeplitz form, in fact,
     is not needed. Everything you need is the firt row of the matrix.

     INPUT  
            RC(n): first row of the autocorrelation matrix (Toeplitz)
            f(n): a filter
            n : size of the vectors to multiply  
            nf: size of the circulant matrix 
            fc(nf): a pointer to an already allocated vector
            gc(nf): a pointer to an already allocated vector

     OUTPUT 
            d(n): the product of A with x

  */
      
  complex czero;
  int i;
  czero.r=czero.i=0.;

  for (i=0;i<n;i++) xc[i]=x[i];      // fill the first half of xc with x
  for (i=n;i<nf;i++) xc[i]=czero;    // pad with zeros the autocorrelation
 
  pfacc(-1,nf,xc);                   // FFT input 
  for (i=0;i<nf;i++) yc[i]=A[i]*xc[i];  // Circular convolution
  pfacc (1,nf,yc);                   // FFT output
  
  for (i=0;i<n;i++) y[i]=yc[i];      // Through away the second half of y
  for (i=0;i<n;i++) y[i]/=nf;        // Scale inverse FFT    

  return;
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
    for (ih=0;ih<nh;ih++) R[iq]+=conjg(L[ih][0])*Cd[ih]*L[ih][iq]/nh; //Top row of LL=LH*L
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


float testadjop(void (*oper) (float *model, float *data, complex ***L,int nt, int nh, int nq, 
			      float dt,complex **m2, complex **d2, int nfft, float fmax, 
			      int adj),complex ***L, int nt, int nh, int nq, float dt, 
		complex **m2, complex **d2, int nfft, float fmax, int adj)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  int iq;
  int ih;
  float test;
  int ny=nt*nh;
  int nx=nt*nq;


  dr1=ealloc1float(ny);
  mr1=alloc1float(nx);
  dr2=alloc1float(ny);
  mr2=alloc1float(nx);
 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();


  oper(mr2,dr1,L,nt,nh,nq,dt,m2,d2,nfft,fmax,1);
  oper(mr1,dr2,L,nt,nh,nq,dt,m2,d2,nfft,fmax,0);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;


  fprintf(stderr,"Test adjoint = %f , dp1=%f, dp2=%f \n",test,dp1,dp2);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}

void conv_corr(float *din, float *dout, int nt, int nh, float *wavelet, int nw, int adj)
{
  int ih, it;
  float *dtemp;

  dtemp=ealloc1float(nt+nw);
  memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 

  if (!adj)
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(0,0,nw,wavelet,nt,&din[ih*nt],dtemp);
      memcpy((void *) &dout[ih*nt],(const void*) dtemp,nt*sizeof(float));
    }
  else
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(1,0,nw,wavelet,nt,dtemp,&din[ih*nt]);
      memcpy((void *) &dout[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
  free1float(dtemp);


  return;
}

 














