#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "radonwtcgls_tfd.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

float wpcgnr(void (*oper)  (float *model, float *data, complex ***L,int nt, int nh, 
			    int nq, float dt,complex **m2, complex **d2, int nf0, 
			    float fmax, int adj),
	     float *x, float *b,  complex ***L, int nt, int nh, 
	     int nq, float dt, complex **m2, complex **d2, int nf0, float fmax, 
	     float *Wm, float *Wd, int itercg, int iter_end, float eps1, float eps2,
	     float eps, float step);

void radonf(float *model, float *data, complex ***L,int nt,int nh,int nq,
	    float dt,complex **m2, complex **d2, int nf0,float fmax, int adj);

void modelweight(float *m, int nx, int norm, float eps1, float *Wm);

int taper(float **data, int nt, int nh, int ntaper,int flag);

float dot(int n, float *a, float *b);

void save_vector(float *d, int n, char* name);

void radonwtcgls_tfd(float *pos, int nh, float **data, float *t, int nt, float dt,
		     float **model, float *q, int nq, float dq, float eps1, float eps2,
		     float eps, float fmax, float **Wd, int itercg, int iter_end,
		     int norm,float step, int testadj, float quantil, int rtmethod, 
		     float depth)
{

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
  complex **d2;
  complex **m2;
  int adj;
  int nf0;
  int nfft;
  int nf;
  float df;
  int maxfreq;
  //float fmax=1./(2*dt)*0.7;
  int nfreq;
  int filtout=0;
  int taperflag=0;
  double pi=PI;
  complex arg;
  float w;
  int mute=0;
  int parmute=0;
  float quantil1=eps1;
  float quantil2=eps2;
  float sigmam;
  float sigmad;
  float *J;
  
  //////////////////////////////////////////////////////////////////////////
  // Set pointer to function for Conjugate Gradients
  void (*radon) (float *model, float *data,  complex ***L,int nt, int nh, int nq, 
		 float dt, complex **m2, complex **d2, int nf0, float fmax, int adj);
  
  // Set Fourier transform parameters  
  nfft = npfaro(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt
  if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
  nf0= nfft + 2;
  nf= nfft/2 + 1;
  df=1/(nf0*dt);
  maxfreq=(int) (fmax/df+0.5);
  nfreq=nf;

  fprintf(stderr,"radon:nfft=%d, nf=%d nf0=%d nf=%d, nfreq=%d, maxfreq=%d, df=%f\n",nfft,
	  nf,nf0,nf,nfreq,maxfreq,df);   

  fprintf(stderr,"step=%f norm=%d\n",step, norm);   

  radon=radonf;

  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc3complex(nq,nh,nfreq);
  g=ealloc1float(nh);
  Wm=ealloc2float(nt,nq);
  J=ealloc1float(iter_end+1);

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
  fprintf(stderr,"++++++++++++++++++2\n");
  // Construct the Radon operator in Frequency
  
  if (rtmethod==1) for (ih=0;ih<nh;ih++) g[ih]=pos[ih];
  else if (rtmethod==2) for (ih=0;ih<nh;ih++) g[ih]=pos[ih]*pos[ih];
  else if (rtmethod==3) for (ih=0;ih<nh;ih++) g[ih]=sqrt(pos[ih]*pos[ih]+depth*depth)-depth;

  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    radon_matrix(L[freq],g,q,nh,nq,w);
  }

  // Compute the adjoint 
  radon(model[0],data[0],L,nt,nh,nq,dt,m2,d2,nf0,fmax,1);

  for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) Wm[iq][it]=1.0;
  fprintf(stderr,"quantil1=%f, quantil2=%f\n",quantil1,quantil2);
  for (j=1;j<=iter_end;j++){
    
    if (j==2) deviations(model[0],nq*nt,data[0],nh*nt,norm,quantil1,quantil2,&sigmam,&sigmad);
    fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
    if (j>1)  weights(model[0],nt*nq,norm,sigmam,Wm[0],j);    
    
    J[j-1]=wpcgnr(radon,model[0],data[0],L,nt,nh,nq,dt,m2,d2,nf0,fmax,Wm[0],Wd[0],itercg,
	   iter_end,eps1,eps2,eps,step);
    
  }
  
  
  if (1) plotcurves(J,iter_end,1,"Cost");    
  if (mute){
    int nmute;
    iq=0; while(fabs(q[iq])<parmute) iq++; 
    //nmute=nq-iq;
    fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
    //taper(m,nt,nq,nmute,2);
    taper(model,nt,nq,iq,1);
  }
  
  //radon(model[0],data[0],L,nt,nh,nq,dt,d2,m2,nf0,fmax,0);
  free1float(J);
  free2float(Wm);
  free1float(g);
  free3complex(L);
  free2complex(d2);
  free2complex(m2);

  return;
  
}

void radonf(float *model, float *data, complex ***L,int nt, int nh, int nq, 
	    float dt,complex **m2, complex **d2, int nf0, float fmax, int adj)
{
  int freq, it, iq, ih; // indexes
  int maxfreq;
  int nfreq;
  float df;
  float w;
  float wa;
  complex czero;
  czero.r=czero.i=0;  

  nfreq=nf0/2;
  df=1/(nf0*dt);
  maxfreq=(int) (fmax/df+0.5);

  //fprintf(stderr,"radonf: nf0=%d df=%f, nfreq=%d maxfreq=%d\n",nf0,df,nfreq,maxfreq);   

  if (adj){
    fftgo1(-1,data,d2,nh,nt,dt,nf0);
    if (0){
      fftback1(1,data,d2,nh,nt,dt,nf0);
      save_gather(data,nh,nt,dt,"dtemp.su");
      system("suxwigb < dtemp.su perc=99 title=\"dtemp\" &");     
      system("sufft < dtemp.su | suamp | suxwigb perc=99 title=\"dtemp\" &");       
    }

    for (freq=1;freq<maxfreq;freq++){
      w=2*PI*freq*df;
      wa=freqweight(freq,df,fmax-10,fmax);
      Atimesx(d2[freq],L[freq],m2[freq],nh,nq,adj);
      if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) m2[freq][iq]*=wa;
    }

    for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nfreq;freq++){
      for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
    }

    fftback1(1,model,m2,nq,nt,dt,nf0);

  }
  else{
    fftgo1(-1,model,m2,nq,nt,dt,nf0);
    for (freq=1;freq<maxfreq;freq++){
      w=2*PI*freq*df;
      wa=freqweight(freq,df,fmax-10,fmax);
      Atimesx(d2[freq],L[freq],m2[freq],nh,nq,adj);
      if ((wa<1)&&(wa>0)) for (ih=0;ih<nh;ih++) m2[freq][ih]*=wa;
    }
    
    for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nfreq;freq++){
      for (ih=0;ih<nh;ih++) d2[freq][ih]=czero;
    }
    
    fftback1(1,data,d2,nh,nt,dt,nf0);
  }

  return;
}


void radonf_new(float *model, float *data, complex ***L,int nt, int nh, int nq, 
	    float dt,complex **m2, complex **d2, int nfft, float fmax, int adj)
{
  int freq, it, iq, ih; // indexes
  int maxfreq;
  int nfreq;
  float df;
  float w;
  float wa;
  complex czero;

  
  int nf=nfft/2+1;

  czero.r=czero.i=0;  
  fft_parameters(nt,dt,&nfft,&nf,&df);

  //nfreq=nf0/2;
  //df=1/(nf0*dt);
  maxfreq=(int) (fmax/df+0.5);

  //fprintf(stderr,"radonf: nf0=%d df=%f, nfreq=%d maxfreq=%d\n",nf0,df,nfreq,maxfreq);   

  if (adj){
    fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);

    for (freq=1;freq<maxfreq;freq++){
      w=2*PI*freq*df;
      wa=freqweight(freq,df,fmax-10,fmax);
      Atimesx(d2[freq],L[freq],m2[freq],nh,nq,adj);
      if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) m2[freq][iq]*=wa;
    }

    for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nfreq;freq++){
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
    for (freq=maxfreq;freq<nfreq;freq++){
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


float wpcgnr(void (*oper)  (float *model, float *data, complex ***L,int nt, int nh, 
			    int nq, float dt,complex **m2, complex **d2, int nf0, 
			    float fmax, int adj),
	     float *x, float *b,  complex ***L, int nt, int nh, 
	     int nq, float dt, complex **m2, complex **d2, int nf0, float fmax, 
	     float *Wm, float *Wd, int itercg, int iter_end, float eps1, float eps2,
	     float eps, float step)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function and square residuals

  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.
  // rc Weighted residuals to compute cost function
  // rc Weighted solution to compute cost function

  float *r2,*r,*g,*s,*z,*w,*rho,*eta,rhold, *rc, *xc;

  int nx=nt*nq;
  int ny=nt*nh;  

  g=ealloc1float(nx);
  s=ealloc1float(nx);
  z=ealloc1float(nx);
  r=ealloc1float(ny);
  r2=ealloc1float(ny);
  w=ealloc1float(ny);
  rho=ealloc1float(itercg+1);
  eta=ealloc1float(itercg+1);
  xc=ealloc1float(nx); 
  rc=ealloc1float(ny);

  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);
  int restart=1;

  if (restart){
    for (i=0;i<nx;i++) x[i]=0.;
    for (i=0;i<ny;i++) r[i]=Wd[i]*b[i];
  }
  else{
    (*oper) (x,r,L,nt,nh,nq,dt,m2,d2,nf0,fmax,0);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }

  for (i=0;i<ny;i++) r2[i]=r[i]*Wd[i]; 
  (*oper) (g,r2,L,nt,nh,nq,dt,m2,d2,nf0,fmax,1);
  for(i=0;i<nx;i++) z[i]=g[i]*Wm[i];


  normb=dot(nx,z,z);
  rho[0]=1;      
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=eps*2;
  k=0;
  while((rhold>eps)&&(k<itercg)){
    k++;
    (*oper) (s,w,L,nt,nh,nq,dt,m2,d2,nf0,fmax,0);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  
 
    for (i=0;i<ny;i++) r2[i]=r[i]*Wd[i];    
    (*oper) (g,r2,L,nt,nh,nq,dt,m2,d2,nf0,fmax,1);

    for(i=0;i<nx;i++) z[i]=g[i]*Wm[i];
    
    rho[k]=dot(nx,z,z)/normb;    
    //fprintf(stderr,"resm[%d]=%e,=======>res[%d]==%e\n",k,rho[k],k,dot(ny,r,r));
    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];
     
    eta[k]=dot(nx,x,x);

  }
  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2
  
  /*****Cost Function *******/

  for (i=0;i<nx;i++) xc[i]=x[i]*Wm[i];

  for (i=0;i<ny;i++) rc[i]=r[i]*Wd[i];

  J=dot(nx,xc,xc)+dot(ny,rc,rc);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);

  free1float(xc);
  free1float(rc);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r2);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);

  return(J);
}






















