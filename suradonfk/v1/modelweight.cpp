/*
  It computes the term Cm for the model weight 
  for the Hessian (L' Cd^{-1} L+ Cm^{-1})
  This term corresponds to the probability model, so that
  the distribution parameters sigma and norm are passed.
  Wm is a vector of dimension nx, but in fact is the diagonal 
  of the Wm matrix of size nx x nx.
  Inqut 
          m: model
          nx: number of model traces
          norm: implemented 1 Huber, 0 Cauchy, else L2
  Output
          Wm 
          eps1: standard deviation of the model

  Daniel Trad- 14 March 2000. UBC- Canada
  Based in Sacchi, 1996. phD thesis. UBC. Canada
*/
#include "su.h"
#include "math.h"
#include "Complex.h"

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=MAX(fabs(m[i]),sigmam);
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}
void weights_cgfft(float *m, int nx, int norm, float sigmam, float *Wm, int iter)

{ 
      int i;

      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
    
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=sqrt(MAX(fabs(m[i]),sigmam));
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=sqrt(sigmam*sigmam+m[i]*m[i]);
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      return;
}

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter)
{ 
      int i;
      
      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
      
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=1./sqrt(abs(m[i])*sigmam);
      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=1./sqrt(sigmam*sigmam+abs(m[i]*m[i]));
      else if(norm==2){ // Mask
	for (i=0;i<nx;i++)
	  Wm[i]=1+200./(1+exp(1*(abs(m[i])-sigmam)+0.5));
      }
      return;
}

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter)
{ 
  /* 
     It computes the model weights using a window in the model space m(f,q).
     For example, Wm[iq] is a function of m(f,iq), m(f-1,iq), ..., m(f-buffer,iq)
     Particularly useful for the dealiased RT (Hermman et al.)

  */     
      int  iq, iw;
      float maveg;

      if (iter==1){ 
	for (iq=0;iq<nq;iq++) Wm[iq]=1;
	return;
      }
      
      if (norm==1){ 
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(fabs(maveg)*sigmam);
	}
      }
      else if(norm==0){
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1./sqrt(sigmam*sigmam+fabs(maveg*maveg));
	}
      }
      else if(norm==2){ // Mask
	for (iq=0;iq<nq;iq++){
	  for (maveg=0, iw=0 ; iw < buffer; iw++) maveg+=abs(m[freq-iw][iq]);
          maveg/=buffer;
	  Wm[iq]=1+200./(1+exp(1*(fabs(maveg)-sigmam)+0.5));
	}
      }
      return;
}

void modelweight(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=fabs(m[isamax(nx,m,1)]);
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=fabs(m[i]);
      }
      else if(norm==0){
	if (maxm>1e-4) 
	  for (i=0;i<nx;i++){
	    // Solved !!!!!!!!!!!!!!!!!
	    // The right Wm from Cauchy is 
	    // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // But if M^-1 ATA x = M^-1 AT b is solved instead
	    // of the satndard form M=WmT *Wm 
	    Wm[i]=(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	    // Actually it works even better with (I don't know why)
	    //Wm[i]=Wm[i]*Wm[i];
	    //if (Wm[i]>2) Wm[i]=2; 
	  }
	else for (i=0;i<nx;i++) Wm[i]=1e-3;
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}

void modelweight_inv(float *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
      float maxm;
      maxm=m[isamax(nx,m,1)];
      //float scale=maxm*maxm*eps2;
      if (norm==1){
	for (i=0;i<nx;i++) Wm[i]=1./MIN(fabs(m[i]),eps1);
      }
      else if(norm==0){
	for (i=0;i<nx;i++){
	  // Solved !!!!!!!!!!!!!!!!!
	  // The right Wm from Cauchy is 
 	  // Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // But if M^-1 ATA x = M^-1 AT b is solved instead
	  // of the satndard form M=WmT *Wm 
	  Wm[i]=1./(eps1*eps1+m[i]*m[i]/(maxm*maxm));
	  // Actually it works even better with (I don't know why)
	  //Wm[i]=Wm[i]*Wm[i];
	  //if (Wm[i]>2) Wm[i]=2; 
	}
      }
      else if(norm==3){
	for (i=0;i<nx;i++) 
	  if (sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm))> 0.1) Wm[i]=1;
	  else Wm[i]=0.1;
      }
      else if(norm==4){
      	for (i=0;i<nx;i++) Wm[i]=(1./(1.+ exp(-1.*(fabs(m[i])-maxm/5.))));
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,eps1,Wm[isamax(nx,Wm,1)]);

      return;
}





void modelweight(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
      int i;
            
      if (norm==1)
	for (i=0;i<nx;i++) Wm[i]=abs(m[i]);

      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=eps1+pow(abs(m[i]),2.0);

      else for (i=0;i<nx;i++) Wm[i]=1.;     


      return;
}

void modelweight_inv(complex *m, int nx, int norm, float eps1, float *Wm)

{ 
  float *mabs;
  int ix;
  float sigma;

  mabs=ealloc1float(nx);
  if (norm==0 || norm ==1){ 
    for (ix=0;ix<nx;ix++) mabs[ix]=abs(m[ix]);
    sigma=quest(eps1,nx,mabs);
    //sigma=MAX(sigma,1e-3);
    //sigma=eps1;
    if (norm==1) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(MAX(mabs[ix],sigma)))+1e-7;
    else if(norm==0) for (ix=0;ix<nx;ix++) Wm[ix]=sqrt(1./(sigma*sigma+pow(mabs[ix],2.0)))+1e-7;
  }
  else for (ix=0;ix<nx;ix++) Wm[ix]=1.;     

  free1float(mabs);
  return;
}



void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=1;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=fabs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=fabs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);
      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);
      free1float(mabs);
      free1float(dabs);

      return;
}

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad)

{ 
      int i;
      float *mabs;
      float *dabs;
      int verbose=0;

      mabs=ealloc1float(nx);
      dabs=ealloc1float(ny);

      for (i=0;i<nx;i++) mabs[i]=abs(m[i]);
      *sigmam=quest(quantil1,nx,mabs);
      for (i=0;i<ny;i++) dabs[i]=abs(d[i]);
      *sigmad=quest(quantil2,ny,dabs);

      if (verbose) fprintf(stderr,"*sigmam=%f, *sigmad=%f\n",*sigmam,*sigmad);

      free1float(mabs);
      free1float(dabs);

      return;
}


void dataweigths(float *pos, int nh, float *Wd, int add)
{
  int ih;
  float *dh;
  float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
  int verbose=0;

  dh=ealloc1float(nh);

  for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   

  dh[0]=pos[1]-pos[0];
  dh[nh-1]=pos[nh-1]-pos[nh-2];
  
  for (ih=0;ih<nh;ih++){
    if (add) Wd[ih]*=MIN(dh[ih],dhav);
    else  Wd[ih]=MIN(dh[ih],dhav);
    if (verbose) fprintf(stderr," Wd[%d]=%f\n",ih,Wd[ih]);
  }

  free1float(dh);
  return;

}









































