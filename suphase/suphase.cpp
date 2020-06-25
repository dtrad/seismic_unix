/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPHASE: $Revision: 1.4 $ ; $Date: 2011/11/16 18:03:07 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
const char *sdoc[] = {
  " 									",
  " SUPHASE - PHASE manipulation by linear transformation			",
  " 									",
  "  suphase  <stdin >sdout      						",
  " 									",
  " Required parameters:							",
  " none									",
  " Optional parameters:							",
  " a=90			constant phase shift              		",
  " b=180/PI              linear phase shift				",
  " c=0.0			phase = a +b*(old_phase)+c*f;			",
  " 									",
  " test=1,2,3,4,5,6  different ways to calculate filters used for migration 2d/3d ",
  " Notes: 								",
  " A program that allows the user to experiment with changes in the phase",
  " spectrum of a signal.							",
  NULL};

/**************** end self doc ********************************/

/* definitions used internally */
#define AMPSP(c) sqrt(c.r*c.r+c.i*c.i)
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

complex mult(complex a, complex b){
  float ar=a.r,ai=a.i,br=b.r,bi=b.i;
  complex c;
  c.r = ar*br-ai*bi; c.i = ar*bi+ai*br;
  return c;
}

complex leakyOmega(int i, float d1, float dt, float rho){
  float omega = (2*PI*(i*d1)*dt);
  complex temp1(cos(omega),sin(omega));
  complex temp2 = 1 - rho*temp1;
  return temp2;
}

complex leakyOmega2(int i, float d1, float dt, float rho, int option){
  float omega = (2*PI*(i*d1)*dt);
  complex temp1(cos(omega),sin(omega));
  complex one(1.,0);
  complex cimag(0.,-1);
  if (option == 0) return(2*PI*(i*d1)*cimag); // i omega (no leaky)
  if (option == 1) return(1- rho*temp1); // 1 - Z (derivative i omega)

  complex temp2 = 0.5* (one+rho*temp1)/(one-rho*temp1); // integration 1/2*(1+z)/(1-z)
  if (option == 2) return(sqrt(one-rho*temp1)*temp2*temp2); // half derivative/double integration
  if (option == 3) return(one-rho*temp1)*temp2*temp2;  // derivative/double integration
  else return(1- rho*temp1); // default 1-Z

}

float normalize(float* a, int nt){
  float sum = 0;
  for (int i=0;i<nt;i++) sum+=(a[i]*a[i]);
  sum = sqrt(sum);
  for (int i=0;i<nt;i++) a[i]/=sum;
  return sum;
}

float freqweight(int j, float df, float f1, float f2){
  //  return weight for each frequency
  float w;
  float f=j*df;
  if(f<=f1) return (1.);
  if(f>=f2) return (0.);
  w = (f2-f)/(f2-f1);
  return (w);
}

/* segy trace */
segy tr;

int
main(int argc, char **argv)
{
  float *rt=NULL;		/* real trace			*/
  float *amp=NULL;	/* amplitude spectra		*/
  float *ph=NULL;		/* phase			*/
  register complex *ct=NULL;	/* complex time trace	*/

  int nt;			/* number of points on input trace	*/
  int nfft;		/* transform length			*/
  int nf;			/* number of frequencies in transform	*/

  float dt;		/* sampling interval in secs		*/
  float d1;		/* output sample interval in Hz		*/
  int count=0;		/* counter				*/

  /* linear phase function */
  float a;		/* bias (intercept) of new phase	*/
  float b;		/* slope of linear phase function	*/
  float c;		/* new phase value			*/
  float onfft;		/* 1/nfft				*/
  int test;
  float rho;
  float fmax1, fmax2;
  float fmin1, fmin2;
  
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);


  /* Get info from first trace */ 
  

  if (!gettr(&tr))  err("can't get first trace");
  nt = tr.ns;


  
  /* get parameters */


  /* dt is used only to set output header value d1 */
  if (!getparfloat("dt", &dt))	dt = ((double) tr.dt)/1000000.0;
  if (!dt) {
    dt = .004;
    warn("dt not set, assumed to be .002");
  }


  /* linear phase paramter values */
  if (!getparfloat("a", &a)) a = 0;
  if (!getparfloat("b", &b)) b = 180/PI;
  if (!getparfloat("c", &c)) c = 0.0;
  if (!getparint("test",&test))  test = 0;
  if (!getparfloat("rho",&rho)) rho = 1;
  if (!getparfloat("fmax1",&fmax1)) fmax1=70;
  if (!getparfloat("fmax2",&fmax2)) fmax2=80;
  if (!getparfloat("fmin1",&fmin1)) fmax1=2;
  if (!getparfloat("fmin2",&fmin2)) fmax2=5;

  a *= PI/180.0;
  b *= PI/180.0;
	
  /* Set up pfa fft */
  int nfft1 = npfaro(nt, LOOKFAC * nt);
  int nfft2 = npfao(nt,LOOKFAC * nt);
  
  if ((test > 0)&&(test<=6)) nfft = nfft2; // use complex to complex
  else if ((!test)||(test >= 7)) nfft = nfft1; // use real to complex.

  if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
  d1 = 1.0/(nfft*dt);
  nf = nfft/2 + 1;
  onfft = 1.0/nfft;
  float fnyquist = d1*(nf-1);

  fprintf(stderr,"nt=%d, dt=%f nf=%d fmin=%f fmax=%f nfft1=%d, nfft2=%d fnyquist=%f \n",
	  nt, dt, nf,fmin1,fmax2, nfft1,nfft2,fnyquist);
  fprintf(stderr,"test = 1; Sam 3d\ntest = 2 multiplying times i omega for pos freq\n");
  fprintf(stderr,"test = 3; Sam 2d\ntest = 4 multiplying times abs(omega) for pos freq\n");
  fprintf(stderr,"test = 5; x i omega for all freq \ntest = 6 x abs(omega) for all freq\n");
  if (test) fprintf(stderr,"performing test %d\n",test);  

  float* weight = ealloc1float(nf);
  for (int iw=0;iw< nf;iw++){
    weight[iw]= freqweight(iw,d1,fmax1,fmax2);
    if ((iw*d1)<fmin2) weight[iw]=1-freqweight(iw,d1,fmin1,fmin2); 
    if (weight[iw] != 1) fprintf(stderr,"weight[%d]=%f freq=%f\n",iw,weight[iw],iw*d1);
  }
  
	 
  checkpars();

  /* Allocate space */
  rt = ealloc1float(nfft);
  if (test)   ct = ealloc1complex(nfft);
  else{
    ct = ealloc1complex(nf);
    amp = ealloc1float(nf);
    ph = ealloc1float(nf);
  }
	
  /* Main loop over traces */
  count=0;
  do {
    register int i;
    count++;
    /* Load trace into rt (zero-padded) */
    memcpy((void *) rt, (const void *) &tr.data, nt*FSIZE);
    memset((void *) (rt + nt), (int) '\0', (nfft-nt)*FSIZE);
    memset((void *) ct,0,(nfft*2)*FSIZE);
    if (test==1){ // Sam 3d derivation
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = nf; i < nfft; ++i){  
	ct[i].r=0;ct[i].i=0;
      }
      for (i = 0; i < nfft; ++i){
	ct[i].r*=i*d1;
	ct[i].i*=i*d1;
      }
      pfacc(-1,nfft,ct);
      for (i = 0; i < nt; ++i) rt[i]=ct[i].i*onfft*(-2);
    }
    else if (test==2){ // multiply positive freqs only by i omega
      complex im;im.r=0;im.i=1;
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = 0; i < nf; ++i){
	ct[i].r*=i*d1;
	ct[i].i*=i*d1;
	ct[i]=mult(im,ct[i]);
      }
      pfacr(-1,nfft,ct,rt);
      for (i = 0; i < nt; ++i) rt[i]=rt[i]*onfft;
    }
    else if (test==3){ // Sam 2d derivation
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = nf; i < nfft; ++i){  
	ct[i].r=0;ct[i].i=0;
      }
      for (i = 0; i < nfft; ++i){
	ct[i].r*=i*d1;
	ct[i].i*=i*d1;
      }
      pfacc(-1,nfft,ct);
      for (i = 0; i < nt; ++i) rt[i]=ct[i].r*onfft*2;
    }
    else if (test==4){ // multiply positive freqs only by abs(omega)
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = 0; i < nf; ++i){
	ct[i].r*=i*d1;
	ct[i].i*=i*d1;
      }
      pfacr(-1,nfft,ct,rt);
      for (i = 0; i < nt; ++i) rt[i]=rt[i]*onfft;
    }

    else if (test==5){ // 3d explicit calculation (x i omega)
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = 0;i < nf;i++){
	int j = nfft-i;
	if (i==0) fprintf(stderr,"DC ct[%d]=(%f,%f) \n",i,ct[i].r,ct[i].i);
	else{
	  fprintf(stderr,"freq=%f Hz ct[%d]=(%f,%f) -- ct[%d]=(%f,%f)\n",
		  i*d1,
		  i,ct[i].r,ct[i].i,
		  j,ct[j].r,ct[j].i);
	}
      }
      complex im;im.r=0;im.i=1;
      for (i = 0;i<nf;i++){
	int j = nfft-i;
	float freq = 2*PI*i*d1;
	ct[i].r*=freq;
	ct[i].i*=freq;
	ct[j].r*=-freq;
	ct[j].i*=-freq;
	ct[i]=mult(ct[i],im);
	ct[j]=mult(ct[j],im);
      }
      pfacc(-1,nfft,ct);
      for (i = 0; i < nt; ++i) rt[i]=ct[i].r*onfft;      
    }
    else if (test==6){ // 2d explicit calculation (x abs(omega))
      for (i = 0; i < nt; ++i)  ct[i].r=rt[i];
      pfacc(1,nfft,ct);
      for (i = 0;i < nf;i++){
	int j = nfft-i;
	if (i==0) fprintf(stderr,"DC ct[%d]=(%f,%f) \n",i,ct[i].r,ct[i].i);
	else{
	  fprintf(stderr,"freq=%f Hz ct[%d]=(%f,%f) -- ct[%d]=(%f,%f)\n",
		  i*d1,
		  i,ct[i].r,ct[i].i,
		  j,ct[j].r,ct[j].i);
	}
      }
      complex im;im.r=0;im.i=1;
      for (i = 0;i<nf;i++){
	int j = nfft-i;
	float freq = 2*PI*i*d1;
	ct[i].r*=freq;
	ct[i].i*=freq;
	ct[j].r*=freq;
	ct[j].i*=freq;
      }
      pfacc(-1,nfft,ct);
      for (i = 0; i < nt; ++i) rt[i]=ct[i].r*onfft;      
    }
    else if (test == 7){ // 3d omega filter with integration and polynomial ratio
      complex czero(0.,0.);
      complex temp1, temp2;
      
      if (count <= 3){
	pfarc(-1, nfft, rt, ct);
	ct[0]=czero;
	for (i = 1; i < nf; ++i) {
	  if (count == 2) ct[i]/= (leakyOmega2(i,d1,dt,rho,0));
	  if (count == 3) ct[i]*= (leakyOmega2(i,d1,dt,rho,3));
	  ct[i] *= weight[i];
	}

	pfacr(1,nfft,ct,rt);
	for (i = 0; i < nt; ++i) rt[i]*=onfft;
	normalize(rt,nfft);
      }
      else if (count == 4){
	pfarc(-1, nfft, rt, ct);
	complex *ca = new complex[nf];
	complex *cb = new complex[nf];
	float *ra = new float[nfft];
	float *rb = new float[nfft];
	ca[0]=cb[0]=czero;
	
	for (i = 1; i < nf; ++i) {
	  ca[i]= ct[i]/(leakyOmega2(i,d1,dt,rho,0))*weight[i];
	  cb[i]= ct[i]*(leakyOmega2(i,d1,dt,rho,3))*weight[i];
	}

	pfacr(1,nfft,ca,ra);
	pfacr(1,nfft,cb,rb);
	for (i = 0; i < nt; ++i) ra[i]*=onfft;
	for (i = 0; i < nt; ++i) rb[i]*=onfft;

	normalize(ra,nfft);
	normalize(rb,nfft);
	for (i=0;i<nfft;i++) rt[i]=ra[i]-rb[i];

	delete [] ca;
	delete [] cb;
	delete [] ra;
	delete [] rb;
      }

    }
    else{ 		/* RtoC FFT */
      pfarc(1, nfft, rt, ct);
      for (i = 0; i < nf; ++i) {
	amp[i] = AMPSP(ct[i]);
	ph[i]  = a+b*atan2(ct[i].i,ct[i].r)+c*i;
      }
      for (i = 0; i < nf; ++i) {
	ct[i].r = amp[i]*cos(ph[i]);
	ct[i].i = amp[i]*sin(ph[i]);
      }
      pfacr(-1,nfft,ct,rt);
      for (i = 0; i < nt; ++i) rt[i]*=onfft;
    }


    memcpy((void *) tr.data, (const void *) rt, nt*FSIZE);
    //tr.ntr =4;
    puttr(&tr);
		
  } while (gettr(&tr));
  
  free1float(weight);
  return(CWP_Exit());
}
