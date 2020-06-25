#include "Complex.h"
#include "su.h"
#include "segy.h"

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void xplotfkspectrum(complex **d, int nh, int nt, float dt, float dx, float fmin, float kmin, const char *s, int num, const char *s2);
void fftshift(int ntfft, int nxfft, float** A, float** Ashift);
void fftshift2(int ntfft, int nxfft, float** A, float** Ashift);
float qamp(float q1, float q2, float q3, float q4);
void qamp2d(int ntfft, int nxfft, float** q1, float** q2, float** q3, float** q4, float** A);
void qfft1d(float sign, int ntfft, int nxfft, float** pfft1, float** pfft2, float** pfft3, float** pfft4, float** qpfft1, float** qpfft2, float** qpfft3, float** qpfft4);
void qfft2d(float sign, int ntfft, int nxfft, float** pfft1, float** pfft2, float** pfft3, float** pfft4, float** qpfft1, float** qpfft2, float** qpfft3, float** qpfft4);
void qfft1d_x_to_k(float sign, int nxfft, float* pfft1, float* pfft2, float* pfft3, float* pfft4);
void xw2kw(complex** cpfft, int nxfft, int nk, int nw, complex* freqslice, int sign);
void save_gather(float **d, int nh, int nt, float dt, const char* name);
void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);

bool add(float** data1, float** data2, int nt, int nh, float value=1);
bool copy(float** data1, float** data2, int nt, int nh, float value=1);
int countNonZero(complex** model, int nt, int nh);

void xplotgather(complex **d, int nh, int nt, float dt, float dx, const char *s, int num, const char *s2);

bool process(float** datain1, float** datain2, float** dataout1, float** dataout2, int nt, int  nx, float dt){
  int it, ix, iw;
  //  for (ih=0;ih<nh;ih++)  for (it=0;it<nt;it++) dataout[ih][it]=datain[ih][it];
  complex czero;czero.r=czero.i=0;
  int ntfft, nxfft, nw, nk;
  float scale;

  // defining sampling for plots.
  float df, dk; // sampling interval in fk
  float dx = 100; // sampling interval in offset
  int Iter = 50; // iterations for thresholding
  float perci = 0.999;
  float percf = 0.001;
  float alpha = 1.0; // denoising parameter, 1 => no noise , ->0 => noisy input traces.
  bool plotInput   = false;
  bool plotInputFK = false;
  bool plotInput2  = false; // after zeroes
  bool plotFXIn    = false;
  bool plotQFKIn   = true; // QFK spectrum after zeroes
  bool plotFKIn    = false; // after zeroes
  bool plotFKOut   = false;
  bool plotOutput  = false;

  int padfactor = 4;      // pad factor for both dimensions

  // plot input
  if (plotInput){
    save_gather(datain1,nx,nt,dt,"input1.su");
    system("suximage < input1.su perc=98 clip=1.0 legend=1 title=\"input\" &");
    save_gather(datain2,nx,nt,dt,"input2.su");
    system("suximage < input2.su perc=98 clip=1.0 legend=1 title=\"input\" &");
  }
  MARK;

  /* copy data from input to FFT array and pad with zeros */
  
  ntfft = npfar(padfactor*nt);
  nxfft = npfa(padfactor*nx);
  nw=ntfft/2+1;
  nk=nxfft;
  scale=nxfft*ntfft;
  df=1./(ntfft*dt);
  dk=1./(nxfft*dx);

  //  complex* freqslice1 = ealloc1complex(nxfft);
  //  complex* freqslice2 = ealloc1complex(nxfft);

  // float** datatrue = 0; //declare outside if you need this variable in another place
  // plot input in FK
  //    if (datatrue){
  //    }

  // zero some traces to test interpolation 
  //  for (ix=0;ix< nx;ix+=3) for (it=0;it<nt;it++) datain[ix][it]=0;
  float rndnum = 0;
  for (ix=0;ix< nx;ix+=1){
    // generate random number distributed between 0 and 1 
    // (== Uniformly distributed == matlab function "rand")
    rndnum = franuni();
    if (rndnum>0.5){
      for (it=0;it<nt;it++){ 
	datain1[ix][it]=0; 
	datain2[ix][it]=0;
      }
    }
  } 
  // plot after zeros
  if (plotInput2){
    save_gather(datain1,nx,nt,dt,"input1_decimated.su");
    system("suximage < input1_decimated.su title=\"after zeros\" legend=1 clip=1.0 &");
    save_gather(datain2,nx,nt,dt,"input2_decimated.su");
    system("suximage < input2_decimated.su title=\"after zeros\" legend=1 clip=1.0 &");
  }
  
  fprintf(stderr,"nt=%d, nx=%d, ntfft=%d, nxfft=%d, nw=%d, nk=%d\n",nt,nx,ntfft,nxfft,nw,nk);

  float**   pfft  = ealloc2float(ntfft,nxfft); // trace oriented (Hale's reversed convention for alloc)
  float**   pfft1  = ealloc2float(ntfft,nxfft); // trace oriented (Hale's reversed convention for alloc)
  float**   pfft2  = ealloc2float(ntfft,nxfft); 
  float**   pfft3  = ealloc2float(ntfft,nxfft); 
  float**   pfft4  = ealloc2float(ntfft,nxfft); 
  float**   qpfft1 = ealloc2float(ntfft,nxfft);     
  float**   qpfft2 = ealloc2float(ntfft,nxfft);
  float**   qpfft3 = ealloc2float(ntfft,nxfft);
  float**   qpfft4 = ealloc2float(ntfft,nxfft);
  //  complex**   cpfft = ealloc2complex(nw,nk);     
  //  complex** cpfft  = ealloc2complex(nk,nw); // freq oriented (not used but may be better)
  
  /* copy data from input to FFT array and pad with zeros */
  for (ix=0;ix<nx;ix++){
    for (it=0; it<nt; it++){ 
      //      pfft[ix][it] = datain[ix][it];
      pfft1[ix][it]= 0; 
      pfft2[ix][it]= datain1[ix][it];  
      pfft3[ix][it]= datain2[ix][it];
      pfft4[ix][it]= 0;
    }
    for (it=nt; it< ntfft;it++){ 
      //      pfft[ix][it] = 0.0;
      pfft1[ix][it] = 0.0;
      pfft2[ix][it] = 0.0;
      pfft3[ix][it] = 0.0;
      pfft4[ix][it] = 0.0;
    }
  }		
  for (ix=nx;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      //      pfft[ix][it] = 0.0;
      pfft1[ix][it] = 0.0;
      pfft2[ix][it] = 0.0;
      pfft3[ix][it] = 0.0;
      pfft4[ix][it] = 0.0;
    }
  }
  save_gather(pfft2,nx,nt,dt,"pfft2.su");
  system("suximage < pfft2.su perc=98 clip=1 legend=1 title=\"component 1 input\" &");	
  save_gather(pfft3,nx,nt,dt,"pfft3.su");
  system("suximage < pfft3.su perc=98 clip=1 legend=1 title=\"component 2 input\" &");	
	
  if (plotQFKIn){
    float** qpfk1 = ealloc2float(ntfft,nxfft);
    float** qpfk2 = ealloc2float(ntfft,nxfft);
    float** qpfk3 = ealloc2float(ntfft,nxfft);
    float** qpfk4 = ealloc2float(ntfft,nxfft);
    /* Quaternion Fourier transform t-x to w-k to plot the FK spectrum */
    qfft2d(1,ntfft,nxfft,pfft1,pfft2,pfft3,pfft4,qpfk1,qpfk2,qpfk3,qpfk4);
    float** Aq = ealloc2float(ntfft,nxfft);
  	  
    qamp2d(ntfft,nxfft,qpfk1,qpfk2,qpfk3,qpfk4,Aq);
    float** Aqshift = ealloc2float(ntfft,nxfft);
    fftshift2(ntfft,nxfft,Aq,Aqshift);
    complex** Aqshift2 = ealloc2complex(ntfft,nxfft);
    for (ix=0;ix<nxfft;ix++){ 
      for(it=0;it<ntfft;it++){ 
	Aqshift2[ix][it].r = Aqshift[ix][it];
      }
    }
  
    float fmin = -(ntfft/2)*df;
    float kmin = -(nxfft/2)*dk;
    xplotfkspectrum(Aqshift2,nxfft,ntfft,df,dk,fmin,kmin,"Aqshift2.su",1,"perc=99 wclip=0 bclip=100 title=\"Quaternion FK spectrum\" &"); // plot [k,w]
  }      
 
  //save_gather(Aqshift2,nxfft,ntfft,df,"Aq.su");
  //system("suximage < Aq.su title=\"Quaternion data in FX\" legend=1 perc=98 clip=1000 &");

  /* Fourier transform t to w */
  //pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]); // cpfft is [x,w]

  /* Quaternion Fourier transform t-x to w-x */
  qfft1d(1,ntfft,nxfft,pfft1,pfft2,pfft3,pfft4,qpfft1,qpfft2,qpfft3,qpfft4);

  MARK;

  //   algorithm starts*/
  float* freqsliceq1= ealloc1float(nxfft);
  float* freqsliceq2= ealloc1float(nxfft);
  float* freqsliceq3= ealloc1float(nxfft);
  float* freqsliceq4= ealloc1float(nxfft);
  float* freqslice2q1= ealloc1float(nxfft);
  float* freqslice2q2= ealloc1float(nxfft);
  float* freqslice2q3= ealloc1float(nxfft);
  float* freqslice2q4= ealloc1float(nxfft);
  
  float* mabs = ealloc1float(nxfft);  
  float* mabsiter = ealloc1float(nxfft);
  float* wd   = ealloc1float(nxfft);
  float sigma;
	
  for (ix=0; ix<nk;ix++){
    float sum = 0;
    float Aqpfft = 0;
    for (iw = 0; iw< ntfft; iw++){
      //		Aqpfft = qpfft1[ix][iw] + qpfft2[ix][iw] + qpfft3[ix][iw] + qpfft4[ix][iw];
      Aqpfft = qamp(qpfft1[ix][iw], qpfft2[ix][iw], qpfft3[ix][iw], qpfft4[ix][iw]);
      //		fprintf(stderr,"Aqpfft[%d][%d]=%f \n",ix,iw,Aqpfft);
      sum += Aqpfft;
    }
    if (sum) wd[ix]=1;
    else{ 
      wd[ix]=0;
      //      fprintf(stderr,"wd[%d]=%f \n",ix,wd[ix]);
    }
    //      fprintf(stderr,"sum[%d]=%f,  ",ix,sum);
    //      fprintf(stderr,"wd[%d]=%f \n",ix,wd[ix]);
	  
  }
    
  // freq by freq algorithm : for the quaternions ensure you process both halves of the spectrum!
  for (iw=0;iw<ntfft;iw++){
    
    // copy freq slice (needs because arrays are trace oriented)
    for (ix=0;ix<nk;ix++){ 
      freqsliceq1[ix]=freqslice2q1[ix]=qpfft1[ix][iw];
      freqsliceq2[ix]=freqslice2q2[ix]=qpfft2[ix][iw];
      freqsliceq3[ix]=freqslice2q3[ix]=qpfft3[ix][iw];
      freqsliceq4[ix]=freqslice2q4[ix]=qpfft4[ix][iw];
    }

    //    pfacc(-1,nxfft,freqslice2);  // FFT x to k
    /* Quaternion Fourier transform w-x to w-k */
    qfft1d_x_to_k(1, nxfft, freqslice2q1, freqslice2q2, freqslice2q3, freqslice2q4);

    // threshold in k
    for (ix=0;ix<nk;ix++) mabs[ix]=qamp(freqslice2q1[ix],freqslice2q2[ix],freqslice2q3[ix],freqslice2q4[ix]);//qamp(freqslice2q1[ix],freqslice2q2[ix],freqslice2q3[ix],freqslice2q4[ix],mabs[ix]);

    //pfacc(1,nxfft,freqslice2);  // iFFT k to x
    qfft1d_x_to_k(-1, nxfft, freqslice2q1, freqslice2q2, freqslice2q3, freqslice2q4);
    for (ix=0;ix<nk;ix++){ 
      freqslice2q1[ix]/=nxfft;
      freqslice2q2[ix]/=nxfft;
      freqslice2q3[ix]/=nxfft;
      freqslice2q4[ix]/=nxfft;
    }
	  
    for (int iter=1;iter<Iter;iter++){  // external thresholding loop
      //      	pfacc(-1,nxfft,freqslice2);  // FFT x to k
      qfft1d_x_to_k(1, nxfft, freqslice2q1, freqslice2q2, freqslice2q3, freqslice2q4);
	
      int count = 0; // count kept
      int kept = 0;  // last index kept
        
      // This is to increase the thresholding within each internal iteration
      sigma=quest(perci - (iter-1)*((perci-percf)/(Iter-1)),nk,mabs); 
      for (ix=0;ix<nk;ix++) mabsiter[ix]=qamp(freqslice2q1[ix],freqslice2q2[ix],freqslice2q3[ix],freqslice2q4[ix]);//qamp(freqslice2q1[ix],freqslice2q2[ix],freqslice2q3[ix],freqslice2q4[ix],mabsiter[ix]);
      for (ix=0;ix<nk;ix++){
	if (mabsiter[ix]<sigma){ 
	  freqslice2q1[ix] = 0;
	  freqslice2q2[ix] = 0;
	  freqslice2q3[ix] = 0;
	  freqslice2q4[ix] = 0;
	}
	else{ kept=ix;count++; }
      }
      //fprintf(stderr,"sigma=%f count = %d kept=%d \n",sigma,count,kept);
      //	pfacc(1,nxfft,freqslice2);  // k to x
      qfft1d_x_to_k(-1, nxfft, freqslice2q1, freqslice2q2, freqslice2q3, freqslice2q4);
      for (ix=0;ix<nk;ix++){ 
	freqslice2q1[ix]/=nxfft;
	freqslice2q2[ix]/=nxfft;
	freqslice2q3[ix]/=nxfft;
	freqslice2q4[ix]/=nxfft;
      }
		
      // add new traces to original freqslice
      for (ix=0;ix<nk;ix++){ 
	freqslice2q1[ix]=(alpha*freqsliceq1[ix]) + (1-alpha*wd[ix])*freqslice2q1[ix]; // x,w
	freqslice2q2[ix]=(alpha*freqsliceq2[ix]) + (1-alpha*wd[ix])*freqslice2q2[ix]; // x,w
	freqslice2q3[ix]=(alpha*freqsliceq3[ix]) + (1-alpha*wd[ix])*freqslice2q3[ix]; // x,w
	freqslice2q4[ix]=(alpha*freqsliceq4[ix]) + (1-alpha*wd[ix])*freqslice2q4[ix]; // x,w
      }
    }
	  
    // put back the slice
    for (ix=0;ix<nk;ix++){ 
      qpfft1[ix][iw]=freqslice2q1[ix];
      qpfft2[ix][iw]=freqslice2q2[ix];
      qpfft3[ix][iw]=freqslice2q3[ix];
      qpfft4[ix][iw]=freqslice2q4[ix];
    }
		
  }
  
  free1float(wd);
  free1float(freqslice2q1);
  free1float(freqslice2q2);
  free1float(freqslice2q3);
  free1float(freqslice2q4);
  free1float(freqsliceq1);
  free1float(freqsliceq2);
  free1float(freqsliceq3);
  free1float(freqsliceq4);
  free1float(mabs);
  // algorithm ends


  //  /* Fourier transform w to t */
  //  pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
  /* Quaternion Fourier transform w-x to t-x */
			
  qfft1d(-1,ntfft,nxfft,qpfft1,qpfft2,qpfft3,qpfft4,pfft1,pfft2,pfft3,pfft4);
  for (ix=0;ix<nx;ix++){ 
    for (it=0; it<nt; it++){
      dataout1[ix][it]=pfft2[ix][it]/ntfft;
      dataout2[ix][it]=pfft3[ix][it]/ntfft;
    }
  }
  
	
  save_gather(dataout1,nx,nt,dt,"output1.su");
  system("suximage < output1.su clip=1 legend=1 title=\"output 1\" &");	
  save_gather(dataout2,nx,nt,dt,"output2.su");
  system("suximage < output2.su clip=1 legend=1 title=\"output 2\" &");	

  free2float(pfft1);
  free2float(pfft2);
  free2float(pfft3);
  free2float(pfft4);
  free2float(qpfft1);
  free2float(qpfft2);
  free2float(qpfft3);
  free2float(qpfft4);
  return true;

}


bool add(float** data1, float** data2, int nt, int nh, float value){
  int it, ih;
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data1[ih][it]+= (value*data2[ih][it]);
  return true;
}

bool copy(float** data1, float** data2, int nt, int nh, float value){
  int it, ih;
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data1[ih][it] = (value*data2[ih][it]);
  return true;
}

int countNonZero(complex** model, int nt, int nh){
  int it, ih, count=0;
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) if (abs(model[it][ih])) count++; 
  return count;
}


void xw2kw(complex** cpfft, int nxfft, int nk, int nw, complex* freqslice, int sign){
  int ix, iw;
  float scale = 1;
  if (sign == 1) scale = 1./nxfft;
  for (iw=0;iw<nw;iw++){
    for (ix=0;ix<nk;ix++) freqslice[ix]=cpfft[ix][iw];
    pfacc(sign,nxfft,freqslice);
    
    for (ix=0;ix<nk;ix++) cpfft[ix][iw]=scale*freqslice[ix];
  }
  return;
}

void wx2wk(complex** cpfft, int nxfft, int nk, int nw, int sign){
  int ix, iw;
  
  for (iw=0;iw<nw;iw++){
    pfacc(sign,nxfft,cpfft[iw]);
    if (sign == 1) for (ix=0;ix<nxfft;ix++) cpfft[iw][ix] /= nxfft;
  }
  return;
}

void qfft1d(float sign, int ntfft, int nxfft, float** pfft1, float** pfft2, float** pfft3, float** pfft4, float** qpfft1, float** qpfft2, float** qpfft3, float** qpfft4){
  // The forward and reverse quaternion Fourier transform along the time axis
  
  int ix, it;
  float Bfw11, Bfw12, Bfw13, Bfw21, Bfw22, Bfw23, Bfw31, Bfw32, Bfw33 ;
  float Brv11, Brv12, Brv13, Brv21, Brv22, Brv23, Brv31, Brv32, Brv33 ;
  float** P1real = ealloc2float(2*ntfft,2*nxfft);
  float** P1imag = ealloc2float(2*ntfft,2*nxfft);
  float** P2real = ealloc2float(2*ntfft,2*nxfft); 
  float** P2imag = ealloc2float(2*ntfft,2*nxfft);
  complex** P1 = ealloc2complex(2*ntfft,2*nxfft+1);
  complex** P2 = ealloc2complex(2*ntfft,2*nxfft+1); 
    
  // float**   example1 = ealloc2float(ntfft,nxfft);
  // complex** example2 = ealloc2complex(ntfft,nxfftq);


  // define the basis for the fast quaternion Fourier transform
  Bfw11 = 1/sqrt(3); // the first vector of the new basis [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]
  Bfw12 = 1/sqrt(3); 
  Bfw13 = 1/sqrt(3);
  Bfw21 = 0;         // the second vector of the new basis [0 1/sqrt(2) -1/sqrt(2)] (orthogonal to the first)
  Bfw22 = 1/sqrt(2);
  Bfw23 = -1/sqrt(2);
  Bfw31 = -1*(Bfw12*Bfw23 - Bfw13*Bfw22); // the third vector of the new basis, -1*(row1Xrow2) = orthogonal to the first two by definition
  Bfw32 = -1*(Bfw13*Bfw21 - Bfw11*Bfw23);
  Bfw33 = -1*(Bfw11*Bfw22 - Bfw12*Bfw21);
  Brv11 = Bfw11;    // B_rv is the transpose of B_fw
  Brv12 = Bfw21;
  Brv13 = Bfw31;
  Brv21 = Bfw12;
  Brv22 = Bfw22;
  Brv23 = Bfw32;
  Brv31 = Bfw13;
  Brv32 = Bfw23;
  Brv33 = Bfw33; 

  // change of basis
  for (ix=0;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      P1real[ix][it] = pfft1[ix][it];
      P1imag[ix][it] = pfft2[ix][it]*Bfw11 + pfft3[ix][it]*Bfw12 + pfft4[ix][it]*Bfw13;
      P2real[ix][it] = pfft2[ix][it]*Bfw21 + pfft3[ix][it]*Bfw22 + pfft4[ix][it]*Bfw23;
      P2imag[ix][it] = pfft2[ix][it]*Bfw31 + pfft3[ix][it]*Bfw32 + pfft4[ix][it]*Bfw33;     
      P1[ix][it].r = P1real[ix][it];  
      P1[ix][it].i = P1imag[ix][it];
      P2[ix][it].r = P2real[ix][it];  
      P2[ix][it].i = P2imag[ix][it];
    }
  }
  
  if (sign == 1){ // forward transform
    pfa2cc(1,1,ntfft,nxfft,P1[0]);  
    pfa2cc(1,1,ntfft,nxfft,P2[0]);  
  }
  else{           // reverse transform
    pfa2cc(-1,1,ntfft,nxfft,P1[0]);  
    pfa2cc(-1,1,ntfft,nxfft,P2[0]);  
  }
  //  if (1)
  for (ix=0;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      P1real[ix][it] = P1[ix][it].r;
      P1imag[ix][it] = P1[ix][it].i;
      P2real[ix][it] = P2[ix][it].r;
      P2imag[ix][it] = P2[ix][it].i;
      qpfft1[ix][it] = P1real[ix][it];
      qpfft2[ix][it] = P1imag[ix][it]*Brv11 + P2real[ix][it]*Brv12 + P2imag[ix][it]*Brv13;
      qpfft3[ix][it] = P1imag[ix][it]*Brv21 + P2real[ix][it]*Brv22 + P2imag[ix][it]*Brv23;
      qpfft4[ix][it] = P1imag[ix][it]*Brv31 + P2real[ix][it]*Brv32 + P2imag[ix][it]*Brv33;
    }
  }
  //MARK; // this tells you where an error is occuring (to solve a segmentation fault)
  return;
}


void qfft1d_x_to_k(float sign, int nxfft, float* pfft1, float* pfft2, float* pfft3, float* pfft4){
  // The forward and reverse quaternion Fourier transform in place
  
  int ix;
  float Bfw11, Bfw12, Bfw13, Bfw21, Bfw22, Bfw23, Bfw31, Bfw32, Bfw33 ;
  float Brv11, Brv12, Brv13, Brv21, Brv22, Brv23, Brv31, Brv32, Brv33 ;
  float* P1real = ealloc1float(nxfft);
  float* P1imag = ealloc1float(nxfft);
  float* P2real = ealloc1float(nxfft); 
  float* P2imag = ealloc1float(nxfft);
  complex* P1 = ealloc1complex(nxfft+1);
  complex* P2 = ealloc1complex(nxfft+1); 
    
  // float**   example1 = ealloc2float(ntfft,nxfft);
  // complex** example2 = ealloc2complex(ntfft,nxfftq);


  // define the basis for the fast quaternion Fourier transform
  Bfw11 = 1/sqrt(3); // the first vector of the new basis [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]
  Bfw12 = 1/sqrt(3); 
  Bfw13 = 1/sqrt(3);
  Bfw21 = 0;         // the second vector of the new basis [0 1/sqrt(2) -1/sqrt(2)] (orthogonal to the first)
  Bfw22 = 1/sqrt(2);
  Bfw23 = -1/sqrt(2);
  Bfw31 = -1*(Bfw12*Bfw23 - Bfw13*Bfw22); // the third vector of the new basis, -1*(row1Xrow2) = orthogonal to the first two by definition
  Bfw32 = -1*(Bfw13*Bfw21 - Bfw11*Bfw23);
  Bfw33 = -1*(Bfw11*Bfw22 - Bfw12*Bfw21);
  Brv11 = Bfw11;    // B_rv is the transpose of B_fw
  Brv12 = Bfw21;
  Brv13 = Bfw31;
  Brv21 = Bfw12;
  Brv22 = Bfw22;
  Brv23 = Bfw32;
  Brv31 = Bfw13;
  Brv32 = Bfw23;
  Brv33 = Bfw33; 

  // change of basis
  for (ix=0;ix<nxfft;ix++){ 
    P1real[ix] = pfft1[ix];
    P1imag[ix] = pfft2[ix]*Bfw11 + pfft3[ix]*Bfw12 + pfft4[ix]*Bfw13;
    P2real[ix] = pfft2[ix]*Bfw21 + pfft3[ix]*Bfw22 + pfft4[ix]*Bfw23;
    P2imag[ix] = pfft2[ix]*Bfw31 + pfft3[ix]*Bfw32 + pfft4[ix]*Bfw33;     
    P1[ix].r = P1real[ix];  
    P1[ix].i = P1imag[ix];
    P2[ix].r = P2real[ix];  
    P2[ix].i = P2imag[ix];
  }
  
  
  if (sign == 1){ // forward transform
    pfacc(-1,nxfft,P1);  // FFT x to k
    pfacc(-1,nxfft,P2);  // FFT x to k
  }
  else{           // reverse transform
    pfacc(1,nxfft,P1);  // FFT k to x
    pfacc(1,nxfft,P2);  // FFT k to x
  }
  for (ix=0;ix<nxfft;ix++){ 
    P1real[ix] = P1[ix].r;
    P1imag[ix] = P1[ix].i;
    P2real[ix] = P2[ix].r;
    P2imag[ix] = P2[ix].i;
    pfft1[ix] = P1real[ix];
    pfft2[ix] = P1imag[ix]*Brv11 + P2real[ix]*Brv12 + P2imag[ix]*Brv13;
    pfft3[ix] = P1imag[ix]*Brv21 + P2real[ix]*Brv22 + P2imag[ix]*Brv23;
    pfft4[ix] = P1imag[ix]*Brv31 + P2real[ix]*Brv32 + P2imag[ix]*Brv33;
  }

  //MARK; // this tells you where an error is occuring (to solve a segmentation fault)
  return;
}


void qfft2d(float sign, int ntfft, int nxfft, float** pfft1, float** pfft2, float** pfft3, float** pfft4, float** qpfft1, float** qpfft2, float** qpfft3, float** qpfft4){
  // The forward and reverse quaternion F-K Fourier transform along the time and space axes.
  
  int ix, it;
  float Bfw11, Bfw12, Bfw13, Bfw21, Bfw22, Bfw23, Bfw31, Bfw32, Bfw33 ;
  float Brv11, Brv12, Brv13, Brv21, Brv22, Brv23, Brv31, Brv32, Brv33 ;
  float** P1real = ealloc2float(ntfft,nxfft);
  float** P1imag = ealloc2float(ntfft,nxfft);
  float** P2real = ealloc2float(ntfft,nxfft); 
  float** P2imag = ealloc2float(ntfft,nxfft);
  complex** P1 = ealloc2complex(ntfft,nxfft+1);
  complex** P2 = ealloc2complex(ntfft,nxfft+1); 
    
  // float**   example1 = ealloc2float(ntfft,nxfft);
  // complex** example2 = ealloc2complex(ntfft,nxfftq);


  // define the basis for the fast quaternion Fourier transform
  Bfw11 = 1/sqrt(3); // the first vector of the new basis [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]
  Bfw12 = 1/sqrt(3); 
  Bfw13 = 1/sqrt(3);
  Bfw21 = 0;         // the second vector of the new basis [0 1/sqrt(2) -1/sqrt(2)] (orthogonal to the first)
  Bfw22 = 1/sqrt(2);
  Bfw23 = -1/sqrt(2);
  Bfw31 = -1*(Bfw12*Bfw23 - Bfw13*Bfw22); // the third vector of the new basis, -1*(row1Xrow2) = orthogonal to the first two by definition
  Bfw32 = -1*(Bfw13*Bfw21 - Bfw11*Bfw23);
  Bfw33 = -1*(Bfw11*Bfw22 - Bfw12*Bfw21);
  Brv11 = Bfw11;    // B_rv is the transpose of B_fw
  Brv12 = Bfw21;
  Brv13 = Bfw31;
  Brv21 = Bfw12;
  Brv22 = Bfw22;
  Brv23 = Bfw32;
  Brv31 = Bfw13;
  Brv32 = Bfw23;
  Brv33 = Bfw33; 

  // change of basis
  for (ix=0;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      P1real[ix][it] = pfft1[ix][it];
      P1imag[ix][it] = pfft2[ix][it]*Bfw11 + pfft3[ix][it]*Bfw12 + pfft4[ix][it]*Bfw13;
      P2real[ix][it] = pfft2[ix][it]*Bfw21 + pfft3[ix][it]*Bfw22 + pfft4[ix][it]*Bfw23;
      P2imag[ix][it] = pfft2[ix][it]*Bfw31 + pfft3[ix][it]*Bfw32 + pfft4[ix][it]*Bfw33;     
      P1[ix][it].r = P1real[ix][it];  
      P1[ix][it].i = P1imag[ix][it];
      P2[ix][it].r = P2real[ix][it];  
      P2[ix][it].i = P2imag[ix][it];
    }
  }

  
  if (sign == 1){ // forward transform
    pfa2cc(1,1,ntfft,nxfft,P1[0]);  /* Fourier transform t to w */
    pfa2cc(-1,2,ntfft,nxfft,P1[0]);  /* Fourier transform wx to wk */
    pfa2cc(1,1,ntfft,nxfft,P2[0]);  /* Fourier transform t to w */
    pfa2cc(-1,2,ntfft,nxfft,P2[0]);  /* Fourier transform wx to wk */
  }
  else{           // reverse transform
    pfa2cc(-1,1,ntfft,nxfft,P1[0]);  /* Fourier transform t to w */
    pfa2cc(1,2,ntfft,nxfft,P1[0]);  /* Fourier transform wx to wk */
    pfa2cc(-1,1,ntfft,nxfft,P2[0]);  /* Fourier transform t to w */
    pfa2cc(1,2,ntfft,nxfft,P2[0]);  /* Fourier transform wx to wk */
  }
  for (ix=0;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      P1real[ix][it] = P1[ix][it].r;
      P1imag[ix][it] = P1[ix][it].i;
      P2real[ix][it] = P2[ix][it].r;
      P2imag[ix][it] = P2[ix][it].i;
      qpfft1[ix][it] = P1real[ix][it];
      qpfft2[ix][it] = P1imag[ix][it]*Brv11 + P2real[ix][it]*Brv12 + P2imag[ix][it]*Brv13;
      qpfft3[ix][it] = P1imag[ix][it]*Brv21 + P2real[ix][it]*Brv22 + P2imag[ix][it]*Brv23;
      qpfft4[ix][it] = P1imag[ix][it]*Brv31 + P2real[ix][it]*Brv32 + P2imag[ix][it]*Brv33;
    }
  }
  return;
}

float qamp(float q1, float q2, float q3, float q4){
  // The quaternion amplitude 
  float Aq;
  Aq = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4);
  return Aq;
}

void qamp2d(int ntfft, int nxfft, float** q1, float** q2, float** q3, float** q4, float** Aq){
  // The quaternion amplitude for a matrix of values 
  int ix, it;
  for (ix=0;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      Aq[ix][it] = sqrt(q1[ix][it]*q1[ix][it] + q2[ix][it]*q2[ix][it] + q3[ix][it]*q3[ix][it] + q4[ix][it]*q4[ix][it]);
    }
  }
  return;
}

void fftshift(int ntfft, int nxfft, float** A, float** Ashift){
  // fftshift the spectrum in the wavenumber direction 
  int ix, it;
  for (ix=0;ix<nxfft/2+1;ix++){ 
    for(it=0;it<ntfft+1;it++){ 
      Ashift[ix][it] = A[nxfft/2 - ix][it];
    }
  }
  for (ix=nxfft/2+1;ix<nxfft;ix++){ 
    for(it=0;it<ntfft;it++){ 
      Ashift[ix][it] = A[nxfft-ix+nxfft/2][it];
    }
  }
  return;
}

void fftshift2(int ntfft, int nxfft, float** A, float** Ashift){
  // fftshift the spectrum in the frequency and wavenumber dimensions 
  int ix, it;

  for (it=0;it<ntfft/2+1;it++){ 
    for(ix=0;ix<nxfft/2+1;ix++){ 
      Ashift[ix][it] = A[nxfft/2 - ix][ntfft/2 - it];
    }
    for(ix=nxfft/2 + 1;ix<nxfft;ix++){ 
      Ashift[ix][it] = A[nxfft - ix + nxfft/2][ntfft/2 - it];
    }
  }
  for (it=ntfft/2 + 1;it<ntfft+1;it++){ 
    for(ix=0;ix<nxfft/2+1;ix++){ 
      Ashift[ix][it] = A[nxfft/2 - ix][ntfft - it + ntfft/2];
    }
    for(ix=nxfft/2 + 1;ix<nxfft;ix++){ 
      Ashift[ix][it] = A[nxfft - ix + nxfft/2][ntfft - it + ntfft/2];
    }
  }
  return;
}

