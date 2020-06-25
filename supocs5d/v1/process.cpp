#include "Complex.h"
#include "su.h"
#include "fftw3.h"

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void xw2kw(complex** cpfft, int nxfft, int nk, int nw, complex* freqslice, int sign);
void save_gather(float **d, int nh, int nt, float dt, const char* name);
void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);

bool add(float** data1, float** data2, int nt, int nh, float value=1);
bool copy(float** data1, float** data2, int nt, int nh, float value=1);
int countNonZero(complex** model, int nt, int nh);

void xplotgather(complex **d, int nh, int nt, float dt, float dx, const char *s, int num, const char *s2);

// bool process(float** datain, float** dataout, int nt, int  nx, float dt){
bool process(float** datain, float** dataout, int nt, int  nx, float dt, int x1_num, int x2_num, int x3_num, int x4_num, int dx1, int dx2, int dx3, int dx4){
  
  MARK;
  // if the sorting is not cdp then offset then the only thing that needs to be changed is x1 and x2 are reversed.
  int nx1, nx2, nx3, nx4;
  nx1 = x1_num;
  nx2 = x2_num;
  nx3 = x3_num;
  nx4 = x4_num;

  int it, ix, ix1, ix2, ix3, ix4, iw;
  //  for (ih=0;ih<nh;ih++)  for (it=0;it<nt;it++) dataout[ih][it]=datain[ih][it];
  complex czero;czero.r=czero.i=0;
  int ntfft, nx1fft, nx2fft, nx3fft, nx4fft, nw, nk, nk1, nk2, nk3, nk4;
  float scale;

  // defining sampling for plots.
  float df, dk1, dk2, dk3, dk4; // sampling interval in fk
  int Iter = 20; // iterations for thresholding
  float perci = 0.999;
  float percf = 0.001;
  float alpha = 1.0; // denoising parameter, 1 => no noise , ->0 => noisy input traces.
  bool plotInput   = true;
  bool plotInput2  = true; // after zeroes
  bool plotOutput  = true;

  int padfactor = 2;      // pad factor for both dimensions


  fprintf(stderr,"Parameters used: \nIter=%d, perci=%6.3f, percf=%6.3f, alpha=%6.3f \n",Iter,perci,percf,alpha);


  // plot input
  if (plotInput){
    save_gather(datain,nx,nt,dt,"input.su");
    system("suximage < input.su perc=98 clip=1.0 legend=1 title=\"input\" &");
  }

  /* copy data from input to FFT array and pad with zeros */
  
  ntfft = npfar(padfactor*nt);
  nx1fft = npfa(padfactor*nx1);
  nx2fft = npfa(padfactor*nx2);
  nx3fft = npfa(padfactor*nx3);
  nx4fft = npfa(padfactor*nx4);
  nw=ntfft/2+1;
  nk=nx1fft*nx2fft*nx3fft*nx4fft;
  nk1=nx1fft;
  nk2=nx2fft;
  nk3=nx3fft;
  nk4=nx4fft;
  scale=nx1fft*nx2fft*nx3fft*nx4fft*ntfft;
  df=1./(ntfft*dt);
  dk1=1./(nx1fft*dx1);
  dk2=1./(nx2fft*dx2);
  dk3=1./(nx3fft*dx3);
  dk4=1./(nx4fft*dx4);

  complex* freqslice = ealloc1complex(nx1fft*nx2fft*nx3fft*nx4fft);

  // zero some traces to test interpolation 
   float rndnum = 0;
   for (ix=0;ix< nx;ix+=1){
     // generate random number distributed between 0 and 1 
     // (== Uniformly distributed == matlab function "rand")
     rndnum = franuni();
     if (rndnum>0.5){
     for (it=0;it<nt;it++) datain[ix][it]=0;
     }
   } 
  // plot after zeros
  if (plotInput2){
    save_gather(datain,nx,nt,dt,"input2.su");
    system("suximage < input2.su title=\"after zeros\" legend=1 perc=98 clip=1.0 &");
  }
  

  fprintf(stderr,"nt=%d, nx=%d, ntfft=%d, nx1fft=%d, nx2fft=%d, nx3fft=%d, nx4fft=%d, nw=%d, nk1=%d, nk2=%d\n",nt,nx,ntfft,nx1fft,nx2fft,nx3fft,nx4fft,nw,nk1,nk2);

  float**   pfft  = ealloc2float(ntfft,nx1fft*nx2fft*nx3fft*nx4fft); // trace oriented (Hale's reversed convention for alloc)
  complex** cpfft = ealloc2complex(nw,nk1*nk2*nk3*nk4);     // trace oriented
  //  complex** cpfft  = ealloc2complex(nk,nw); // freq oriented (not used but may be better)

  
  /* copy data from input to FFT array and pad with zeros in time dimension*/
  for (ix=0;ix<nx;ix++){
    for (it=0; it<nt; it++) pfft[ix][it]=datain[ix][it];
    for (it=nt; it< ntfft;it++) pfft[ix][it] = 0.0;
  }

//******************************************************************************************** TX to FX
  // transform 3d data from t-x to w-x using FFTW
  int N = ntfft; 
  complex* out = ealloc1complex(nw);
  fftwf_plan p1;
  float* in = ealloc1float(N);
  p1 = fftwf_plan_dft_r2c_1d(N, in, (fftwf_complex*)out, FFTW_ESTIMATE);
  complex** datain_wx = ealloc2complex(nw,nx);
  for (ix=0;ix<nx;ix++){
    for(it=0;it<ntfft;it++){
      in[it] = pfft[ix][it];
    }
    fftwf_execute(p1); /* take the FFT along the time dimension */
    for(iw=0;iw<nw;iw++){
      cpfft[ix][iw] = out[iw]; 
    }
  }
  fftwf_destroy_plan(p1);
  fftwf_free(in); fftwf_free(out);
  //  xplotgather(cpfft,nx,nw,df,dx,"cpfft.su",1,"perc=99 wclip=0 bclip=50 title=\"input FX done with pfa2rc!\" &"); //plot [x,w]
//********************************************************************************************


  /* pad along the 2 spatial dimensions */
  ix = 0;
  for (ix1=0;ix1<nx1fft;ix1++){ 
    for (ix2=0;ix2<nx2fft;ix2++){
      for (ix3=0;ix3<nx3fft;ix3++){
	for (ix4=0;ix4<nx4fft;ix4++){
	  if (ix1>=nx1 && ix2>=nx2 && ix3>=nx3  && ix4>=nx4){
	    for (iw=0;iw< nw;iw++) cpfft[ix][iw] = 0.0;
	  }
	  ix++;
	}
      }
    }
  }

  /* algorithm starts*/
  complex* freqslice2= ealloc1complex(nx1fft*nx2fft*nx3fft*nx4fft);
  
  float* mabs = ealloc1float(nx1fft*nx2fft*nx3fft*nx4fft);  
  float* mabsiter = ealloc1float(nx1fft*nx2fft*nx3fft*nx4fft);
  float* wd   = ealloc1float(nx1fft*nx2fft*nx3fft*nx4fft);
  float sigma;

  for (ix=0; ix<nk;ix++){
    float sum = 0;
    for (iw = 0; iw< nw; iw++) 
      sum += abs(cpfft[ix][iw]);
    if (sum) wd[ix]=1;
    else{ 
      wd[ix]=0;
      //fprintf(stderr,"wd[%d]=%f \n",ix,wd[ix]);
    }
  }
    
//******************************************************************************************** FX1X2 to FK1K2
  // make the plan that will be used for each frequency slice
  // written as a 4D transform with length =1 for two of the dimensions. 
  // This is to make it easier to upgrade to reconstruction of 4 spatial dimensions. 
  int rank = 4;
  int n[4];
  n[0] = nx1fft;
  n[1] = nx2fft;
  n[2] = nx3fft;
  n[3] = nx4fft;
  N = nx1fft*nx2fft*nx3fft*nx4fft;
  fftwf_plan p2;
  p2 = fftwf_plan_dft(rank, n, (fftwf_complex*)freqslice2, (fftwf_complex*)freqslice2, FFTW_FORWARD, FFTW_ESTIMATE);
//********************************************************************************************

//******************************************************************************************** FK1K2 to FX1X2
  // make the plan that will be used for each frequency slice
  // written as a 4D transform with length =1 for two of the dimensions. 
  // This is to make it easier to upgrade to reconstruction of 4 spatial dimensions. 
  fftwf_plan p3;
  p3 = fftwf_plan_dft(rank, n, (fftwf_complex*)freqslice2, (fftwf_complex*)freqslice2, FFTW_BACKWARD, FFTW_ESTIMATE);
//********************************************************************************************

  // freq by freq algorithm
  for (iw=0;iw<nw;iw++){
    
    // copy freq slice (needs because arrays are trace oriented)
    for (ix=0;ix<nk;ix++) freqslice[ix]=freqslice2[ix]=cpfft[ix][iw];

    fftwf_execute(p2); // FFT x to k
    //pfacc(-1,nxfft,freqslice2);  // FFT x to k


    // threshold in k
    for (ix=0;ix<nk;ix++) mabs[ix]=abs(freqslice2[ix]);
    fftwf_execute(p3); // FFT k to x
    //pfacc(1,nxfft,freqslice2);  // iFFT k to x
    for (ix=0;ix<nk;ix++) freqslice2[ix]/=nx1fft*nx2fft*nx3fft*nx4fft;
    for (int iter=1;iter<Iter;iter++){  // external thresholding loop
      	fftwf_execute(p2); // FFT x to k
	//pfacc(-1,nxfft,freqslice2);  // FFT x to k
	
	int count = 0; // count kept
	int kept = 0;  // last index kept
        
        // This is to increase the thresholding within each internal iteration
        sigma=quest(perci - (iter-1)*((perci-percf)/(Iter-1)),nk,mabs); 
	for (ix=0;ix<nk;ix++) mabsiter[ix]=abs(freqslice2[ix]);
        for (ix=0;ix<nk;ix++){
	  if (mabsiter[ix]<sigma) freqslice2[ix] = czero;
	  else{ kept=ix;count++; }
	}

	//fprintf(stderr,"sigma=%f count = %d kept=%d \n",sigma,count,kept);
        fftwf_execute(p3); // FFT k to x
	//pfacc(1,nxfft,freqslice2);  // k to x
	for (ix=0;ix<nk;ix++) freqslice2[ix]/=nx1fft*nx2fft*nx3fft*nx4fft;

	// add new traces to original freqslice
	//if (iter2 != 5) 
	for (ix=0;ix<nk;ix++) freqslice2[ix]=(alpha*freqslice[ix]) + (1-alpha*wd[ix])*freqslice2[ix]; // x,w

    }

//     // put back the slice
//     for (ix=0;ix<nk;ix++) cpfft[ix][iw]=freqslice2[ix];
 

    /* put back the slice and remove padding along the 2 spatial dimensions */
    ix = 0;
    for (ix1=0;ix1<nx1fft;ix1++){ 
      for (ix2=0;ix2<nx2fft;ix2++){
	for (ix3=0;ix3<nx3fft;ix3++){
	  for (ix4=0;ix4<nx4fft;ix4++){
	    if (ix1<nx1 && ix2<nx2 && ix3<nx3 && ix4<nx4){
	      cpfft[ix][iw] = freqslice2[ix];
	      ix++;
	    }
	  }
	}
      }
    }


  
  }
  
  fftwf_destroy_plan(p2);
  fftwf_destroy_plan(p3);

  free1float(wd);
  free1complex(freqslice2);
  free1complex(freqslice);
  free1float(mabs);
  // algorithm ends



//******************************************************************************************** FX to TX
  // transform 2d data from w-x to t-x using IFFTW
  N = ntfft; 
  float* out2 = ealloc1float(ntfft);
  fftwf_plan p4;
  complex* in2 = ealloc1complex(N);
  p4 = fftwf_plan_dft_c2r_1d(N, (fftwf_complex*)in2, out2, FFTW_ESTIMATE);
  for (ix=0;ix<nx;ix++){
    for(iw=0;iw<nw;iw++){
      in2[iw] = cpfft[ix][iw];
    }
    fftwf_execute(p4); /* take the FFT along the time dimension */
    for(it=0;it<nt;it++){
      pfft[ix][it] = out2[it]; 
    }
  }
  fftwf_destroy_plan(p4);
  fftwf_free(in2); fftwf_free(out2);
//********************************************************************************************

//   /* Fourier transform w to t */
//   pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);

  for (ix=0;ix<nx;ix++) for (it=0; it<nt; it++) dataout[ix][it]=pfft[ix][it]/ntfft;
  

  save_gather(dataout,nx,nt,dt,"output.su");

  // plot after pocs
  if (plotOutput){
    save_gather(dataout,nx,nt,dt,"output.su");
    system("suximage < output.su title=\"after pocs\" legend=1 perc=98 clip=1.0 &");
  }

  free2float(pfft);
  free2complex(cpfft);
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

