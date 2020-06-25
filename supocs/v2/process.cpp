#include "Complex.h"
#include "su.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void xw2kw(complex** cpfft, int nxfft, int nk, int nw, complex* freqslice, float sign);
void save_gather(float **d, int nh, int nt, float dt, const char* name);
void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);

bool add(float** data1, float** data2, int nt, int nh, float value=1);
bool copy(float** data1, float** data2, int nt, int nh, float value=1);
int countNonZero(complex** model, int nt, int nh);

void xplotgather(complex **d, int nh, int nt, float dt, float dx, const char *s, int num, const char *s2);

// bool process(float** datain, float** dataout, int nt, int  nx, float dt){
bool process(float** datain, float** dataout, int nt, int  nx, float dt){
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
  float alpha = 0.8; // denoising parameter, 1 => no noise , ->0 => noisy input traces.
  bool plotInput   = true;
  bool plotInputFK = false;
  bool plotInput2  = true; // after zeroes
  bool plotFXIn    = false;
  bool plotFKIn    = false; // after zeroes
  bool plotFKOut   = false;
  bool plotOutput  = true;

  int padfactor = 2;      // pad factor for both dimensions


  // plot input
  if (plotInput){
    save_gather(datain,nx,nt,dt,"input.su");
    system("suximage < input.su perc=98 clip=1.0 legend=1 title=\"input\" &");
  }

  /* copy data from input to FFT array and pad with zeros */
  
  ntfft = npfar(padfactor*nt);
  nxfft = npfa(padfactor*nx);
  nw=ntfft/2+1;
  nk=nxfft;
  scale=nxfft*ntfft;
  df=1./(ntfft*dt);
  dk=1./(nxfft*dx);

  complex* freqslice = ealloc1complex(nxfft);

  // float** datatrue = 0; //declare outside if you need this variable in another place
  // plot input in FK
  if (plotInputFK){
    float** datatrue  = ealloc2float(ntfft,nxfft); // trace oriented (Hale's reversed convention for alloc)
    for (ix=0;ix<nx;ix++) for(it=0;it<nt;it++) datatrue[ix][it] = datain[ix][it];
    for (ix=nx;ix<nxfft;ix++) for(it=nt;it<ntfft;it++) datatrue[ix][it] = 0;
    complex** datatrueFK = ealloc2complex(nw,nxfft);
    pfa2rc(1,1,ntfft,nx,datatrue[0],datatrueFK[0]); // now it is in FX
    xw2kw(datatrueFK,nxfft,nk,nw,freqslice,-1); // now it is in FK
    xplotgather(datatrueFK,nk,nw,df,dk,"datatruefk.su",1,"perc=99 wclip=0 bclip=1000 title=\"input FK\" &"); // plot [k,w]
    //dealloc2float(datatrue); // if you declare this variable in the if statement, then it doesn't exist outside and you should deallocate the memory
  }


  // zero some traces to test interpolation 
  //  for (ix=0;ix< nx;ix+=3) for (it=0;it<nt;it++) datain[ix][it]=0;
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
  

  fprintf(stderr,"nt=%d, nx=%d, ntfft=%d, nxfft=%d, nw=%d, nk=%d\n",nt,nx,ntfft,nxfft,nw,nk);

  float**   pfft  = ealloc2float(ntfft,nxfft); // trace oriented (Hale's reversed convention for alloc)
  complex** cpfft = ealloc2complex(nw,nk);     // trace oriented
  //  complex** cpfft  = ealloc2complex(nk,nw); // freq oriented (not used but may be better)

  
  /* copy data from input to FFT array and pad with zeros */
  for (ix=0;ix<nx;ix++){
    for (it=0; it<nt; it++) pfft[ix][it]=datain[ix][it];
    for (it=nt; it< ntfft;it++) pfft[ix][it] = 0.0;
  }
  for (ix=nx;ix<nxfft;ix++) for(it=0;it<ntfft;it++) pfft[ix][it] = 0.0;

  /* Fourier transform t to w */
  pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]); // cpfft is [x,w]


  /* algorithm starts*/
  complex* freqslice2= ealloc1complex(nxfft);
  
  float* mabs = ealloc1float(nxfft);  
  float* mabsiter = ealloc1float(nxfft);
  float* wd   = ealloc1float(nxfft);
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
    
  if (plotFXIn){
    xplotgather(cpfft,nk,nw,df,dx,"datafxin.su",1,"perc=99 wclip=0 bclip=50 &"); // plot [x,w]
  }      
  if (plotFKIn){
    xw2kw(cpfft,nxfft,nk,nw,freqslice,-1);                      // x to k
    xplotgather(cpfft,nk,nw,df,dk,"datafkin.su",1,"perc=99 wclip=0 bclip=1000 title=\"after zeros FK\" &"); // plot [k,w]
    xw2kw(cpfft,nxfft,nk,nw,freqslice,1);                       // k to x
  }      

  // freq by freq algorithm
  for (iw=0;iw<nw;iw++){
    
    // copy freq slice (needs because arrays are trace oriented)
    for (ix=0;ix<nk;ix++) freqslice[ix]=freqslice2[ix]=cpfft[ix][iw];

    pfacc(-1,nxfft,freqslice2);  // FFT x to k
    // threshold in k
    for (ix=0;ix<nk;ix++) mabs[ix]=abs(freqslice2[ix]);
    pfacc(1,nxfft,freqslice2);  // iFFT k to x
    for (ix=0;ix<nk;ix++) freqslice2[ix]/=nxfft;
    for (int iter=1;iter<Iter;iter++){  // external thresholding loop
      pfacc(-1,nxfft,freqslice2);  // FFT x to k
	
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
      pfacc(1,nxfft,freqslice2);  // k to x
      for (ix=0;ix<nk;ix++) freqslice2[ix]/=nxfft;

      // add new traces to original freqslice
      //if (iter2 != 5) 
      for (ix=0;ix<nk;ix++) freqslice2[ix]=(alpha*freqslice[ix]) + (1-alpha*wd[ix])*freqslice2[ix]; // x,w

    }

    // put back the slice
    for (ix=0;ix<nk;ix++) cpfft[ix][iw]=freqslice2[ix];
      
  }
  
  if (plotFKOut){   // plot final spectra
    xw2kw(cpfft,nxfft,nk,nw,freqslice,-1); // x to k
    xplotgather(cpfft,nk,nw,df,dk,"datafkout.su",1,"perc=99 wclip=0 bclip=1000 title=\"pocs FK\" &");
    xw2kw(cpfft,nxfft,nk,nw,freqslice,1);                       // k to x
  }

  free1float(wd);
  free1complex(freqslice2);
  free1complex(freqslice);
  free1float(mabs);
  // algorithm ends

  /* Fourier transform w to t */
  pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
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


void xw2kw(complex** cpfft, int nxfft, int nk, int nw, complex* freqslice, float sign){
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

void wx2wk(complex** cpfft, int nxfft, int nk, int nw, float sign){
  int ix, iw;
  
  for (iw=0;iw<nw;iw++){
    pfacc(sign,nxfft,cpfft[iw]);
    if (sign == 1) for (ix=0;ix<nxfft;ix++) cpfft[iw][ix] /= nxfft;
  }
  return;
}

