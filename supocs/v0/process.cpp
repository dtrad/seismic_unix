#include "Complex.h"
#include "su.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);

bool add(float** data1, float** data2, int nt, int nh, float value=1);
bool copy(float** data1, float** data2, int nt, int nh, float value=1);
int countNonZero(complex** model, int nt, int nh);

void xplotgather(complex **d, int nh, int nt, float dt, float dx, char *s, int num, char *s2);

bool process(float** datain, float** dataout, int nt, int  nx, float dt){
  int it, ix, iw, sign;
  //  for (ih=0;ih<nh;ih++)  for (it=0;it<nt;it++) dataout[ih][it]=datain[ih][it];
  complex czero; 
  int ntfft, nxfft, nw, nk;
  float scale;

  // defining sampling for plots.
  float df, dk; // sampling interval in fk
  float dx = 100; // sampling interval in offset
  
  ntfft = npfar(nt);
  nxfft = npfa(nx);
  nw=ntfft/2+1;
  nk=nxfft;
  scale=nxfft*ntfft;
  df=1./(ntfft*dt);
  dk=1./(nxfft*dx);

  fprintf(stderr,"nt=%d, nx=%d, ntfft=%d, nxfft=%d, nw=%d, nk=%d\n",nt,nx,ntfft,nxfft,nw,nk);

  float**   pfft  = ealloc2float(ntfft,nxfft);
  complex** cpfft = ealloc2complex(nw,nk);
  complex phaseShift;phaseShift.r=phaseShift.i=0;
  
  /* copy data from input to FFT array and pad with zeros */
  for (ix=0;ix<nx;ix++){
    /* if ix odd, negate to center transform of dimension 2 */
    if (ISODD(ix)) sign=-1; else sign=1;
    for (it=0; it<nt; it++) pfft[ix][it]=sign*datain[ix][it];
    for (it=nt; it< ntfft;it++) pfft[ix][it] = 0.0;
  }
  for (ix=nx;ix<nxfft;ix++) for(it=0;it<ntfft;it++) pfft[ix][it] = 0.0;

  /* Fourier transform t to w */
  pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]);

  if (0){ // one pass 2D FK
    fprintf(stderr,"ONE pass inv fft\n");
    pfa2cc(-1,2,nw,nxfft,cpfft[0]);
    xplotgather(cpfft,nk,nw,df,dk,"datafk.su",1,"perc=99 &");  

    for (ix=0;ix<nk;ix++){
      for (iw=0;iw<nw;iw++){
	cpfft[ix][iw].r+=phaseShift.r;
	cpfft[ix][iw].i+=ix;
      }
    }
    xplotgather(cpfft,nk,nw,df,dk,"datafkshifted.su",1,"perc=99 &");  
    pfa2cc(1,2,nw,nxfft,cpfft[0]);
  }
  else if (0){
    complex* freqslice = ealloc1complex(nxfft);
    for (iw=0;iw<nw;iw++){
      for (ix=0;ix<nk;ix++) freqslice[ix]=cpfft[ix][iw];
      pfacc(-1,nxfft,freqslice);
      for (ix=0;ix<nk;ix++) cpfft[ix][iw]=freqslice[ix];
    }
    xplotgather(cpfft,nk,nw,df,dk,"datafk.su",1,"perc=99 &");
    for (iw=0;iw<nw;iw++){
      for (ix=0;ix<nk;ix++) freqslice[ix]=cpfft[ix][iw];
      pfacc(1,nxfft,freqslice);
      for (ix=0;ix<nk;ix++) cpfft[ix][iw]=freqslice[ix];
    }
    free1complex(freqslice);
  }
  else{ // one call per frequency
    /* algorithm starts*/
    complex* freqslice = ealloc1complex(nxfft);
    complex* freqslice2= ealloc1complex(nxfft);
    
    float* mabs = ealloc1float(nxfft);
    float* wd   = ealloc1float(nxfft);
    for (ix=0;ix<nxfft;ix++) wd[ix] = 1;
    //for (ix=0;ix<nxfft;ix+=3) wd[ix] = 0; // these traces are zero.
    
    float sigma;
    complex czero; czero.r=czero.i=0;
    for (iw=0;iw<nw;iw++){
      
      // copy freq slice
      for (ix=0;ix<nk;ix++) freqslice[ix]=freqslice2[ix]=cpfft[ix][iw];
      
      // iteration inside one freq slice
      for (int iter=0;iter<0;iter++){
	/* Fourier transform x to k */
	pfacc(-1,nxfft,freqslice2);
	
	// threshold
	for (ix=0;ix<nk;ix++) mabs[ix]=abs(freqslice2[ix]);
	sigma=quest(0.9,nk,mabs);
	for (ix=0;ix<nk;ix++) if (mabs[ix]<sigma) freqslice2[ix]=czero;
	
	/* Fourier transform k to x */
	pfacc(1,nxfft,freqslice2);
	// subtract
	for (ix=0;ix<nk;ix++) freqslice2[ix]=wd[ix]*freqslice[ix]-freqslice2[ix];
	
	
      }
      
      /* Fourier transform x to k */
      pfacc(-1,nxfft,freqslice2);

      // put back the slice
      for (ix=0;ix<nk;ix++) if (wd[ix]) cpfft[ix][iw]=freqslice2[ix];
      
    }

    
    /* thresholding and subtraction */
    xplotgather(cpfft,nk,nw,df,dk,"datafk.su",1,"perc=99 &");
    
    for (iw=0;iw<nw;iw++){
      for (ix=0;ix<nk;ix++) freqslice[ix]=cpfft[ix][iw];
      /* Fourier transform k to x */
      pfacc(1,nxfft,freqslice);
      for (ix=0;ix<nk;ix++) cpfft[ix][iw]=freqslice[ix];
    }
    free1float(wd);
    free1complex(freqslice2);
    free1complex(freqslice);
    free1float(mabs);
  }
  // algorithm ends

  /* Fourier transform w to t */
  pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
  

  for (ix=0;ix<nx;ix++){
    /* if ix odd, pfft has been negated to center transform of dimension 2 */
    if (ISODD(ix)) sign=-1; else sign=1; 
    for (it=0; it<nt; it++) dataout[ix][it]=sign*pfft[ix][it]/scale;
  }


  free2float(pfft);
  free2complex(cpfft);
  
  return true;

}

bool process2(float** datain, float** dataout, int nt, int  nh, float dt){
  int ih;

  int ifreq, maxfreq, nfreq, nf0;
  float threshold;

  float fmax=0, df;
  complex czero; 
  complex** model=ealloc2complex(nt,nh);
  
  czero.r = czero.i=0;
  // first fft get back nf0
  fftgo0(-1,datain ,model,nh,nt,dt,&nf0);
  
// use nf0 to calculate number of freqs (first time).  
  nfreq=nf0/2;
  df=1/(nf0*dt);
  if (fmax==0) maxfreq=nfreq;
  else         maxfreq=(int) (fmax/df);

  //calculate threshold
  threshold=0;
  for (ifreq=10;ifreq<maxfreq;ifreq++) 
    for (ih=0;ih<nh;ih++)
      if (abs(model[ifreq][ih]) > threshold)  threshold=abs(model[ifreq][ih]);
  
  if (1){
    threshold*=1;
    // thresholding
    for (ifreq=0;ifreq<maxfreq;ifreq++) 
      for (ih=0;ih<nh;ih++)
	if (abs(model[ifreq][ih]) < threshold)  model[ifreq][ih]= czero;
  }

  fprintf(stderr,"non zero = %d\n",countNonZero(model,maxfreq,nh));

  // back to time
  fftback0(1,dataout ,model,nh,nt,dt,nf0);
  
  // subtract from datain
  add(datain,dataout,nt,nh,-1.);

  // put datain into output
  copy(dataout,datain,nt,nh,1.);
  
  free2complex(model);

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
