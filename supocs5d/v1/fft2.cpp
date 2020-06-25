#include "su.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#include "Complex.h"

void fksize(int nt, int nx, int *pnw, int *pnk)
{
  *pnw = npfar(nt)/2+1;
  *pnk = npfa(nx);
}


void fft2_parameters(int nt, float dt, int nx, float dx, int *pnw, float *pdw,
		     int *pnk, float *pdk, float *w, float *k)
{
  /********** fft_parameters *************
    Inputs: Time and space parameters
       nt
       dt
       nx
       dx
     Outputs: Freq parameters
     *ntfft  (address of number of freq) 
     *dw  (address of dw)
     *nxfft  (address of number of freq) 
     *dk  (address of dw)
     w[nw]  freq axis
     k[nk]  wavenumber axis

  Daniel Trad - Dec 2000   
  ***************************************/
  int ntfft;
  int nxfft;
  int nw;
  int nk;
  float fw;
  float fk;
  float dw;
  float dk;
  int iw;
  int ik;
  /* determine lengths and scale factors for prime-factor FFTs */
	
  ntfft = npfar(nt);
  nxfft = npfa(nx);
	
  /* determine frequency and wavenumber sampling */
  *pnw = nw = ntfft/2+1;
  *pdw = dw = 1.0/(ntfft*dt);

  fw = 0.0;

  *pnk = nk = nxfft;
  *pdk = dk = 1.0/(nxfft*dx);

  fk = -1.0/(2*dx);

  for (k[0]=fk,ik=1; ik<nk; ik++) k[ik]=k[ik-1]+dk;
  for (w[0]=fw,iw=1; iw<nw; iw++) w[iw]=w[iw-1]+dw;


  return;
}


void fft2_xt2kf(float **dataxt, complex **datakf, int nt, int nx, int nw, int nk)
{
  /****************************************************************************
    void fft2_xt2fk(float **dataxt, complex **datafk, int nt, int nx,float *w, float *k,   
		int nw, float dw, float fw, int nk, float dk, float fk)
 
		Input :
		int nw;			  number of frequencies 
		float dw;		  frequency sampling interval 
		float fw;		  first frequency 
		int nk;			  number of wavenumbers 
		float dk;		  wavenumber sampling interval 
		float fk;                 first wavenumber 
  ***************************************************************************/
		
	int ntfft=2*(nw-1);		 /* nt after padding for FFT */
	int nxfft=nk;		         /* nx after padding for FFT */
	int it,ix,iw,ik;	/* sample indices */
	complex **cpfft;	/* complex FFT workspace */
	float **pfft;		/* float FFT workspace */
	float sign=1;

	/* allocate real and complex workspace for FFTs */
	cpfft = ealloc2complex(nw,nk);
	pfft = ealloc2float(ntfft,nxfft);

	/* copy data from input to FFT array and pad with zeros */
	for (ix=0; ix<nx; ix++) {
	  /* if ix odd, negate to center transform of dimension 2 */
	  if (ISODD(ix)) sign=-1; else sign=1; 
	  for (it=0; it<nt; it++) pfft[ix][it]=sign*dataxt[ix][it];
	  for (it=nt; it<ntfft; it++) pfft[ix][it] = 0.0;
	}
	for (ix=nx; ix<nxfft; ix++) for (it=0; it<ntfft; it++) pfft[ix][it] = 0.0;
	
	/* Fourier transform t to w */
	pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]);
	/* Fourier transform x to k */
	pfa2cc(-1,2,nw,nxfft,cpfft[0]);
	/* FK output */
	for (ik=0; ik<nxfft; ik++) 
	  for (iw=0; iw<nw; iw++) 
	    datakf[ik][iw]=cpfft[ik][iw];

	free2complex(cpfft);
	free2float(pfft);
	return;

}
void fft2_kf2xt(float **dataxt, complex **datakf, int nt, int nx, int nw, int nk)
{
  /****************************************************************************
    void fft2_xt2fk(float **dataxt, complex **datafk, int nt, int nx,float *w, float *k,   
		int nw, float dw, float fw, int nk, float dk, float fk)
 
		Input :
		int nt;			  number of time samples 
		int nk;			  number of offsets 
		int nw;			  number of frequencies 
		int nk;			  number of wavenumbers 
  ***************************************************************************/
		
	int ntfft=2*(nw-1);		 /* nt after padding for FFT */
	int nxfft=nk;		         /* nx after padding for FFT */
	//	float sfft= 1.0/(ntfft*nxfft);	 /* scale factor for FFT */
	int it,ix,iw,ik;	/* sample indices */
	complex **cpfft;	/* complex FFT workspace */
	float **pfft;		/* float FFT workspace */
	float sign=1;

	/* allocate real and complex workspace for FFTs */
	cpfft = ealloc2complex(nw,nk);
	pfft = ealloc2float(ntfft,nxfft);

	for (ik=0; ik<nk; ik++) 
	  for (iw=0; iw<nw; iw++) 
	    cpfft[ik][iw]=datakf[ik][iw];
	
	/* Fourier transform k to x */
	pfa2cc(1,2,nw,nxfft,cpfft[0]);
	/* Fourier transform w to t */
	pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
	/* copy filtered data from FFT array to output */
	for (ix=0; ix<nx; ix++){
	  /* if ix odd, pfft has been negated to center transform of dimension 2 */
	  if (ISODD(ix)) sign=-1; else sign=1; 
	  for (it=0; it<nt; it++) dataxt[ix][it]=sign*pfft[ix][it];
	}
	/* free workspace */
	free2complex(cpfft);
	free2float(pfft);

}







