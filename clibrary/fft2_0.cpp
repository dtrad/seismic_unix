#include "fft2.h"
//#include "clibrary.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

/*

Daniel Trad - June 9- 2000
*/


void fft2_0(float **data, complex **model, float *h, int nh, float *t, int nt, 
		    float *k, int nk, float *w, int nw, float *vel, float eps1, 
		    float smute, float nmofactor, float fmax, float **data2, float *h2, 
		    int nh2)
{
  int it, ih;
  float *dtemp;
  float dt=t[1]-t[0];
  float *htemp;

  dtemp=alloc1float(nt);
  htemp=ealloc1float(nh);

  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT
  // All methods are the same with the only difference that the offsets are
  // transformed to pseudo offsets.

  //////////////////////////////////////////////////////////////////////

  if (nmofactor)
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  
  fft2(data,model,t,h,nt,nh,w,k,&nw,&nk,data2);
  fprintf(stderr,"In hrfft2_interface pnw=%d, pnk=%d\n",nw,nk);

  /////////////////////////////
  if (nmofactor)
    for (ih=0;ih<nh2;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data2[ih][it];    
      nmo(data2[ih],dtemp,t,nmofactor*h2[ih],vel,1,nt,dt,smute);  
    }    
  

  TRACE;

  free1float(htemp);
  free1float(dtemp);

  return;

}

void fft2(float **dataxt, complex **datafk, float *t, float *x, int nt, int nx, 
	  float *w, float *k, int *pnw, int *pnk, float **data2xt)
{
	int ntfft;		/* nt after padding for FFT */
	int nxfft;		/* nx after padding for FFT */
	float sfft;		/* scale factor for FFT */
	int nw;			/* number of frequencies */
	float dw;		/* frequency sampling interval */
	float fw;		/* first frequency */
	int nk;			/* number of wavenumbers */
	float dk;		/* wavenumber sampling interval */
	int it,ix,iw,ik;	/* sample indices */
	complex **cpfft;	/* complex FFT workspace */
	float **pfft;		/* float FFT workspace */
	float dt=t[1]-t[0];
	float dx=x[1]-x[0];
	float fk;
	float sign=1;

	/* determine lengths and scale factors for prime-factor FFTs */
	
	ntfft = npfar(nt);
	nxfft = npfa(nx);
	sfft = 1.0/(ntfft*nxfft);
	
	/* determine frequency and wavenumber sampling */
	nw = ntfft/2+1;
	dw = 1.0/(ntfft*dt);
	fw = 0.0;

	nk= nxfft;
	*pnk=nk;
	*pnw=nw;
	dk = 1.0/(nxfft*dx);
	fk = -1.0/(2*dx);


	for (k[0]=fk,ik=1; ik<nk; ik++) k[ik]=k[ik-1]+dk;
	for (w[0]=fw,iw=1; iw<nw; iw++) w[iw]=w[iw-1]+dw;
	fprintf(stderr,"dt=%f,dx=%f,dw=%f,fk=%f,k[nk-1]=%f,fw=%f,w[nk-1]=%f\n",
		dt,dx,dw,fk,k[nk-1],fw,w[nw-1]);

	fprintf(stderr,"nt=%d,nx=%d,nw=%d,ntfftk=%d\n",nt,nx,nw,ntfft);

	/* allocate real and complex workspace for FFTs */
	cpfft = alloc2complex(nw,nk);
	pfft = alloc2float(ntfft,nxfft);

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
	for (ik=0; ik<nxfft; ik++)  for (iw=0; iw<nw; iw++) datafk[ik][iw]=cpfft[ik][iw];
	/* Fourier transform k to x */
	pfa2cc(1,2,nw,nxfft,cpfft[0]);
	/* Fourier transform w to t */
	pfa2cr(-1,1,ntfft,nx,cpfft[0],pfft[0]);
	/* copy filtered data from FFT array to output */
	for (ix=0; ix<nx; ix++){
	  /* if ix odd, pfft has been negated to center transform of dimension 2 */
	  if (ISODD(ix)) sign=-1; else sign=1; 
	  for (it=0; it<nt; it++) data2xt[ix][it]=sign*pfft[ix][it];
	}
	/* free workspace */
	free2complex(cpfft);
	free2float(pfft);
}

void fksize(int nt, int nx, int *pnw, int *pnk)
{
  *pnw = npfar(nt)/2+1;
  *pnk = npfa(nx);
}




















