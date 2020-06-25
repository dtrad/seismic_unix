#include "interpfk.h"
#include "segy.h"
#define testpad 1

void stoltzopinv2(float **data, float **model, int nt, int nx, float *t, 
		  float *x, float vel, complex **F)
{
  
  float **px;
  complex **pk;
  int ix,it,ik;
  float dt=t[1]-t[0];
  float fmax=1./(2*dt);
  float dx = (x[nx-1]-x[0])/(nx-1);
  float vmax = vel;
  /* wavenumber (k) sampling */
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  if (testpad) nxpad=nx;
  int nxfft = npfar(nx+nxpad);
  //int nxfft=nx+nxpad; // For the DFT we do not need npfar
  int nk = (int) (nxfft/2+1);
  float dk = 2.0*PI/(nxfft*dx);
  float *k,kmin=0; // Wavenumber axis
  complex *d,*m; // Temporal arrays
  int verbose1=1;
  float scale=1.0/nxfft,kx;

  if (verbose1) fprintf(stderr," dx=%f, nxpad=%d, nxfft=%d \n", dx, nxpad, nxfft);

  /* allocate and zero common-offset gather p(t,x) */
  pk = alloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = alloc2float(nt,nxfft); // px requires zero padding for pfarc

  // Zeroed arrays //
  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;

  // Generate k axis
  k=ealloc1float(nk);

  // Create temporal complex arrays
  m=ealloc1complex(nk);
  d=ealloc1complex(nx);

  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  
  for (ix=0; ix<nx;ix++) for (it=0; it<nt; ++it) px[ix][it] = model[ix][it];
    
  /* Fourier transform p(u,x) to p(u,k) */
  pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 

  if (0) plotAmpSpec(pk,nk,nt,dt);

  /* migrate each wavenumber */
  for (ik=1,kx=dk; ik<nk; ++ik,kx+=dk)
    stolt1k(kx,vel,1,fmax,nt,dt,pk[ik],pk[ik]);

  /* Fourier transform p(k) to p(u,x) and scale */
  if (0){
    for (it=0;it<nt;it++){
      // We need to respect the symmetry along 
      for (ik=0;ik<nk;ik++) m[ik]=pk[ik][it];
      Atimesx_DFT(d,F,m,nx,nk,0);
      for (ix=0;ix<nx;ix++) px[ix][it]=d[ix].r;
    }
  }
  else pfa2cr(1,2,nt,nxfft,pk[0],px[0]);

  for (ix=0; ix<nx; ++ix)
    for (it=0; it<nt; ++it)
      data[ix][it] = scale*px[ix][it];
  

  free1complex(d);
  free1complex(m);
  free1float(k);
  free2float(px);
  free2complex(pk);
  
  return;
}




