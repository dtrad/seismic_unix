#include "interpfk.h"
#include "segy.h"


void stoltzop(float **datain, float **dataout, int nt, int nh, float *t, 
	       float *h, float vel, complex **F, complex **FLS, int adj)
{
  
  float **px,**qq;
  complex **pk;
  int ix,it,ik,ih;
  float dt=t[1]-t[0];
  float fmax=1./(2*dt);
  int nx=nh;
  float dx = (h[nh-1]-h[0])/(nh-1);
  float vmax = vel;
  /* wavenumber (k) sampling */
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  //  int nxfft = npfar(nx+nxpad);
  int nxfft=nx+nxpad; // For the DFT we do not need npfar
  int nk = (int) (nxfft/2+1);
  float dk = 2.0*PI/(nxfft*dx);
  float *k,kmin=0; // Wavenumber axis
  complex *d,*m; // Temporal arrays
  int verbose1=1;


  if (verbose1) fprintf(stderr," dx=%f, nxpad=%d, nxfft=%d \n", dx, nxpad, nxfft);

  /* allocate and zero common-offset gather p(t,x) */
  pk = alloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = alloc2float(nt,nxfft); // px requires zero padding for pfarc
  qq = ealloc2float(nt,nx);  // The original t-h space

  // Zeroed arrays //
  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;
  for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) qq[ix][it] = 0;

  // Generate k axis
  k=ealloc1float(nk);
 
  // Create temporal complex arrays
  m=ealloc1complex(nk);
  d=ealloc1complex(nh);

  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  fprintf(stderr,"Inside stoltzop nk=%d,dk=%f\n",nk,dk); 
  for (ix=0; ix<nx; ++ix)    
    for (it=0; it<nt; ++it) qq[ix][it] = 0;
  

  for (ih=0; ih<nh; ++ih) for (it=0; it<nt; ++it) px[ih][it] = datain[ih][it];
  /* local variables */
  float scale=1.0/nxfft,kx;
    
  /* Fourier transform p(u,x) to p(u,k) */
  // Substitute the regular sampled FFT with DFT
  // the input data has a size nh x nt, and the offset is h

  for (it=0;it<nt;it++){ 
    for (ih=0;ih<nh;ih++) {
      d[ih].r=px[ih][it];
      d[ih].i=0;
    }
    Atimesx_DFT(d,F,m,nh,nk,1);
    //if (adj) Atimesx_DFT(d,FLS,m,nh,nk,1);
    //else Atimesx_DFT(d,FLS,m,nh,nk,1);
    for (ik=0;ik<nk;ik++) pk[ik][it]=m[ik];
  }

  if (0) plotAmpSpec(pk,nk,nt,dt);


  /* migrate each wavenumber */
  if (1)
  for (ik=1,kx=dk; ik<nk; ++ik,kx+=dk){
    if (0) fprintf(stderr,"ik=%d.....\n",ik);
    //TRACE;
    if (adj) stolt1kadj(kx,vel,1,fmax,nt,dt,pk[ik],pk[ik]);
    else stolt1k(kx,vel,1,fmax,nt,dt,pk[ik],pk[ik]);
  }  

  /* Fourier transform p(k) to p(u,x) and scale */
  for (it=0;it<nt;it++){
    // We need to respect the symmetry along 
    for (ik=0;ik<nk;ik++) m[ik]=pk[ik][it];
    Atimesx_DFT(d,F,m,nh,nk,0);
    //    if (adj) Atimesx_DFT(d,F,m,nh,nk,0);
    //    else Atimesx_DFT(d,F,m,nh,nk,0);
    for (ih=0;ih<nh;ih++) px[ih][it]=d[ih].r;
  }

  for (ix=0; ix<nx; ++ix)
    for (it=0; it<nt; ++it)
      px[ix][it] *= scale;

  /* add migrated traces to mix */
  for (ix=0; ix<nx; ++ix)
    for (it=0; it<nt; ++it)
      qq[ix][it] = px[ix][it];
  
  /* zero common-offset gather */
  for (ix=0; ix<nxfft; ++ix)
    for (it=0; it<nt; ++it)
      px[ix][it] = 0.0;

  if (0){
    save_gather(qq,nx,nt,0.004,"migrated.su");
    system("suxwigb < migrated.su perc=100 title=migrated ");
  }	
  
  for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) dataout[ix][it]=qq[ix][it];

  free1complex(d);
  free1complex(m);
  free1float(k);
  free2float(qq);
  free2float(px);
  free2complex(pk);


  
  return;
}


