#include "interpfk.h"
#include "segy.h"
#define DFT 1
#define testpad 1

void stoltzop2(float **data, float **model, int nt, int nh, int nx, float *t, 
	       float *h, float *x, float vel, complex **F, complex **F2, float *wavelet, 
	       int nw, int adj)
{
  /* This operator performs Stoltz migration on an irregularly sampled data set
     such that the output is a regularly sampled data set (called model).
     The time axis are the same for data and model, the offset axis are,
     in general, different. The offset axis for data is h, the offset axis for
     model is x. The irregularly sampled data requires use of the least squares
     Fourier transform (x-->k).
     Hence the operator is 

     d= W F M FFT m
     m= FFTI M F d 
      
  */

  float **px;
  complex **pk;
  int ix,it,ik,ih;
  float dt=t[1]-t[0];
  float fmax=1./(2*dt);
  float dx = (x[nx-1]-x[0])/(nx-1);
  float vmax = vel;
  /* wavenumber (k) sampling */
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  if (testpad) nxpad*=2;
  int nxfft;
  if (DFT) nxfft=nx+nxpad; // For the DFT we do not need npfar
  else nxfft = npfar(nx+nxpad);
  int nk = (int) (nxfft/2+1);
  float dk = 2.0*PI/(nxfft*dx);
  float *k,kmin=0; // Wavenumber axis
  complex *d,*d2,*m; // Temporal arrays
  complex czero;czero.r=czero.i=0;
  int verbose1=1;
  float scale=1.0/nxfft,kx;
  float **model2;

  if (verbose1) fprintf(stderr," dx=%f, nxpad=%d, nxfft=%d \n", dx, nxpad, nxfft);

  /* allocate and zero common-offset gather p(t,x) */
  pk = ealloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = ealloc2float(nt,nxfft); // px requires zero padding for pfarc

  // Zeroed arrays //

  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;

  // Generate k axis
  k=ealloc1float(nk);
 
  // Create temporal complex arrays
  model2=ealloc2float(nt,nx);
  m=ealloc1complex(nk);
  d=ealloc1complex(nh);
  d2=ealloc1complex(nx);
  
  zero_vector(m,nk);
  zero_vector(d,nh);
  zero_vector(d2,nx);

  for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) model2[ix][it] = 0;
  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  fprintf(stderr,"Inside stoltzop nk=%d,dk=%f\n",nk,dk);
  
  if (adj){
    for (ih=0; ih<nh; ++ih) for (it=0; it<nt; ++it) px[ih][it] = data[ih][it];
    /* Fourier transform p(u,x) to p(u,k) */
    // Substitute the regular sampled FFT with DFT
    // the input data has a size nh x nt, and the offset is h
    if (DFT){   
      for (it=0;it<nt;it++){ 
	for (ih=0;ih<nh;ih++) {
	  d[ih].r=px[ih][it];
	  d[ih].i=0;
	}
	Atimesx_DFT(d,F,m,nh,nk,1);
	for (ik=0;ik<nk;ik++) pk[ik][it]=m[ik];
      }
    }
    else{
      fprintf(stderr,"using fft\n");      
      pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 
    }
    if (0) plotAmpSpec(pk,nk,nt,dt);
    
    /* migrate each wavenumber */
        
    if (1)
      for (ik=1,kx=dk; ik<nk; ++ik,kx+=dk)
	stolt1kadj(kx,vel,0.6,fmax,nt,dt,pk[ik],pk[ik]);
    
    /* Fourier transform p(u,k) to p(u,x) and scale */
    /* Fourier transform p(k) to p(u,x) and scale */

    if (DFT){
      for (it=0;it<nt;it++){
	// We need to respect the symmetry along 
	for (ik=0;ik<nk;ik++) m[ik]=pk[ik][it];
	Atimesx_DFT(d2,F2,m,nx,nk,0);
	//    if (adj) Atimesx_DFT(d,F,m,nh,nk,0);
	//    else Atimesx_DFT(d,F,m,nh,nk,0);
	for (ix=0;ix<nx;ix++) px[ix][it]=d2[ix].r;
      }
    }
    else{
      fprintf(stderr,"using fft\n");
      pfa2cr(1,2,nt,nxfft,pk[0],px[0]);
    }
    for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) model2[ix][it] = scale*px[ix][it];
    if(0){ // see the wavelet
      char buf[80];
      sprintf(buf,"xgraph < wavelet n=%d pairs=2 d1=1 style=normal title=\"wavelet\" & ",nw);
      plotvector(wavelet,nw,buf);
    }

    int plotmodel=0;
    if (plotmodel){
      fprintf(stderr,"nx=%d,nt=%d\n",nx,nt);
      save_gather(model2,nx,nt,0.004,"migrated_bef.su");
      system("suxwigb < migrated_bef.su perc=100 title=migrated_bef ");
    }	
    convgather(model2[0],model[0],nx,nt,wavelet,nw,1);
    
    if (plotmodel){
      fprintf(stderr,"nx=%d,nt=%d\n",nx,nt);
      save_gather(model,nx,nt,0.004,"migrated_aft.su");
      system("suxwigb < migrated_aft.su perc=100 title=migrated_aft ");
    }	
  }
  else{ // (!adj)
    convgather(model2[0],model[0],nx,nt,wavelet,nw,adj);
    for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) px[ix][it] = scale*model2[ix][it];
    /* Fourier transform p(u,x) to p(u,k) */
    if (DFT){
      for (it=0;it<nt;it++){
	// We need to respect the symmetry along 
	for (ix=0;ix<nx;ix++){ 
	  d2[ix].r=px[ix][it];
	  d2[ix].i=0;
	}
	Atimesx_DFT(d2,F2,m,nx,nk,1);
	//    if (adj) Atimesx_DFT(d,F,m,nh,nk,0);
	//    else Atimesx_DFT(d,F,m,nh,nk,0);
	for (ik=0;ik<nk;ik++) pk[ik][it]=m[ik];
      }
    }
    else{
      fprintf(stderr,"using fft\n");
      pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 
    }

    if (0) plotAmpSpec(pk,nk,nt,dt);
    
    /* migrate each wavenumber */
    for (ik=1,kx=dk; ik<nk; ++ik,kx+=dk){
      if (0) fprintf(stderr,"ik=%d.....\n",ik);
      //stolt1kadj(kx,vel,1,fmax,nt,dt,pk[ik],pk[ik]);
      stolt1k(kx,vel,0.6,fmax,nt,dt,pk[ik],pk[ik]);
    }  
    
    /* Fourier transform p(u,k) to p(u,x) and scale */
    if (DFT){
      for (it=0;it<nt;it++){ 
	for (ik=0;ik<nk;ik++) {
	  m[ik]=pk[ik][it];
	}
	Atimesx_DFT(d,F,m,nh,nk,0);
	for (ih=0;ih<nh;ih++) data[ih][it]=d[ih].r;
      }
    }
    else{
      fprintf(stderr,"using fft\n");
      pfa2cr(1,2,nt,nxfft,pk[0],px[0]);
      for (it=0;it<nt;it++)
	for (ih=0;ih<nh;ih++)
	  data[ih][it]=px[ih][it];
    }

  }
  free1complex(d2);
  free1complex(d);
  free1complex(m);
  free1float(k);
  free2float(px);
  free2complex(pk);
  free2float(model2);
  

  
  return;
}

void kaxis(int nx, float vmax, float dt, int nt, float dx, int *pnk, float *pdk)
{
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  if (testpad) nxpad*=2;
  int nxfft;
  if (DFT) nxfft=nx+nxpad; // For the DFT we do not need npfar
  else nxfft = npfar(nx+nxpad);

  *pnk = (int) (nxfft/2+1);
  *pdk = 2.0*PI/(nxfft*dx);
  
  return;
}



void stoltzopinv2(float **data, float **model, int nt, int nx, float *t, 
		  float *x, float vel, complex **F, float *wavelet, int nw)
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
  if (testpad) nxpad*=2;
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
  pk = ealloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = ealloc2float(nt,nxfft); // px requires zero padding for pfarc

  // Zeroed arrays //
  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;

  // Generate k axis
  k=ealloc1float(nk);

  // Create temporal complex arrays
  m=ealloc1complex(nk);
  d=ealloc1complex(nx);

  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  convgather(model[0],model[0],nx,nt,wavelet,nw,0);  
  for (ix=0; ix<nx;ix++) for (it=0; it<nt; ++it) px[ix][it] = model[ix][it];
    
  /* Fourier transform p(u,x) to p(u,k) */
  pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 

  if (0) plotAmpSpec(pk,nk,nt,dt);

  /* migrate each wavenumber */
  for (ik=1,kx=dk; ik<nk; ++ik,kx+=dk)
    stolt1k(kx,vel,0.6,fmax,nt,dt,pk[ik],pk[ik]);

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






