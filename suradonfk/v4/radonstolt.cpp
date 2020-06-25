#include "radonfk.h"
#define STRETCH 0.6 // use 0.6 for other variable velocity and 1 for constant
#define testpad 0 /* Pad with extra zeros in offset */
#define PLOTSPEC 0 /* plot amplitude spectrum (1 w symm., 2 k symm.)*/
#define standard 0 /* type amplitude spectrum (1 w symm., 2 k symm.)*/

static int callCounter = 0;
void stoltzop2(float **data, float **model, int nt, int nx, float *t, 
	       float *x, float vel, int adj)
{
  /* This operator performs Stoltz migration on data set:
     such that the output is a regularly sampled data set (called model).
 
     If the data are irregularly sampled it uses DFT 
     or fft otherwise.

     The time axis are the same for data and model, the offset axis are,
     in general, different. The offset axis for data is h, the offset axis for
     model is x. 

     For regularly sampled data 

     d =  FFT2I M FFT2 m
     m =  FFT2I M FFT2 d
     
  */
  fprintf(stderr,"call number = %d \n",callCounter);
  callCounter++;
  float **px;
  complex **pk;
  int ix,it,ik,ih;
  float dt=t[1]-t[0];
  float fmax=1./(2*dt);
  fprintf(stderr,"nx = %d , x [%d] = %f, x[%d]= %f vel=%f \n", nx, nx-1, x[nx-1], 0, x[0], vel);
  float dx = (x[nx-1]-x[0])/(nx-1);
  float vmax = vel;
  /* wavenumber (k) sampling */
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  if (testpad) nxpad*=2;
  int nxfft;

  nxfft = npfar(nx+nxpad);
  int nk = (int) (nxfft/2+1);
  float dk = 2.0*PI/(nxfft*dx);
  float *k,kmin=0; // Wavenumber axis
  complex czero;czero.r=czero.i=0;
  int verbose1=1;
  float scale=1.0/nxfft,kx;
 
  // Test: limit the maximum kx to migrate
  int nkmax=(int) (nk/1);

  if (verbose1) fprintf(stderr," dx=%f, nxpad=%d, nxfft=%d \n", dx, nxpad, nxfft);

  /* allocate and zero common-offset gather p(t,x) */
  pk = ealloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = ealloc2float(nt,nxfft); // px requires zero padding for pfarc

  // Zeroed arrays //

  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;



  // Generate k axis
  k=ealloc1float(nk);
 
  
  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  fprintf(stderr,"Inside stoltzop nk=%d,dk=%f\n",nk,dk);
  
  if (adj){
    for (ih=0; ih<nx; ++ih) for (it=0; it<nt; ++it) px[ih][it] = data[ih][it];

    save_gather(px,nxfft,nt,0.004,"ppk.su");
      
    /* Fourier transform p(u,x) to p(u,k) */
    // Substitute the regular sampled FFT with DFT
    // the input data has a size nx x nt, and the offset is h
 
    pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 
    
    if (PLOTSPEC){
      xplotgather(data,nx,nt,dt,"datatx",callCounter," perc=100 xbox=500 &");
      if (standard) plotAmpSpec(data,nx,nt,dt,"datafk");
      else          plotAmpSpec(pk,nk,nt,dt,"datafk");

    }
    
    /* migrate each wavenumber */
    for (ik=1,kx=dk; ik<nkmax; ++ik,kx+=dk)
      stolt1k(kx,vel,STRETCH,fmax,nt,dt,pk[ik],pk[ik],1);
    
    /* Fourier transform p(u,k) to p(u,x) and scale */
    /* Fourier transform p(k) to p(u,x) and scale */

    pfa2cr(1,2,nt,nxfft,pk[0],px[0]);
    for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) model[ix][it] = scale*px[ix][it];

    /* The following two lines apply a mute to the migrated domain to get rid of aliasing */

    //for (ix=0; ix<120; ++ix) for (it=0; it<nt; ++it) model[ix][it] = model[nx-ix-1][it] = 0;
    //for (ix=0; ix<nx; ++ix) for (it=0; it<185; ++it) model[ix][it] = model[nx-ix-1][it] = 0;
    int plotmodel=0;
    
    if (plotmodel){
      fprintf(stderr,"nx=%d,nt=%d\n",nx,nt);
      save_gather(model,nx,nt,0.004,"migrated_aft.su");
      system("suxwigb < migrated_aft.su perc=100 title=migrated_aft ");
    }	
  }
  else{ // (!adj)
    for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) px[ix][it] = scale*model[ix][it];
    /* Fourier transform p(u,x) to p(u,k) */
    pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 

    if (PLOTSPEC){
      xplotgather(model,nx,nt,dt,"modeltx",callCounter," perc=100 xbox=500 &");
      if (standard) plotAmpSpec(model,nx,nt,dt,"modelfk");
      else          plotAmpSpec(pk,nk,nt,dt,"modelfk");

    }
    
    /* migrate each wavenumber */
    for (ik=1,kx=dk; ik<nkmax; ++ik,kx+=dk){
      if (0) fprintf(stderr,"ik=%d.....\n",ik);
      //stolt1kadj(kx,vel,1,fmax,nt,dt,pk[ik],pk[ik]);
      stolt1k(kx,vel,STRETCH,fmax,nt,dt,pk[ik],pk[ik],0);
    }  
    
    /* Fourier transform p(u,k) to p(u,x) and scale */
    pfa2cr(1,2,nt,nxfft,pk[0],px[0]);
    for (it=0;it<nt;it++) for (ih=0;ih<nx;ih++) data[ih][it]=px[ih][it];

  }

  free1float(k);
  free2float(px);
  free2complex(pk);
 
  

  
  return;
}

void kaxis(int nx, float vmax, float dt, int nt, float dx, int *pnk, float *pdk, int dft)
{
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  if (testpad) nxpad*=2;
  int nxfft;
  if (dft) nxfft=nx+nxpad; // For the DFT we do not need npfar
  else nxfft = npfar(nx+nxpad);

  *pnk = (int) (nxfft/2+1);
  *pdk = 2.0*PI/(nxfft*dx);
  
  return;
}



void stoltzopinv2(float **data, float **model, int nt, int nx, float *t, float *x, float vel)
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
  int nkmax=(int) (nk/1);
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

  for (ik=0;ik<nk;ik++) k[ik]=kmin+ik*dk;
  for (ix=0; ix<nx;ix++) for (it=0; it<nt; ++it) px[ix][it] = model[ix][it];
    
  /* Fourier transform p(u,x) to p(u,k) */
  pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 

  if (0) plotAmpSpec(pk,nk,nt,dt,"");

  /* migrate each wavenumber */
  for (ik=1,kx=dk; ik<nkmax; ++ik,kx+=dk)
    stolt1k(kx,vel,STRETCH,fmax,nt,dt,pk[ik],pk[ik],0);

  /* Fourier transform p(k) to p(u,x) and scale */
  pfa2cr(1,2,nt,nxfft,pk[0],px[0]);

  for (ix=0; ix<nx; ++ix)
    for (it=0; it<nt; ++it)
      data[ix][it] = scale*px[ix][it];

  free1float(k);
  free2float(px);
  free2complex(pk);
  
  return;
}






