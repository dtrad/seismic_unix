#include "su.h"
#include "dan.h"
#include "segy.h"
#include "Complex.h"
#define DFT 0

/* Interpolates spatial direction by padding zeros in the f-k domain 
   returns 
   nk2 the number of (new) k elements 
   k2 the new wavenumber axis
   nx the number of (new) k elements 
   x  the new offset axis
   dataout(nx,nt);
   
*/


int fft2_zeropad(float **datain, float **dataout, float **Wd, int nt, int nx, float interpfact)
{
  int ix,it,ik;
  complex czero;czero.r=czero.i=0;
  int nxpad;
  int nxfft;
  int nxfft2;
  int nx2;
  
  nxfft = npfar(nx);
  nxpad=(int) (interpfact-1)*nxfft;  
  nxfft2 = npfar(nx+nxpad);
  int nk = (int) (nxfft/2+1);
  int nk2 = (int) (nxfft2/2+1);

  nx2=(int) (interpfact*nx);

  float **px;
  complex **pk;
  float **px2;
  complex **pk2;

  int verbose1=1;
  float scale=1.0/nxfft;

  /* for data weights ony */
  int count;
  int iz,iy; 

  /*********** Data interpolation *****************/
  if (verbose1) 
    fprintf(stderr," nx=%d, nx2=%d, nk=%d, nk2=%d, \n nxpad=%d, nxfft=%d, nxfft2=%d \n", nx, nx2, nk, nk2, nxpad, nxfft, nxfft2);

  /* allocate and zero common-offset gather p(t,x) */
  pk = ealloc2complex(nt,nk); // k axis is symmetric, hence nk=nxfft/2+1
  px = ealloc2float(nt,nxfft); // px requires zero padding for pfarc
  pk2 = ealloc2complex(nt,nk2); // k axis is symmetric, hence nk=nxfft/2+1
  px2 = ealloc2float(nt,nxfft2); // px requires zero padding for pfarc

  // Zeroed arrays //
 
  for (ix=0; ix<nxfft; ++ix) for (it=0; it<nt; ++it) px[ix][it] = 0;
  for (ix=0; ix<nk; ++ix) for (it=0; it<nt; ++it) pk[ix][it].r = pk[ix][it].i=0;
 
  for (ix=0; ix<nxfft2; ++ix) for (it=0; it<nt; ++it) px2[ix][it] = 0;
  for (ix=0; ix<nk2; ++ix) for (it=0; it<nt; ++it) pk2[ix][it].r = pk2[ix][it].i=0;
 

  for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) px[ix][it] = datain[ix][it];
 
  pfa2rc(-1,2,nt,nxfft,px[0],pk[0]); 
 
  for (ik=0; ik<nk; ++ik) for (it=0; it<nt; ++it) pk2[ik][it] = pk[ik][it];
  for (ik=nk; ik<nk2; ++ik) for (it=0; it<nt; ++it) pk2[ik][it] = czero;
 
  pfa2cr(1,2,nt,nxfft2,pk2[0],px2[0]);
 
  for (ix=0; ix<nx2; ++ix) for (it=0; it<nt; ++it)  dataout[ix][it] = scale*px2[ix][it] ;

  free2float(px);
  free2complex(pk);
  free2float(px2);
  free2complex(pk2);
  
  /****************** Data weights **************/
  nx2=((int) interpfact)*nx;
  count=0;
  fprintf(stderr,"setting data weights \n");
  for (ix=0; ix<nx; ++ix){ 
    iz=ix*(int) interpfact;
    for (it=0; it<nt; ++it) Wd[iz][it]=1;
    for (iy=1; iy<interpfact;++iy){ 
      for (it=0; it<nt; ++it) Wd[iz+iy][it]=0;
      count++;
    }
  }
 
  return(nx2);


}

int add_zerotraces(float **datain, float **dataout, float **Wd, int nt, int nx, float interpfact)
{
  int ix,it,iy,iz;
  int count;
  int nx2;
  nx2=((int) interpfact)*nx;
  count=0;
  fprintf(stderr,"Interpolating with zero traces (in between traces)\n");
  for (ix=0; ix<nx; ++ix){ 
    iz=ix*(int) interpfact;
    for (it=0; it<nt; ++it){
      dataout[iz][it] = datain[ix][it] ;
      Wd[iz][it]=1;
    }
    for (iy=1; iy<interpfact;++iy){ 
      for (it=0; it<nt; ++it){
	dataout[iz+iy][it] = 0 ;
	Wd[iz+iy][it]=0;
      }
      count++;
    }
  }
  return(nx2);


}
 
int add_zerotraces(float **datain, float **dataout, float **Wd, int nt, int nx, int nx2, float *h, float *h2)
{
  int ix,it,ix2;
  fprintf(stderr,"--Interpolating with zero traces inside gaps\n");
  fprintf(stderr,"nx=%d,nx2=%d\n",nx,nx2);
  for (ix=0,ix2=0; (ix<nx||ix2<nx2); ++ix,++ix2){ 
    if (h[ix]==h2[ix2]){ 
      fprintf(stderr,"h[%d]=%f---->h2[%d]=%f\n",ix,h[ix],ix2,h2[ix2]);
      for (it=0; it<nt; ++it){
	dataout[ix2][it] = datain[ix][it] ;
	Wd[ix2][it]=1;
      }
    }
    else{
      fprintf(stderr,"h[%d]=%f-- h2[%d]=%f<----- zero trace \n",ix,h[ix],ix2,h2[ix2]);
      for (it=0; it<nt; ++it){
	dataout[ix2][it] = 0 ;
	Wd[ix2][it]=0;
      }
      ix--;
    }
  }
  fprintf(stderr,"ix=%d,ix2=%d\n",ix,ix2);
  return(ix2++);
}

















