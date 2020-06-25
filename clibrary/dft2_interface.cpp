#include "dft2.h"

/*

Daniel Trad - December - 2000
*/

void dft2_op(float **d, complex **m2, float*t, float* pos, float*k, int nt, int nh, int nk, float fmax, float& df);
void dft2_interface(float **data, complex **model, float *h, int nh, float *t, int nt, float *k, 
		      int nk, float *vel, float smute, float nmofactor, float fmax, 
		    float **data2, float *h2, int nh2, float **Wd, inv_par inv, int method, float& df)
{
  int ik, it, ih;
  float *dtemp;
  float dt=t[1]-t[0];
  float dk=k[1]-k[0];
  float kmin=k[0];
  float *htemp;
  int nf0;

  dtemp=alloc1float(nt);
  htemp=ealloc1float(nh);

  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT
  // All methods are the same with the only difference that the offsets are
  // transformed to pseudo offsets.

  //////////////////////////////////////////////////////////////////////

    
  //////////////////////////////////////////////////////////////////////
  if (nmofactor)
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  if (0){ 
    save_gather(data,nh,nt,t[1]-t[0],"data");
    system("suxwigb < data  title=\"plot data\" &");
  } 

  memset( (void *) model[0], (int) '\0', nk * 2  * nt *FSIZE);  
  
  dft2_op(data,model,t,h,k,nt,nh,nk,fmax,df);
  return;

  if (method==0) dft2_toep(data,model,t,h,k,nt,nh,nk,inv.eps1,fmax,&nf0);
  if (method==1) dft2_wtcgls(data,model,t,h,k,nt,nh,nk,Wd,inv,fmax,&nf0);

  fprintf(stderr,"In hrfft2_interface nf0=%d\n",nf0);

  dft2_inv(data2,model,t,h2,k,nt,nh2,nk,fmax,nf0);
  
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
