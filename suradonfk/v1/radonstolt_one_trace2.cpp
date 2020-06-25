/* Stolt operator */
/* Taken from program SUSTOLT: $Revision: 1.12 $ ; $Date: 1997/07/28 22:36:46 $		*/


/*****************************************************************************
Stolt's migration for one wavenumber k
******************************************************************************
Input:
k		wavenumber
v		velocity
s		Stolt stretch factor (0<s<2; use s=1 for constant velocity)
fmax		maximum frequency (in cycles per unit time)
nt		number of time samples
dt		time sampling interval (first time assumed to be zero)
p		array[nt] containing data to be migrated

Output
q		array[nt] containing migrated data (may be equivalenced to p)
*****************************************************************************/

void plotimage(float **d, int n1, int n2, const char *s, float perc);
#include "radonfk.h"
#include "segy.h"
#define testpad 1
#define npad 2

void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
	      complex *p, complex *q, int  adj)
{



  /* arguments that should be added to this operator */
  float wmax;
  float dw;
  float fw;
  float wtau;
  float wtauh;
  float wtaul;
  float scale;
  float fftscl;  
  int nw;
  complex *pp, *qq;
  /*****************************************************/  
  int it,nwtau,iwtau,ntau,itau,iwtaul,iwtauh;
  float dtau, dwtau, fwtau;
  float a,b;
  complex czero;
  float *wwtau, *wwtau2;
  float vko2s;
  cmplx(0.0,0.0,czero);

  /* modify stolt stretch factor to simplify calculations below */
  if (s!=1.0) s = 2.0-s;

  /* frequency sampling - must pad to avoid interpolation error;
   * pad by factor of 2 because time axis is not centered;
   * pad by factor of 1/0.6 because 8-point sinc is valid
   * only to about 0.6 Nyquist
   */
	
  if (testpad){
    nw = (int) (nt*npad/0.6);
    nw=npfao(nw,nw*2) ;
  }	
  else{
    nw = (int) (nt*2/0.6);
    nw = npfao(nw,nw*2);
  }
  dw = 2.0*PI/(nw*dt);
  fw = -PI/dt;

      
  /* migrated time */
  ntau = nt;
  dtau = dt;

  /* migrated frequency - no need to pad since no interpolation */
  nwtau=nw;
  dwtau = dw;
  fwtau = fw;
  
  /* tweak first migrated frequency to avoid wtau==0.0 below */
  fwtau += 0.001*dwtau;
  
  /* maximum frequency to migrate in radians per unit time */
  wmax = 2.0*PI*MIN(fmax,0.5/dt);
  
  /* (v*k/2)^2 */
  vko2s = 0.25*v*v*k*k;

  /* high and low migrated frequencies - don't migrate evanescent */
  wtauh = sqrt(MAX(0.0,wmax*wmax-s*vko2s));
  iwtauh = MAX(0,MIN(nwtau-1,NINT((wtauh-fwtau)/dwtau)));
  iwtaul = MAX(0,MIN(nwtau-1,NINT((-wtauh-fwtau)/dwtau)));
  wtauh = fwtau+iwtauh*dwtau;
  wtaul = fwtau+iwtaul*dwtau;
  if (0) fprintf(stderr,"iwtaul=%d,iwtauh=%d,nw=%d,nwtau=%d\n",iwtaul,iwtauh,nw,nwtau);

  /* workspace */
  pp = ealloc1complex(nw);
  qq = ealloc1complex(nwtau);
  wwtau = ealloc1float(nwtau);
  wwtau2 = ealloc1float(nwtau);
  
  /* pad with zeros and Fourier transform t to w, with w centered */
  for (it=0; it<nt; it+=2)
    pp[it] = p[it];
  for (it=1; it<nt; it+=2) {
    pp[it].r = -p[it].r;
    pp[it].i = -p[it].i;
  }
  for (it=nt; it<nw; it++)
    pp[it].r = pp[it].i = 0.0;
  pfacc(1,nw,pp);
  

  /* zero -Nyquist frequency for symmetry */
  pp[0] = czero;
  
  /* frequencies at which to interpolate */
  if (s==1.0) {
    for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
      wwtau[iwtau] = -sqrt(wtau*wtau+vko2s);
    for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
      wwtau[iwtau] = sqrt(wtau*wtau+vko2s);
  } else {
    a = 1.0/s;
    b = 1.0-a;
    for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
      wwtau[iwtau] = b*wtau-a*sqrt(wtau*wtau+s*vko2s);
    for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
      wwtau[iwtau] = b*wtau+a*sqrt(wtau*wtau+s*vko2s);
  }
  if (0) fprintf(stderr,"iwtau=%d\n",iwtau);
  
  
  for (iwtau=0;iwtau<nwtau;iwtau++) qq[iwtau]=czero;
  
  for (iwtau=iwtaul,wtau=wtaul;iwtau<iwtauh; ++iwtau,wtau+=dwtau){
    if ((NINT((wwtau[iwtau]-fw)/dw) < nwtau) && (iwtau < nw)){
      if (!adj){
	qq[NINT((wwtau[iwtau]-fwtau)/dwtau)]=pp[iwtau];
      }else
	qq[iwtau]=pp[NINT((wwtau[iwtau]-fw)/dw)];	
    }
  }
  
  /* fft scaling and obliquity factor */
  fftscl = 1.0/nwtau;
  scale = fftscl;
  
  if (s==1.0){ 
    for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh;++iwtau,wtau+=dwtau) {
      qq[iwtau].r *= scale;
      qq[iwtau].i *= scale;
    }
  } else {
    a = 1.0/(s*s);
    b = 1.0-1.0/s;
    for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh;++iwtau,wtau+=dwtau) {
      qq[iwtau].r *= scale;
      qq[iwtau].i *= scale;
    }
  }
  
  /* zero evanescent frequencies */
  for (iwtau=0; iwtau<iwtaul; ++iwtau)
    qq[iwtau] = czero;
  for (iwtau=iwtauh+1; iwtau<nwtau; ++iwtau)
    qq[iwtau] = czero;
  
  /* Fourier transform wtau to tau, accounting for centered wtau */
  pfacc(-1,nwtau,qq);
  for (itau=0; itau<ntau; itau+=2)
    q[itau] = qq[itau];
  for (itau=1; itau<ntau; itau+=2) {
    q[itau].r = -qq[itau].r;
    q[itau].i = -qq[itau].i;
  }  

  /* free workspace */
  free1complex(pp);
  free1complex(qq);
  free1float(wwtau2);
  free1float(wwtau);
  
}


void taper (int lxtaper, int lbtaper, 
	    int nx, int ix, int nt, float *trace, int inv)
/*****************************************************************************
Taper traces near left and right sides of trace window
******************************************************************************
Input:
lxtaper		length (in traces) of side taper
lbtaper		length (in samples) of bottom taper
nx		number of traces in window
ix		index of this trace (0 <= ix <= nx-1)
nt		number of time samples
trace		array[nt] containing trace

Output:
trace		array[nt] containing trace, tapered if within lxtaper of side
*****************************************************************************/
{
  int it;
  float xtaper;
  int debug = 1;
  //inv=0;
  /* if near left side */
  if (ix<lxtaper) {
    xtaper = 0.54+0.46*cos(PI*(lxtaper-ix)/lxtaper);
    
    /* else if near right side */
  } else if (ix>=nx-lxtaper) {
    xtaper = 0.54+0.46*cos(PI*(lxtaper+ix+1-nx)/lxtaper);
    
    /* else x tapering is unnecessary */
  } else {
    xtaper = 1.0;
  }
  if (debug > 2) fprintf(stderr,"xtaper=%f\n",xtaper);
  
  /* if x tapering is necessary, apply it */
  if (xtaper!=1.0){
    if (!inv)
      for (it=0; it<nt; ++it)
	trace[it] *= xtaper;
    else
      if (xtaper>0.001)
	for (it=0; it<nt; ++it)
	  trace[it] /= xtaper;
  }

  /* if requested, apply t tapering */
  for (it=MAX(0,nt-lbtaper); it<nt; ++it)
    trace[it] *= (0.54+0.46*cos(PI*(lbtaper+it+1-nt)/lbtaper));
}


void plotAmpSpec(complex **p, int nk, int nt, float dt, char *c)
{
  int nw;
  float dw;
  float fw;
  complex *pp;
  int it,ik;
  complex czero;czero.r=czero.i=0;
  
  float **aspec; // Amplitude spectrum
  char buf[12];
  
  nw = (int) (nt*2/0.6);
  nw = npfao(nw,nw*2);
  dw = 2.0*PI/(nw*dt);
  fw = -PI/dt;
  aspec=ealloc2float(nw,nk);
  pp = ealloc1complex(nw);
  
  for (ik=0;ik<nk;ik++){
    for (it=0; it<nt; it+=2)
      pp[it] = p[ik][it];
    for (it=1; it<nt; it+=2) {
      pp[it].r = -p[ik][it].r;
      pp[it].i = -p[ik][it].i;
    }
    for (it=nt; it<nw; ++it) pp[it].r = pp[it].i = 0.0;
    pfacc(1,nw,pp);
    /* zero -Nyquist frequency for symmetry */
    pp[0] = czero;
    for (it=0; it<nw; it++) aspec[ik][it]=abs(pp[it]);
      
  }
  fprintf(stderr,"Inside plotAmpSpec nw=%d nk=%d\n",nw,nk);
  if (0){
    save2d(aspec,nk,nw,dw,"temp");
    system("suxwigb < temp perc=100 title=temp &"); 
  }
  sprintf(buf,"fkspec%s",c);
  if (1) plotimage(aspec,nw,nk,buf,99.0);
  free1complex(pp);
  free2float(aspec);
  
}

void plotimage(float **d, int n1, int n2, const char *s, float perc)
{
  FILE* fp;

  char buf[80];
  int i2;
  fp=efopen(s,"w");

  for (i2=0;i2<n2;i2++)  
    efwrite(d[i2],sizeof(float),n1,fp);
  
  fclose(fp);
  sprintf(buf,"ximage < %s n1=%d n2=%d  perc=%f title=%s legend=1 & \n",s,n1,n2,perc,s);
  system(buf);

  return;
  
}


void save2d(float **d, int nh, int nt, float dt, char* name)
{
  segy tr2;
  
  int  itr;
  FILE* fp;
  
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr2.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr2.tracl=itr+1;
      tr2.dt=(int) (dt*1e6);
      tr2.ns=nt;
      tr2.ntr=nh;
      fprintf(stderr,"==>itr=%d\n",itr);
      fputtr(fp,&tr2);
      
      //for (int it=0;it<nt;it++) fprintf(stderr,"tr2.data[%d]=%f\n",it,tr2.data[it]); 
  }    

  fflush(fp);
  fclose(fp);
  
  return;
  
}




