/* SUSTOLT: $Revision: 1.12 $ ; $Date: 1997/07/28 22:36:46 $		*/

#include "dan.h"
#include "inversion_par.h"
#include "segy.h"
void save_gather(float **d, int nh, int nt, float dt, char* name);

static void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
		     complex *p, complex *q);
static void stolt1kadj (float k, float v, float s, float fmax, int nt, float dt,
			complex *p, complex *q);

static void taper (int lxtaper, int lbtaper, 
		   int nx, int ix, int nt, float *trace);

int radonhypfk(float *t, float *q, float *h, float **model, float **data,
		int nt, int nh, int nq, float dt, float t0, inv_par inv, 
	       int plot,float dx0, int adj)
{
  float **px,**qq;
  complex **pk;
  int ix,it,iq;
  float fmax=1./(2*dt);

  int nx = nh;
  int apex;
  if (h[0]<0) apex=NINT(nx/2);
  else apex=0;
  float dx = (h[nh-1]-h[0])/(nh-1);
  dx=dx0; // Test for finding the right dx
  float vmax = q[nq-1];
  /* wavenumber (k) sampling */
  int nxpad = (int) (0.5*vmax*nt*dt/dx);
  int nxfft = npfar(nx+nxpad);
  int nk = nxfft/2+1;
  float dk = 2.0*PI/(nxfft*dx);
  int lstaper=0;
  int lbtaper=0;
  
  fprintf(stderr," dx=%f, nxpad=%d, nxfft=%d \n", dx, nxpad, nxfft);
  /* allocate and zero common-offset gather p(t,x) */
  pk = alloc2complex(nt,nk);
  px = alloc2float(nt,nxfft);
  qq = ealloc2float(nt,nx);
  for (ix=0; ix<nx; ++ix)
    for (it=0; it<nt; ++it)
      qq[ix][it] = 0;
  
  for (iq=0;iq<nq;iq++){
    for (ix=0; ix<nx; ++ix) for (it=0; it<nt; ++it) px[ix][it] = data[ix][it];
    /* local variables */
    int ik,ix,it;
    float scale=1.0/nxfft,k;
    
    /* apply side and bottom tapers */
    for (ix=0; ix<nx; ++ix)
      taper(lstaper,lbtaper,nx,ix,nt,px[ix]);
    
    /* Fourier transform p(u,x) to p(u,k) */
    pfa2rc(-1,2,nt,nxfft,px[0],pk[0]);
    
    /* migrate each wavenumber */
    fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
    for (ik=1,k=dk; ik<nk; ++ik,k+=dk){
      if (0) fprintf(stderr,"ik=%d.....\n",ik);
      //TRACE;
      if (adj) stolt1kadj(k,q[iq],1,fmax,nt,dt,pk[ik],pk[ik]);
      else stolt1k(k,q[iq],1,fmax,nt,dt,pk[ik],pk[ik]);
    }  
    /* Fourier transform p(u,k) to p(u,x) and scale */
    pfa2cr(1,2,nt,nxfft,pk[0],px[0]);
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
    
    if (plot==1){
      save_gather(qq,nx,nt,0.004,"migrated.su");
      system("suxwigb < migrated.su perc=100 title=migrated ");
    }	

    if (adj) for (it=0; it<nt; ++it) model[iq][it] = qq[apex][it];
    
  }
  free2float(qq);
  free2float(px);
  free2complex(pk);

  return EXIT_SUCCESS;
}
	

static void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
	complex *p, complex *q)
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
{
	int nw,it,nwtau,iwtau,ntau,itau,iwtaul,iwtauh;
	float vko2s,wmax,dw,fw,dwtau,fwtau,wtau,dtau,
		wtauh,wtaul,scale,fftscl,a,b,*wwtau;
	complex czero,*pp,*qq;
        float *wwtau2;
	cmplx(0.0,0.0,czero);
	/* modify stolt stretch factor to simplify calculations below */
	if (s!=1.0) s = 2.0-s;

	/* (v*k/2)^2 */
	vko2s = 0.25*v*v*k*k;

	/* maximum frequency to migrate in radians per unit time */
	wmax = 2.0*PI*MIN(fmax,0.5/dt);

	/* frequency sampling - must pad to avoid interpolation error;
	 * pad by factor of 2 because time axis is not centered;
	 * pad by factor of 1/0.6 because 8-point sinc is valid
	 * only to about 0.6 Nyquist
	 */
       	nw = (int) (nt*2/0.6);
	//nw=nt;
	nw = npfao(nw,nw*2);
	dw = 2.0*PI/(nw*dt);
	fw = -PI/dt;

	/* migrated time */
	ntau = nt;
	dtau = dt;

	/* migrated frequency - no need to pad since no interpolation */
	nwtau = npfao(ntau,ntau*2);
        nwtau=nw;
	dwtau = 2.0*PI/(nwtau*dtau);
	fwtau = -PI/dtau;

	/* tweak first migrated frequency to avoid wtau==0.0 below */
	fwtau += 0.001*dwtau;

	/* high and low migrated frequencies - don't migrate evanescent */
	wtauh = sqrt(MAX(0.0,wmax*wmax-s*vko2s));
	iwtauh = MAX(0,MIN(nwtau-1,NINT((wtauh-fwtau)/dwtau)));
	iwtaul = MAX(0,MIN(nwtau-1,NINT((-wtauh-fwtau)/dwtau)));
	wtauh = fwtau+iwtauh*dwtau;
	wtaul = fwtau+iwtaul*dwtau;
	fprintf(stderr,"iwtaul=%d,iwtauh=%d,nw=%d,nwtau=%d\n",iwtaul,iwtauh,nw,nwtau);

	/* workspace */
	pp = alloc1complex(nw);
	qq = alloc1complex(nwtau);
	wwtau = alloc1float(nwtau);
	wwtau2 = alloc1float(nwtau);

	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it+=2)
		pp[it] = p[it];
	for (it=1; it<nt; it+=2) {
		pp[it].r = -p[it].r;
		pp[it].i = -p[it].i;
	}
	for (it=nt; it<nw; ++it)
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
	fprintf(stderr,"iwtau=%d\n",iwtau);
	
	/* interpolate */
	if (0){ 
	  yxtoxy(iwtau,1.0,iwtaul,wwtau,iwtau,1.0,iwtaul,-iwtaul,iwtauh,wwtau2);
	  ints8c(nw,dw,fw,pp,czero,czero,iwtauh-iwtaul+1,wwtau2+iwtaul,qq+iwtaul);
	}
	else{
	  for (iwtau=0;iwtau<nwtau;iwtau++) qq[iwtau]=czero;
	  for (iwtau=iwtaul,wtau=wtaul;iwtau<iwtauh; ++iwtau,wtau+=dwtau){
	    if ((NINT((wwtau[iwtau]-fw)/dw) < nwtau) && (iwtau < nw)){
	      if (0) fprintf(stderr,"iwtau=%d,iwtau2=%d\n",iwtau,NINT((wtau-fw)/dwtau));
	      if (0) fprintf(stderr,"iwwtau[%d]=%d,nwtau=%d,nw=%d,fw=%f\n",iwtau,NINT((wwtau[iwtau]-fw)/dw),
			     nwtau,nw,fw);
	      ////TRACE;
	      qq[NINT((wwtau[iwtau]-fwtau)/dwtau)]=pp[iwtau];
	      //qq[iwtau]=pp[iwtau];
	      ////TRACE;
	    }
	  }
	}

	/* Do the same but without interpolation */
	if (0) {
		for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
			qq[iwtau+iwtaul]=pp[NINT(-sqrt(wtau*wtau+vko2s)/dw)];
		for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
			qq[iwtau+iwtaul]=pp[NINT(sqrt(wtau*wtau+vko2s)/dw)];
	}	
        
	/* fft scaling and obliquity factor */
	fftscl = 1.0/nwtau;
	if (s==1.0) {
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh; 
			++iwtau,wtau+=dwtau) {
		        scale = fftscl*wtau/wwtau[iwtau];
			qq[iwtau].r *= scale;
			qq[iwtau].i *= scale;
		}
	} else {
		a = 1.0/(s*s);
		b = 1.0-1.0/s;
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh; 
			++iwtau,wtau+=dwtau) {
			scale = fftscl*(b+a*wtau/(wwtau[iwtau]-b*wtau));
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

	//TRACE;
	/* free workspace */
	free1complex(pp);
	//TRACE;
	free1complex(qq);
	//TRACE;
	free1float(wwtau2);
	//TRACE;
	free1float(wwtau);
	//TRACE;
}

static void stolt1kadj (float k, float v, float s, float fmax, int nt, float dt,
	complex *p, complex *q)
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
{
	int nw,it,nwtau,iwtau,ntau,itau,iwtaul,iwtauh;
	float vko2s,wmax,dw,fw,dwtau,fwtau,wtau,dtau,
		wtauh,wtaul,scale,fftscl,a,b,*wwtau;
	complex czero,*pp,*qq;
        float *wwtau2;
	cmplx(0.0,0.0,czero);

	/* modify stolt stretch factor to simplify calculations below */
	if (s!=1.0) s = 2.0-s;

	/* (v*k/2)^2 */
	vko2s = 0.25*v*v*k*k;

	/* maximum frequency to migrate in radians per unit time */
	wmax = 2.0*PI*MIN(fmax,0.5/dt);

	/* frequency sampling - must pad to avoid interpolation error;
	 * pad by factor of 2 because time axis is not centered;
	 * pad by factor of 1/0.6 because 8-point sinc is valid
	 * only to about 0.6 Nyquist
	 */
	nw = (int) (nt*2/0.6);
	nw = npfao(nw,nw*2);
	dw = 2.0*PI/(nw*dt);
	fw = -PI/dt;

	/* migrated time */
	ntau = nt;
	dtau = dt;

	/* migrated frequency - no need to pad since no interpolation */
        nwtau=nw;
	//nwtau = npfao(ntau,ntau*2);
	dwtau = 2.0*PI/(nwtau*dtau);
	fwtau = -PI/dtau;

	/* tweak first migrated frequency to avoid wtau==0.0 below */
	fwtau += 0.001*dwtau;

	/* high and low migrated frequencies - don't migrate evanescent */
	wtauh = sqrt(MAX(0.0,wmax*wmax-s*vko2s));
	iwtauh = MAX(0,MIN(nwtau-1,NINT((wtauh-fwtau)/dwtau)));
	iwtaul = MAX(0,MIN(nwtau-1,NINT((-wtauh-fwtau)/dwtau)));
	wtauh = fwtau+iwtauh*dwtau;
	wtaul = fwtau+iwtaul*dwtau;
	if (0) fprintf(stderr,"iwtaul=%d,iwtauh=%d,nw=%d,nwtau=%d\n",iwtaul,iwtauh,nw,nwtau);

	/* workspace */
	pp = alloc1complex(nw);
	qq = alloc1complex(nwtau);
	wwtau = alloc1float(nwtau);
	wwtau2 = alloc1float(nwtau);

	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it+=2)
		pp[it] = p[it];
	for (it=1; it<nt; it+=2) {
		pp[it].r = -p[it].r;
		pp[it].i = -p[it].i;
	}
	for (it=nt; it<nw; ++it)
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
	
	/* interpolate */
	if (0){
	  if (1) fprintf(stderr,"iwtaul=%d,dw=%f,dwtau=%f,fw=%f,nw=%d\n",iwtaul,dw,dwtau,fw,nw);
	  ints8c(nw,dw,fw,pp,czero,czero,iwtauh-iwtaul+1,wwtau+iwtaul,qq+iwtaul);
	}
	else{
	  for (iwtau=iwtaul,wtau=wtaul; iwtau<iwtauh; ++iwtau,wtau+=dwtau)
	    qq[iwtau]=pp[NINT((wwtau[iwtau]-fw)/dw)];
	}
	
	/* Do the same but without interpolation */
	if (0) {
		for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
			qq[iwtau+iwtaul]=pp[NINT(-sqrt(wtau*wtau+vko2s)/dw)];
		for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
			qq[iwtau+iwtaul]=pp[NINT(sqrt(wtau*wtau+vko2s)/dw)];
	}	
        
	/* fft scaling and obliquity factor */
	fftscl = 1.0/nwtau;
	if (s==1.0) {
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh; 
			++iwtau,wtau+=dwtau) {
		        scale = fftscl*wtau/wwtau[iwtau];
			qq[iwtau].r *= scale;
			qq[iwtau].i *= scale;
		}
	} else {
		a = 1.0/(s*s);
		b = 1.0-1.0/s;
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh; 
			++iwtau,wtau+=dwtau) {
			scale = fftscl*(b+a*wtau/(wwtau[iwtau]-b*wtau));
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

	//TRACE;
	/* free workspace */
	free1complex(pp);
	//TRACE;
	free1complex(qq);
	//TRACE;
	free1float(wwtau2);
	//TRACE;
	free1float(wwtau);
	//TRACE;
}


static void taper (int lxtaper, int lbtaper, 
	int nx, int ix, int nt, float *trace)
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

	/* if x tapering is necessary, apply it */
	if (xtaper!=1.0)
		for (it=0; it<nt; ++it)
			trace[it] *= xtaper;
	
	/* if requested, apply t tapering */
	for (it=MAX(0,nt-lbtaper); it<nt; ++it)
		trace[it] *= (0.54+0.46*cos(PI*(lbtaper+it+1-nt)/lbtaper));
}


void save_gather(float **d, int nh, int nt, float dt, char* name)
{
  segy tr;
  
  int  itr;
  FILE* fp;
  
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr.tracl=itr+1;
      tr.dt=(int) (dt*1e6);
      tr.ns=nt;
      tr.ntr=nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  //fflush(fp);
  fclose(fp);

  return;
  
}
