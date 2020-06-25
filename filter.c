#include "su.h"
#include "segy.h"
char *sdoc[] = {
" 									",
"SUKDMIG2D - Kirchhoff Depth Migration of 2D poststack/prestack data	",
" 									",
"    off0 bin.",
"									",
NULL};
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);
/* segy trace */
segy tr, tro;
int
main (int argc, char **argv)
{
  char *infile="stdin",*outfile="stdout",*ttfile,*jpfile,*tvfile,
    *csfile,*outfile1;
  FILE *infp,*outfp,*ttfp,*jpfp,*tvfp=NULL,*out1fp=NULL,*csfp=NULL;
  int ntr = 10000;
  float dt = 1;
  int nt = 1;
  float fmax = 0;
  int jtr;
  initargs(argc, argv);
  requestdoc(1);
  
  /* open input and output files	*/
  if( !getparstring("infile",&infile)) {
    infp = stdin;
  } else  
    if ((infp=fopen(infile,"r"))==NULL)
      err("cannot open infile=%s\n",infile);
  if( !getparstring("outfile",&outfile)) {
    outfp = stdout;
  } else  
    outfp = fopen(outfile,"w");
  efseeko(infp,(off_t) 0,SEEK_CUR);
  efseeko(outfp,(off_t) 0,SEEK_CUR);
  jtr=1;
  fgettr(infp,&tr);
  if (!getparfloat("dt",&dt)) dt = ((double) tr.dt)/1000000.0; 
  if (!getparint("ntr",&ntr)) ntr = (int) tr.ntr;
  if (!getparfloat("fmax",&fmax)) fmax = 1./(2*dt)*0.6;
  nt = tr.ns;
  do{
    jtr++;
    tro=tr;
    //tro.tracl = tr.tracl;
    filt(tr.data, nt, dt,fmax,1,nt, tro.data);
    fputtr(outfp,&tr); 
    
  } while (fgettr(infp,&tr) && jtr<=ntr);
  for (int j = 0; j< ntr; j++)

  return 0;
}
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf)
/* Low-pass filter, integration and phase shift for input data	 
   input: 
    trace(nt)	single seismic trace
   fmax	high cut frequency
    ls		ls=1, line source; ls=0, point source
  output:
    trace(nt) 	filtered and phase-shifted seismic trace 
    tracei(nt) 	filtered, integrated and phase-shifted seismic trace 
 */
{
	static int nfft=0, itaper, nw, nwf;
	static float *taper, *amp, *ampi, dw;
	int  it, iw, itemp;
	float temp, ftaper, const2, *rt;
	complex *ct;

	fmax *= 2.0*PI;
	ftaper = 0.1*fmax;
	const2 = sqrt(2.0);

	if(nfft==0) {
	  	/* Set up FFT parameters */
	  	nfft = npfaro(nt+m, 2 * (nt+m));
	  	if (nfft >= SU_NFLTS || nfft >= 720720)
		    	err("Padded nt=%d -- too big", nfft);

	  	nw = nfft/2 + 1;
		dw = 2.0*PI/(nfft*dt);

		itaper = 0.5+ftaper/dw;
		taper = ealloc1float(2*itaper+1);
		for(iw=-itaper; iw<=itaper; ++iw){
			temp = (float)iw/(1.0+itaper); 
			taper[iw+itaper] = (1-temp)*(1-temp)*(temp+2)/4;
		}

		nwf = 0.5+fmax/dw;
		if(nwf>nw-itaper-1) nwf = nw-itaper-1;
		amp = ealloc1float(nwf+itaper+1);
		ampi = ealloc1float(nwf+itaper+1);
		amp[0] = ampi[0] = 0.;
		for(iw=1; iw<=nwf+itaper; ++iw){
			amp[iw] = sqrt(dw*iw)/nfft;
			ampi[iw] = 0.5/(1-cos(iw*dw*dt));
		}
	}

	  /* Allocate fft arrays */
	  rt   = ealloc1float(nfft);
	  ct   = ealloc1complex(nw);

	  memcpy(rt, trace, nt*FSIZE);
	  memset((void *) (rt + nt), (int) '\0', (nfft-nt)*FSIZE); 
	  pfarc(1, nfft, rt, ct);

	for(iw=nwf-itaper;iw<=nwf+itaper;++iw){
		itemp = iw-(nwf-itaper);
		ct[iw].r = taper[itemp]*ct[iw].r; 
		ct[iw].i = taper[itemp]*ct[iw].i; 
	}
	for(iw=nwf+itaper+1;iw<nw;++iw){
		ct[iw].r = 0.; 
		ct[iw].i = 0.; 
	}

	 	if(!ls){
		for(iw=0; iw<=nwf+itaper; ++iw){
			/* phase shifts PI/4 	*/
			temp = (ct[iw].r-ct[iw].i)*amp[iw]*const2;
			ct[iw].i = (ct[iw].r+ct[iw].i)*amp[iw]*const2;
			ct[iw].r = temp;
		    }
	} else {
		for(iw=0; iw<=nwf+itaper; ++iw){
			ct[iw].i = ct[iw].i*amp[iw];
			ct[iw].r = ct[iw].r*amp[iw];
		}
	}		  
	  pfacr(-1, nfft, ct, rt);
		
	  /* Load traces back in */
	for (it=0; it<nt; ++it) trace[it] = rt[it];

	  /* Integrate traces   */
	for(iw=0; iw<=nwf+itaper; ++iw){
		ct[iw].i = ct[iw].i*ampi[iw];
		ct[iw].r = ct[iw].r*ampi[iw];
	}
	  pfacr(-1, nfft, ct, rt);
	  for (it=0; it<m; ++it)  trf[it] = rt[nfft-m+it];
	  for (it=0; it<nt+m; ++it)  trf[it+m] = rt[it];

	free1float(rt);
	free1complex(ct);
}
