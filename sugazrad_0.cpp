/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* SUGAZMIG: $Revision: 1.16 $ ; $Date: 1997/07/28 22:36:46 $	    */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
#include "gazmig.h"
#include "complex.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" SUGAZMIG - SU version of Jeno GAZDAG's phase-shift migration 		",
"	     for zero-offset data.					",
"									",
" sugazmig <infile >outfile vfile= [optional parameters]		",
"									",
" Optional Parameters:							",
" dt=from header(dt) or	.004	time sampling interval			",
" dx=from header(d2) or 1.0	midpoint sampling interval		",
" ft=0.0			first time sample			",
" ntau=nt(from data)	number of migrated time samples			",
" dtau=dt(from header)	migrated time sampling interval			",
" ftau=ft		first migrated time sample			",
" tmig=0.0		times corresponding to interval velocities in vmig",
" vmig=1500.0		interval velocities corresponding to times in tmig",
" vfile=		name of file containing velocities		",
"									",
" verbose=0	verbose = 1 echoes information				",
"									",
" tmpdir= 	 if non-empty, use the value as a directory path	",
"		 prefix for storing temporary files; else if the	",
"	         the CWP_TMPDIR environment variable is set use		",
"	         its value for the path; else use tmpfile()		",
"									",
" Note: ray bending effects not accounted for in this version.		",
"									",
" The tmig and vmig arrays specify an interval velocity function of time.",
" Linear interpolation and constant extrapolation is used to determine	",
" interval velocities at times not specified.  Values specified in tmig	",
" must increase monotonically.						",
"									",
" Alternatively, interval velocities may be stored in a binary file	",
" containing one velocity for every time sample in the data that is to be",
" migrated.  If vfile is specified, then the tmig and vmig arrays are ignored.",
"									",
NULL};

/* 
 * Credits: CWP John Stockwell 12 Oct 1992
 * 	Based on a constant v version by Dave Hale.
 * 
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */ 
/**************** end self doc ********************************/

/* prototypes for functions defined and used below */
void gazdagvt (float k,
	int nt, float dt, float ft,
	int ntau, float dtau, float ftau,
	float *vt, complex *p, complex *q);
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

segy tr;

int
main (int argc, char **argv)
{
	int nt;			/* number of time samples */
	int ntau;		/* number of migrated time samples */
	int nx;			/* number of midpoints 	*/
	int ik,ix,it,itau,itmig;/* loop counters 	*/
	int nxfft;		/* fft size		*/
	int nk;			/* number of wave numbers */	

	int ntmig,nvmig;	

	float dt;		/* time sampling interval 	*/
	float ft;		/* first time sample		*/
	float dtau;		/* migrated time sampling interval */
	float ftau;		/* first migrated time value 	*/
	float dk;		/* wave number sampling interval */
	float fk;		/* first wave number 		*/
	float t,k;		/* time,wave number		*/
	float *tmig, *vmig;	/* arrays of time, interval velocities */
	float dx;		/* spatial sampling interval	*/
	float *vt;		/* velocity v(t)		*/
	float **p,**q;		/* input, output data		*/

	complex **cp,**cq;	/* complex input,output		*/

	char *vfile="";		/* name of file containing velocities */
	int verbose;		/* flag for echoing info		*/
	char *tmpdir;		/* directory path for tmp files		*/
	cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/
	/* Changes for radon */
        float **q2;
        int nrad;
        int irad;
	float dv;

	/* hook up getpar to handle the parameters */
	initargs(argc,argv);
	requestdoc(1);

	/* get info from first trace */
	if (!gettr(&tr))  err("can't get first trace");
	nt = tr.ns;

	/* let user give dt and/or dx from command line */
	if (!getparfloat("dt", &dt)) {
		if (tr.dt) { /* is dt field set? */
			dt = ((double) tr.dt)/1000000.0;
		} else { /* dt not set, assume 4 ms */
			dt = 0.004;
			warn("tr.dt not set, assuming dt=0.004");
		}
	}
	if (!getparfloat("dx",&dx)) {
		if (tr.d2) { /* is d2 field set? */
			dx = tr.d2;
		} else {
			dx = 50.0;
			warn("tr.d2 not set, assuming dx=50.0 (offset sampling)");
		}
	}


	/* get optional parameters */
	if (!getparfloat("ft",&ft)) ft = 0.0;
	if (!getparint("ntau",&ntau)) ntau = nt; CHECK_NT("ntau",ntau);
	if (!getparfloat("dtau",&dtau)) dtau = dt;
	if (!getparfloat("ftau",&ftau)) ftau = ft;
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (!getparint("nrad", &nrad)) nrad = 20;
	if (!getparfloat("dv", &dv)) dv = 200; /* velocity sampling */

	/* Look for user-supplied tmpdir */
	if (!getparstring("tmpdir",&tmpdir) &&
	    !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
	if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
		err("you can't write in %s (or it doesn't exist)", tmpdir);


	/* store traces and headers in tempfiles while getting a count */
	if (STREQ(tmpdir,"")) {
		tracefp = etmpfile();
		headerfp = etmpfile();
		if (verbose) warn("using tmpfile() call");
	} else { /* user-supplied tmpdir */
		char directory[BUFSIZ];
		strcpy(directory, tmpdir);
		strcpy(tracefile, temporary_filename(directory));
		strcpy(headerfile, temporary_filename(directory));
		/* Trap signals so can remove temp files */
		signal(SIGINT,  (void (*) (int)) closefiles);
		signal(SIGQUIT, (void (*) (int)) closefiles);
		signal(SIGHUP,  (void (*) (int)) closefiles);
		signal(SIGTERM, (void (*) (int)) closefiles);
		tracefp = efopen(tracefile, "w+");
		headerfp = efopen(headerfile, "w+");
      		istmpdir=cwp_true;		
		if (verbose) warn("putting temporary files in %s", directory);
	}

	nx = 0;
	do {
		 ++nx;
		efwrite(&tr,HDRBYTES,1,headerfp);
		efwrite(tr.data, FSIZE, nt, tracefp);
	} while (gettr(&tr));
	erewind(tracefp);
	erewind(headerfp);
	
	/* determine wavenumber sampling (for real to complex FFT) */
	nxfft = npfar(nx);
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);
	fk = 0.0;
	
	/* allocate space */
	p = alloc2float(nt,nxfft);
	q = alloc2float(ntau,nxfft);
	q2 = alloc2float(ntau,nrad);
	cp = alloc2complex(nt,nk);
	cq = alloc2complex(ntau,nk);

	/* load traces into the zero-offset array and close tmpfile */
	efread(*p, FSIZE, nt*nx, tracefp);
	efclose(tracefp);

	/* determine velocity function v(t) */
	vt = ealloc1float(ntau);
	if (!getparstring("vfile",&vfile)) {
		ntmig = countparval("tmig");
		if (ntmig==0) ntmig = 1;
		tmig = ealloc1float(ntmig);
		if (!getparfloat("tmig",tmig)) tmig[0] = 0.0;
		nvmig = countparval("vmig");
		if (nvmig==0) nvmig = 1;
		if (nvmig!=ntmig) err("number of tmig and vmig must be equal");
		vmig = ealloc1float(nvmig);
		if (!getparfloat("vmig",vmig)) vmig[0] = 1500.0;
		for (itmig=1; itmig<ntmig; ++itmig)
			if (tmig[itmig]<=tmig[itmig-1])
				err("tmig must increase monotonically");
		for (it=0,t=0.0; it<ntau; ++it,t+=dt)
			intlin(ntmig,tmig,vmig,vmig[0],vmig[ntmig-1],
				1,&t,&vt[it]);
	} else {
		if (fread(vt,sizeof(float),nt,fopen(vfile,"r"))!=nt)
			err("cannot read %d velocities from file %s",nt,vfile);
	}
	/* set to zero Radon space */
	for (irad=0; irad<nrad; irad++)
	  for (itau=0; itau < ntau; itau++)
	    q2[irad][itau]=0;

	fprintf(stderr,"irad=%d,nrad=%d,dv=%f\n",irad,nrad,dv);

	/* pad with zeros and Fourier transform x to k */
	for (ix=nx; ix<nxfft; ix++)
		for (it=0; it<nt; it++)
			p[ix][it] = 0.0;
	pfa2rc(-1,2,nt,nxfft,p[0],cp[0]); /* cp is data(kx,t) */
	
	for (irad=0; irad<nrad; irad++){

	  /* migrate each wavenumber */
	  for (ik=0,k=fk; ik<nk; ik++,k+=dk)
	    gazdagvt(k,nt,dt,ft,ntau,dtau,ftau,vt,cp[ik],cq[ik]);
  
	  /* Fourier transform k to x (including FFT scaling) */
	  pfa2cr(1,2,ntau,nxfft,cq[0],q[0]);
	  for (ix=0; ix<nx; ix++)
	    for (itau=0; itau<ntau; itau++)
	      q[ix][itau] /= nxfft;

	  /* 
	     This part of the code will be modified to stack the 
	     migrated traces  and output a single trace with one velocity
	  */
	  
	  if (1){
	    save_gather(q,nx,ntau,0.004,"migrated.su");
	    system("suxwigb < migrated.su perc=100 key=offset  title=migrated &");
	  }
	  
	  /* Stack migrated traces in trace irad */
	  
	  for (itau=0; itau < ntau; itau++)
	    for (ix=(nx-1)/2; ix<(nx+1)/2; ix++)                
	      q2[irad][itau]+=q[ix][itau];

	  /* Increment velocity for next migration */
	  for (itau=0; itau < ntau; itau++){
	    vt[itau]+=dv;
	    fprintf(stderr,"vt[%d]=%f\n",itau,vt[itau]);
	  }	
	  fprintf(stderr,"irad=%d,nrad=%d\n",irad,nrad);	
	}     
	/* restore header fields and write output */
	for (irad=0; irad<nrad; ++irad) {
                tr.f2=(int) irad;
                tr.ntr=(int) nrad;     
		/*efread(&tr,HDRBYTES,1,headerfp);*/
		tr.ns = ntau ;
		tr.dt = dtau * 1000000.0 ;
		tr.delrt = ftau * 1000.0 ;
		memcpy( (void *) tr.data, (const void *) q2[irad],ntau*FSIZE);
		puttr(&tr);
	}
	
	/* Clean up */
	efclose(headerfp);
	if (istmpdir) eremove(headerfile);
	if (istmpdir) eremove(tracefile);
        free2float(q2);
	return EXIT_SUCCESS;	
}

void gazdagvt (float k, 
	int nt, float dt, float ft, 
	int ntau, float dtau, float ftau,
	float *vt, complex *p, complex *q)
/*****************************************************************************
Gazdag's phase-shift zero-offset migration for one wavenumber
adapted to v(tau) velocity profile
******************************************************************************
Input:
k		wavenumber
nt		number of time samples
dt		time sampling interval
ft		first time sample
ntau		number of migrated time samples
dtau		migrated time sampling interval
ftau		first migrated time sample
vt		velocity v[tau]
p		array[nt] containing data to be migrated

Output:
q		array[ntau] containing migrated data
******************************************************************************/
{
	int ntfft,nw,it,itau,iw;
	float dw,fw,tmax,w,tau,phase,coss;
	complex cshift,*pp;

	/* determine frequency sampling */
	ntfft = npfa(nt);
	nw = ntfft;
	dw = 2.0*PI/(ntfft*dt);
	fw = -PI/dt;
	
	/* determine maximum time */
	tmax = ft+(nt-1)*dt;

	/* allocate workspace */
	pp = alloc1complex(nw);
	
	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it++)
		pp[it] = (it%2 ? -(p[it]) : p[it]);
	for (it=nt; it<ntfft; it++)
		pp[it] = cmplx(0.0,0.0);
	pfacc(1,ntfft,pp);
	
	/* account for non-zero ft and non-zero ftau */
	for (itau=0 ; itau < ftau ; itau++){
		for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
			if (w==0.0) w = 1e-10/dt;
			coss = 1.0-pow(0.5 * vt[itau] * k/w,2.0);
			if (coss>=pow(ftau/tmax,2.0)) {
				phase = w*(ft-ftau*sqrt(coss));
				cshift = cmplx(cos(phase),sin(phase));
				pp[iw] = cmul(pp[iw],cshift);
			} else {
				pp[iw] = cmplx(0.0,0.0);
			}
		}
	}
	
	/* loop over migrated times tau */
	for (itau=0,tau=ftau; itau<ntau; itau++,tau+=dtau) {
		
		/* initialize migrated sample */
		q[itau] = cmplx(0.0,0.0);
		
		/* loop over frequencies w */
		for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
			
			/* accumulate image (summed over frequency) */
			q[itau] = cadd(q[itau],pp[iw]);
			
			/* compute cosine squared of propagation angle */
			if (w==0.0) w = 1e-10/dt;
			coss = 1.0-pow(0.5 * vt[itau] * k/w,2.0);
			
			/* if wave could have been recorded in time */
			if (coss>=pow(tau/tmax,2.0)) {
			
				/* extrapolate down one migrated time step */
				phase = -w*dtau*(sqrt(coss));
				cshift = cmplx(cos(phase),sin(phase));
				pp[iw] = cmul(pp[iw],cshift);
			
			/* else, if wave couldn't have been recorded in time */
			} else {
				
				/* zero the wave */
				pp[iw] = cmplx(0.0,0.0);
			}
		}
		
		/* scale accumulated image just as we would for an FFT */
		q[itau] = crmul(q[itau],1.0/nw);
	}
		
	/* free workspace */
	free1complex(pp);	
	
}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}











