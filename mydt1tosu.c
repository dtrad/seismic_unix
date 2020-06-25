/* Copyright (c) Colorado School of Mines, 2000.*/
/* All rights reserved.                       */

/* DT1TOSU: $Revision: 1.13 $ ; $Date: 2000/09/22 16:15:31 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									 ",
" DT1TOSU - Convert ground-penetrating radar data in the Sensors & Software",
"	    X.dt1 GPR format to SU format.				",
"									",
" dt1tosu < gpr_data_in_dt1_format  > stdout				",
"									",
" Optional parameters:							",
" ns=from header	number of samples per trace			",
" dt=.8			time sample interval (see below)		",
" swap=1		perform byte swapping (big endian machines)	",
"			=0 don't swap bytes (little endian machines)	",
" verbose=0		silent						",
"			=1 S & S header values from first trace		",
"				sent to outpar				",
"			=2 S & S header values from all traces		",
"				sent to outpar				",
" outpar=/dev/tty	output parameter file				",
" list=0		silent						",
"			=1 list explaining labels used in verbose	",
"			     is printed to stderr			",
"									",
" Caution: An incorrect ns field will munge subsequent processing.	",
"									",
" Notes:								",
" For compatiblity with SEGY header, apparent dt is set to		",
" .8 ms (800 microsecs).  Actual dt is .8 nanosecs.			",
" Using TRUE DISTANCES, this scales velocity				",
" and frequency by a factor of 1 million.				",
"	Example: v_air = 9.83X10^8 ft/s	 (real)				",
"		 v_air = 983 ft/s	(apparent for su)		",
"		Example: fnyquist = 625 MHz	(real)			",
"			fnyquist = 625 Hz	(apparent for su)	",
"									",
" IBM RS6000, NeXT, SUN are examples of big endian machines		",
" PC's and DEC are examples of little endian machines			",
"									",
" Caveat:								",
" This program has not been tested on DEC, some modification of the	",
" byte swapping routines may be required.				",
"									",
NULL};

/* Credits:
 *	CWP: John Stockwell, Jan 1994   Based on a code "sugpr" by
 *	UTULSA: Chris Liner & Bill Underwood  (Dec93)
 * modifications permit S & S dt1 header information to be transferred
 * directly to SU header
 *
 * Trace header fields set: ns, tracl, tracr, dt, delrt, trid,
 *			    hour, minute, second
 */
/**************** end self doc *******************************************/

/* define hed structure */
#define SSHDRBYTES 128 /* size of ssdt1 */
typedef struct {
	float  tracr;	/* trace number */
	float  posit;	/* position */
	float  ns;	/* number of points per trace */
	float  topog;	/* Topographic data if any */
	float  nua;	/* not used "a" */ 
	float  bpp;	/* bytes/point (always 2 for Rev 3 firmware) */ 
	float  tracl;	/* second trace number */ 
	float  nstks;	/* number of stacks */ 
	float  twind;	/* time window */ 
	float  nub[11];	/* not used "b" */ 
	float  tza;	/* time zero adjustment */ 
				/* where point(x) = point(x + adjustment) */
	float  zflag;	/* zero flag: 0=data ok ,  1 = zero data */
	float  nuc;	/* not used "c" */ 
	float  tod;	/* time of day data collected */
				/* in secs past midnight */ 
	float  cflag;	/* comment flag =1 if comment is attached */
	char comm[28];	/* comment 28 characters */
	short data[SU_NFLTS];  /* use SU maximum number data values */
} ssdt1;

ssdt1 sstr;
segy tr;

/* list explaining the ssdt1 convention */
char *list[] = {
" float	 tracr;	first trace counter				",
" float	 posit;	position					",
" float	 ns;		number of points per trace		",
" float	 topog;	Topographic data if any				",
" float	 nua;		not used \"a\"				",
" float	 bpp;		bytes/point (always 2 for Rev 3 firmware)",
" float	 tracl;	second trace counter				",
" float	 nstks;	number of stacks				",
" float	 twind;	time window					",
" float	 nub[11];	not used \"b\"				",
" float	 tza;		time zero adjustment			",
"			where point(x) = point(x + adjustment)	",
" float	 zflag;	zero flag: 0=data ok ,	1 = zero data		",
" float	 nuc;		not used \"c\"				",
" float	 tod;		time of day data are collected		",
"				in secs past midnight		",
" float	 cflag;	comment flag =1 if comment is attached		",
" char comm[28];	comment 28 characters			",
" short data[SU_NFLTS];	 data part of the ssdt1 structure	",
NULL};

/* pointer to list */
char **listptr=list;

/* function prototypes for internally defined subroutines */
void swap_ss_sstr(void);
void fprintf_sstr(FILE *outparfp);


int
main(int argc, char **argv)
{
	int i,ns;		/* counter, number of samples */
	float dt;		/* time sample interval */	
	size_t databytes;
	size_t ssdatabytes;
	size_t ssbytes;
	int swap,verbose,list;	/* flags */
	float *data;
	float hour,minute;	/* hour and minute */
	char *outpar;		/* outpar filename */
	FILE *outparfp;		/* outpar file pointer */


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Make sure stdout is a file or pipe */
	switch(filestat(STDOUT)) {
	case TTY:
		err("stdout can't be tty");
	break;
	case DIRECTORY:
		err("stdout must be a file, not a directory");
	break;
	case BADFILETYPE:
		err("stdout is illegal filetype");
	break;
	default: /* rest are OK */
	break;
	}

	/* Read header information from first header */
	efread(&sstr, 1, SSHDRBYTES, stdin);
	if (!getparint("swap", &swap))		swap = 1;
	if (swap) swap_ss_sstr();

	/* Get parameters */
	if (!getparint("ns", &ns)) 	 ns = sstr.ns ;
	if (!ns) err ("ns not set!");

	if (!getparfloat("dt", &dt))		dt= .8 ;
	dt *= 1000.;
	if (!getparint("list", &list))		list= 0;
	if (!getparint("verbose", &verbose))	verbose= 0;
	if (!getparstring("outpar", &outpar))	outpar = "/dev/tty" ;

	outparfp = efopen(outpar, "w");

	/* print out ss header field for first trace */
	if (verbose==1) fprintf_sstr(outparfp);
	if (list) while (*listptr) fprintf(stderr,"%s\n", *listptr++);

	/* assign tr header values */
	tr.tracr = sstr.tracr;
	tr.ns = ns;	
	tr.tracl = sstr.tracl;
	tr.dt = NINT(dt); 
	tr.delrt = sstr.tza;
	tr.trid = sstr.zflag + 1;

	/* compute hour and minute on first trace */
	hour = sstr.tod/3600;
	minute = (hour - (int) hour) * 60;

	/* set hour, minute, and second on first trace */
	tr.hour = hour;
	tr.minute = minute;
	tr.sec = ( minute - (int) minute ) * 60;

	/* assign size of data portion on first trace */
	databytes = ns*FSIZE;
	ssdatabytes = ns*sizeof(short);
	
	/* allocate space for the data array */
	data = alloc1float(databytes);

	/* read in and swap the data portion of the first trace */
	/* loop over data values in X.dt1 data */
	for (i=0 ; i < ns ; i++) {
		short temp;	/* temporary variable to store data vals in */
		efread(&temp, 1, sizeof(short), stdin);

		if (swap) swap_short_2(&temp);

		data[i] = (float) temp;
	}

	/* write out the first trace */
	efwrite(&tr, 1, HDRBYTES, stdout);
	efwrite(data, 1, databytes, stdout);

	/* total size of an ss trace */
	ssbytes = SSHDRBYTES + ssdatabytes;

	/* load selected sstr values and data into tr */
	while (efread(&sstr,1,ssbytes,stdin)) {

		/* swap header values */
		if (swap) swap_ss_sstr();

		/* set trace header fields */
		tr.tracr = sstr.tracr;
		tr.ns = ns;	
		tr.tracl = sstr.tracl;
		tr.dt = NINT(dt); 
		tr.delrt = sstr.tza;
		tr.trid = sstr.zflag + 1;

		/* compute hour and minute */
		hour = sstr.tod/3600;
		minute = (hour - (int) hour) * 60;

		/* set hour, minute, and second */
		tr.hour = hour;
		tr.minute = minute;
		tr.sec = ( minute - (int) minute ) * 60;

		/* swap values in ss data */
		if (swap) for (i=0 ; i < ns ; i++) swap_short_2(&sstr.data[i]);

		/* assign data values to tr.data[] */
		for (i=0; i<ns; i++) tr.data[i] = (float) sstr.data[i];
		if (verbose==2) fprintf_sstr(outparfp);
		fprintf(stderr,"trace # %d, ns=%d\n ",tr.tracr,ns);
		puttr(&tr);
	}
	if (verbose) efclose(outparfp);
	return EXIT_SUCCESS;

}


/* swap */
void swap_ss_nub(float *nub);

void swap_ss_sstr(void)
{
	swap_float_4(&sstr.tracr);
	swap_float_4(&sstr.posit);
	swap_float_4(&sstr.ns);
	swap_float_4(&sstr.topog);
	swap_float_4(&sstr.nua);
	swap_float_4(&sstr.bpp);
	swap_float_4(&sstr.tracl);
	swap_float_4(&sstr.nstks);
	swap_float_4(&sstr.twind);
	swap_ss_nub(sstr.nub);
	swap_float_4(&sstr.tza);
	swap_float_4(&sstr.zflag);
	swap_float_4(&sstr.nuc);
	swap_float_4(&sstr.tod);
}

void swap_ss_nub(float *nub)
{
	int i;
	for ( i=0 ; i < 11; i++)
		/* swap_float_4(&sstr.nub[i]);*/
		swap_float_4(nub+i);
}


void fprintf_sstr(FILE *outparfp)
{ /* send ssdt1 values to outpar */
	int i;		/* counter */
	int tracr=sstr.tracr; /* first trace counter */
	float  posit=sstr.posit; /* position */
	int ns=sstr.ns;	 /* number of samples */
	float  topog=sstr.topog; /* topographic data */
	float nua=sstr.nua;	/* unused "a" */
	int bpp=sstr.bpp;	/* bytes per point */ 
	int tracl=sstr.tracl;	/* second trace counter */
	int nstks=sstr.nstks;	/* number of stacks */
	float  twind=sstr.twind; /* time window */
	float nub;		/* unused "b" */
	float  tza=sstr.tza;	/* time zero adjustment */
	int  zflag=sstr.zflag;	/* zero flag */
	float nuc=sstr.nuc; /* unused "c" */
	int  cflag=sstr.cflag; /* comment flag */


	fprintf(outparfp, "tracr = %d\n", tracr);
	fprintf(outparfp, "posit = %.3f\n", posit);
	fprintf(outparfp, "ns = %d\n", ns);
	fprintf(outparfp, "topog = %.3f\n", topog);
	fprintf(outparfp, "nua = %.3f\n", nua);
	fprintf(outparfp, "bpp = %d\n", bpp);
	fprintf(outparfp, "tracl = %d\n", tracl);
	fprintf(outparfp, "nstks = %d\n", nstks);
	fprintf(outparfp, "twind = %.3f\n", twind);
	for(i=0; i < 11 ; i++){
		nub = sstr.nub[i];		
		fprintf(outparfp, "nub[%d] = %.3f\n",i,nub);
	}
	fprintf(outparfp, "tza = %.3f\n", tza);
	fprintf(outparfp, "zflag = %d\n", zflag);
	fprintf(outparfp, "nuc = %.3f\n", nuc);
	fprintf(outparfp, "cflag = %d\n", cflag);
	if(cflag) fprintf(outparfp, "comm [%d] = %s\n",i,sstr.comm);
	else      fprintf(outparfp, "no comm\n");
}


