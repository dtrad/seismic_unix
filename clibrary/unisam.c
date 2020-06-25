/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* UNISAM: $Revision: 1.8 $ ; $Date: 1997/02/14 17:16:06 $	*/

#include "par.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 									",
" UNISAM - UNIformly SAMple a function y(x) specified as x,y pairs	",
" 									",
"   unisam xin= yin= nout= [optional parameters] >binaryfile		",
"    ... or ...								",
"   unisam xfile= yfile= npairs= nout= [optional parameters] >binaryfile",
"    ... or ...								",
"   unisam xyfile= npairs= nout= [optional parameters] >binaryfile	",
" 									",
" Required Parameters:							",
" xin=,,,	array of x values (number of xin = number of yin)	",
" yin=,,,	array of y values (number of yin = number of xin)	",
"  ... or								",
" xfile=	binary file of x values					",
" yfile=	binary file of y values					",
"  ... or								",
" xyfile=	binary file of x,y pairs				",
" npairs=	number of pairs input (active only if xfile= and yfile=	",
" 		or xyfile= are set)					",
" 									",
" nout=		 number of y values output to binary file		",
" 									",
" Optional Parameters:							",
" dxout=1.0	 output x sampling interval				",
" fxout=0.0	 output first x						",
" method=linear  =linear for linear interpolation (continuous y)	",
"		 =mono for monotonic cubic interpolation (continuous y')",
"		 =akima for Akima's cubic interpolation (continuous y') ",
"		 =spline for cubic spline interpolation (continuous y'')",
" isint=,,,	 where these sine interpolations to apply		",
" amp=,,,	 amplitude of sine interpolations			",
" phase0=,,,	 starting phase (defaults: 0,0,0,...,0)			",
" totalphase=,,, total phase (default pi,pi,pi,...,pi.)			",
" nwidth=0       apply window smoothing if nwidth>0                     ",
" 									",
NULL};
/*
 * AUTHOR:  Dave Hale, Colorado School of Mines, 07/07/89
 *          Zhaobo Meng, Colorado School of Mines, 
 * 	    added sine interpolation and window smoothing, 09/16/96 
 *          CWP: John Stockwell,  added file input options, 24 Nov 1997
 *
 * Remarks: In interpolation, suppose you need 2 pieces of 
 * 	    sine interpolation before index 3 to 4, and index 20 to 21
 *	    then set: isint=3,20. The sine interpolations use a sine
 *	    function with starting phase being phase0, total phase 
 *	    being totalphase (i.e. ending phase being phase0+totalphase
 *	    for each interpolation).
 * 	    
 */
/**************** end self doc ********************************/

int
main(int argc, char **argv)
{
	int iout;		/* counters for values output		*/
	int nin=0;		/* number of values read in		*/
	int nout;		/* number of y values output		*/
	float dxout;		/* x interval in output			*/
	float fxout;		/* first x value output			*/
	float *xin=NULL,*yin=NULL;/* pointer to input x and y values	*/
	float *xyin;		/* pointer to input x,y pairs		*/
	float (*yind)[4];	/* array used for cubic splines		*/
	float *xout,*yout;	/* arrays of output x and y values	*/
	char *method="linear";	/* method of interpolation		*/
	FILE *outfp=stdout;	/* output file pointer			*/

	/* added by Z. Meng */
	int nsint;		/* number of sine interpolations to apply */
	int *isint=NULL;	/* the indices where sine interpolations to apply */
	float *amp=NULL;	/* the amplitude of sine interpolations */   
	float *phase0=NULL;	/* initial phase for each interpolation */
	float *totalphase=NULL;	/* total phase for sine interpolation */
	int ix0,ix1;    	/* indices for left and right sine int points */
	float denorm;   	 /* denominator used		*/
	int ix,is,ip;   	 /* work indices		*/
        int nwidth;     	 /* number of samples in input	*/
        float *yout0;		 /* pointer to temporary array */
	float sum;

	/* xfile, yfile, and xyfile  stuff */
	char *xfile="";		/* name of file of x values 	*/
	char *yfile="";		/* name of file of y values 	*/
	char *xyfile="";	/* name of file of x,y pairs	*/

	FILE *xfilefp=NULL;	/* pointer to xfile */
	FILE *yfilefp=NULL;	/* pointer to yfile */
	FILE *xyfilefp=NULL;	/* pointer to xyfile */

	int npairs;		/* number of x,y pairs read from files */

	/* hook up getpar */
	initargs(argc,argv);
	requestdoc(0);


	/* get parameters */
	getparstring("xfile",&xfile);
	getparstring("yfile",&yfile);
	getparstring("xyfile",&xyfile);

	if ( (*xfile=='\0') && (*yfile=='\0') && (*xyfile=='\0')) {
		xin = alloc1float(countparval("xin"));
		yin = alloc1float(countparval("xin"));
		if ((nin=getparfloat("xin",xin))==0)
			err("Must specify xin!");
		if (getparfloat("yin",yin)!=nin) 
			err("Number of yins must equal number of xins!");
	} else if (*xfile!='\0' && *yfile!='\0') {
		MUSTGETPARINT("npairs",&npairs);
		nin=npairs;

		/* open x file */
		if((xfilefp=efopen(xfile,"r"))==NULL)
                        err("cannot open xfile=%s\n",xfile); 

		/* allocate space for x values */
		xin = ealloc1float(npairs);

		/* read x values */
		if (fread(xin,sizeof(float),npairs,xfilefp)!=npairs)
                        err("error reading xfile=%s\n",xfile);

		/* open y file */
		if((yfilefp=efopen(yfile,"r"))==NULL)
                        err("cannot open yfile=%s\n",yfile); 

		/* allocate space for y values */
		yin = ealloc1float(npairs);

		/* read y values */
		if (fread(yin,sizeof(float),npairs,yfilefp)!=npairs)
                        err("error reading yfile=%s\n",yfile);


	} else if (*xyfile!='\0') {
		MUSTGETPARINT("npairs",&npairs);
		nin=npairs;

		/* open xyfile */
		if((xyfilefp=fopen(xyfile,"r"))==NULL)
                        err("cannot open xyfile=%s\n",xyfile); 

		/* allocate space for x,y pairs */
		xin = ealloc1float(npairs);
		yin = ealloc1float(npairs);
		xyin = ealloc1float(2*npairs);

		/* read x,y pairs */
		if (fread(xyin,sizeof(float),2*npairs,xyfilefp)!=2*npairs)
                        err("error reading xyfile=%s\n",xyfile);

		/* read x and y values into separate arrays */
		for (ix=0; ix<2*npairs ; ++ix) {

			if( !(ISODD(ix)))
				xin[ix/2]=xyin[ix];
			else
				yin[ix/2]=xyin[ix];
		}
	}


	if (!getparint("nout",&nout))
		err("Must specify nout!");
	dxout = 1.0;  getparfloat("dxout",&dxout);
	fxout = 0.0;  getparfloat("fxout",&fxout);
	getparstring("method",&method);

	/* allocate space for output */
	xout = ealloc1float(nout);
	yout = ealloc1float(nout);

	/* compute uniformly sampled xs */
	for (iout=0; iout<nout; iout++)
		xout[iout] = fxout+iout*dxout;

	/* if linear interpolation or only one input sample */
	if (method[0]=='l' || nin==1) {
		intlin(nin,xin,yin,yin[0],yin[nin-1],nout,xout,yout);

	/* else, if monotonic interpolation */
	} else if (method[0]=='m') {
		yind = (float (*)[4])ealloc1float(nin*4);
		cmonot(nin,xin,yin,yind);
		intcub(0,nin,xin,yind,nout,xout,yout);

	/* else, if Akima interpolation */
	} else if (method[0]=='a') {
		yind = (float (*)[4])ealloc1float(nin*4);;
		cakima(nin,xin,yin,yind);
		intcub(0,nin,xin,yind,nout,xout,yout);

	/* else, if cubic spline interpolation */
	} else if (method[0]=='s') {
		yind = (float (*)[4])ealloc1float(nin*4);;
		csplin(nin,xin,yin,yind);
		intcub(0,nin,xin,yind,nout,xout,yout);

	/* else, if unknown method specified */
	} else {
		err("%s is an unknown interpolation method!\n",method);
	}

	/* added by Z. Meng */
	if (!(nsint=countparval("isint"))) nsint=0;
        if (nsint>0 && nsint<nin-1) {
		isint=alloc1int(nsint);
		amp=alloc1float(nsint);
		phase0=alloc1float(nsint);
		totalphase=alloc1float(nsint);
	
		if (!getparint("isint",isint))
			err("Must specify %d isint values",nsint);
		if (!getparfloat("amp",amp))
			err("Must specity %d amp values",nsint);

		if (!getparfloat("phase0",phase0))
			for (is=0;is<nsint;is++)
				phase0[is]=0.0;
	
		if (!getparfloat("totalphase",totalphase))
			for (is=0;is<nsint;is++)
				totalphase[is]=PI;
		for (is=0;is<nsint;is++)
			if (totalphase[is]==0.0) 
				totalphase[is]=PI;	
	}

        for (is=0;is<nsint;is++) {
                if (isint[is]+1>=nin-1) break;
		ix0=(xin[isint[is]]-fxout)/dxout;
		ix1=(xin[isint[is]+1]-fxout)/dxout;
                denorm=(xin[isint[is]+1]-xin[isint[is]])/totalphase[is];
                for (ix=ix0+1;ix<ix1;ix++)
                        yout[ix]+=amp[is]*
                        sin(phase0[is]+(dxout*ix-xin[isint[is]]+fxout)/denorm);
        }

        if (!getparint("nwidth",&nwidth)) nwidth=0;
        nwidth=nwidth/2*2+1;  /*make it odd*/
        if (nwidth>1) {
                yout0=alloc1float(nout);
                for (ix=nwidth/2;ix<nout-nwidth/2;ix++) {
                        sum=0;
                        for (ip=-nwidth/2;ip<=nwidth/2;ip++)
                                sum+=yout[ix+ip];
                        yout0[ix]=sum/nwidth;
                }
                for (ix=nwidth/2;ix<nout-nwidth/2;ix++)
                        yout[ix]=yout0[ix];
        }

	/* output */
	efwrite(yout,sizeof(float),nout,outfp);

	return EXIT_SUCCESS;
}
