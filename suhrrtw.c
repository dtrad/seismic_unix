/* Copyright (c) Colorado School of Mines, 1998.*/
/* All rights reserved.                       */

/* SUMUTE: $Revision: 1.25 $ ; $Date: 1997/05/05 17:40:30 $	*/
#define NNX 228
#define NT 2048
#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUMUTE - mute above (or below) a user-defined polygonal curve with	", 
"	   the distance along the curve specified by key header word 	",
" 	   								",
" sumute <stdin >stdout xmute= tmute= [optional parameters]		",
" 									",
" Required parameters:							",
" xmute=		array of position values as specified by	",
" 			the `key' parameter				",
" tmute=		array of corresponding time values (sec)	",
"  ... or input via files:						",
" nmute=		number of x,t values defining mute		",
" xfile=		file containing position values as specified by	",
" 			the `key' parameter				",
" tfile=		file containing corresponding time values (sec)	",
" 									",
" Optional parameters:							",
" key=offset		Key header word specifying trace offset 	",
" 				=tracl  use trace number instead	",
" ntaper=0		number of points to taper before hard		",
"			mute (sine squared taper)			",
" below=0		=1 to zero BELOW the polygonal curve		",
" 									",
" Notes: 								",
" The tmute interpolant is extrapolated to the left by the smallest time",
" sample on the trace and to the right by the last value given in the	",
" tmute array.								",
"									",
" The files tfile and xfile are files of binary (C-style) floats.	",
"									",
NULL};

/* Credits:
 *
 *	SEP: Shuki Ronen
 *	CWP: Jack K. Cohen, Dave Hale, John Stockwell
 *	DELPHI: Alexander Koek
 *
 * Trace header fields accessed: ns, dt, delrt, key=keyword
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/

segy tr;
void hrrt_(double *pos,float *data,int *nh, int *nt,double *dt,int *method);

int main(int argc, char **argv)
{

	char *key;	/* header key word from segy.h		*/	
	char *type;	/* ... its type				*/
	int index;	/* ... its index			*/
	Value val;	/* ... its value			*/
	float fval;	/* ... its value cast to float		*/
	float *xmute;	/* array of key mute curve values	*/
	float *tmute;	/* ...		mute curve time values 	*/
	float *taper=NULL;/* ...		taper values	*/
	int nxmute;	/* number of key mute values		*/
	int ntmute;	/* ...		mute time values 	*/
 	int ntaper;	/* ...		taper values		*/
	int below;	/* mute below curve			*/
	int nxtmute;	/* number of mute values 		*/
	cwp_String xfile="";	/* file containing positions by key	*/
	FILE *xfilep;		/* ... its file pointer			*/
	cwp_String tfile="";	/* file containing times	 	*/
	FILE *tfilep;		/* ... its file pointer			*/
	/***************************************************************/
	/*For RADON*/
	
	FILE *myfilep;
	int ii,i;
	float *data;
        double *pos, dtf;
	int ntf, nhf, method;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	pos=ealloc1double(NNX);
	data=ealloc1float(NNX*NT);

	/* Get parameters */
	if (!(getparstring("xfile",&xfile) && getparstring("tfile",&tfile))) {
		if (!(nxmute = countparval("xmute")))
			err("must give xmute= vector");
		if (!(ntmute = countparval("tmute")))
			err("must give tmute= vector");
		if (nxmute != ntmute)
			err("lengths of xmute, tmute must be the same");
		xmute = ealloc1float(nxmute);	getparfloat("xmute", xmute);
		tmute = ealloc1float(nxmute);	getparfloat("tmute", tmute);
	
	} else {
		MUSTGETPARINT("nmute",&nxtmute);
		nxmute = nxtmute;
		xmute = ealloc1float(nxtmute);
		tmute = ealloc1float(nxtmute);

                if((xfilep=fopen(xfile,"r"))==NULL)
                        err("cannot open xfile=%s\n",xfile);
                if (fread(xmute,sizeof(float),nxtmute,xfilep)!=nxtmute)
                        err("error reading xfile=%s\n",xfile);
                fclose(xfilep);

                if((tfilep=fopen(tfile,"r"))==NULL)
                        err("cannot open tfile=%s\n",tfile);
                if (fread(tmute,sizeof(float),nxtmute,tfilep)!=nxtmute)
                        err("error reading tfile=%s\n",tfile);
                fclose(tfilep);
	}

	if (!getparint("ntaper", &ntaper))	ntaper = 0;
	if (!getparint("below", &below))	below = 0;
	if (!getparstring("key", &key))		key = "offset";

	/* get key type and index */
	type = hdtype(key);
	index = getindex(key);

	/* Set up taper weights if tapering requested */
	if (ntaper) {
		register int k;
		taper = ealloc1float(ntaper);
		for (k = 0; k < ntaper; ++k) {
			float s = sin((k+1)*PI/(2*ntaper));
			taper[k] = s*s;
		}
	}

						
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr) err("ntr header field must be set");
	
	for (ii=0;ii<228;ii++)
	                  pos[ii]=0;

	nhf=tr.ntr;
	ntf=tr.dt;
        if((myfilep=fopen("myfile","r"))==NULL)
                        err("cannot open myfile=%s\n","myfile");
	for (ii=0;ii<nhf;ii++){
		for (i=0;i<ntf;i++){
			fscanf(myfilep,"%f \n",&data[i+ii*ntf]);
		}
	}
	fclose(myfilep);

	ii=0;
	/* Loop over traces */
	do {
		int nt     = (int) tr.ns;
		float tmin = tr.delrt/1000.0;
		float dt   = ((double) tr.dt)/1000000.0;
		float t;
		int nmute;
		register int i;
		if (ii==0) {
		    ntf=nt;
		    dtf=dt;
		    method=1;
		}

		/* get value of key and convert to float */
		gethval(&tr, index, &val);
		fval = vtof(type,val);
		

		/* linearly interpolate between (xmute,tmute) values */
		intlin(nxmute,xmute,tmute,tmin,tmute[nxmute-1],1,&fval,&t); 

		/* do the mute */
		if (!below) {	/* mute above */
			nmute = NINT((t - tmin)/dt);
			memset( (void *) tr.data, (int) '\0', nmute*FSIZE);
			for (i = 0; i < ntaper; ++i)
				tr.data[i+nmute] *= taper[i];
			tr.muts = NINT(t*1000);
		} else {	/* mute below */
			nmute = NINT((tmin + nt*dt - t)/dt);
			memset( (void *) (tr.data+nt-nmute),
					(int) '\0', nmute*FSIZE);
			for (i = 0; i < ntaper; ++i)
				tr.data[nt-nmute-1-i] *= taper[i];
			tr.mute = NINT(t*1000);
		}
		
		pos[ii]=tr.offset;
		for (i=0;i<nt;i++){
		       tr.data[i]=data[i+ii*nt];
		}
		ii++;
       		puttr(&tr);
	} while (gettr(&tr));

	return EXIT_SUCCESS;
}
