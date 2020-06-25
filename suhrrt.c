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
void hrrt_(double *pos,float *data,int *nh, int *nt,double *dt,int *method,
double *eps1,int *iter_end, double *qmin);

int main(int argc, char **argv)
{
	
	FILE *myfilep;
	int ii,i;
	float *data;
        double *pos, dtf, eps1, qmin;
	int ntf, nhf, method, iter_end;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile");

	pos=ealloc1double(NNX);
	data=ealloc1float(NNX*NT);

	/* Get parameters */


	if (!getparint("method", &method))  method = 1;
	if (!getpardouble("eps1", &eps1))  eps1 = 3;
	if (!getparint("iter_end", &iter_end))  iter_end = 5;
	if (!getpardouble("qmin", &qmin))  qmin = 0;
						
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	
	for (ii=0;ii<228;ii++)
	                  pos[ii]=0;
	ii=0;
	/* Loop over traces */
	do {
		int nt     = (int) tr.ns;
		float dt   = ((double) tr.dt)/1000000.0;
   		register int i;

		if (ii==0) {
		    ntf=nt;
		    dtf=dt;
		}

				
		pos[ii]=tr.offset;
		for (i=0;i<nt;i++){
		        data[i+ii*nt]=tr.data[i];
		}
		ii++;
       		/*puttr(&tr);*/
	} while (gettr(&tr));
	nhf=ii;	
        hrrt_(pos,data,&nhf,&ntf,&dtf,&method,&eps1,&iter_end,&qmin);
	for(ii=0;ii<nhf;ii++){
	    for(i=0;i<ntf;i++){
			fprintf(myfilep,"%f \n",data[i+ii*ntf]);
			
	    }
	}
	
	fclose(myfilep);
       
	return EXIT_SUCCESS;
}
