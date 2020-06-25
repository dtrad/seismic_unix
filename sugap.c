/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUGAP -Forward  High Resolution Parabolic Radon transform           ", 
"	   Program in development                                 	",
" 	   								",
" suhrrt2 < stdin > stdout [optional parameters]          		",
" 									",
" 									",
" Optional parameters:							",
" methdo=1		2 High resolution (Cauchy Gaus method)  	",
" 			1 Normal RT (Gauss- Gauss)			",
" eps1 =3		nosie level                             	",
"                               					",
" iter =5		number of iterations for Cauchy Gauss	        ",
" qmin =0               minimum Radon parameter in sec/offset^2         ",
" 		                           				",
" nq=150		Number of Radon Traces                         	",
" 									",
" Required parameters:		[None]					",
"                                                                	",
" Output : traces with header in the Radon domain.                     	",
" Input : sudata file                                    		",
"		                                        		",
" Example: 		                                                ",
" #Forward  Radon transform                                             ",
" suhrrt2 method=2 eps1=3 iter=5 qmin=-1.e-8 nq=$NP <sudata > sudatarad ", 
"                                                                       ",
" key=offset contains the radon parameter *1e12                         ",
" Other word contain similar values than trace 1         	        ",
"                                                                       ",
"									",
NULL};

/* Credits:
 *
 *      Mauricio Sacchi
 *	Daniel Trad
 *      Tad Ulrych
 *
 * Trace header fields accessed: ns, dt, delrt, key=keyword
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/

segy tr,*trr;
void gapfill_(float *data,int *nt, int *nx,float *dt,float *dx,
int *n1,int *n2, int *ip, float *perc);



int main(int argc, char **argv)
{
	
	FILE *myfilep;
	int ii,i;
	float *data, dtf, dx, perc ; 
	int ntf, nx, nq, n1, n2, ip;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile");

	data=ealloc1float(NNX*NT);

	/* Get parameters */

	if (!getparint("n1", &n1))  n1 = 0;
	if (!getparint("n2", &n2))  n2 = 0;
	if (!getparint("ip", &ip))  ip = 0;
	if (!getparfloat("perc",&perc)) perc=3;

						
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr) err("ntr header field must be set");
	nx=(int) tr.ntr;
	fprintf(stderr,"nx=%d\n",nx);
	 

	if ((trr=malloc(nx*sizeof(segy)))==NULL)
	  fprintf(stderr,"Sorry, space for traces could not be allocated\n");	

	ii=0;

	/* Loop over traces */
	do {

		int nt     = (int) tr.ns;
		float dt   = ((float) tr.dt)/1000000.0;
   		register int i;
		trr[ii]=tr;
		if (ii==0) {
		    ntf=nt;
		    dtf=dt;
		}

		if (ii==2) dx=trr[ii].offset-trr[ii-1].offset;			
		for (i=0;i<nt;i++){
		        data[i+ii*nt]=(float) tr.data[i];
		}
		ii++;
	} while (gettr(&tr));
	nx=ii;	
	fprintf(stderr,"ntf=%d, nx= %d,dtf=%f, dx=%f, n1=%d, n2=%d, ip=%d",
ntf,nx,dtf,dx,n1,n2,ip);
 
	gapfill_(data,&ntf,&nx,&dtf,&dx,&n1,&n2,&ip,&perc);

	
	ii=0;
	do{
	 
	  for (i=0;i<ntf;i++)
	  trr[ii].data[i]=data[i+ntf*ii];
	  puttr(&trr[ii]);
	  ii++;
	 } while(ii<nx);

	 
	free(trr);       
	return EXIT_SUCCESS;
}
