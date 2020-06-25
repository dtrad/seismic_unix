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
" SUHRRTF -Forward  High Resolution Radon transform                     ", 
"	   Program in development                                 	",
" 	   								",
" suhrrt < stdin > stdout [optional parameters]          		",
" 									",
" 									",
" Optional parameters:							",
" methdo=1		2 High resolution (Cauchy Gaus method)  	",
" 			1 Normal RT (Gauss- Gauss)			",
" eps1 =3		nosie level                             	",
"                               					",
" iter =5		number of iterations for Cauchy Gauss		",
" qmin =0               minimum Radon parameter in sec/offset^2        	",
" 		                           				",
" nq=150		Number of Radon Traces                         	",
" 									",
" Required parameters:		[None]					",
"                                                                	",
" Output : Bare traces in the Radon domain.                         	",
" Values of Radon paramtere are in binary (C) file myfile		",
"		                                        		",
" Example: 		                                                ",
" #Forward  Radon transform                                             ",
" suhrrtf method=2 eps1=3 iter=5 qmin=-1.e-8 nq=$NP < sudata |          ", 
" suaddhead ns=$NT | sushw key=dt,ntr a=$DT,$NP > temp                  ",
" sushw key=offset infile=myfile < temp  > sudata1rad 			",
"          								",
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

segy tr,trr;
void hrrtf_(double *pos,float *data,int *nh, int *nt,double *dt,int *method,
double *eps1,int *iter_end, double *qmin, float *model,float *q,double *dq,
int *nq);

int main(int argc, char **argv)
{
	
	FILE *myfilep;
	int ii,i;
	float *data, *model, *q; 
        double *pos, dtf, eps1, qmin, dq;
	int ntf, nhf, nq, method, iter_end;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile");

	pos=ealloc1double(NNX);	
	q=ealloc1float(NNX);
	data=ealloc1float(NNX*NT);
	model=ealloc1float(NNX*NT);
	/* Get parameters */


	if (!getparint("method", &method))  method = 1;
	if (!getpardouble("eps1", &eps1))  eps1 = 3;
	if (!getparint("iter_end", &iter_end))  iter_end = 5;
	if (!getpardouble("qmin", &qmin))  qmin = 0;
	if (!getparint("nq", &nq))  nq = 200;
						
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

				
		pos[ii]=(double) tr.offset;
		for (i=0;i<nt;i++){
		        data[i+ii*nt]=(float) tr.data[i];
		}
		ii++;
       		/*puttr(&tr);*/
	} while (gettr(&tr));
	nhf=ii;	
        hrrtf_(pos,data,&nhf,&ntf,&dtf,&method,&eps1,&iter_end,&qmin,
	       model,q,&dq,&nq);
	
	fwrite(model,sizeof(float),ntf*nq,stdout);
	
       	for (ii=0;ii<nq;ii++) q[ii]=1e12*q[ii];
	 
      	fwrite(q,sizeof(float),nq,myfilep);

	fclose(myfilep);



       
	return EXIT_SUCCESS;
}
