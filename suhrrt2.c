/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 512
#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUHRRT2 -Forward  High Resolution Parabolic Radon transform           ", 
"	   Program in development                                 	",
" 	   								",
" suhrrt2 < stdin > stdout [optional parameters]          		",
" 									",
" 									",
" Optional parameters:							",
" method=1                                                      	",
" 		1 Normal RT (Gauss- Gauss)			        ",
"		2 Cholesky        (Cauchy Gaus method)                  ",
"               3 Stepest descent (Cauchy Gauss method)                 ",
"               4 Conjugate gradient (Cauchy Gauss method)              ",
"                                                                       ",
" eps1 =3		noise level  for Cholesky                       ",
"                                                                       ",
" eps=1e-7              small number for conjugate gradient             ",
"                               					",
" iter =5		number of iterations               	        ",
" qmin =0               minimum Radon parameter in sec/offset^2         ",
" 		                           				",
" nq=150		Number of Radon Traces                         	",
" rtmethod              1-PRT 2-LRT                                     ",
" freq                  defines the freq max value used to compute      ",
"                       qmax according to the formula                   ",
"                  PRT:    qmax=1.5* {1/(2*fmax*(xmax-xmin)*dx)}        ",
"                  LRT:    qmax=1.5*(1/fmax*dx)                         ",
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
" key=f2 contains the radon parameter                                   ",
" If nh < np ( as usual) header words of first nh traces are preserved  ",
" and radon parameter is kept in f2 header word. Hence, when suhrrti2   ",
" is used to recover the data the header is preserved.                  ",
" eps1 is only for method 1 and 2, and eps for method 3 ad 4            ",
" qmax is computed from fmax, qmin, and nq	      		    	",
" and dq results from qmin, qmax and nq. Check in the messages that     ",
" the used dq is smaller than the critical. Otherwise increase np.      ",
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
void hrrtf_(double *pos,float *data,int *nh, int *nt,double *dt,int *method,
double *eps1,int *iter_end, double *qmin, float *model,float *q,double *dq,
int *nq, int *freq, int *rtmethod);

int main(int argc, char **argv)
{
	
  //FILE *myfilep;
	int ii,i;
	float *data, *model, *q; 
        double *pos;
        /*float data[466944],model[466944],q[228];
	  double pos[228];*/
        double  dtf, eps1, eps, qmin, dq;
	int ntf, nhf, nq, method, iter_end, maxtr, freq, rtmethod;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
	/* if((myfilep=fopen("myfile","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile");
	*/
	
	pos=alloc1double(NNX);	
	q=alloc1float(NNX);
	data=alloc1float(NNX*NT);
	model=alloc1float(NNX*NT);
	
	/* Get parameters */

	//fprintf(stderr,"**************\n");
	if (!getparint("method", &method))  method = 1;
	if (!getpardouble("eps1", &eps1))  eps1 = 3;
	if (!getpardouble("eps", &eps))  eps = 1e-7;
	if (!getparint("iter_end", &iter_end))  iter_end = 5;
	if (!getpardouble("qmin", &qmin))  qmin = 0;
	if (!getparint("nq", &nq))  nq = 220;
	if (!getparint("freq", &freq))  freq =0; 
	                                  /* freq=0 use Fnyquist form data*/
	if (!getparint("rtmethod", &rtmethod))  rtmethod =1; /*PRT default*/
	if (method==3) eps1=eps;
        if (method==4) eps1=eps;

						
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ns) err("ns header field must be set");
	ntf = (int) tr.ns;
	if (!tr.ntr) {
	            tr.ntr=NNX;
		    fprintf(stderr,"ntr is not set, I will take nh=228\n");
	}
	nhf=(int) tr.ntr;

       	maxtr= (nq>nhf) ? nq : nhf; 
		fprintf(stderr,"maxtr=%d\n",maxtr); 

	/* To decide nq we need dx, xmax and xmin, 
	   so we take the maximum possible for now 
           maxtr=NNX;
	*/

        if (maxtr>NNX) 
	      fprintf(stderr,"max number of x or p traces=%d\n",NNX);

	if (ntf>NT) 
	      fprintf(stderr,"max number of samples per traces=%d\n",NT);
 

	if ((trr=malloc(maxtr*sizeof(segy)))==NULL)
	  fprintf(stderr,"Sorry, space for traces could not be allocated\n");	      
         
	for (ii=0;ii<NNX;ii++)
	                  pos[ii]=0;
	ii=0;

	/* Loop over traces */
	do {

	        int nt     = (int) tr.ns;
		float dt   = ((double) tr.dt)/1000000.0;
   		register int i;
		trr[ii]=tr;
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
              model,q,&dq,&nq,&freq,&rtmethod);
	fprintf(stderr,"nq=%d, ntf=%d, nhf=%d\n",nq,ntf,nhf);

       	

	ii=0;
	do{

	  if (ii<nhf)  trr[ii]=trr[ii];  /*keep the header if nh=np*/
	  else trr[ii]=trr[0];    /* otherwise replace it with first one */     

	  trr[ii].f2=q[ii];
       	  trr[ii].ntr=nq;
	  for (i=0;i<ntf;i++)
	  trr[ii].data[i]=model[i+ntf*ii];
	  puttr(&trr[ii]);
	  ii++;
	 } while(ii<nq);

	 
      	//fwrite(q,sizeof(float),nq,myfilep);

	//fclose(myfilep);
        
	free(trr);
        free1float(data);
	free1float(model);
        free1double(pos);
        free1float(q);
	
	return EXIT_SUCCESS;
}













