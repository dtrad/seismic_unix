/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUHRRTFPI:  $Date: March 1999  */
void hrrtfpi(float *pos, float **d ,float dt,float eps1,float eps2,float qmin,
	   float **m,float *q,float dq,float freq, float eps);

#define NNX 228
#define NT 2048
#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
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
 *	Daniel Trad
 *       Mauricio Sacchi
 *      Tad Ulrych
 *
 * Trace header fields accessed: ns, dt, delrt, key=keyword
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/

segy tr,*trr; // for malloc use *trr
complex czero=(0.0);
int nt, nh, nq, method, iter_end, rtmethod, norm, itercg, freqflag, costflag,
  nfpm;
complex **pm;
int main(int argc, char **argv)
{
	
	FILE *myfilep, *radonfilep;
	int j,i, maxtr;
	float **d, **m;
        float *q, eps1,eps2, eps, qmin, dq; 
        float *pos, dt, freq;
	extern int nt, nh, nq, method, iter_end, rtmethod, norm, 
	               itercg, freqflag, costflag, nfpm;
        extern complex **pm; 
	complex *pmtemp;
	cwp_String radonfile=""; /* file containing prior model */	
        fprintf(stderr,"**************************\n");
	// Initialize 
	initargs(argc, argv);
	requestdoc(1);

        if((myfilep=fopen("radonpar","w"))==NULL)
                        err("cannot open file=%s\n","radonpar");

	// Get parameters 
	if (!getparint("costflag", &costflag))  costflag = 1; //=1 goes to Nyq
	if (!getparint("freqflag", &freqflag))  freqflag = 0; //=1 goes to Nyq
	if (!getparint("method", &method))  method = 1;
	if (!getparfloat("eps1", &eps1))  eps1 = 3;
	if (!getparfloat("eps2", &eps2))  eps2 = 1-7;
	if (!getparfloat("eps", &eps))  eps = 1e-7;
	if (!getparint("iter_end", &iter_end))  iter_end = 5;
	if (!getparfloat("qmin", &qmin))  qmin = 0;
	if (!getparint("nq", &nq))  nq = 220;
	if (!getparint("itercg", &itercg))  itercg = nq/4;
	if (!getparfloat("freq", &freq))  freq =0; 
	                                       // freq=0 use Fnyquist form data
	if (!getparint("rtmethod", &rtmethod))  rtmethod =1; // PRT default
	if (!getparint("norm", &norm))  norm =10; // PRT default
	if (!getparint("nfpm", &nfpm))  nfpm =100; // NF in prior model	    
					
	// Get info from first trace 

	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ns) err("ns header field must be set");
	if (!tr.ntr) err("ntr header field must be set");
	dt   = ((float) tr.dt)/1000000.0;
	nt = (int) tr.ns;
        nh= (int) tr.ntr;
        fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt); 

	// Allocate memory for data and model

	if ((d=ealloc2float(nh,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
 	
        //fprintf(stderr,"&d=%p,dd=%g\n",&d,d[15][20]);
	if ((m=ealloc2float(nq,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");

	if ((q=ealloc1float(nq))==NULL)
	  fprintf(stderr,"***Sorry, space for q could not be allocated\n");

	if ((pos=ealloc1float(nh))==NULL)
	  fprintf(stderr,"***Sorry, space for pos could not be allocated\n");
      
	// Because we want to use same struct array for data and model 
	// the maximun number of traces between them will be taken.
	///////////////////////////////////////////////////////////////
	if (getparstring("radonfile",&radonfile)){
	  if ((pm=ealloc2complex(nq,nfpm))==NULL)
	    fprintf(stderr,"***Space for pm could not be allocated\n");
 	  if ((pmtemp=ealloc1complex(nfpm))==NULL)
	    fprintf(stderr,"***Space for pmtemp could not be allocated\n");          
          
	  fprintf(stderr,"Prior information given by %s\n",radonfile); 
	  if ((norm==11)&&((method==2)||(method==7))) { 
	  if((radonfilep=fopen(radonfile,"r"))==NULL)
	    err("cannot open offset file=%s\n",radonfile);
	    for (i=0;i<nq;i++) {
	      fprintf(stderr,"i=%d\n",i);  
	      fread(pmtemp,sizeof(complex),nfpm,radonfilep);             
	      for (j=0;j<nfpm;j++) pm[j][i]=pmtemp[j];
	    }           
	    fclose(radonfilep);
	  }
	}
	/////////////////////////////////////////
       	maxtr= (nq>nh) ? nq : nh; 
		fprintf(stderr,"maxtr=%d\n",maxtr); 

        if ((trr=(segy*) malloc(maxtr*sizeof(segy)))==NULL)
          fprintf(stderr,"**Sorry, space for traces could not be allocated\n");

   //  if ((trr=(segy*) ealloc1(maxtr,sizeof(segy)))==NULL)
   //  fprintf(stderr,"**Sorry, space for traces could not be allocated\n");	
	for (j=0;j<nh-1;++j) 
	         pos[j]=0;

	j=0;
	// Loop over traces 
	do {
   		register int i;

		trr[j]=tr;
      		pos[j]=(float) tr.offset;
                //fprintf(stderr,"pos[j]=%f\n",pos[j]);
		for (i=0;i<nt;i++){
		        d[i][j]=(float) tr.data[i];
			//	fprintf(stderr,"i=%d\n",i);
		}
      		j++;
	} while (gettr(&tr));

       nh=j;	
       hrrtfpi(pos,d,dt,eps1,eps2,qmin,m,q,dq,freq,eps);
       fprintf(stderr,"After hrrtf nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
       if ((norm==1)&&((method==2)||(method==7))) {
	 if((radonfilep=fopen(radonfile,"w+"))==NULL)
	   err("cannot open offset file=%s\n",radonfile);
	 for (i=0;i<nq;i++) {
	   fprintf(stderr,"i=%d\n",i);
 	   for (j=0;j<nfpm;j++) pmtemp[j]=pm[j][i]; 
	   fwrite(pmtemp,sizeof(complex),nfpm,radonfilep);             
	 }           
	 fclose(radonfilep);
       }
       /*fwrite(model,sizeof(float),ntf*nq,stdout);*/

	j=0;
	do{

	  if (j>=nh)  trr[j]=trr[0];  // otherwise replace it with first one
      	  trr[j].f2=q[j];  // copy radon parameter in tr.f2
       	  trr[j].ntr=nq;
 
	  for (i=0;i<nt;i++)
	    trr[j].data[i]=m[i][j];

	  puttr(&trr[j]);
	  j++;
	 } while(j<nq);
 
	// Binary file with radon parameter for other non su programs

       	fwrite(q,sizeof(float),nq,myfilep);
        fprintf(stderr," nq=%d, nt=%d, nh=%d, j=%d\n",nq,nt,nh,j);
	fclose(myfilep);

	free(trr);
	free1float(pos);
	//             free1float(q);
	//    	free2float(m);
	free2float(d);
        free2complex(pm);
        free1complex(pmtemp);
       
	return EXIT_SUCCESS;
}




















