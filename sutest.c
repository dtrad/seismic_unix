/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
/*#include "/usr/local/cwp/include/Complex.h"*/
/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUTEST -                                                              ", 
"	                                                          	",
" 	   								",

"          								",
"                                                                       ",
"									",
NULL};

/* 
 * Trace header fields accessed: ns, dt, delrt, key=keyword
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/

segy tr,*trr;
complex  x,y,z,**Z;
float **R,**RR;

/*void testcomplex(complex *x, complex *y, complex *z);*/

int main(int argc, char **argv)
{
	
	FILE *myfilep;
	int ii,jj,i;
	float *data, *model, *q; 
        double *pos, dtf, eps1, qmin, dq;
	int ntf, nhf, nq, method, iter_end, maxtr;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile");

	pos=ealloc1double(NNX);	
	q=ealloc1float(NNX);
	data=ealloc1float(NNX*NT);
	model=ealloc1float(NNX*NT);
	Z=alloc2complex(3,3);
	R=alloc2float(3,3);
	RR=alloc2float(3,3);
	/* Get parameters */


	if (!getparint("method", &method))  method = 1;
	if (!getpardouble("eps1", &eps1))  eps1 = 3;
	if (!getparint("iter_end", &iter_end))  iter_end = 5;
	if (!getpardouble("qmin", &qmin))  qmin = 0;
	if (!getparint("nq", &nq))  nq = 200;
						
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr) tr.ntr=64;
	nhf=(int) tr.ntr;
	maxtr= (nq>nhf) ? nq : nhf ;
	  fprintf(stderr,"maxtr=%d\n",maxtr);
	 

	if ((trr=malloc(maxtr*sizeof(segy)))==NULL)
	  fprintf(stderr,"Sorry, space for traces could not be allocated\n");	

	for (ii=0;ii<228;ii++)
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
		  /*data[i+ii*nt]=(float) tr.data[i];*/
		}
		ii++;

       		/*puttr(&tr);*/
	} while (gettr(&tr));
	nhf=ii;	
	x.r=2;x.i=1;
	y.r=3;y.i=2;
	for (jj=0;jj<3;jj++){
	  for (ii=0;ii<3;ii++){
	    Z[jj][ii].r=ii+1+jj*jj;
	    Z[jj][ii].i=ii*(jj+1);
	    R[jj][ii]=ii+1+jj*jj;
	    if (jj==ii) R[jj][ii]=R[jj][ii]+3;
	  }
	}

	for (jj=0;jj<3;jj++){  
	  for (ii=0;ii<3;ii++)
	    /* fprintf(stderr,"%f %f ",Z[jj][ii].r,Z[jj][ii].i);*/
	    fprintf(stderr,"%f ",R[jj][ii]);
	  fprintf(stderr,"\n");
	}
	/*	inverse_matrix(3,R);*/
	inverse_matrix_multiply(3,R,3,3,R,RR);
	for (jj=0;jj<3;jj++){  
	  for (ii=0;ii<3;ii++)
	    /* fprintf(stderr,"%f %f ",Z[jj][ii].r,Z[jj][ii].i);*/
	    fprintf(stderr,"%f ",RR[jj][ii]);
	  fprintf(stderr,"\n");
	}
 	/*sutestcomplex(&x,&y,&z);*/
	z=cadd(cexp(x),cexp(y));
	fprintf(stderr,"r=%f,i=%f",z.r,z.i);

	ii=0;
	do{
    	  puttr(&trr[ii]);
	  ii++;
	 } while(ii<nhf);

       	for (ii=0;ii<nq;ii++) q[ii]=1e12*q[ii];
	 
      	fwrite(q,sizeof(float),nq,myfilep);

	fclose(myfilep);
	free(trr);
       
	return EXIT_SUCCESS;
}






