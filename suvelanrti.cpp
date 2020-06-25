/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
void velop(float *m, float *t, float *h, float *q, float *d);
/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SURADTDI -Inverse High Resolution Hyperbolic Radon transform          ", 
"	                                                          	",
" 	   								",
" suradtd < stdin > stdout offsetfile=  [no optional parameters]	",
" 									",
" Required parameters:							",
" offsetfile=	ascii file of new offset values. If there are no        ",
" 		changes use                                             ",
"               sugethw output=geom key=offset < sudata>offsetfile	",
" stdin: Radon model as it is produced by suhrrt2.                      ", 
"        It is a sufile with header.                                    ",  
"        This sufile must have the radon parameter value in             ",
"        key=f2                                                         ",
"                                                                       ",
" stdout: is a  sufile with the header and traces. The offset is        ",
"         given by offset file. Except offset, all other header words   ",
"         copied from original traces. As the field geometry do not     ",
"         exist for interpolated traces the header must corrected if    ",
"         interpolation is preformed. If original and final nh          ",
"         are the same, for example for multiple removal, the header    ",
"         is preserved.                                                 ",
" Optional parameters:              					",
"                                                                       ",
"									",
" Example : # Inverse Radon Transform                                   ",
"  suradtdi offsetfile=sudata.off   < sudatarad > sudatarec             ",
"                                                                       ",
"                                                                       ",
"									",
NULL};

/* Credits:
 *      Daniel Trad
 *
 * Trace header fields accessed: ns, dt, key=f2
 *
 */
/**************** end self doc ***********************************/

segy tr,*trr;
int nt,nh,nq,nx,ny,rtmethod,itercg;
float dt,dh,dq;


int main(int argc, char **argv)
{
	FILE *offsetfilep;
	FILE *myfilep;
	int i,j, nn, maxtr,flag;
	float *d, *m, temp; 
        float *q, *t, *h, *tempos, qmin, qmax, t0=0.;
	extern int nt,nh,nq,rtmethod;
	extern float dt,dh,dq;
	cwp_String offsetfile=""; /* file containing positions */

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile2","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile2");

       	

	/* Get parameters */
	if (!getparint("rtmethod", &rtmethod))  rtmethod =2; /*PRT default*/

	
	if (getparstring("offsetfile",&offsetfile)){
	  if ((tempos=ealloc1float(500))==NULL)
	   fprintf(stderr,"***Space for temppos could not be allocated\n");
          
          flag=0;  
	  fprintf(stderr,"New offset given by %s\n",offsetfile);  
	  if((offsetfilep=fopen(offsetfile,"r"))==NULL)
		err("cannot open offset file=%s\n",offsetfile);
	  j=0;
	  do{
	     nn=fscanf(offsetfilep,"%f",&temp); 
             tempos[j]=temp;
             j++;
	  }while(nn==1);
	  nh=j-1;
          fprintf(stderr,"nh= %d\n",nh);
	  fclose(offsetfilep);
	}
        else {
	  flag=1;
          if (!getparint("nh", &nh))  
	    err("I need either offsetfile or original nh");
        }

	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr){ 
             nq=NNX; 
             fprintf(stderr,"ntr header field is not set");
             fprintf(stderr,"nq set to %d\n",nq);
        } 
	else nq=(int) tr.ntr;
	
        nt     = (int) tr.ns;
	dt   = ((float) tr.dt)/1000000.0;
        if ((!nt)||(!nh)||(!nq)) err("Error reading nt,nh or nq");

	
	// Allocate memory for data and model
  
	if ((d=ealloc1float(nh*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
	
	//fprintf(stderr,"&d=%p,dd=%g\n",&d,d[15][20]);
	if ((m=ealloc1float(nq*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");
	
	if ((q=ealloc1float(nq))==NULL)
	  fprintf(stderr,"***Sorry, space for q could not be allocated\n");
  
	if ((h=ealloc1float(nh))==NULL)
	  fprintf(stderr,"***Sorry, space for h could not be allocated\n");

	if ((t=ealloc1float(nt))==NULL)
	  fprintf(stderr,"***Sorry, space for t could not be allocated\n");         
	// Because we want to use same struct array for data and model 
	// the maximun number of traces between them will be taken.

       	maxtr= (nq>nh) ? nq : nh; 
		fprintf(stderr,"maxtr=%d\n",maxtr); 
        if ((trr=(segy*) malloc(maxtr*sizeof(segy)))==NULL)
          fprintf(stderr,"**Sorry, space for traces could not be allocated\n");
	j=0;
       	/* Loop over traces */
	do {
       		register int i;
		q[j]= (float) tr.f2;
	        trr[j]=tr;
                if ((flag==1)&&(j<nh)) h[j]=(float) tr.offset; 					
		for (i=0;i<nt;i++){
		        m[i+nt*j]=(float) tr.data[i];
		}
		j++;
       		/*puttr(&tr);*/
                if (j > maxtr) err("Number of traces > %d\n",maxtr);     
	} while (gettr(&tr));

        nq=j;
	fprintf(stderr,"nq=%d\n",nq);
        fprintf(myfilep,"nt%d,dt %f\n ",nt,dt);

        if (flag==0) {
	  for (j=0;j<nh;j++) h[j]=tempos[j];
          free1float(tempos);
        }
	for (i=0;i<nt;i++) t[i]=t0+i*dt;
        nx=nq*nt;
        ny=nh*nt;
	//for (j=0;j<nt;j++)  fprintf(stderr,"t[%d]=%f\n",j,t[j]);
	velop(m,t,h,q,d);
         	
	j=0;
	do{

	  if (j>=nq) trr[j]=trr[0];  //replace it with first trace header 
	  trr[j].offset=(int) h[j];
       	  trr[j].ntr=nh;
	  for (i=0;i<nt;i++)
	    trr[j].data[i]=d[i+j*nt];
	  puttr(&trr[j]);
	  j++;
	 } while(j<nh);
	free(trr);
	/*for (ii=0;ii<64;ii++)*/
	/*fwrite(data,sizeof(float),ntf*nh,stdout);*/
	/*       	fwrite(q,sizeof(float),nq,myfilep);*/
        fprintf(stderr,"nt%d, nh %d \n ",nt,nh);
	fclose(myfilep); 
	free1float(t);
	free1float(h);
	free1float(q);
	free1float(m);
	free1float(d);            
	return EXIT_SUCCESS;
}











