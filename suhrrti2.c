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
" SUHRRTI -Inverse High Resolution Radon transform                      ", 
"	                                                          	",
" 	   								",
" suhrrt < stdin > stdout offsetfile=  [no optional parameters]		",
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
" rtmethod=1       1- PRT, 2LRT                                         ",
"									",
" Example : # Inverse Radon Transform                                   ",
"  suhrrti offsetfile=sudata.off  rtmethdo=1 < sudatarad > sudatarec    ",
"                                                                       ",
"                                                                       ",
"									",
NULL};

/* Credits:
 *
 *      Mauricio Sacchi
 *	Daniel Trad
 *      Tad Ulrych
 *
 * Trace header fields accessed: ns, dt, key=f2
 *
 */
/**************** end self doc ***********************************/

segy tr,*trr;

void hrrti_(double *pos,float *data,int *nhf,int *ntf,double *dtf,
float *model,double *qmin,double *qmax,int *nq, int *rtmethod);

int main(int argc, char **argv)
{
	FILE *offsetfilep;
	FILE *myfilep;
	int ii,i, nn;
	float *data, *model, temp; 
        double *pos, dtf, qmin,qmax, dq;
	int nt,nh,ntf, nq, maxtr, rtmethod;
	cwp_String offsetfile=""; /* file containing positions */

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile2","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile2");

       

	pos=ealloc1double(NNX);	
	data=ealloc1float(NNX*NT);
	model=ealloc1float(NNX*NT);

	/* Get parameters */
	
	if (!getparstring("offsetfile",&offsetfile)) 
	    err("must give pos=offset file");
       	fprintf(myfilep,"%s\n",offsetfile);
	if((offsetfilep=fopen(offsetfile,"r"))==NULL)
		err("cannot open offset file=%s\n",offsetfile);
	if (!getparint("rtmethod", &rtmethod))  rtmethod =1; /*PRT default*/ 

        ii=0;
        do{
	     nn=fscanf(offsetfilep,"%f",&temp); 
             pos[ii]=temp;
             ii++;
	}while(nn==1);

        nh=ii-1;
               fprintf(stderr,"nh= %d\n",nh);

	/*
	for (ii=0;ii<nh;ii++){
	    fprintf(myfilep,"%f\n",pos[ii]);
        }
	*/

	fclose(offsetfilep);

	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr) err("dt header field must be set");
	nq=(int) tr.ntr;
	maxtr= (nq>nh) ? nq : nh ;
	  fprintf(stderr,"maxtr=%d\n",maxtr);

	if ((trr=malloc(maxtr*sizeof(segy)))==NULL)
	  fprintf(stderr,"Sorry, space for traces could not be allocated\n");	
	
	ii=0;
       
	/* Loop over traces */
	do {
		int nt     = (int) tr.ns;
		float dt   = ((double) tr.dt)/1000000.0;
   		register int i;
		double f2= (float) tr.f2;
	        trr[ii]=tr;

		if (ii==0){  
		    qmin=f2;
                    ntf=nt;
		    dtf=dt;
		}
		qmax=f2;
					
		for (i=0;i<nt;i++){
		        model[i+ii*nt]=(float) tr.data[i];
		}
		ii++;
       		/*puttr(&tr);*/
   
	} while (gettr(&tr));

        nq=ii;
	fprintf(stderr,"nq=%d\n",nq);

        qmin=qmin;
	qmax=qmax;

        fprintf(myfilep,"ntf %d, dtf %e, qmin %e, qmax %e\n ",ntf,dtf,qmin,qmax);

        hrrti_(pos,data,&nh,&ntf,&dtf,model,&qmin,&qmax,&nq,&rtmethod);
 	
	ii=0;
	do{

	  if (ii<nq)  trr[ii]=trr[ii];  /*keep the header if nh=np*/
	  else trr[ii]=trr[0];    /* otherwise replace it with first one */   
	  trr[ii].offset=pos[ii];
       	  trr[ii].ntr=nh;
	  for (i=0;i<ntf;i++)
	  trr[ii].data[i]=data[i+ntf*ii];
	  puttr(&trr[ii]);

	  ii++;
	 } while(ii<nh);
	free(trr);
	/*for (ii=0;ii<64;ii++)*/
	/*fwrite(data,sizeof(float),ntf*nh,stdout);*/
	/*       	fwrite(q,sizeof(float),nq,myfilep);*/
        fprintf(myfilep,"ntf %d, nh %d \n ",ntf,nh);

	fclose(myfilep);



       
	return EXIT_SUCCESS;
}
