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
"	   Program in development                                 	",
" 	   								",
" suhrrt < stdin > stdout offsetfile=  [no optional parameters]		",
" 									",
" Required parameters:							",
" offsetfile=	ascii file of new offset values. If there are no        ",
" 		changes use                                             ",
"               sugethw output=geom key=offset < sudata>offsetfile	",
" stdin: output from suhrrtf or any other sudata with teh radon model.  ",  
"        This sufile must have the radon parameter value in key=offset  ",
"                                                                       ",
" stdout: is a binary file with the bare traces. To put the same        ",
" original header use sustrip to get the header and supaste to put it   ",
" back.                                                              	",
"               							",
"									",
" Example : # Inverse Radon Transform                                   ",
"      suhrrti offsetfile=sudata.off  < sudatarad > temp                ",
"      supaste ns=$NT head=headers < temp > sudata                      ",
"									",
NULL};

/* Credits:
 *
 *      Mauricio Sacchi
 *	Daniel Trad
 *      Tad Ulrych
 *
 * Trace header fields accessed: ns, dt, key=offset
 *
 */
/**************** end self doc ***********************************/

segy tr;

void hrrti_(double *pos,float *data,int *nhf,int *ntf,double *dtf,
float *model,double *qmin,double *qmax,int *nq);

int main(int argc, char **argv)
{
	FILE *offsetfilep;
	FILE *myfilep;
	int ii,i, nn;
	float *data, *model, temp; 
        double *pos, dtf, qmin,qmax, dq;
	int nt,nh,ntf, nq;
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
	 

        ii=0;
        do{
	     nn=fscanf(offsetfilep,"%f",&temp); 
             pos[ii]=temp;
             ii++;
	}while(nn==1);

        nh=ii-1;
	nh=64;
        fprintf(myfilep,"%d\n",nh);

	/*
	for (ii=0;ii<nh;ii++){
	    fprintf(myfilep,"%f\n",pos[ii]);
        }
	*/

	fclose(offsetfilep);

	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	
	ii=0;
	/* Loop over traces */
	do {
		int nt     = (int) tr.ns;
		float dt   = ((double) tr.dt)/1000000.0;
   		register int i;
		double offset= (int) tr.offset;

		if (ii==0){  
		    qmin=offset;
                    ntf=nt;
		    dtf=dt;
		}
		qmax=offset;
					
		for (i=0;i<nt;i++){
		        model[i+ii*nt]=(float) tr.data[i];
		}
		ii++;
       		/*puttr(&tr);*/
 
   
	} while (gettr(&tr));

        nq=ii;

        qmin=qmin*1e-12;
	qmax=qmax*1e-12;
        fprintf(myfilep,"ntf %d, dtf %e, qmin %e, qmax %e\n ",ntf,dtf,qmin,qmax);

        hrrti_(pos,data,&nh,&ntf,&dtf,model,&qmin,&qmax,&nq);
 	/*for (ii=0;ii<64;ii++)*/
	fwrite(data,sizeof(float),ntf*nh,stdout);
	/*       	fwrite(q,sizeof(float),nq,myfilep);*/
        fprintf(myfilep,"ntf %d, nh %d \n ",ntf,nh);

	fclose(myfilep);



       
	return EXIT_SUCCESS;
}
