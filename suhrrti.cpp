/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "su.h"
#include "segy.h"
#include "header.h"
#include "Complex.h"
#include "clibrary.h"
void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax);
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

segy tr;
int nt, nh, nq, rtmethod;
float depth;
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	FILE *offsetfilep;
	FILE *myfilep;
	int i,ih,iq, nn,flag;
	float **data, **model, temp, fmax; 
        float *pos, dt, *q, *tempos;
	extern int nt,nh,nq,rtmethod;
	extern float depth;

	cwp_String offsetfile=""; /* file containing positions */

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        if((myfilep=fopen("myfile2","w"))==NULL)
                        err("cannot open myfile=%s\n","myfile2");

       	

	/* Get parameters */
	if (!getparint("rtmethod", &rtmethod))  rtmethod =2; /*PRT default*/
	if (!getparfloat("fmax", &fmax))  fmax=70; /* default*/
	if (!getparfloat("depth", &depth))  depth=500; 
	
	if (getparstring("offsetfile",&offsetfile)){
	  if ((tempos=ealloc1float(500))==NULL)
	   fprintf(stderr,"***Space for temppos could not be allocated\n");
          
          flag=0;  
	  fprintf(stderr,"New offset given by %s\n",offsetfile);  
	  if((offsetfilep=fopen(offsetfile,"r"))==NULL)
		err("cannot open offset file=%s\n",offsetfile);
	  ih=0;
	  do{
	     nn=fscanf(offsetfilep,"%f",&temp); 
             tempos[ih]=temp;
             ih++;
	  }while(nn==1);
	  nh=ih-1;
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

	if ((data=ealloc2float(nh,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
 
	if ((model=ealloc2float(nq,nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");

	if ((q=ealloc1float(nq))==NULL)
	  fprintf(stderr,"***Sorry, space for q could not be allocated\n");
               
	if ((pos=ealloc1float(nh))==NULL)
	  fprintf(stderr,"***Sorry, space for pos could not be allocated\n");
        
	// Because we want to use same struct array for data and model 
	// the maximun number of traces between them will be taken.

	headerfp = etmpfile();
	if (verbose) warn("using tmpfile() call");  
	
	iq=0;
       	/* Loop over traces */
	do {
       		register int i;
		efwrite(&tr,HDRBYTES,1,headerfp);    
		q[iq]= (float) tr.f2;
                if ((flag==1)&&(iq<nh)) pos[iq]=(float) tr.offset;    
		for (i=0;i<nt;i++){
		  model[i][iq]=(float) tr.data[i];
		}
		iq++;   
	} while (gettr(&tr));
	erewind(headerfp); 
        nq=iq;
	if (verbose) fprintf(stderr,"nq=%d\n",nq);
        fprintf(myfilep,"nt%d,dt %f\n ",nt,dt);

        if (flag==0) {
	  for (ih=0;ih<nh;ih++) pos[ih]=tempos[ih];
          free1float(tempos);
        }

        hrrti(data,pos,dt,model,q,fmax);
 	
	ih=0;
	do{
	  if (ih<nq)  efread(&tr,HDRBYTES,1,headerfp);
	  tr.offset=(int) pos[ih];
       	  tr.ntr=nh;
	  for (i=0;i<nt;i++)
	    tr.data[i]=data[i][ih];
	  puttr(&tr);
	  ih++;
	 } while(ih<nh);
       
        if (verbose) fprintf(stderr,"nt%d, nh %d \n ",nt,nh);
	fclose(myfilep);
 	free1float(pos);
	free1float(q);
	free2float(model);
	free2float(data);
	efclose(headerfp);
              
	return EXIT_SUCCESS;
}











