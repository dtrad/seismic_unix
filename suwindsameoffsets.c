/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUWINDSAMEOFFSETS:  $Date: June 2000  */

#include "su.h"
#include "segy.h"

#include "header.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUWINDSAMEOFFSETS - input two gathers and it let pass the traces in ",
  " first gather that are also in the second one                        ",
  " It is useful if gathers with different sampling have to be compared ",
  "                                                                     ",
  " suwindsameoffsets file1 file2 > stdout [optional parameters]        ",
  " It takes no parameteres                                             ",
  "                                                                     ",
  " IMPORTANT: Both files have to be sorted by offset                   ",
  " Mixes two gathers keeping only the traces of the first file         ",
  " if the offset is the same. The purpose is to substitute only        ",
  " traces non existing in file1 by traces interpolated store in file2. ", 
  " Example. If file1 is original data file and file 2 is obtained by   ",
  " resampling with Radon transform, then the output contains original  ",
  " traces with gaps filled                                             ",
  "                                                                     ", 
  " 									",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr,tr2; 
FILE *fp1;		/* file pointer for first file		*/
FILE *fp2;		/* file pointer for second file		*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

int main(int argc, char **argv)
{
  int j,ih,k;
  register int it;
  int nt;
  int ntr;
  int nh;
  int flag;
  float dt;
  float *h;
  float scale=1;
  int scaling;
  float closeunder=0.999;
  float closeupper=1.001;

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);

   

  /* Open two files given as arguments for reading */
  fp1 = efopen(argv[1], "r");
  fp2 = efopen(argv[2], "r");  
  tracefp = etmpfile();
  headerfp = etmpfile();
  if (!getparint("scaling",&scaling)) scaling=0;  
  if (!fgettr(fp1,&tr)) err("can't read first trace");

  nt = (int) tr.ns;  
  if (!(ntr=(unsigned long int) tr.ntr)) err("***ntr must be set\n");
  h=ealloc1float(ntr);

  j=0;
  do{   // Loop over traces    
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);  	   
    h[j]=tr.offset;
    j++;
  }while(fgettr(fp1, &tr));

  erewind(tracefp);
  erewind(headerfp);
  nh=j;
  fprintf(stderr,"nh=%d\n",nh);
  scale=1;
  /* compare one by one the headers of the two files
     if they are equal or almost equal set flag=0
     which means that trace of file 2 is not going to be 
     in the output
  */
 
  while(fgettr(fp2, &tr2)){
    /* in principle every trace should be output */
    /* in principle every trace should be output */
    flag=1;
    for (ih=0;ih<nh;ih++){
      if (h[ih]>0) 
	if ((tr2.offset>closeunder*h[ih])&&(tr2.offset<closeupper*h[ih])) flag=0;
      if (h[ih]<0) 
	if ((tr2.offset<closeunder*h[ih])&&(tr2.offset>closeupper*h[ih])) flag=0;
      if ((h[ih]==0)&&(tr2.offset==0)) flag=0; /* special case */
    }
    if (flag==0){
      if (scaling && fabs(tr2.offset) > 0){
	scale=(1+0.03*(fabs(tr2.offset)/1000.));
	fprintf(stderr,"tr2.offset=%d,scale=%f\n",tr2.offset,scale);
	for (it=0;it<nt;it++) tr2.data[it]*=scale;
      }
      puttr(&tr2);
    }

  }
  

  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);
  free1float(h);

  return EXIT_SUCCESS;
}




















