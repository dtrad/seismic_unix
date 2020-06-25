/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTEMPLATE1:  $Date: June 2000  */


#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " sumixgathers file1 file2 > stdout [optional parameters]   		",
  "        								",
  "        								",
  " 	This template is for reading a model file generated from some   ",
  "     and from the model generate the original data                   ",
  "     The original data file is read to put back the original offset  ",
  "     but no computation is done with the tr.data                     ",
  "        								",
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
  int nq;
  int flag;
  float dt;
  float *h;
  float scale;
  int plot;
  float **d;
  float **m;

  
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;
  
  /* Open two files given as arguments for reading */
  fp1 = efopen(argv[1], "r");
  fp2 = efopen(argv[2], "r");  
  tracefp = etmpfile();
  headerfp = etmpfile();


  // Read first file and save it to a temporal file 
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
  nh=j;
  erewind(tracefp);
  erewind(headerfp);
  fprintf(stderr,"nh=%d,nt=%d\n",nh,nt);

  d=ealloc2float(nt,nh);


  // Read second file and save the data to a 2D array m
 
  if (!fgettr(fp2,&tr)) err("can't read first trace");
  if (!(nq=(unsigned long int) tr.ntr)) err("***ntr must be set\n");

  m=ealloc2float(nt,nq);
  
  ih=0;
  do{   // Loop over traces    
    memcpy((void *) m[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
  }while(fgettr(fp2, &tr));
  nq=ih;
  fprintf(stderr,"nh2=%d\n",nq);

  
  // Here we do something with m[iq][it] and put the result to d[ih][it]
  // Just for testing we make data1 = data2 
  
  for (ih=0;ih<nh;ih++){
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data, (const void *) m[ih], nt*sizeof(float));
    puttr(&tr);
  }

  
  free2float(d);
  free2float(m);

  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);
  free1float(h);

  return EXIT_SUCCESS;
}




















