/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUTAPERINV:  $Date: October 2002  */


#include "su.h"
#include "segy.h"
#include "header.h"

float dot(int n, float *a, float *b);
void save_gather(float **d, int nh, int nt, float dt, char* name);

/*********************** self documentation **********************/
char *sdoc[] = {
  " sutaperinv file1 file2 > output   					",
  "                                      				",
  " given file1 and file2, where file2 is a processed version of file1  ",
  " that has suffered some tapering, this program calculates a scale    ",
  " factor in a least squares sense to eliminate the effect of the taper",
  " in file 2.                                                          ",
  "                                                                     ",
  " Optional parameter                                                  ",  
  " ntaper=5   the number of traces that had suffered tapering at each  ",
  " side of the gather                                                  ",
  " 	   								",
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
  int nh2;
  int flag;
  float dt;
  float *h;
  float scale=1;
  int ntaper;
  int plot;
  float dpd; /* dot product dprec * data */
  float dpdp;/* dot product dprec * dpred */

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;

  //if (!getparfloat("scale",&scale)) scale=0;
  if (!getparint("ntaper",&ntaper)) ntaper=5;

  /* Open two files given as arguments for reading */
  fp1 = efopen(argv[1], "r");
  fp2 = efopen(argv[2], "r");  
  tracefp = etmpfile();
  headerfp = etmpfile();
  
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

  fprintf(stderr,"nh=%d,nt=%d\n",nh,nt);

  float **d=ealloc2float(nt,nh);
  float **d2=ealloc2float(nt,nh);


  if (!fgettr(fp2,&tr)) err("can't read first trace");
  
  ih=0;
  do{   // Loop over traces    
    memcpy((void *) d2[ih],(const void *) tr.data,nt*sizeof(float));
    //for (it=0;it<nt;it++) d2[ih][it]=tr.data[ih];
    ih++;
  }while(fgettr(fp2, &tr));
  nh2=ih;

  if (plot){
    save_gather(d2,nh,nt,0.004,"d2");
    system("suxwigb  < d2 title=\"d2\" perc=99 & ");
  }

  fprintf(stderr,"nh2=%d\n",nh2);

  if (nh!=nh2) err("nh is different from nh2\n");

  erewind(tracefp);
  erewind(headerfp);
  for (ih=0;ih<nh;ih++){
    efread(tr.data,FSIZE, nt, tracefp);   
    memcpy((void *) d[ih],(const void *) tr.data,nt*sizeof(float));
  }

  if (plot){
    save_gather(d,nh,nt,0.004,"d");
    system("suxwigb  < d title=\"d\" perc=99 & ");
  }




  erewind(tracefp);
  erewind(headerfp);

  for (ih=0;ih<nh;ih++){


    efread(&tr,HDRBYTES,1,headerfp);
    efread(tr.data,FSIZE, nt, tracefp);


    dpd=dot(nt,d[ih],d2[ih]);
    dpdp=dot(nt,d2[ih],d2[ih]);

    /* if scale=0 (default) then calculate the least square scale */
    if ((ih<ntaper)||(ih>=nh-ntaper)) 
      scale=dpd/dpdp;
    else
      scale=1;

    fprintf(stderr,"scale===>%f\n",scale);
    
    for (it=0;it<nt;it++) tr.data[it]=scale*d2[ih][it];
    puttr(&tr);
  }

  
  free2float(d);
  free2float(d2);

  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);
  free1float(h);

  return EXIT_SUCCESS;
}




















