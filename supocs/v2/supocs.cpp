/* Copyright (c) University of British Columbia, 2002.*/
/* All rights reserved.                       */
/* suradonfk0  :  $Date: April    2002- Last version July    2002  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>
void save_gather(float **d, int nh, int nt, float dt, const char* name);
bool process(float** datain, float** dataout, int nt, int  nh, float dt);
bool process2(float** datain, float** dataout, int nt, int  nh, float dt);
/*********************** self documentation **********************/
const char *sdoc[] = {
  " 	   								",
  " SUPOCS  template for su programs                                ",
  "                                                                     ",
 NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 * Last changes: July : 2012 
 */
/**************** end self doc ***********************************/

segy tr, tr2;
FILE* fp1;      // file pointer to first file
FILE* fp2;      // file pointer to second file
FILE* fpout1;      // file pointer to first  file
FILE* fpout2;      // file pointer to second file


int main(int argc, char **argv)
{
  int verbose;
  segy tr; 
  cwp_String modelfile=""; /* output sufile for the model */ 	
  //FILE *modelfilep;

  time_t start,finish;
  double elapsed_time;
  int it,ih;
  float t0=0;

  float **datain    = 0;
  float **datain2   = 0;
  float **dataout   = 0;
  float **dataout2  = 0;

  float **Wd        = 0;
  float *t, *h;

  int nt, nh; 
  int method;
  int plot; 
  float dt,dh;
  //float kmin;
  float fmax;
  //char buf[80];
  ////////////////
    
  fprintf(stderr,"*******SUPOCS*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  for (int ia=0; ia< argc; ia++) fprintf(stderr,"arg[%d]=%s\n",ia,argv[ia]);
  

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =1;
  if (!getparint("plot",&plot)) plot = 0;

  // if (!gettr(&tr)) err("can't read first trace"); // if reading from stdinput.
  fp1 = efopen(argv[1], "r");
  fp2 = efopen(argv[2], "r");
  fpout1 = efopen(argv[3],"w");
  fpout2 = efopen(argv[4],"w");

  
  if(!fgettr(fp1,&tr)) err("can't read first trace of first file");

  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");
  

  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh = (int) tr.ntr;

  if (!getparfloat("fmax",&fmax)) fmax = 0.5/dt;
  fmax = MIN(fmax,0.5/dt);

  // Allocate memory for data and model

  datain   = ealloc2float(nt,nh);
  datain2  = ealloc2float(nt,nh);
  dataout  = ealloc2float(nt,nh);
  dataout2 = ealloc2float(nt,nh);


  Wd=ealloc2float(nt,nh);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);

  // Loop over traces from first file.
  // Dead traces are marked with trid = 2. They are weighted to zero
  ih=0;
  do {
    if (tr.trid==2) for(it=0;it<nt;it++) Wd[ih][it]=0;
    else  for(it=0;it<nt;it++) Wd[ih][it]=1;
    h[ih]=(float) tr.offset;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (fgettr(fp1,&tr));
  //erewind(stdin);
  nh=ih;

  
  // Loop over traces from second file
  if (!fgettr(fp2,&tr)) err("can't read first trace of second file");
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) datain2[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh);
  } while (fgettr(fp2,&tr));

  if (verbose) fprintf(stderr,"processing %d traces \n", nh);
  for (it=0;it<nt;it++) t[it]=t0+it*dt;  /* Not implemented for t0 != 0  */
  dh=(h[nh-1]-h[0])/(nh-1);
  
  if (plot){ /* additional plots only for degugging  */  
    save_gather(datain,nh,nt,dt,"datain.su");
    system("suximage < datain.su perc=98 key=offset title=datain xbox=%d curve=curve1 npair=5 &");  
  }
  
  //process(datain,datain2,dataout,dataout2,nt,nh,dt);
  //for (ih=0;ih<nh;ih++)  for (it=0;it<nt;it++) dataout[ih][it]=datain[ih][it];
  
  rewind(fp1);
  for (ih=0;ih<nh;ih++){ 
    fgettr(fp1,&tr);
    memcpy((void *) tr.data,(const void *) datain[ih],nt*sizeof(float));
    tr.offset=(int) h[ih];
    tr.ntr=nh;
    if (Wd[ih][0]==0) tr.trid=2; // dead trace
    fputtr(fpout1,&tr);
  }
  
  rewind(fp2);
  for (ih=0;ih<nh;ih++){ 
    fgettr(fp2,&tr);
    memcpy((void *) tr.data,(const void *) datain2[ih],nt*sizeof(float));
    tr.offset=(int) h[ih];
    tr.ntr=nh;
    if (Wd[ih][0]==0) tr.trid=2; // dead trace
    fputtr(fpout2,&tr);
  }




  /******** End of interpolated output **********/


  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);
  
  free2float(datain);
  free2float(datain2);
  free2float(dataout);
  free2float(dataout2);
  
  efclose(fp1);
  efclose(fp2);
  return 0;
  efclose(fpout1);
  efclose(fpout2);

  return EXIT_SUCCESS;
}










