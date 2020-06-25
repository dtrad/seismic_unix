/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suinterpfk4 :  $Date: November 2001- Last version January 2000  */

#include "interpfk.h"
#include "segy.h"
#include "header.h"
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUFILL       Interpolation with LS square time migration            ", 
  "	         by using Stolt fk operator                       	",
  " 	   								",
  " sufill  < stdin > stdout [optional parameters]               	",
  " 									",
  " Given a data set with gaps or poor sampling this program adds       ",
  " new traces into the gaps or between traces.                         ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Standard input : data file  (offset time domain)          		",
  " Standard output : Interpolated data file                            ",
  " Optional parameters:		       				",
  " option=1        1 sinc interpolation with (interpfact-1)            ",
  "                   non-zero traces                                   ",
  "                 2 interpolate first with (interpfact-1)             ", 
  "                   zero traces and then fill them with LS migration  ",
  "                 3 interpolate gaps first with zero traces           ",
  "                   as specified by offsetfile                        ",
  "                   and then fill them with LS migration              ",
  "                 other.. do not interpolate                          ",
  " offsetfile=     ascii file with new offset (1 columm)               ",
  " interpfact=2    Factor for increasing sampling      		",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 * Last changes: January 18: 2002 
 */
/**************** end self doc ***********************************/



int main(int argc, char **argv)
{
  int verbose;
  segy tr; 
  FILE *offsetfile; 
  time_t start,finish;
  double elapsed_time;
  int it,ih;
  float t0=0;
  float **datain=0;
  float **dataout=0;
  float *t, *h, *h2;
  int nt, nh, nh2; 
  int plot; 
  int plot2;
  float dt,dh2;
  int option;

  float interpfact;
  float dh=0;

  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  //////////////////////////////////////////////
  fprintf(stderr,"******* SUFILLOFFSETS *********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("verbose", &verbose))  verbose =1;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparint("plot2",&plot2)) plot2 = 0;
  if (!getparint("testadj",&testadj)) testadj=0; 
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  if (!getparfloat("interpfact",&interpfact)) interpfact=2;
  if (!getparint("option",&option)) option=2;

  if (offsetfile==0) err("Requires to input a new offset");

  // Velocity law
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  nh2=10*nh;

  ft = tr.delrt/1000.0;
  if (ft!=0.0) err("cannot handle non-zero time of first sample");

  /* for mute in the migrated space we need to set the geometry of the mask */
  if (verbose){
    fprintf(stderr,"New time axis after stretching: nt=%d ,dt=%f \n",nu,du);  
    fprintf(stderr,"Option=%d is selected\n",option);  
    if (option==1) 
      fprintf(stderr,"Means sinc interpolation before migration \n");
    else if (option==2)  
      fprintf(stderr,"Means zero traces interpolation before migration \n");
    else if (option==3)  
      fprintf(stderr,"Means fill in gaps according to offsetfile \n");
    if ((option==1)||(option==2)) 
      fprintf(stderr,"%f traces will be inserted between original traces\n",(interpfact-1));
  }

  // Allocate memory for data
  datain=ealloc2float(nt,nh);
  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);
  /* If offetfile name is given read it */
  h2=ealloc1float(nh2); // allocate more to play safe
  memset( (void *) h2, (int) '\0', nh2 * FSIZE);


  if (option==3) if (offsetfile) nh2=read_ascii_file(offsetfile,h2); 
  else nh2=nh;

  dataout0=ealloc2float(nt,nh);

  // Loop over traces 
  if (verbose) fprintf(stderr,"Original traces are:\n");
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
    if (verbose) fprintf(stderr,"ih=%d\n",ih);   
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;

  /* Time axis */
  for (it=0;it<nt;it++) t[it]=t0+it*dt;  /* Not implemented for t0 != 0  */

  /***************************************************************** 
     From this point the program has to create two new data spaces
     first: original sampling in offset to new offset sampling (interpolation) 
     second: original time sampling to new time sampling (stretching)
  *******************************************************************/


  /* Create the interpolated data space */
  if ((option==1)||(option==2)){
    dataout1=ealloc2float(nt,200);
    Wd=ealloc2float(nt,nh*(int) interpfact);
  }
  else if (option==3){
    dataout1=ealloc2float(nt,nh2);
    Wd=ealloc2float(nt,nh2);
  }
  else{
    dataout1=ealloc2float(nt,nh);
    Wd=ealloc2float(nt,nh);
  } 

  /***** Initial offset interpolation *****/

  if (plot){ /* additional plots only for degugging  */  
    save_gather(datain,nh,h,nt,dt,"datain.su");
    system("suxwigb < datain.su perc=100 key=offset title=datain &");  
  }
  if (option==1)   /* initial interpolation by zero padding in f-k */
    nh2=fft2_zeropad(datain,dataout1,Wd,nt,nh,interpfact);
  else if (option==2)   /* initial interpolation by zero padding in t-x */
    nh2=add_zerotraces(datain,dataout1,Wd,nt,nh,interpfact);
  else if (option==3)   /* initial interpolation by filling gaps by zero padding in t-x */
    nh2=add_zerotraces(datain,dataout1,Wd,nt,nh,nh2,h,h2);
  else for (ih=0;ih<nh2;ih++) h2[ih]=h[ih]; /* for wrong option  */

  /************************************************/
  /* Create new offset axis */  
  if ((option==1)||(option==2)){ 
    dh=(h[nh-1]-h[0])/(nh-1);
    dh2=dh/interpfact;
    for (ih=0;ih<nh2;ih++) h2[ih]=h[0]+ih*dh2;
  }
  else if (option==3){
    dh=(h[nh-1]-h[0])/(nh-1);
    dh2=(h2[nh2-1]-h2[0])/(nh2-1);
  }
  else dh2=(h2[nh2-1]-h2[0])/(nh2-1);
  
  /* Replace original offset axis with the new one */
  h=realloc1float(h,nh2);
  for (ih=0;ih<nh2;ih++) h[ih]=h2[ih];
  nh=nh2;

  if (verbose) fprintf(stderr,"New offset axis has dh=%f, nh=%d\n",dh,nh);
  if (plot){
    save_gather(dataout1,nh2,h2,nt,dt,"dataout1.su");
    system("suxwigb < dataout1.su perc=100 key=offset title=dataout2 &");    
  }

  /************************************************/
  
  /* Replace  original data array with the new one*/
  free2float(datain);
  datain=ealloc2float(nt,nh2);
  for (ih=0;ih<nh2;ih++) for (it=0;it<nt;it++) datain[ih][it]=dataout1[ih][it]; 
  /* Temporal array dataout1 is no longer needed */
  free2float(dataout1);
 
  dataout=ealloc2float(nt,nh2);    

  erewind(stdin);
  for (ih=0;ih<nh2;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    for (it=0;it<nt;it++) tr.data[it]*=ascale;
    tr.offset=(int) h2[ih];
    tr.ntr=nh2;
    fputtr(stdout,&tr);
  }

  /******** End of interpolated output **********/

  free1float(h2);
  free1float(t);
  free1float(h);
  free2float(datain);
  free2float(dataout);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}


















