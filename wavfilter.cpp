/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

#include <math.h>
#include <stdio.h>
#include "wavfilter.h"
#include "su.h"
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " wavefilter wavelet filering for time series                         ", 
  "	   Program in development                                 	",
  " 	   								",
  " wavefilter   < stdin > stdout [optional parameters]          	",
  " 									",
  " 									",
  " Optional parameters:		       				",
  " method=0                                                      	",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Input : sudata file  (offset time domain)              		",
  " Output : adjoint model                                              ",
  "		                                        		",
  " Example: 		                                                ",
  " #Forward  Radon transform                                           ",
  " suradon1 pervmin=10 dperv=0.2 nq=30 itercg=5 < sudata  > sudata     ", 
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad.
 */
/**************** end self doc ***********************************/



int main(int argc, char* argv[])
{
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
	
  int maxscale;   //Level Wavelet 
  int nx;         //Lenght of signal
  float *x; // input
  float *y; // output

  if (!getparint("maxscale", &maxscale))  maxscale = 0;
  
  nx=100000;
  
  x=ealloc1float(nx);	
  nx=read_ascii_stdin(x);
  fprintf(stderr,"nx=%d",nx);
  realloc1float(x,nx);
  y=ealloc1float(nx);
  
  

  write_ascii_stdout(x,nx);
  free1float(y);
  free1float(x);
  return 0;
}

int read_ascii_stdin(float *x)
{
  int nn;
  int ix;
  int nx;

  ix=0;
  do{
    nn=fscanf(stdin,"%f",&x[ix]);
    ix++;
  }while(nn==1);
  nx=ix-1;
  
  return(nx);
}

void write_ascii_stdout(float *x, int nx)
{
  int ix;
  for(ix=0;ix<nx;ix++)
    fprintf(stdout,"%f\n",x[ix]);
  return;
}
