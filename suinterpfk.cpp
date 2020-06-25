/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADONHYPFK:  $Date: June 1999 - Last version October 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>
#include "dan.h"
#include "inversion_par.h"
//#include "radonhypfk.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONHYPFK Forward  Hyperbolic Radon transform  in the f-k domain ", 
  "	   Program in development                                 	",
  " 	   								",
  " suradonhypfk < stdin > stdout [optional parameters]          	",
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
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
//inv_par inv;




int verbose;
int interpfk(float *t, float *h, float *h2, float vel, float **datain, float **dataout,
	       int nt, int nh, int nh2, float dt, float t0, inv_par inv, int plot, 
	       float dx, int adj, int dft);

int main(int argc, char **argv)
{
  segy tr; 
  inv_par inv;
  cwp_String modelfile=""; /* output sufile for the model */ 	
  FILE *modelfilep;
  //  FILE *offsetfile; 
  time_t start,finish;
  double elapsed_time;
  int it,ih;
  float vel,dv;
  float t0=0;
  float **datain;
  float **dataout;
  float *t, *h, *h2;
  int nt, nh, nh2, rtmethod; 
  int method;
  int plot;
  float dt,dx;
  float factor; // multiply dq by factor
  int adj;
  int dft;
  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  //////////////////////////////////////////////
  fprintf(stderr,"*******SURADONHYPFK*********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("method", &method))  method = 0;
  if (!getparfloat("vel", &vel))  vel = 1500;
  if (!getparfloat("dv", &dv))  dv = 100;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparfloat("step", &inv.step))  inv.step =0.9;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("rtmethod", &rtmethod))  rtmethod =3; // HRT default
  if (!getparint("norm", &inv.norm))  inv.norm =1; 
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";    
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparfloat("dx",&dx)) dx=10; 
  if (!getparint("adj",&adj)) adj=1; 
  if (!getparint("dft",&dft)) dft=0; 
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;

  // Get info from first trace 

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  nh2=2*nh;
  // Allocate memory for data and model

  datain=ealloc2float(nt,nh);

  h=ealloc1float(nh);
  t=ealloc1float(nt);

  memset( (void *) h, (int) '\0', nh * FSIZE);
  /* If offetfile name is given read it */
  h2=ealloc1float(nh2); // allocate more to play safe
  memset( (void *) h2, (int) '\0', nh2 * FSIZE);
  if (offsetfile) nh2=read_ascii_file(offsetfile,h2);  

  dataout=ealloc2float(nt,nh2);

  // Loop over traces 
  ih=0;
  do {
    h[ih]=(float) tr.offset;
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  
  for (it=0;it<nt;it++) t[it]=0+it*dt;
  if (plot==1){
    save_gather(datain,nh,nt,0.004,"datain.su");
    system("suxwigb < datain.su perc=100 title=datain &");
  }	
  interpfk(t,h,h2,vel,datain,dataout,nt,nh,nh2,dt,t0,inv,plot,dx,1,dft);
 
  if (plot==1){
    save_gather(dataout,nh2,nt,0.004,"migrated.su");
    system("suxwigb < migrated.su perc=100 title=migrated &");
  }	

  interpfk(t,h2,h2,vel,dataout,dataout,nt,nh2,nh2,dt,t0,inv,plot,dx,0,dft);

  if (plot==1){
    save_gather(dataout,nh2,nt,0.004,"demigrated.su");
    system("suxwigb < demigrated.su perc=100 title=demigrated &");
  }	

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);

  rewind(stdin);
  for (ih=0;ih<nh2;ih++){ 
    fgettr(stdin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    tr.offset=h2[ih];
    fputtr(stdout,&tr);
  }

  free1float(h2);
  free1float(t);
  free1float(h);
  free2float(datain);
  free2float(dataout);

  efclose(modelfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}













