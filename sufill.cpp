/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suinterpfk4 :  $Date: November 2001- Last version January 2000  */
#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUFILL       fill gaps                                              ", 
  " 	   								",
  " sufill  < stdin (from disk)  > stdout [optional parameters]         ",
  " 									",
  " Given a data set with gaps or poor sampling this program adds       ",
  " new traces into the gaps or between traces.                         ",
  " Required parameters:		[None]		       		",
  "                                                                	",
  " Standard input : data file  (offset time domain)          		",
  " Standard output : file filled with zero traces                      ",
  " Optional parameters:		       				",
  " dhnew=          if not given is calculated from the data            ",
  " option=2        2 fill with zero traces                             ",
  "                 1 read ascii file with offset                       ", 
  " offsetfile=     ascii file with new offset (1 columm)               ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt, offset
 * Last changes: January 18: 2002 
 */
/**************** end self doc ***********************************/

void maxmin(float *x,int lx, float *pmin, float *pmax);
int findgaps(float *h,float *h2,int nh,float dh0, float hmin, float hmax);
int read_ascii_file(const char *name,float *x);

int main(int argc, char **argv)
{
  int verbose;
  segy tr, trz; 
  //FILE *offsetfile; 
  time_t start,finish;
  double elapsed_time;
  int it,ih=0, ih2=0;
  float *t, *h, *h2;
  int nt, nh, nh2; 
  int plot; 
  int plot2;
  int option;
  //  float dh0;
  float maxdh; 
  float mindh;
  float dhnew;


  cwp_String offsetfile=NULL; /*input ascii file for offset if interpolation is desired */
  //////////////////////////////////////////////
  fprintf(stderr,"******* SUFILLOFFSETS *********\n");
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  // Get parameters 
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("plot",&plot)) plot = 0;
  if (!getparint("plot2",&plot2)) plot2 = 0;
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  if (!getparint("option",&option)) option=2;

  if (offsetfile==0) option=2;

  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  nt = (int) tr.ns;
  nh= (int) tr.ntr;
  nh2=10*nh;

  /* for mute in the migrated space we need to set the geometry of the mask */
  if (verbose){
    fprintf(stderr,"Option=%d is selected\n",option);  
    if (option==1)  
      fprintf(stderr,"Read offset file \n");
    else if (option==2)  
      fprintf(stderr,"Use findgaps to create the new offset \n");
  }

  // Allocate memory for offset
  h=ealloc1float(nh);
  t=ealloc1float(nt);
  h2=ealloc1float(nh2); // allocate more to play safe
  memset( (void *) h, (int) '\0', nh * FSIZE);
  memset( (void *) h2, (int) '\0', nh2 * FSIZE);

  /* If offetfile name is given read it */
  if ((option==1)&&(offsetfile)) nh2=read_ascii_file(offsetfile,h2);
  else{
    fprintf(stderr,"Using findgaps\n");
    option=2;
  }
  // Loop over traces 

  if (option==2){
    ih=0;
    do{   // Loop over traces    
      h[ih]=tr.offset;
      ih++;
    }while(gettr(&tr));
    nh=ih;

    maxmin(h,nh,&mindh,&maxdh);
    fprintf(stderr,"min=%f, max=%f, \n",mindh,maxdh);
    if (!getparfloat("dhnew",&dhnew)) dhnew=mindh;
    nh2=findgaps(h,h2,nh,dhnew,h[0],h[nh-1]);

    if (verbose) for (ih=0;ih<nh2;ih++)  fprintf(stderr,"%f\n",h2[ih]);
  }
  fprintf(stderr,"nh=%d,nh2=%d\n",nh,nh2);
  if (verbose) fprintf(stderr,"Original traces are:\n");
  if (option==2){
    erewind(stdin);
    gettr(&tr);
  }

  ih=0;ih2=0;
  do {

    h[ih]=(float) tr.offset;
    if (verbose)  fprintf(stderr,"h[%d]=%f---->h2[%d]=%f\n",ih,h[ih],ih2,h2[ih2]);

    if (h[ih]==h2[ih2]){ 
      tr.ntr=nh2;
      fputtr(stdout,&tr);
      gettr(&tr);
      ih++;ih2++;
    }
    else if (h[ih] > h2[ih2]){
      fprintf(stderr,"h[%d]=%f-- h2[%d]=%f<----- zero trace \n",ih,h[ih],ih2,h2[ih2]);
      //trz=tr; // change later to copy the trace properly
      memcpy((void *) &trz,(const void *) &tr,HDRBYTES);
      for (it=0; it<nt; ++it) trz.data[it] = 0 ;
      trz.trid=2;
      trz.offset=(int) h2[ih2];
      trz.ntr=nh2;
      fputtr(stdout,&trz);
      ih2++;
    }
    // There is a problem here, needs to be fixed.
    // temporal solution (trial and error)
    else{ // if h[ih] < h2[ih2] output h[ih]
      fprintf(stderr,"h2 > h1 ??? \n"); 
      fprintf(stderr,"h[%d]=%f-- h2[%d]=%f<----- zero trace \n",ih,h[ih],ih2,h2[ih2]);
      //trz=tr; // change later to copy the trace properly
      memcpy((void *) &trz,(const void *) &tr,HDRBYTES);
      for (it=0; it<nt; ++it) trz.data[it] = 0 ;
      trz.trid=2;
      trz.offset=(int) h2[ih2];
      trz.ntr=nh2;
      fputtr(stdout,&trz);
      ih2++;



      fprintf(stderr,"ih=%d,nh=%d\n",ih,nh);
      break;
      //gettr(&tr);
      //ih++;
      
    }
  } while ((ih<nh) || (ih2 < nh2));


  nh=ih;


  free1float(h2);
  free1float(t);
  free1float(h);
  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}


int read_ascii_file(const char *name,float *x)
{
  int nn;
  int ix;
  int nx;
  FILE *fp;
  fp=efopen(name,"r");
  
  ix=0;
  do{
    nn=fscanf(fp,"%f",&x[ix]);
    ix++;
    //fprintf(stderr,"ix=%d\n",ix);
  }while(nn==1);
  nx=ix-1;
  
  efclose(fp);
  
  return(nx);
}




int findgaps(float *h,float *h2,int nh,float dh0, float hmin, float hmax)
{
  /* count nh2 */
  //  int nh2l;
  //  int nh2r;
  int ih;
  int ig;
  int ih2;
  int ngaps;
  float dh;
  int *gaps;
  int nh2=5*nh; // 5 times nh maximum for nh2 but can be increased as needed
  int nhtemp;
  // This routine is a mesh, but it works
  // Perhaps it should be implemented as a linked list instead.
  // The purpose is, given some aproximate regular offset with gaps,  
  // 1- Find gaps
  // 2- Create a new offset axis with gaps infilled
  gaps=ealloc1int(10*nh);

  // find the location of the gaps and save them in gaps
  ig=0;
  for (ih=1;ih<nh;ih++){
    dh=h[ih]-h[ih-1];
    if (dh<0) err("Input data must be sort by increasing offset  \n");
    if (dh>(1.9*dh0)){
      gaps[ig]=ih-1;
      fprintf(stderr,"gaps[%d]=%d,ih=%d,dh0=%f\n",ig,gaps[ig],ih,dh0);
      ig++;
    }
  }
  ngaps=ig;
  gaps[ngaps]=nh-1;

  // find the number of traces that should be in every gap (nhtemp)
  // and create the new offset locations by assuming interval dh0
   
  //if (ngaps==0) for (ih=1;ih<nh;ih++) { h2[ih]=h[ih];return(nh);}
  for (ih=0;ih<=gaps[0];ih++) h2[ih]=h[ih];
  ih2=ih;
  for (ig=0;ig<ngaps;ig++){
    nhtemp=(int) (0.5+(h[gaps[ig]+1]-h[gaps[ig]])/dh0);
    fprintf(stderr,"nhtemp=%d, %f, %f, ih2=%d\n",nhtemp,h[gaps[ig]+1],h[gaps[ig]],ih2);
    for (ih=1;ih<nhtemp;ih++){
      h2[ih2]=h[gaps[ig]]+ih*dh0;
      //fprintf(stderr,"h2[%d]=%f,h[%d]=%f\n",ih2,h2[ih2],gaps[ig],h[gaps[ig]]);
      ih2++;
    }
    
    for (ih=1;h2[ih2-1]<h[gaps[ig+1]];ih++){
      if ((gaps[ig]+ih)<nh){
	h2[ih2]=h[gaps[ig]+ih];
	//fprintf(stderr,"--->h2[%d]=%f,h[%d]=%f\n",ih2,h2[ih2],gaps[ig]+ih,h[gaps[ig]+ih]);
	ih2++;
      }      
    }
  }
  nh2=ih2;

  //fprintf(stderr,"dh0=%f,nh2=%d,ngaps=%d,h2[0]=%f,h2[%d]=%f\n",dh0,nh2,ngaps,h2[0],nh2-1,h2[nh2-1]);

  free1int(gaps);

  return(nh2);  
}


void maxmin(float *x,int lx, float *pmin, float *pmax)
{
/********************************************************************  
return the  min and max of increments
*********************************************************************/

  /* 
     Given a vector x of length lx, with irregular increments 
     find the max and min of the increments.
  */
   int i;
   float dx, max, min;		
   max=0;
   min=1e10;	
   for (i=1;i<lx;i++){
     dx=x[i]-x[i-1];
     if (dx>max) max=dx;
     if (dx<min) min=dx;
   }
   *pmax=max;
   *pmin=min;
   return;
}










