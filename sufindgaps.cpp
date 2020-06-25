/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUFINDGAPS:  Module to search for gaps in a data file and 
   generate an ascii file with new offset 
   $Date: June 2000  
*/


#include "su.h"
#include "segy.h"
#include "header.h"
#include "dan.h"

int findgaps(float *h,float *h2,int nh,float dh0, float hmin, float hmax);


/*********************** self documentation **********************/
char *sdoc[] = {
  " sufindgaps < data > offsetfile (ASCII)   dhnew= (optional)          ",
  "        								",
  "        								",
  " 	This program reads a  file with gaps and computes a new offset  ",
  "  axis without gaps to use with any of the interpolation codes.      ",
  "  The new increment is given by the min of the offset increments.    ",
  "  It preserves the original locations, but if a gap is found it      ",
  "  fills it. The gap is found when the delta offset is greater than  	",
  "  1.9 times the minimum delta offset.                               	",
  "  A min delta offset can be given by the parameter dhnew.           	",
  "                                                                    	",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt                          */
/**************** end self doc ***********************************/

segy tr; 

int main(int argc, char **argv)
{
  int j,ih;
  int ntr;
  int nh;
  int nh2;
  float *h;
  float *h2;
  float dh0;
  float maxdh; 
  float mindh;
  float dhnew;
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!(ntr=(unsigned long int) tr.ntr)) err("***ntr must be set\n");
  nh = ntr;

  h=ealloc1float(nh); // allocate larger for now

  j=0;
  do{   // Loop over traces    
    h[j]=tr.offset;
    j++;
  }while(gettr(&tr));
  nh=j;

  h2=ealloc1float(5*nh);
  fprintf(stderr,"nh=%d\n",nh);

  //dh0=MIN(h[1]-h[0],h[nh-1]-h[nh-2]);
  maxmin(h,nh,&mindh,&maxdh);
  fprintf(stderr,"min=%f, max=%f, \n",mindh,maxdh);

  if (!getparfloat("dhnew",&dhnew)) dhnew=mindh;

  nh2=findgaps(h,h2,nh,dhnew,h[0],h[nh-1]);

  for (ih=0;ih<nh2;ih++)  fprintf(stdout,"%f\n",h2[ih]);

  free1float(h);
  fprintf(stderr,"%s Finished successfully\n",__FILE__);

  return EXIT_SUCCESS;
}



int findgaps(float *h,float *h2,int nh,float dh0, float hmin, float hmax)
{
  /* count nh2 */
  int nh2l;
  int nh2r;
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
  if (ngaps==0) for (ih=1;ih<nh;ih++) { h2[ih]=h[ih];return(nh);}
  for (ih=0;ih<=gaps[0];ih++) h2[ih]=h[ih];
  ih2=ih;
  for (ig=0;ig<ngaps;ig++){
    nhtemp=(int) ((h[gaps[ig]+1]-h[gaps[ig]])/dh0);
    fprintf(stderr,"nhtemp=%d, %f, %f, ih2=%d\n",nhtemp,h[gaps[ig]+1],h[gaps[ig]],ih2);
    for (ih=1;ih<=nhtemp;ih++){
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
























