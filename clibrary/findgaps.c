#include "su.h"

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
  // The purpose is given some aproximate regular offset with gaps 
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
  nh2=ih2-2;

  //fprintf(stderr,"dh0=%f,nh2=%d,ngaps=%d,h2[0]=%f,h2[%d]=%f\n",dh0,nh2,ngaps,h2[0],nh2-1,h2[nh2-1]);

  free1int(gaps);

  return(nh2);  
}


