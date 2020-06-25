#include "su.h"

void eomig_old(float *m, float *t, float *h, float *q,  float *d, float *vel,
	   int adj, int nt, int nh, int nq) 
{
  /* 
     Hyperbolic operator  of a csp to Radon or viceversa

     if (adj) data to model
     else model to data


     Input 
        csp[ih][it]  or m[iq*nt+it] 
	he[ih]  equivalent offset axis
	t[it]   time axis
	vel[it] Velocity axis
     Output 
	m[iq*nt+it]  or   csp[ih][it]

     Daniel Trad - UBC- January 2000. 
	
  */
  
  int it;
  int ih,iq;
  int j;
  float ftime;
  float dint;
  float time,hxh,moveout,veltot,slowness2;
  int iqxnt,ihxnt,itint;
  int itime;
  float a1;
  float a2;
  int nx=nt*nq;
  int ny=nt*nh;
  float dt=t[1]-t[0];

  if (adj) for (it=0;it<nx;it++) m[it]=0.;
  else for (it=0;it<ny;it++) d[it]=0.;
  
  for (it=0;it<nt;it++){
    for (iq=0;iq<nq;iq++){
      veltot=q[iq]+vel[it];
      slowness2=1./(veltot*veltot);
      iqxnt=iq*nt;
      for (ih=0;ih<nh;ih++){
	hxh=h[ih]*h[ih];
	ihxnt=ih*nt;
	moveout=hxh*slowness2;
	time=2.*sqrt(t[it]*t[it]/4.+moveout);
	ftime=time/dt;      
	if (adj){
	  if (ftime <nt){
	    ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,1,&ftime,&dint);
	    m[iqxnt+it]+=dint;
	  }
        }
	else{	
	  itime=(int) floor(ftime);
          a2=ftime-itime;
	  a1=1.-a2;	  
	  if (itime<nt)  d[ihxnt+itime]+=a1*m[iqxnt+it];
	  if ((itime+1) < nt) d[ihxnt+itime+1]+=a2*m[iqxnt+it];
	}
      }
    }      
  }  
  return;
}

