#include <math.h>
#include "su.h"
/* Kirchhoff Migration for prestack data */
void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **vel, int nt, int nx, float dt) 
{
  int i,k,j,ih,ix;
  register int it;
  float *itime,*dint,*obliq,time,sx,gx,t02;
  float s2,g2,s2v,g2v,xm, temp, temp2;
  float tsource;
  float trec;
  float cossour;
  float cosrec;
  float angs;
  float angr;    

  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");

  //fprintf(stderr,"cdpspace=%f\n",cdpspace);
  for (ix=0;ix<nx;ix++){
    xm=fabs(cdp-x[ix]);
    if (1){
      sx=xm-h;
      gx=xm+h;
      s2=sx*sx;
      g2=gx*gx;
      // DSR equation
      for (it=0;it<nt;it++){
	t02=(t[it]*t[it]/4);
	s2v=s2/(vel[ix][it]*vel[ix][it]);
	g2v=g2/(vel[ix][it]*vel[ix][it]);

	tsource=sqrt(t02+s2v);
	trec=sqrt(t02+g2v);
	time=tsource+trec;
	
	itime[it]=time/dt;
      
	/* The following lines are different obliquity factor
	   computations. I kept the the simpler one (t0/t)^(3/2)
	   
	   //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
	   //if (tsource>1e-3) cossour=t[it]/2/tsource;else cossour=0;
	   //if (trec>1e-3)	cosrec=t[it]/2/trec;else cosrec=0;
	   
      //if (cossour>1) cossour=0;
      //if (cosrec>1)  cosrec=0;
      
      //angs=acos(cossour);
      //angr=acos(cosrec);  
      //obliq[it]=MAX(fabs((cossour+cosrec)*cos(0.5*(angs-angr))),1);

      // obliq[it]=(cossour+cosrec)*0.5;
	*/
	if (0){ /* speed up */
	  temp=t[it]/time;
	  if (time>1e-2) obliq[it]=sqrt(temp*temp*temp);
	  else obliq[it]=0;
	  if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
	}
      }  
      ints8r(nt,1.0,0.,d,0.0,0.0,nt,itime,dint);
      
      for (it=0;it<nt;it++) m[ix][it]+=dint[it];
      //m[ix][it]+=(dint[it]*obliq[it]);
    }                  
  }
  free1float(itime);
  free1float(dint);
  free1float(obliq);
  return;
}

void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **ivel2, int nt, int nx, float dt, float cdpspace, float aper) 
{

  /* ivel2 is  1/vel^2 */

  int ix;
  register int it;
  float *itime,*dint,*obliq,time,sx,gx,t02;
  float s2,g2,s2v,g2v,xm, temp;
  float tsource;
  float trec;
  //float cossour;
  //float cosrec;
  //float angs;
  //float angr;    
  
   
  itime=ealloc1float(nt);
  dint=ealloc1float(nt);
  obliq=ealloc1float(nt);
  


  //fprintf(stderr,"cdpspace=%f\n",cdpspace);
  for (ix=0;ix<nx;ix++){
    xm=fabs(cdp-x[ix]);
    if (xm < aper){
      xm*=cdpspace;
      sx=xm-h;
      gx=xm+h;
      s2=sx*sx;
      g2=gx*gx;
      // DSR equation
      for (it=0;it<nt;it++){
	t02=(t[it]*t[it])*0.25;

	s2v=s2*ivel2[ix][it];
	g2v=g2*ivel2[ix][it];

	tsource=sqrt(t02+s2v);
	trec=sqrt(t02+g2v);
	time=tsource+trec;
	
	itime[it]=time/dt;
      
	/* The following lines are different obliquity factor
	   computations. I kept the the simpler one (t0/t)^(3/2)
	   
	   //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
	   //if (tsource>1e-3) cossour=t[it]/2/tsource;else cossour=0;
	   //if (trec>1e-3)	cosrec=t[it]/2/trec;else cosrec=0;
	   
      //if (cossour>1) cossour=0;
      //if (cosrec>1)  cosrec=0;
      
      //angs=acos(cossour);
      //angr=acos(cosrec);  
      //obliq[it]=MAX(fabs((cossour+cosrec)*cos(0.5*(angs-angr))),1);

      // obliq[it]=(cossour+cosrec)*0.5;
	*/
	if (0){ /* speed up */
	  temp=t[it]/time;
	  if (time>1e-2) obliq[it]=sqrt(temp*temp*temp);
	  else obliq[it]=0;
	}
	if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
      }  
      ints8r(nt,1.0,0.,d,0.0,0.0,nt,itime,dint);
      
      for (it=0;it<nt;it++) m[ix][it]+=dint[it];
      //m[ix][it]+=(dint[it]*obliq[it]);
    }                  
  }
  free1float(itime);
  free1float(dint);
  free1float(obliq);
  return;
}























