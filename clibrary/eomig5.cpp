#include <math.h>
#include "su.h"

void eomig1(float **csp, float *he,float *m, float *t, float *vel,float t0) 
{
  /* 
     NMO-Kirchoff stack of a csp (includes obliquity factor)

     Input 
        csp[ih][it]  Common Scatter Point 
	he[ih]  equivalent offset axis
	t[it]   time axis
	vel[it] Velocity axis
     Output 
	m[it]   Migrated trace        
     Daniel Trad - UBC- January 2000. 	
  */

  int i,k,j,ih,ix;
  register int it;
  float itime,*dint,*dtemp,*obliq,time,sx,gx,t02;
  float temp, moveout;
  extern int nt,nh,nx,nhcsp;
  extern float dt,dh,dx;
  int it0=(int) (t0/dt+0.5);
  
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");

  for (it=0;it<nt;it++){
    m[it]=0;
    dint[it]=0;
    dtemp[it]=0;
  } 
  for (ih=0;ih<nhcsp;ih++){
    // DSR equation
    for (it=it0;it<nt;it++){
      moveout=he[ih]*he[ih]/(vel[it]*vel[it]);
      dtemp[it]=csp[ih][it];
      t02=(t[it]*t[it]/4.0);
      time=2*sqrt(fabs(t02+moveout));
      itime=time/dt;
      //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
      if (time>1e-2){
	temp=t[it]/time;//temp=pow(temp,1.5);
        temp=sqrt(temp*temp*temp);
      	obliq[it]=temp;
      }
      else obliq[it]=0;
      if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
      
      if (itime<nt) ints8r(nt,1.0,0.,dtemp,0.0,0.0,1,&itime,&dint[it]);
      else dint[it]=0;

    }
    //    ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    for (it=it0;it<nt;it++)
      m[it]+=(dint[it]*obliq[it]);                  
  }
 
  free1float(dtemp);
  free1float(dint);
  free1float(obliq);
  return;
}






void eomig2(float **csp, float *he,float *m, float *t, float *vel,float t0, int nt, int nhcsp) 
{
  /* 
     NMO-Kirchoff stack of a csp (includes obliquity factor)

     Input 
        csp[ih][it]  Common Scatter Point 
	he[ih]  equivalent offset axis
	t[it]   time axis
	vel[it] Velocity axis
     Output 
	m[it]   Migrated trace        
     Daniel Trad - UBC- January 2000. 	
  */

  int i,k,j,ih,ix;
  register int it;
  float itime,*dint,*dtemp,*obliq,time,sx,gx,t02;
  float temp, moveout;

  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");

  for (it=0;it<nt;it++){
    m[it]=0;
    dint[it]=0;
    dtemp[it]=0;
  } 
  for (ih=0;ih<nhcsp;ih++){
    // DSR equation
    for (it=it0;it<nt;it++){
      moveout=he[ih]*he[ih]/(vel[it]*vel[it]);
      dtemp[it]=csp[ih][it];
      t02=(t[it]*t[it]/4.0);
      time=2*sqrt(fabs(t02+moveout));
      itime=time/dt;
      //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
      if (time>1e-2){
	temp=t[it]/time;//temp=pow(temp,1.5);
        temp=sqrt(temp*temp*temp);
      	obliq[it]=temp;
      }
      else obliq[it]=0;
      if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
      
      if (itime<nt) ints8r(nt,1.0,0.,dtemp,0.0,0.0,1,&itime,&dint[it]);
      else dint[it]=0;

    }
    //    ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    for (it=it0;it<nt;it++)
      m[it]+=(dint[it]*obliq[it]);                  
  }
 
  free1float(dtemp);
  free1float(dint);
  free1float(obliq);
  return;
}













