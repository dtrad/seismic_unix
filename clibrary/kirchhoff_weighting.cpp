#include <math.h>
#include "su.h"

void kirchhoff_weighting(float *trace, float *t, float *vel,float t0, int nt) 
{
  /* 
     NMO-Kirchoff weighting for the stack of a csp (includes obliquity factor)

     Input 
        trace[it]  trace after weighting 
      
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
  
  for (it=it0;it<nt;it++){
    if (time>1e-2){
      temp=t[it]/time;//temp=pow(temp,1.5);
      temp=sqrt(temp*temp*temp);
      obliq[it]=temp;
    }
    else obliq[it]=0;
    if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
  }

  for (it=it0;it<nt;it++) trace[it]*=obliq[it];                  

  return;
}













