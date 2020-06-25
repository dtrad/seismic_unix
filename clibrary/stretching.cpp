#include "su.h"
#include "interpfk.h"

void stretch(float **data, float **datas, int nt, int nu, int nh, float *t, float *tu, float *ut,
	     float dt, float du, int str)
{
  int ih, it, iu;

  if (str==1){/* stretch */
    if (nu!=nt && tu!=NULL) {
      float *temp=ealloc1float(nu);
      for (ih=0; ih<nh; ++ih) {
	ints8r(nt,dt,0.0,data[ih],0.0,0.0,nu,tu,temp);
	for (iu=0; iu<nu; ++iu) datas[ih][iu] = temp[iu];
      }
      free1float(temp);
    }
  }
  else{   /* unstretch */
    if (nu!=nt && ut!=NULL) {
      float *temp=ealloc1float(nt);
      for (ih=0; ih<nh; ++ih) {
	ints8r(nu,du,0.0,datas[ih],0.0,0.0,nt,ut,temp);
	for (it=0; it<nt; ++it) data[ih][it] = temp[it];
      }
      free1float(temp);
    }
  }
  
  return;
}


void makev (int nmig, float *tmig, float *vmig, float vscale,
	    int nt, float dt, float ft, float **v, float *vmin, float *vmax)
/*****************************************************************************
make uniformly sampled rms velocity function v(t) for migration
******************************************************************************
Input:
nmig		number of tmig,vmig pairs
tmig		array[nmig] of times
vmig		array[nmig] of rms velocities
vscale		velocity scale factor
nt		number of time samples
dt		time sampling interval
ft		first time sample

Output:
v		array[nt] of rms velocities
vmin		minimum velocity
vmax		maximum velocity
******************************************************************************
Author:	 Dave Hale, Colorado School of Mines, 10/22/91
*****************************************************************************/
{
	int it;
	float t,*vel,velmin,velmax,(*vmigd)[4];
	
	vmigd = (float(*)[4])ealloc1float(nmig*4);
	cmonot(nmig,tmig,vmig,vmigd);
	vel = ealloc1float(nt);
	for (it=0,t=ft; it<nt; ++it,t+=dt)
		intcub(0,nmig,tmig,vmigd,1,&t,&vel[it]);
	for (it=0; it<nt; ++it)
		vel[it] *= vscale;
	for (it=1,velmin=velmax=vel[0]; it<nt; ++it) {
		velmin = MIN(velmin,vel[it]);
		velmax = MAX(velmax,vel[it]);
	}
	free1float((float*)vmigd);
	*v = vel;
	*vmin = velmin;
	*vmax = velmax;
}

void makeut (float vstolt, float fmax, float *vrms,
	     int nt, float dt, float **ut, int *nu, float *du, float **tu)
/*****************************************************************************
Compute u(t) and t(u) for Stolt stretch
******************************************************************************
Input:
vstolt		Stolt migration velocity
fmax		maximum frequency
vrms		array[nt] of RMS velocities
nt		number of t samples
dt		t sampling interval (first t assumed to be zero)

Output
ut		array[nt] of u(t); NULL if constant velocity
nu		number of u samples
du		u sampling interval (first u assumed to be zero)
tu		array[nu] of t(u); NULL if constant velocity
*****************************************************************************/
{
	int it;
	
	/* check for constant velocity */
	for (it=1; it<nt; ++it)
		if (vrms[it]!=vrms[0]) break;
		
	/* if constant velocity */
	if (it==nt) {
		*ut = NULL;
		*tu = NULL;
		*nu = nt;
		*du = dt;

	/* else if velocity not constant */
	} else {
		int it,nuu;
		float duu,delu,umax,*u,*t;
		
		/* u(t) */
		u = alloc1float(nt);
		makeu(vstolt,vrms,nt,dt,u);

		/* smallest du and maximum u */
		duu = FLT_MAX;
		for (it=1; it<nt; ++it) {
			delu =	u[it]-u[it-1];
			if (delu<duu) duu = delu;
		}
		umax = u[nt-1];

		/* u sampling */
		duu = duu/(2.0*fmax*dt);
		nuu = 1+NINT(umax/duu);

		/* t(u) */
		t = alloc1float(nuu);
		yxtoxy(nt,dt,0.0,u,nuu,duu,0.0,0.0,(nt-1)*dt,t);

		/* set output parameters before returning */
		*ut = u;
		*tu = t;
		*nu = nuu;
		*du = duu;
	}
}

void makeu (float vstolt, float *v, int nt, float dt, float *u)
/*****************************************************************************
Compute				     t
	u(t) = sqrt( 2/vstolt^2 * Integral ds*s*(v(s)^2) )
				     0
via the trapezoidal rule.
******************************************************************************
Input:
vstolt		Stolt migration velocity
v		array[nt] of RMS velocities
nt		number of t samples
dt		t sampling interval

Output
u		array[nt] of u(t)
*****************************************************************************/
{
	int it;
	float t,scale,sum;

	scale = 2.0/(vstolt*vstolt);
	u[0] = sum = 0.0;
	for (it=1,t=dt; it<nt; ++it,t+=dt) {
		sum += 0.5*dt*(t*v[it]*v[it]+(t-dt)*v[it-1]*v[it-1]);
		u[it] = sqrt(scale*sum);
	}
}







