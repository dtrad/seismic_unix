#include <math.h>
#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"
#include "eomig.h"


void eomigls_sparse_sinc(float *m,float *mm,float *t, float *h, float *q, 
			  float *vel,float **data, int nt, int nh, int nq, 
			  float fmax, float step, float eps, float eps1, 
			  float eps2, int itercg, int iter_end, int smooth, 
			  int testadj, int buffer, float pervmin, float dperv,
			  int norm, float parmute, unsigned short mute,
			  float t0, float *Wm, float *Wd, int centralq)
{
/*
    function: emigls.cpp
    Least Squares Equivalent Offset Migration

    Migrates the hyperbolic data[h][t] (a csp gather) to 
    a single trace m[t] located at midpoint zero, using
    a velocity function vel[t]. 
    1 - The adjoint model m[iq*nt+it] is computed , 
    2 - An approximation to the Covariance matrix
    is estimated from the adjoint. 
    3 - WTCGLS solves the system
     (LH WdT Wd L + WmT Wm) m = LH WdT Wd d
    4 - The model m[iq*nt+iq] is Kirchhoff-nmo stacked to obtain
        the migrated trace m[it]
 

    Input:
         data[ih][it]  Common Scatter Point
	 h[ih]  equivalent offset axis
	 t[it]   time axis
	 q[iq]   Radon axis
	 vel[it] velocity function for a the particular CMP
	 nt size of t
	 nh size of h
	 nq size of q
	 fmax maximum frequency fo interest
         Parameters for inversion: eps1, eps2, itercg, iter_end, eps
	 testadj (0 or 1)
	 smooth  (0 or 1)
	 buffer  (length/2 of weighted summation) 

    Output:
 	 mm[iq*nt+it] Radon model
	 m[it]        migrated trace
	 
    
    Daniel Trad - UBC- January 2000. 	

    E-mail: dtrad@geop.ubc.ca
*/ 
 
  register int it;
  int  j, iter;
  int ih;
  int iq;
  float *vel2;          // 1/vel^2
  int ny=nt*nh;         // length of data
  int nx=nt*nq;         // length of model
  int nsparse;          // Number of non zero elements in index matrix
  float **vgrid;        // Velocity grid
  float **index;        // Migration matrix operator 
  void (*oper) (float *,float *,float **, int ,int , int, int, int);
  float *J; // Cost function array
  oper=eomig;
  int restart=1;
  extern int FLAG, KEEP_WM;
  int nmute;
  // Filter
  int nl=3;
  int nr=3;
  int flag=2;
  float *filter;
  float test;
  //////////////
  fprintf(stderr,"Inside eomigls_sparse\n*************");

  if ((Wd=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");
  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  
  if ((filter=alloc1float(nl+nr+1))==NULL)
    err("cannot allocate memory for filter \n");
  if ((J=alloc1float(iter_end+1))==NULL)
    err("cannot allocate memory for J \n");

  for (it=0;it<ny;it++) Wd[it]=1;
  int qh=nq/2+1;

  // Irregular velocity grid
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  nsparse=nt*nh*nq;
  // allocate the big monster
  index=ealloc2float(nsparse+1,2);
  // Assign the elements of index
  build_index_inv(t,h,vgrid,nt,nh,nq,index,nsparse,t0);
  if (testadj) test=testadjop(eomig,index,nt,nh,nq,nsparse);
  // Adjoint     
  eomig(mm,data[0],index,1,nt,nh,nq,nsparse);
  if ( FLAG && KEEP_WM ){
    fprintf(stderr,"Using Wm for the first time from adjoint\n");
    modelweight(mm,nx,norm,eps1,Wm);
    FLAG=0;
  }
  
  for (j=1;j<=iter_end;j++){
    //if (0) smoothing(mm,nt,nq,nl,nr,flag);
    // norm==1 ==> L1 , ==0  Cauchy else l2
    if (!KEEP_WM) modelweight(mm,nx,norm,eps1,Wm);
    // if (0) smoothing(Wm,nt,nq,nl,nr,flag);
    J[j]=wpcgnr(eomig,nt,nh,nq,nsparse,mm,data[0],Wd,Wm,index,eps,step,
		       itercg,restart);
    if (KEEP_WM) modelweight(mm,nx,norm,eps1,Wm);
    //if (0) smoothing(mm,nt,nq,nl,nr,flag); 
  }

  float qaux=0;
  
  if (mute){
    iq=0; while(fabs(q[iq])<parmute) iq++; 
    nmute=nq-iq;
    fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
    taper(mm,nt,nq,nmute,2);
  }

  eomig(mm,data[0],index,0,nt,nh,nq,nsparse);

  if (smooth) smoothing(data[0],nt,nh,nl,nr,flag);

  //if (1) fprintf(stderr,"nt=%d,nh=%d\n",nt,nh);

  free2float(index);
  free2float(vgrid);
  free1float(J);
  free1float(filter);
  free1float(vel2);
  //free1float(d);

  return;
}


















