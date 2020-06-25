#include <math.h>
#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"
#include "eomig.h"


void eomigls_sparse_large(float *m,float *mm,float *t, float *h, float *q, 
			  float *vel,float **data, int nt, int nh, int nq, 
			  float fmax, float step, float eps, float eps1, 
			  float eps2, int itercg, int iter_end, int smooth, 
			  int testadj, int buffer, float pervmin, float dperv,
			  int norm, unsigned short nmute, unsigned short mute,
			  float t0, float *Wm, float *Wd)
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
  float *d;             // temporal vector containing all data
  float *vel2;          // 1/vel^2
  int ny=nt*nh;         // length of data
  int nx=nt*nq;         // length of model
  int nsparse;          // Number of non zero elements in index matrix
  float **vgrid;        // Velocity grid
  unsigned int **index;        // Migration matrix operator 
  void (*oper) (float *,float *,unsigned int **, int ,int , int, int, int);
  float *J; // Cost function array
  oper=eomig;
  int restart=1;
  extern int FLAG, KEEP_WM;
  // Filter
  int nl=3;
  int nr=3;
  int flag=2;
  float *filter;

  //////////////
  fprintf(stderr,"Inside eomigls_sparse\n*************");

  if ((d=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  if ((Wd=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");
  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  
  if ((filter=alloc1float(nl+nr+1))==NULL)
    err("cannot allocate memory for filter \n");
  if ((J=alloc1float(iter_end+1))==NULL)
    err("cannot allocate memory for J \n");

  //for (it=0;it<nt;it++) m[it]=0.;
  //for (it=0;it<nx;it++) mm[it]=0.;
 
  //for (it=0;it<nt;it++) vel2[it]=1/(vel[it]*vel[it]);
  
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) d[ih*nt+it]=data[ih][it];
  for (it=0;it<ny;it++) Wd[it]=1;
  int qh=nq/2+1;

  // Irregular velocity grid
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);

  // count the number of elements we need for index
  // nsparse=count_nsparse(t,h,q,vgrid,nt,nh,nq,t0);
  nsparse=nt*nh*nq;

  // allocate the big monster
  size_t size=sizeof(unsigned int);
  if ((index=(unsigned int **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");
  // Assign the elements of index
  build_index(t,h,q,vgrid,nt,nh,nq,index,nsparse,t0);

  // Adjoint     
  eomig(mm,d,index,1,nt,nh,nq,nsparse);
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
    J[j]=wpcgnr(eomig,nt,nh,nq,nsparse,mm,d,Wd,Wm,index,eps,step,
		       itercg,restart);
    if (KEEP_WM) modelweight(mm,nx,norm,eps1,Wm);
    //if (0) smoothing(mm,nt,nq,nl,nr,flag); 
  }

  float qaux=0;
  
  if (mute){
    fprintf(stderr,"MUTING ********************************\n");
    for (iq=0;iq<nmute;iq++) for (it=0;it<nt;it++)  mm[iq*nt+it]=0.;
    for (iq=nq-nmute;iq<nq;iq++) for (it=0;it<nt;it++)  mm[iq*nt+it]=0.;
  }
  

  eomig(mm,d,index,0,nt,nh,nq,nsparse);

  if (smooth) smoothing(d,nt,nh,nl,nr,flag);

  //if (1) fprintf(stderr,"nt=%d,nh=%d\n",nt,nh);
  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
       data[ih][it]=d[ih*nt+it];
    
  //  eomig1(data,h,m,t,vel);
  free2((void **) index);
  free2float(vgrid);
  free1float(J);
  free1float(filter);
  free1float(vel2);
  free1float(d);

  return;
}

void eomgetvel_large(float *m,float *mm,float *t, float *h, float *q, 
		     float *vel,float **data, int nt, int nh, int nq, 
		     float fmax, float step, float eps, float eps1, 
		     float eps2, int itercg, int iter_end, int smooth, 
		     int testadj, int buffer, float pervmin,
		     float dperv, int norm, int itervel, float t0, 
		     float *Wm, float *Wd)
{
  /*
    function: eomgetvel.cpp

    Velocity correction for Least Squares Equivalent Offset Migration

    Migrates the hyperbolic data[h][t] (a csp gather) to 
    a single trace m[t] located at midpoint zero, using
    a velocity function vel[t]. 
    1 - The adjoint model m[iq*nt+it] is computed , 
    2 - An approximation to the Covariance matrix
    is estimated from the adjoint. 
    3 - WTCGLS solves the system
     (LH WdT Wd L + WmT Wm) m = LH WdT Wd d
    4 - A weighted sum of the model m[iq*nt+iq] is 
    used to find the velocity corrections 
    5 - New velocity = Input velocity + correction
 

    Input:
         1D arrays  m[nt] and mm[nt*nq] (working space) 
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
         vel2[it] corrected velocity function for a CMP location
    
    Daniel Trad - UBC- January 2000. 	

    E-mail: dtrad@geop.ubc.ca
  */  
 
  register int it;
  int  j, iter;
  int ih;
  int iq;
  float *d;             // temporal vector containing all data
  float *vel2;          // 1/vel^2
  int ny=nt*nh;         // length of data
  int nx=nt*nq;         // length of model
  int flag=2;
  int nsparse;          // Number of non zero elements in index matrix
  float **vgrid;        // Velocity grid
  unsigned int **index;        // Migration matrix operator 
  void (*oper) (float *,float *,unsigned int **, int ,int , int, int, int);
  float *J; // Cost function array
  oper=eomig;
  int restart=1;        // Every cg iteration restart the solution

  if ((d=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  

  for (it=0;it<nt;it++) m[it]=0.;
  
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) d[ih*nt+it]=data[ih][it];
  for (it=0;it<ny;it++) Wd[it]=1.;
  for (it=0;it<nx;it++) Wm[it]=1;

  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  // count the number of elements we need for index
  //nsparse=count_nsparse(t,h,q,vgrid,nt,nh,nq,t0);
  nsparse=nt*nq*nh;
  // allocate the big monster
  size_t size=sizeof(unsigned int);
  // allocate the big monster
  if ((index=(unsigned int **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");

  for (int jj=1;jj<=itervel;jj++){  
    // Irregular velocity grid
    irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
    // Assign the elements of index
    build_index(t,h,q,vgrid,nt,nh,nq,index,nsparse,t0);
    // Adjoint     
    eomig(mm,d,index,1,nt,nh,nq,nsparse);
    for (j=1;j<=iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(mm,nx,norm,eps1,Wm);
      wpcgnr(eomig,nt,nh,nq,nsparse,mm,d,Wd,Wm,index,eps,step,
		  itercg,restart);
    }
    get_vel(mm,t,q,vel,vel2,nt,nq,buffer);
    for (it=0;it<nt;it++) vel[it]=vel2[it];
  }

  free2((void **) index);
  free2float(vgrid);
  smoothing(vel,nt,1,10,10,flag);

  free1float(vel2);
  free1float(d);

  return;
}

















