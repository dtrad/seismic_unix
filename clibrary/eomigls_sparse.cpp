#include <math.h>
#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"
#include "eomig.h"

// This file contains 
// eomigls_sparse: to perform (sort of) least squares migration
// eomgetvel: to iteratively apply velocity corrections
// There  2 versions of each one: one for unsigned short and other 
// for unsigned int. Except of this, they are equal.
// In the future will be substitute for templates.
// Then there is a function called 
// getvel: to compute the correction in velocities

void eomigls_sparse(float *m,float *mm,float *t, float *h, float *q, 
		    float *vel,float **data, int nt, int nh, int nq, 
		    float fmax, float step, float eps, float eps1, 
		    float eps2, int itercg, int iter_end, int smooth, 
		    int testadj, int buffer, float pervmin, float dperv,
		    int norm, float parmute, unsigned short mute, 
		    float t0, float *Wm, float *Wd,int centralq)
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
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  float *d;             // temporal vector containing all data

  //float *Wm;            // Factoriztion of inv model covariance matrix  
  float *vel2;          // 1/vel^2
  int ny=nt*nh;         // length of data
  int nx=nt*nq;         // length of model
  int nsparse;          // Number of non zero elements in index matrix
  float **vgrid;        // Velocity grid
  unsigned short **index;        // Migration matrix operator 
  void (*oper) (float *,float *,unsigned short **, int ,int , int, int, int);
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

  //////////////
  fprintf(stderr,"Inside eomigls_sparse\n*************");

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
  //irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
  irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  nsparse=nt*nq*nh;
  // allocate the big monster
  size_t size=sizeof(unsigned short);
  if ((index=(unsigned short **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");
  // Assign the elements of index
  build_index_inv(t,h,vgrid,nt,nh,nq,index,nsparse,t0);
  // Adjoint     
  eomig(mm,data[0],index,1,nt,nh,nq,nsparse);  
  // This is a filter for outliers or dominant bad data
  if (1){
    float datum;
    //eomig(mm,Wd,index,0,nt,nh,nq,nsparse);
    float q99=quest(0.90,nh*nt,data[0]);
    float q50=quest(0.50,nh*nt,data[0]);
    fprintf(stderr,"q99=%f, q50=%f \n",q99,q50);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++){
	datum=fabs(data[ih][it]);
	if (datum>q99) Wd[ih*nt+it]=q50/datum;
	else Wd[ih*nt+it]=1.0;
      }
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih*nt+it]=1.0;

  if ( FLAG && KEEP_WM ){
    fprintf(stderr,"Using Wm for the first time from adjoint\n");
    modelweight(mm,nx,norm,eps1,Wm);
    FLAG=0;
  }

  for (j=1;j<=iter_end;j++){
    //if (0) smoothing(mm,nt,nq,nl,nr,flag);
    // norm==1 ==> L1 , ==0  Cauchy else l2
    if (!KEEP_WM) modelweight(mm,nx,norm,eps1,Wm);
    //if (0) smoothing(Wm,nt,nq,nl,nr,flag);
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
  free2((void **) index);
  free2float(vgrid);
  free1float(J);
  free1float(filter);
  free1float(vel2);
  return;
}

void eomgetvel(float *m,float *mm,float *t, float *h, float *q, float *vel,
	     float **data, int nt, int nh, int nq, float fmax, 
	     float step, float eps, float eps1, float eps2, int itercg, 
	     int iter_end, int smooth, int testadj, int buffer, float pervmin,
	     float dperv, int norm, int itervel, float t0,float *Wm,float *Wd)
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
  unsigned short **index;        // Migration matrix operator 
  void (*oper) (float *,float *,unsigned short **, int ,int , int, int, int);
  float *J; // Cost function array
  oper=eomig;
  int restart=1;        // Every cg iteration restart the solution
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);

  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  

  for (it=0;it<nt;it++) m[it]=0.;
  for (it=0;it<ny;it++) Wd[it]=1.;
  for (it=0;it<nx;it++) Wm[it]=1;

  // Irregular velocity grid
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  nsparse=nt*nq*nh;

  // allocate the big monster
  size_t size=sizeof(unsigned short);
  if ((index=(unsigned short **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");

  for (int jj=1;jj<=itervel;jj++){  
    // Irregular velocity grid
    irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
    // Assign the elements of index
    build_index(t,h,vgrid,nt,nh,nq,index,nsparse,t0);
    // Adjoint     
    eomig(mm,data[0],index,1,nt,nh,nq,nsparse);
    
    for (j=1;j<=iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(mm,nx,norm,eps1,Wm);
      wpcgnr(eomig,nt,nh,nq,nsparse,mm,data[0],Wd,Wm,index,eps,step,
	     itercg,restart);
    }
    get_vel(mm,t,q,vel,vel2,nt,nq,buffer);
    for (it=0;it<nt;it++) vel[it]=vel2[it];
  }
  free2((void **) index);
  free2float(vgrid);
  smoothing(vel,nt,1,10,10,flag);
  free1float(vel2);
  return;
}



void get_vel(float *m, float *t, float *q, float *vel, float *vel2, 
	   int nt, int nq, int buffer) 
{
  /* Given an initial velocity vel and the radon model space 
  it computes the new axis vel2=vel+correction
  The correction is estimated using a weighted summation of the radon
  model.
  Input:
       m[iq*nt+it] Radon model
       t[it]       time axis
       q[iq]       Radon axis
       vel[it]     Initial velocity axis
       nt 
       nq 
       buffer      length/2 of the weighted summation
  Output:
       vel2[it]    corrected axis 

  Daniel Trad - UBC- January 2000. 
  */

  int it;
  int it2;
  int iq;
  int iv;
  int iqxnt;
  int nx=nt*nq;
  int qh=nq/2;
  float qtemp;
  float qtot;
  float dt=t[2]-t[1];  

  for(it=0;it<nt;it++) vel2[it]=vel[it];
  
  for(it=buffer;it<nt-buffer;it++){
    qtot=qtemp=0;
    for (it2=it-buffer;it2<it+buffer;it2++) 
      for(iq=0;iq<nq;iq++) qtot+=fabs(m[iq*nt+it2]);
    // If there is enough energy in the Radon space apply velocity correction
    if (qtot > 1e-2){
      for (it2=it-buffer;it2<it+buffer;it2++)
	for(iq=0;iq<nq;iq++)
	  qtemp+=(fabs(m[iq*nt+it2])*q[iq]);    
      vel2[it]=vel[it]+qtemp/qtot;
    }
    
    if (0)  
      fprintf(stderr,"vel[%d]=%f,vel2[%d]=%f,t=%f,qtemp=%f\n",it,vel[it],it,vel2[it],it*dt,qtemp);

  } 
  return;
}


void eomigls_sparse_large(float *m,float *mm,float *t, float *h, float *q, 
			  float *vel,float **data, int nt, int nh, int nq, 
			  float fmax, float step, float eps, float eps1, 
			  float eps2, int itercg, int iter_end, int smooth, 
			  int testadj, int buffer, float pervmin, float dperv,
			  int norm, float parmute, unsigned short mute,
			  float t0, float *Wm, float *Wd,  int centralq)
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
  int nmute;
  // Filter
  int nl=3;
  int nr=3;
  int flag=2;
  float *filter;

  //////////////
  fprintf(stderr,"Inside eomigls_sparse_large\n*************");

  vel2=ealloc1float(nt);
  filter=alloc1float(nl+nr+1);
  J=alloc1float(iter_end+1);
  
  for (it=0;it<ny;it++) Wd[it]=1;
  int qh=nq/2+1;

  // Irregular velocity grid
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  //irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
  irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  nsparse=nt*nh*nq;

  // allocate the big monster
  size_t size=sizeof(unsigned int);
  if ((index=(unsigned int **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");
  // Assign the elements of index
  build_index_inv(t,h,vgrid,nt,nh,nq,index,nsparse,t0);

  // Adjoint     
  eomig(mm,data[0],index,1,nt,nh,nq,nsparse);
  //eomig(mm,data[0],index,0,nt,nh,nq,nsparse);
  // This is a filter for outliers or dominant bad data

  if (1){
    //eomig(mm,Wd,index,0,nt,nh,nq,nsparse);
    float q99=quest(0.98,nh*nt,data[0]);
    float q50=quest(0.50,nh*nt,data[0]);
    fprintf(stderr,"q99=%f, q50=%f \n",q99,q50);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++) 
	if (fabs(data[ih][it])>q99) Wd[ih*nt+it]=0;//q50/q99;
    else Wd[ih*nt+it]=1.0;
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih*nt+it]=1.0;

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
  int side=0; // side 2 --> right mute ; side 1 left mute
  float **datarec=0;
  float dpd;
  float dpdp;
  float scale;
  float parmute1=2.5e-8;
  float parmute2=parmute;
  float mute1=mute;
  float mute2=mute;
  float *mmtemp=0;
  int nqleft=0; // Test variable, set in irregaxis.cpp to the same value
  mmtemp=ealloc1float(nt*nq);
  datarec=ealloc2float(nt,nh);
  
  if (mute1){
    side=2; 
    iq=0; while(fabs(q[iq])<parmute1) iq++; 
    if (side==2) nmute=nq-iq;
    else if (side==1) nmute=iq;
    fprintf(stderr,"MUTING 1 at nmute=%d************************\n",nmute);
    for (iq=0;iq<nq*nt;iq++) mmtemp[iq]=mm[iq];

    //taper(mmtemp,nt,nq,nmute,side);
    if (0){
      for (iq=nqleft;iq<nq;iq++) for (it=0;it<nt;it++) mmtemp[iq*nt+it]=0;
      if (0) plotgather_pipe(mmtemp,nq,nt,"mute1");
      eomig(mmtemp,datarec[0],index,0,nt,nh,nq,nsparse);
      plotgather_pipe(datarec,nh,nt,"mute1");
      fprintf(stderr,"substracting artifacts \n");
      dpd=dot(nh*nt,data[0],datarec[0]);
      dpdp=dot(nh*nt,datarec[0],datarec[0]);
      scale=dpd/dpdp;

      fprintf(stderr,"scale 1 ===>%f\n",scale);
      for (ih=0;ih<nh;ih++) 
	for (it=0;it<nt;it++)
	  data[ih][it]=data[ih][it]-scale*datarec[ih][it];
    }
    else{
      for (iq=0;iq<nqleft;iq++) for (it=0;it<nt;it++) mmtemp[iq*nt+it]=0;
      if (0) plotgather_pipe(mmtemp,nq,nt,"mute1");
      eomig(mmtemp,data[0],index,0,nt,nh,nq,nsparse);
      if (0) plotgather_pipe(data,nh,nt,"mute1");
      fprintf(stderr,"Muting artifacts \n");
    }
  }

  if (mute2){
    side=1;
    iq=0; while(fabs(q[iq])<parmute2) iq++; 
    if (side==2) nmute=nq-iq;
    else if (side==1) nmute=iq;
    fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);
    for (iq=0;iq<nq*nt;iq++) mmtemp[iq]=mm[iq];

    taper(mmtemp,nt,nq,nmute,side);
    if (0) plotgather_pipe(mmtemp,nq,nt,"mute 2");
    eomig(mmtemp,datarec[0],index,0,nt,nh,nq,nsparse);
    if (0) plotgather_pipe(datarec,nh,nt,"mute 2");
    fprintf(stderr,"substracting multiples \n");
    dpd=dot(nh*nt,data[0],datarec[0]);
    dpdp=dot(nh*nt,datarec[0],datarec[0]);
    scale=dpd/dpdp;
    fprintf(stderr,"scale 2===>%f\n",scale);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++)
	data[ih][it]=data[ih][it]-scale*datarec[ih][it];


  }

  if ((mute1==0)&&(mute2==0)) eomig(mm,data[0],index,0,nt,nh,nq,nsparse);
  if (smooth) smoothing(data[0],nt,nh,nl,nr,flag);



  //if (1) fprintf(stderr,"nt=%d,nh=%d\n",nt,nh);
  //  eomig1(data,h,m,t,vel);
  free1float(mmtemp);
  free2float(datarec);
  free2((void **) index);
  free2float(vgrid);
  free1float(J);
  free1float(filter);
  free1float(vel2);

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
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);

  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  

  for (it=0;it<nt;it++) m[it]=0.;
  for (it=0;it<ny;it++) Wd[it]=1.;
  for (it=0;it<nx;it++) Wm[it]=1;

  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");

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
    build_index(t,h,vgrid,nt,nh,nq,index,nsparse,t0);
    // Adjoint     
    eomig(mm,data[0],index,1,nt,nh,nq,nsparse);
    for (j=1;j<=iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(mm,nx,norm,eps1,Wm);
      wpcgnr(eomig,nt,nh,nq,nsparse,mm,data[0],Wd,Wm,index,eps,step,
		  itercg,restart);
    }
    get_vel(mm,t,q,vel,vel2,nt,nq,buffer);
    for (it=0;it<nt;it++) vel[it]=vel2[it];
  }

  free2((void **) index);
  free2float(vgrid);
  smoothing(vel,nt,1,10,10,flag);

  free1float(vel2);

  return;
}















