#include "su.h"
#include "clibrarytd.h"

void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin, float alum, float blum, float clum, int norm, float t0, int mute, float parmute)

  /*
    subroutine: radtd.cpp
    RADON TRANSFORM IN TIME DOMAIN
    Three possible curves, linear, parabolic and hyperbolic
    are implemented in the file radon.cpp
    Inside this file the functions are called radonlin, radonpar, radonhyp
    At the beginning a pointer to function is declared and 
    associated with the chosen function. 
    This pointer is used as a reference to the operator 
    in WTCGLS, TESTADJ or any other function. 
    


    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
  */
{
  extern int nt,nh,nq, nx, ny, method, iter_end, reorth, rtmethod, nh2;
  extern float dt,dh,dq,eps1,eps2,thres,theta;
  register int i;
  float *J;
  int  j, iter, ih, iq;
  float *Wm; // Model weights
  float *Wd;// Data weights
  float *M; // preconditioner: prec=1 acting on model, =2 on data  
  float dx_max;
  float dx_av;
  float qmaxt;
  float test;
  extern float factor; // Used to compute dq, i.e., dq=factor*dq;
  extern float step;
  extern int itercg;
  extern int testadj;
  extern int taperflag;
  unsigned int **index; // sparse matrix L
  int nsparse;
  int restart=1;
  float *wdvec;
  int it0=(int) (t0/dt+0.5);
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////////////////////////////////////////////////////////////////
  void (*radon) (float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon2) (float *,float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon3) (float *,float *,float *,float *,float *,float **,int ,int ,int,int);
  void (*radon3s) (float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);
  void (*radon4s) (float *m, float *d, float **index, int adj, int nt, int nh, int nq); 
 
  if (rtmethod==3) {
    radon=radonhyp;
    radon2=radonhyp;
    radon3=radonhyp;
    radon3s=radonhyp;
    radon4s=radonhyp;
  }
  else if (rtmethod==2) radon=radonpar;
  else if (rtmethod==1) radon=radonlin;
  
  /////////////////////////////////////////////////////////////////////////
  if ((J=ealloc1float(iter_end+2))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");      
  if ((Wm=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n"); 
  if ((Wd=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");
  if ((M=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for M could not be allocated\n");
  for (i=0;i<ny;i++) M[i]=1.;
  // Lumley_precond(M,t,h,nt,nh,alum,blum,clum);

  ///////////////////////////////////////////////////////////////////
  // Compute inv(CD)=WdT Wd       

  for (i=0;i<ny;i++) Wd[i]=1.;
  if (taperflag==1)
    for (ih=0;ih<5;ih++) 
      for (i=0;i<nt;i++){
	Wd[ih*nt+i]=1-exp(-(ih+.3));
	Wd[(nh-ih-1)*nt+i]=Wd[ih*nt+i];
      }
  if ((wdvec=alloc1float(nh+1))==NULL)
    err("cannot allocate wdvec\n");
  for (ih=0;ih<nh;ih++) wdvec[ih]=Wd[ih*nt];
  save_vector(&wdvec[0],nh,"wd");
  free1float(wdvec);


  //   Velocity axis //////////////////////////////////////////////////
  
  float **vgrid;
  
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid);

  for (i=0;i<nx;i++) m[i]=0;

  if (method == 0){ 
    nsparse=(nt-it0)*nq*nh; 
    fprintf(stderr,"nsparse=%d\n",nsparse);
    // allocate the big monster
    size_t size=sizeof(unsigned int);
    if ((index=(unsigned int **) alloc2(nsparse+1,2,size))==NULL) 
      err("Cannot allocate index\n");
    // Assign the elements of index
    fprintf(stderr,"++++++++++++++++++1\n");
    build_index_rad(t,h,q,vgrid,nt,nh,nq,index,it0);
    fprintf(stderr,"++++++++++++++++++2\n");
    // Test the adjoint
    if (testadj) test=testadjop(radon3s,index,nt,nh,nq,nsparse);
    fprintf(stderr,"++++++++++++++++++3\n");
    // Adjoint   
    radonhyp(m,d,index,1,nt,nh,nq,nsparse);
    fprintf(stderr,"++++++++++++++++++4\n");    
    for (i=0;i<ny;i++) Wd[i]=1.;
    for (i=0;i<nx;i++) Wm[i]=1.;
    // WTCGLS

    for (j=1;j<=iter_end;j++){
      //if (0) smoothing(m,nt,nq,nl,nr,flag);
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m,nx,norm,eps1,Wm);
      //if (0) for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      J[j]=wpcgnr(radon3s,nt,nh,nq,nsparse,m,d,Wd,Wm,index,eps,step,itercg,restart);
      //if (0) wtcgls(radon3s,nt,nh,nq,nsparse,m,d,Wm,Wd,index,eps,step,eps1,eps2,itercg);
    }
    fprintf(stderr,"++++++++++++++++++5\n");
    //if (0) smoothing(m,nt,nq,nl,nr,flag);    
    // Make true the following condition to save  a file with INDEX
    if (0){
      FILE *indexp;
      if ((indexp=fopen("index.bin","w"))==NULL) 
	err("Can't open index.bin\n");
      efwrite(index[0],2*nsparse,sizeof(unsigned int),indexp);
      efclose(indexp);
      if ((indexp=fopen("weight.bin","w"))==NULL)
	err("Can't open index.bin\n");
      efwrite(Wm,nx,sizeof(float),indexp);
      fclose(indexp);
    }

    // Let us kill the monster
    free2((void **) index);

    if (mute){
      int nmute;
      iq=0; while(q[iq]<parmute) iq++; nmute=iq;
      fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
      taper(m,nt,nq,nmute,5);
    }

    // Now again operators    
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq); 
  }
  else if (method==1){
    radon3(m,t,h,q,d,vgrid,1,nt,nh,nq);
    //modelweight(m,nx,norm,eps1,Wm);    
    if (testadj) test=testadjop(radon3,t,h,q,vgrid,nt,nh,nq);
    //for (i=0;i<ny;i++) Wd[i]=1.;
    for (j=1;j<=iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m,nx,norm,eps1,Wm);    
      J[j]=wpcgnr(radon3,nt,nh,nq,t,h,q,m,d,Wd,Wm,vgrid,eps,step,itercg,restart); 
      //for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      //for (i=0;i<ny;i++) Wd[i]=1.;
      //wtcgls(radon3,nt,nh,nq,t,h,q,m,d,Wm,Wd,vgrid,eps,step,eps1,eps2,itercg); 
    }
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq);
  }
  if (method == 2){ 
    // allocate the big monster
    nsparse=nt*nh*nq;
    float **indexf;
    indexf=ealloc2float(nsparse,2);
    // Assign the elements of index
    radhypsp(t,h,q,vgrid,nt,nh,nq,indexf);
    if (testadj) test=testadjop(radon4s,indexf,nt,nh,nq);
    // Adjoint   
    radonhyp(m,d,indexf,1,nt,nh,nq);
    for (i=0;i<ny;i++) Wd[i]=1.;
    for (i=0;i<nx;i++) Wm[i]=1.;
    // WTCGLS

    for (j=1;j<=iter_end;j++){
      //if (0) smoothing(m,nt,nq,nl,nr,flag);
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m,nx,norm,eps1,Wm);
      J[j]=wpcgnr(radon4s,nt,nh,nq,m,d,Wd,Wm,indexf,eps,step,itercg,restart);
    }
    // Let us kill the monster
    free2float(indexf);

    // Now again operators    
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq); 
  }

  free2float(vgrid);  


  if (0) save_vector(&J[1],iter_end,"costfunc");

  free1float(M);
  free1float(Wm);
  free1float(Wd);  
  free1float(J);
  return;
  
}























