#include "su.h"
#include "clibrarytd.h"

void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin)

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
  int  j, iter, ih, iq;
  float *Wm;
  float *Wd;
  float dx_max;
  float dx_av;
  float qmaxt;
  float test;
  extern float factor; // Used to compute dq, i.e., dq=factor*dq;
  extern float step;
  extern int itercg;
  extern int testadj;
  extern int taper;
  unsigned short **index; // sparse matrix L
  int nsparse;
  //////////////////////////////////////////////////////////////////////////
  void (*radon) (float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon2) (float *,float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon3) (float *,float *,float *,float *,float *,float **,int ,int ,int,int);
  void (*radon3s) (float *m, float *d, unsigned short **index, int adj, int nt, int nh, int nq, int nsparse);
  
  if (rtmethod==3) {
    radon=radonhyp;
    radon2=radonhyp;
    radon3=radonhyp;
    radon3s=radonhyp;
  }
  else if (rtmethod==2) radon=radonpar;
  else if (rtmethod==1) radon=radonlin;
  
  /////////////////////////////////////////////////////////////////////////
 
  if ((Wm=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n"); 
  if ((Wd=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n"); 
  ///////////////////////////////////////////////////////////////////
  // Compute inv(CD)=WdT Wd       
  for (i=0;i<ny;i++) Wd[i]=1.;
  if (taper==1)
    for (ih=0;ih<5;ih++) 
      for (i=0;i<nt;i++){
	Wd[ih*nt+i]=1-exp(-(ih+.3));
	Wd[(nh-ih-1)*nt+i]=Wd[ih*nt+i];
      }
  
  //   Velocity axis //////////////////////////////////////////////////
  
  float **vgrid;
  float *perv;

  if (method==0 || method==1){ // Irregular velocity grid
    // Define irregular spacing along the horizontal
    // Here q  is perturbartion around the central velocity law
    // dperq is a parameter used to generate the increasing space
    // perqmin defines the minimum distance to the perturbation
    int nqh=nq/2;
       
    unsigned int it;
    if ((perv=alloc1float(nq))==NULL) err(" Cannot allocate perv");
    if ((vgrid=alloc2float(nt,nq))==NULL) err("Cannot allocate vgrid");
    perv[nqh]=0;
    for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
    for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
    for (iq=0;iq<nq;iq++) fprintf(stderr,"perv[%d]=%f\n",iq,perv[iq]);
    // Define irregualr qgrid as a sumation of 
    for (iq=0;iq<nq;iq++){
      for (it=0;it<nt;it++){
	vgrid[iq][it]=vel[it]*vel[it]+perv[iq]*perv[iq]+2*vel[it]*perv[iq];
        q[iq]=perv[iq];
      }
    } 
 

    ///////////////////////////////////////////////////////////////////////  

    fprintf(stderr,"Inside radtd nq=%d, nt=%d, nh=%d dq=%f\n",nq,nt,nh,dq);
  
    for (i=0;i<nx;i++) m[i]=0.;
 
  }
  if (method == 0){ 
    // count the number of elements we need for index
    nsparse=radhypsp(t,h,q,vgrid,nt,nh,nq);
    // allocate the big monster
    size_t size=sizeof(unsigned short);
    if ((index=(unsigned short **) alloc2(nsparse,2,size))==NULL) 
      err("Cannot allocate index\n");
    // Assign the elements of index
    radhypsp(t,h,q,vgrid,nt,nh,nq,index,nsparse);
    // Test the adjoint
    if (testadj) test=testadjop(radon3s,index,nt,nh,nq,nsparse);
    // Adjoint   
    radonhyp(m,d,index,1,nt,nh,nq,nsparse);
    // WTCGLS
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon3s,nt,nh,nq,nsparse,m,d,Wm,Wd,index,eps,step,eps1,eps2,itercg);
    }
    if (1){
      FILE *indexp;
      if ((indexp=fopen("index.bin","w"))==NULL) 
	err("Can't open index.bin\n");
      efwrite(index[0],2*nsparse,sizeof(unsigned short),indexp);
      //efwrite(index[0],nsparse,sizeof(unsigned short),indexp);
      //efwrite(index[1],nsparse,sizeof(unsigned short),indexp);
      efclose(indexp);
      if ((indexp=fopen("weight.bin","w"))==NULL)
	err("Can't open index.bin\n");
      efwrite(Wm,nx,sizeof(float),indexp);
      fclose(indexp);
    }
        



    // Let us kill the monster
    free2((void **) index);
    // Now again operators    
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq); 
  }
  else if (method==1){
    radon3(m,t,h,q,d,vgrid,1,nt,nh,nq);
    
    if (testadj) test=testadjop(radon3,t,h,q,vgrid,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon3,nt,nh,nq,t,h,q,m,d,Wm,Wd,vgrid,eps,step,eps1,eps2,itercg); 
    }
    radon3(m,t,h2,q,dint,vgrid,0,nt,nh2,nq);
  }

  if (method==0 || method==1){
    free2float(vgrid);  
    free1float(perv);
  }


  free1float(Wm);
  free1float(Wd);  
  return;
  
}






















