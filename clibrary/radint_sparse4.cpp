#include "su.h"
#include "clibrarytd.h"
#include "inversion_par.h"


void radint_sparse4(float *t, float *q, float *h, float *h2, float **m,float **d,float **dint,int nt, int nh, int nq, int nh2, float dt, float *vel, float dperv, float pervmin, float t0, options_par opt, inv_par inv, int flagvel, int centralq, int filtout)

  /*
    RADINT
    RADON TRANSFORM IN TIME DOMAIN
    
    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
  */
{
  register int it;
  int  j, ih, iq;
  float **Wm; // Model weights
  float **Wd;// Data weights
  unsigned int **index; // sparse matrix L
  int nsparse;
  int ntaper=5; // Number of lateral traces to taper 
  int it0=(int) (t0/dt+0.5);
  float test;
  int nx=nt*nq;
  int ny=nh*nt;
  float fpeak=20;
  float *wavelet;
  int nw=30;
  char buf[80];
  //////////////////////////////////////////////////////////////////////////
  void (*radon3s) (float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw);

  void (*radon3) (float *,float *,float *,float *,float *,float **,int ,int ,int,int);
  radon3s=radonhyp;
  //radon3=radonhyp;  
  void (*convop)  (int, int, int, float *,int, float *,float *);
  convop=contran;
  // Generate a wavelet for the forward and adjoint operator
  wavelet=ealloc1float(nw);
  if (1){ 
    ricker1_wavelet (nw,dt,fpeak,wavelet);
    sprintf(buf,"xgraph < wavelet n=%d pairs=2 d1=1 style=normal title=\"wavelet\"",nw);
    plotvector(wavelet,nw,buf);  
  }
  fprintf(stderr,"Test dot product for contran\n");
  test=testadjop_conv(contran,nw,wavelet,nt,nt+nw);

/* This option requires velocity constant with time because the
   requirement that ttn[it] must increase monotonically with time (see
   yxtoxy.c ). It is like nmo with crossing events that truncates the
   offset. */

  if (flagvel==1) radon3=radonhyp_crude; //Use this instead until I can solved it. 
  else radon3=radonhyp_crude_inv;  
  Wm=ealloc2float(nt,nq);
  Wd=ealloc2float(nt,nh);

  ///////////////////////////////////////////////////////////////////
  // Compute inv(CD)=WdT Wd       

  //for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1;
  
  if (opt.taperflag==1) taper(d,nt,nh,ntaper,0);
    
  //   Velocity axis //////////////////////////////////////////////////
  float **vgrid;
  if ((vgrid=alloc2float(nt,nq+1))==NULL) err("Cannot allocate vgrid");
  if (flagvel==1) irreg_velaxis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  else if (flagvel==2) irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  else reg_slowness_axis(nq,nt,dperv,t,vel,q,vgrid);
  //  else irreg_vel2inv_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
  for (iq=0;iq<nq;iq++) memset((void *)m[iq],(int)'\0',nt*FSIZE);

  nsparse=(nt-it0)*nq*nh; 
  fprintf(stderr,"nsparse=%d\n",nsparse);
  // allocate the big monster
  size_t size=sizeof(unsigned int);
  if ((index=(unsigned int **) alloc2(nsparse+1,2,size))==NULL) 
    err("Cannot allocate index\n");

  // Assign the elements of index
  fprintf(stderr,"++++++++++++++++++1\n");
  if (flagvel==1) build_index_rad(t,h,q,vgrid,nt,nh,nq,index,it0);
  else build_index_rad_inv(t,h,q,vgrid,nt,nh,nq,index,it0);
  fprintf(stderr,"++++++++++++++++++2\n");
  if (0){
    for (it=0;it<nw;it++) wavelet[it]=0;
    wavelet[0]=0.5;
    wavelet[1]=1;
    wavelet[2]=0;
  }
  // Test the adjoint
  if (opt.testadj) test=testadjop(radon3s,index,nt,nh,nq,nsparse,wavelet,nw);
  fprintf(stderr,"++++++++++++++++++3\n");
  // Adjoint   
  radonhyp(m[0],d[0],index,1,nt,nh,nq,nsparse,wavelet,nw);
  //radonhyp(m[0],dint[0],index,0,nt,nh,nq,nsparse);
  //plotgather(m,nq,nt,dt,"suxwigb perc=90 title=\"plotgather\"");  
  save_gather(m,nq,nt,dt,"model");
  system("suxwigb < model perc=90 title=\"plotgather\"");

  float *dtemp=ealloc1float(nt+nw);
  

  if (0){
    for (iq=0;iq<nq;iq++){
      xcor(nt,0,m[iq],nw,0,wavelet,nt+nw-1,0,dtemp);
      //sprintf(buf,"xgraph < wavelet n=%d pairs=2 d1=1 style=normal title=\"wavelet\"",nt);
      //plotvector(m[iq],nt,buf);
      memcpy((void *) m[iq],(const void *) &dtemp[0],nt*sizeof(float));
    }
  }
  save_gather(m,nq,nt,dt,"model2");
  system("suxwigb < model perc=90 title=\"plotgather\"");
  //plotgather(m,nq,nt,dt,"suxwigb perc=90 title=\"plotgather\"");  
  // This is a filter for outliers or dominant bad data
  if (filtout){
    //radonhyp(m[0],Wd[0],index,0,nt,nh,nq,nsparse);
    float qup=quest(0.99,nh*nt,d[0]);
    float qmean=quest(0.50,nh*nt,d[0]);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++) 
	if (fabs(d[ih][it])>qup) Wd[ih][it]=qmean/qup;
	else Wd[ih][it]=1.0;
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1.0;


  fprintf(stderr,"inv.norm=%d,inv.eps1=%f,inv.eps2=%f,inv.itercg=%d,inv.iter_end=%d,inv.eps=%f,inv.restart=%d,inv.step=%f\n",inv.norm,inv.eps1,inv.eps2,inv.itercg,inv.iter_end,inv.eps,inv.restart,inv.step);    

  // WTCGLS
  if(0){
    int temp=inv.itercg;
    inv.itercg=5;
    for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) Wm[iq][it]=1;
    wpcgnr(radon3s,nt,nh,nq,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv,wavelet,nw);
    inv.itercg=temp;
  }

  for (j=1;j<=inv.iter_end;j++){
    // norm==1 ==> L1 , ==0  Cauchy else l2
    modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
    wpcgnr(radon3s,nt,nh,nq,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv,wavelet,nw);
  }
  fprintf(stderr,"++++++++++++++++++5\n");
  
  // Let us kill the monster
  free2((void **) index);
  //free(index[0]);
  //free(index);
  if (opt.mute){
    int nmute;
    if (flagvel==1 ){ 
      iq=0; while(fabs(q[iq])>opt.parmute) iq++;
      nmute=iq;
      fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
      taper(m,nt,nq,nmute,4);
    }
    else{
      iq=0; while(fabs(q[iq])<opt.parmute) iq++; 
      //nmute=nq-iq;
      fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
      //taper(m,nt,nq,nmute,2);
      taper(m,nt,nq,iq,1);
    }
  }
  //radonhyp(m[0],dint[0],index,0,nt,nh,nq,nsparse);
  
  for (iq=0;iq<nq;iq++){
    contran(0,0,nw,wavelet,nt,m[iq],dtemp);
    //conv(nt,0,m[iq],nw,0,wavelet,nt,0,dtemp);
    //contruc(0,0,nw,nt,m[iq],nw,wavelet,nt,dtemp);
    memcpy((void *) m[iq],(const void *) dtemp,nt*sizeof(float));
  } 
  free1float(dtemp);

  // Now again operators    
  radon3(m[0],t,h2,q,dint[0],vgrid,0,nt,nh2,nq); 
  
  free2float(vgrid);  
  free2float(Wm);
  free2float(Wd);  
  free1float(wavelet);

  return;
  
}





































