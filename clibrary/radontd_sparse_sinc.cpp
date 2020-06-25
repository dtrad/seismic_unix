#include "su.h"
#include "clibrarytd.h"
#include "radontd_sparse.h"


void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0,  inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI)

  /*
    RADONTD_SPARSE
    RADON TRANSFORM IN TIME DOMAIN
    
    Daniel Trad - UBC October 2000
    E-mail: dtrad@geop.ubc.ca
  */
{
  register int it;
  int  j, ih, iq, iw;
  float **Wm; // Model weights
  float **Wd;// Data weights
  float **index; // sparse matrix L
  int nsparse;
  int ntaper=5; // Number of lateral traces to taper 
  int it0=(int) (t0/dt+0.5);
  float test;
  char buf[80];
  int testadj=1;

  //////////////////////////////////////////////////////////////////////////
  // Define pointers to functions to be used with operators
  void (*radon3s) (float *m, float *d, float **index, int adj, int nt, 
		   int nh, int nq, int it0);

  // The actual functions to be used are
  radon3s=radonhyp;

  Wm=ealloc2float(nt,nq);
  Wd=ealloc2float(nt,nh);

  ///////////////////////////////////////////////////////////////////
  
  //if (inv.taperflag==1) taper(d,nt,nh,ntaper,0);
    
  //   Velocity axis //////////////////////////////////////////////////
  float **vgrid;
  if ((vgrid=alloc2float(nt,nq))==NULL) err("Cannot allocate vgrid");
  fprintf(stderr,"nq=%d,nt=%d,pervmin=%f,dperv=%f,centralq=%d\n",
	  nq,nt,pervmin,dperv,centralq);
  irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  save2dfile(vgrid,nq,nt,dt,"vgrid");
    //  else irreg_vel2inv_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid);
  for (iq=0;iq<nq;iq++) memset((void *)m[iq],(int)'\0',nt*FSIZE);
  
  nsparse=(nt-it0)*nq*nh; 
  fprintf(stderr,"nsparse=%d\n",nsparse);
  // allocate the big monster
  index=ealloc2float(nsparse,2);
  
  // Assign the elements of index
  build_index_slowness(t,h,q,vgrid,nt,nh,nq,index,it0);

  // Test the adjoint
  if (testadj) test=testadjop(radon3s,index,nt,nh,nq,it0);

  // Adjoint   
  radon3s(m[0],d[0],index,1,nt,nh,nq,it0);

  if (0){ 
    save_gather(m,nq,nt,dt,"model");
    system("suxwigb < model  title=\"plotgather\"");
  }

  // This is a filter for outliers or dominant bad data
  if (filtout){
    //radonhyp(m[0],Wd[0],index,0,nt,nh,nq,nsparse);
    float qup=quest(0.999,nh*nt,d[0]);
    float qmean=quest(0.50,nh*nt,d[0]);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++) 
	if (fabs(d[ih][it])>qup) Wd[ih][it]=qmean/qup;
	else Wd[ih][it]=1.0;
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1.0;

  fprintf(stderr,"inv.norm=%d,inv.eps1=%f,inv.eps2=%f,inv.itercg=%d,inv.iter_end=%d,inv.eps=%f,inv.restart=%d,inv.step=%f\n",inv.norm,inv.eps1,inv.eps2,inv.itercg,inv.iter_end,inv.eps,inv.restart,inv.step);    

  for (j=1;j<=inv.iter_end;j++){
    // norm==1 ==> L1 , ==0  Cauchy else l2
    modelweight(m[0],nq*nt,inv.norm,inv.eps1,Wm[0]);
    wpcgnr(radon3s,nt,nh,nq,it0,m[0],d[0],Wd[0],Wm[0],index,inv);
  }

  radon3s(m[0],d[0],index,0,nt,nh,nq,it0);

  free2float(index);
  free2float(vgrid);  
  free2float(Wm);
  free2float(Wd);  

  return;
  
}
void build_index_slowness(float *t, float *h, float *q, float **slow2, int nt, int nh, int nq,
			  float **index, int it0)
{
  register int it;
  int ih,iq;
  int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  int iqxnt,ihxnt;
  int nx=nt*nq;
  int ny=nt*nh;
  int nt0=nt-it0;
  int ns=nt0*nh*nq;


  /* The psarse matrix has times only from it0 to nt, i.e., nt0 elements,
     however teh data and model vectors have times from it to nt, hence
     the pointer for index uses nt0, the pointers for the data and model
     use nt. 
  */
  
  float dt=t[1]-t[0];
  float dt2=dt*dt;

  ttn=ealloc1float(nt);
  tnt=ealloc1float(nt);

  memset((void *) index[0],(int)'\0',ns*FSIZE);
  memset((void *) index[1],(int)'\0',ns*FSIZE);

  //for (it=0;it<ns;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	pxhxh=hxh*slow2[iq][it]; 
        ttn[it]=sqrt(t[it]*t[it]+pxhxh)/dt;
	index[0][ih*nq*nt0+iq*nt0+it-it0]=ttn[it];
      }
      // For sinc interpolation the times must be regularly sampled
      yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
      for (it=it0;it<nt;it++) index[1][ih*nq*nt0+iq*nt0+it-it0]=tnt[it];
    }
  }

  free1float(ttn);
  free1float(tnt);
  return;
}

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int it0)
{
  int  j;
  int it,ih,iq,iqxnt,ihxnt;
  
  int ny=nh*nt;
  int nx=nq*nt;
  int nt0=nt-it0;
  int ns=nt0*nq*nt;
  float *dint;
  float nqnh=sqrt(nq*nh);
  int nqxnt0=nq*nt0;

  dint=ealloc1float(nt);
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (ih=0;ih<nh;ih++){
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){
	iqxnt=iq*nt;
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt0,&index[1][ih*nqxnt0+iq*nt0],dint);
	for (it=it0;it<nt;it++) d[ihxnt+it]+=dint[it-it0];
      }
    }
    //for (j=0;j<ny;j++) d[j]/=nqnh;
  }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);   
    for (iq=0;iq<nq;iq++){
      iqxnt=iq*nt;
      for (ih=0;ih<nh;ih++){
	ihxnt=ih*nt;
	ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt0,&index[0][ih*nqxnt0+iq*nt0],dint);
	for (it=it0;it<nt;it++) m[iqxnt+it]+=dint[it-it0];
      }
    }
    //for (j=0;j<nx;j++) m[j]/=nqnh;
  }
  free1float(dint);
  return;
}



float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int it0)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  int iq;
  int ih;
  float test;
  int ny=nt*nh;
  int nx=nt*nq;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

  for (it=0;it<it0;it++)  for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=0;
  for (it=0;it<it0;it++)  for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=0;

  for (it=it0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=it0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  oper(mr2,dr1,index,1,nt,nh,nq,it0);
  oper(mr1,dr2,index,0,nt,nh,nq,it0);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"Test adjoint = %f \n",test);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}


void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here q  is perturbartion around the central velocity law
   dperq is a parameter used to generate the increasing space
   perqmin defines the minimum distance to the perturbation */
     
  float *perv;
  int it;
  int iq;
  int nqh=centralq;
  float vaux;
  fprintf(stderr,"pervmin=%f,dperv=%f\n",pervmin,dperv);
  if ((perv=alloc1float(nq+1))==NULL) err(" Cannot allocate perv");
  perv[nqh]=0;
  for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
  for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
  for (iq=0;iq<nq;iq++){
    //fprintf(stderr,"perv[%d]=%f\n",iq,perv[iq]);
    for (it=0;it<nt;it++){
      vaux=perv[iq];
      vgrid[iq][it]=1./(vel[it]*vel[it])+vaux;
      //if (vgrid[iq][it]<0) vgrid[iq][it]=0;//vgrid[iq][MAX(it-1,0)];
    }
    //    q[iq]=perv[iq];
    q[iq]=(vgrid[iq][0]);
    //fprintf(stderr,"vtop[%d]=%6.0f<======>,vbot[%d]=%6.0f\n",iq,sqrt(1./q[iq]),iq,sqrt(1./vgrid[iq][nt-1]));
  }
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}

void modelweight(float *m, int nx, int norm, float eps1, float *Wm)
  /*
  The right Wm from Cauchy is 
  Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
  But if M^-1 ATA x = M^-1 AT b is solved instead
  of the satndard form M=WmT *Wm 
  Actually it works even better with (I don't know why)
  Wm[i]=Wm[i]*Wm[i];
  if (Wm[i]>2) Wm[i]=2; 
  */

{ 
      int i;
      float maxm;
      maxm=fabs(m[isamax(nx,m,1)]);

      if (norm==1) for (i=0;i<nx;i++) Wm[i]=fabs(m[i])+eps1;
      else if(norm==0){
	if (maxm>1e-4) 
	  for (i=0;i<nx;i++) Wm[i]=(eps1*eps1+m[i]*m[i]/(maxm*maxm));
 	else for (i=0;i<nx;i++) Wm[i]=1e-3;
      }
      else for (i=0;i<nx;i++) Wm[i]=1.;     
      fprintf(stderr,"+++++++++++norm=%d,maxm=%f,eps1=%f,Wmmax=%f\n",norm,maxm,
	      eps1,Wm[isamax(nx,Wm,1)]);

      return;
}





























