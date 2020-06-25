#include <math.h>
#include "su.h"
#include "clibrarytd.h"

float dot(int n, float *a, float *b);


float testadjop(void (*oper) (float *,float *,float *,float *,float *,int ,int ,int,int), float *t, float *h, float *q, int nt, int nh, int nq)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  float test;
  int ny=nt*nh;
  int nx=nt*nq;

  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
 
  for (it=0;it<ny;it++) dr1[it]=frannor();
  for (it=0;it<nx;it++) mr1[it]=frannor();

  
  oper(mr2,t,h,q,dr1,1,nt,nh,nq);

  oper(mr1,t,h,q,dr2,0,nt,nh,nq); 

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"++++Test adjoint = %f \n",test);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}



float testadjop(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), float *t, float *h, float *q, float *vel,int nt, int nh, int nq)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */
  oper(mr2,t,h,q,dr1,vel,1,nt,nh,nq);

  oper(mr1,t,h,q,dr2,vel,0,nt,nh,nq); 

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

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), float *t, float *h, float *q, float **vel,int nt, int nh, int nq)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */
  oper(mr2,t,h,q,dr1,vel,1,nt,nh,nq);

  oper(mr1,t,h,q,dr2,vel,0,nt,nh,nq); 

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



float testadjop(void (*oper) (float *,float *,unsigned short **,int ,int ,int, int, int),unsigned short **index,int nt, int nh, int nq, int nsparse)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */

  oper(mr2,dr1,index,1,nt,nh,nq,nsparse);
  oper(mr1,dr2,index,0,nt,nh,nq,nsparse);

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

float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int),unsigned int **index,int nt, int nh, int nq, int nsparse)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */

  oper(mr2,dr1,index,1,nt,nh,nq,nsparse);
  oper(mr1,dr2,index,0,nt,nh,nq,nsparse);

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

float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int, 
			      float *wavelet, int nw),unsigned int **index,int nt, int nh, 
		int nq, int nsparse, float *wavelet, int nw)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */

  oper(mr2,dr1,index,1,nt,nh,nq,nsparse,wavelet,nw);
  oper(mr1,dr2,index,0,nt,nh,nq,nsparse,wavelet,nw);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"Test adjoint = %f dp1=%f, dp2=%f \n",test,dp1,dp2);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}


float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int),float **index,int nt, int nh, int nq)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  oper(mr2,dr1,index,1,nt,nh,nq);
  oper(mr1,dr2,index,0,nt,nh,nq);

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

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int ns)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */

  oper(mr2,dr1,index,1,nt,nh,nq,ns);
  oper(mr1,dr2,index,0,nt,nh,nq,ns);

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



float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int, float *wavelet, int nw),float **index,int nt, int nh, int nq, int ns, float *wavelet, int nw)
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
  // Filter
  int nl=5;
  int nr=5;
  int flag=1;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  //smoothing(dr1,nt,nh,nl,nr,flag);
  //smoothing(mr1,nt,nq,nl,nr,flag);
  
  /*
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) 
    fprintf(stderr,"dr1[%d]=%f\n",ih*nt+it,dr1[ih*nt+it]);
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) 
    fprintf(stderr,"mr1[%d]=%f\n",iq*nt+it,mr1[iq*nt+it]); 
  */

  oper(mr2,dr1,index,1,nt,nh,nq,ns,wavelet,nw);
  oper(mr1,dr2,index,0,nt,nh,nq,ns,wavelet,nw);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;


  fprintf(stderr,"Test adjoint = %f , dp1=%f, dp2=%f \n",test,dp1,dp2);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}






