#include <math.h>
#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"

float testadj_mig(float *t, float *h, float *vel);

float dot(int n, float *a, float *b);

void migration(float *m, float *t, float *h, float *vel, float *d); 

void modelling(float *m, float *t, float *h, float *vel, float *d);
 
void wtcglsmig(float *t, float *vel, float *h,float *x,float *b,float *Qp, float tol, float step, float eps1, float eps2, int itercg);

void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

/*
    subroutine: emigls.cpp
    Least Squares Equivalent Offset Migration

    Migrates the hyperbolic data[t][h] (a csp gather) to 
    a single trace m[t] located at midpoint zero, using
    a velocity function vel[t]. 
    First it is computed the adjoint model, and  a Covariance matrix
    is computed from it. WTCGLS solves the system 
    (LH*L + Cm) m = LH d

    Sparse approach (Thompson, 1984) 
    
    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
*/

void eomigls(float *t, float *vel, float *h, float *m,float **data,
	   float fmax, int nvel, float step, float eps1, float eps2, 
	   int itercg, int iter_end)
{
 
  
  register int it;
  int  j, iter;
  int ih;
  float t0=0;
  float *mtemp;
  int lfilter=50;
  float *d;             // temporal vector containing all data
  float *Wm;            // Factoriztion of inv model covariance matrix
  extern float dt;   
  extern float dhe;     // offset interval for equivalent offset
  extern int nt;        
  extern int nhcsp;     // size of CSP
  extern float eps;     // small number
  extern int testadj;   // =1 apply test adjoint
  extern int smooth; // =1 apply smmthing filter
  float test;
  float *vel2; // 1/vel^2
  int ny=nt*nhcsp;  // length of data
  int nx=nt;        // length of model

  // Filter
  int nl=3;
  int nr=3;
  int flag=2;
  float *filter;
  //////////////

  for (it=0;it<nt;it++) m[it]=0.;


  if ((d=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  if ((Wm=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");
  if ((mtemp=alloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for mtemp could not be allocated\n");   
  if ((vel2=alloc1float(nt))==NULL)
    err("cannot allocate memory for vel2 \n");  
  if ((filter=alloc1float(nl+nr+1))==NULL)
    err("cannot allocate memory for filter \n"); 

  for (it=0;it<nt;it++) vel2[it]=1/(vel[it]*vel[it]);
  
  for (ih=0;ih<nhcsp;ih++) for (it=0;it<nt;it++) d[ih*nt+it]=data[it][ih];

  for (it=0;it<nt;it++) Wm[it]=1;

  if (testadj) testadj_mig(t,h,vel2);

  migration(m,t,h,vel2,d);
  
  for (j=1;j<=iter_end;j++){
    for (it=0;it<nt;it++) Wm[it]=1./(eps2+fabs(m[it]))+eps1;
    wtcglsmig(t,vel2,h,m,d,Wm,eps,step,eps1,eps2,itercg);  
  }
  
  modelling(m,t,h,vel2,d);

  for (ih=0;ih<nhcsp;ih++) for (it=0;it<nt;it++) data[it][ih]=d[ih*nt+it];  
  
  if (smooth){
    smoothing(d,nt,nhcsp,nl,nr,flag);
    smoothing(m,nt,1,nl,nr,flag);
  }
  
  /*
  if (smooth){
    rwa_smoothing_filter (flag,nl,nr,filter);
    conv (nl+nr+1,-nl,filter,nt,0,m,nt,0,mtemp);
    memcpy(m,mtemp,nt*FSIZE);
  }
  */

  free1float(filter);
  free1float(vel2);
  free1float(Wm);
  free1float(mtemp);
  free1float(d);

  return;
}


void migration(float *m, float *t, float *h, float *velin2, float *d) 
{
  register int it;
  int ih,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh;
  int itint;
  int ihxnt;
  extern float dt;
  extern int nt;
  extern int nhcsp;
  float moveout;

  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");  
  for (it=0;it<nt;it++) m[it]=0;
  
  for (ih=0;ih<nhcsp;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (itau=0;itau<nt;itau++){
	moveout=hxh*velin2[itau];
	time=2*sqrt(t[itau]*t[itau]/4+moveout);
	itime[itau]=time/dt;
	dtemp[itau]=d[ihxnt+itau];
    }
    ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    for (itau=0;itau<nt;itau++) m[itau]=m[itau]+dint[itau];            	      
  }
  
 
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}

void modelling_out_dr(float *m, float *t, float *h, float *velin2, float *d) 
{
  /* This function models data according to output driven scheme
     i.e., every output sample is scanned to find the corresponding model.
     The reason to do this is for using intsinc8r routine, that interpolates
     from regular sample (for example m[tn]) to irregular m[tni]
     The input driven routine instead, computes a irregular 
     grid that must be interpolated to regular.
     I keep the negative times only to make the code simpler when using 
     intsinc but they are not used. 
     This way every itime[itau] corresponds to d[itau]. 
     For some reason this operator does not pass the test for adjoints. 
   */ 
  int it;
  int ih;
  int itau, ntau;
  float itfloat;
  int itint;
  int ihxnt;
  unsigned int j;
  float time,hxh,moveout,a1,a2;
  extern float dt;
  extern int nhcsp;
  extern int nt;
  int ny=nt*nhcsp;
  int nx=nt;
  float *dtemp;
  float *itime2;
  float *itime;


  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dtemp=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for dtemp could not be allocated\n");  
  if ((itime2=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for itime2 could not be allocated\n");
 
  for (it=0;it<ny;it++) d[it]=0;


  for (ih=0;ih<nhcsp;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;

    for (it=0;it<nt;it++) dtemp[it]=0;
    for (it=0;it<nt;it++) itime2[it]=0;
    j=0;
    for (itau=0;itau<nt;itau++){
      moveout=hxh*velin2[itau];
      time=t[itau]*t[itau]/4-moveout;
      if (time>=0){
	time=2*sqrt(time);
	itime[itau]=time/dt;
      }
      else{  // skip sqrt of negative time
        time=fabs(time);
	time=2*sqrt(time);
	itime[itau]=-1*time/dt;;
        
      }
      //fprintf(stderr,"itime[%d]=%f\n",itau,itime[itau]);
    }     
    ints8r(nt,1.0,0.,m,0.0,0.0,nt,itime,dtemp);

    for (itau=0;itau<nt;itau++)
      if (itau < nt && itime[itau] >= 0) d[ihxnt+itau]+=dtemp[itau];      
    
  }
  
  free1float(dtemp);
  free1float(itime2);
  free1float(itime);


  return;
}


void modelling(float *m, float *t, float *h, float *velin2, float *d) 
{
  int it;
  int ih,itau;
  float itfloat;
  int itint;
  int ihxnt;
  unsigned int j;
  float time,hxh,moveout,a1,a2;
  extern float dt;
  extern int nhcsp;
  extern int nt;
  int ny=nt*nhcsp;
  int nx=nt;
  float *dtemp;
  float *dtemp2;

  float fpass=0.20;
  float apass=0.9;
  float fstop=0.30;
  float astop=0.1;
  int npoles;
  float f3db;

  if ((dtemp=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for dtemp could not be allocated\n");  
  if ((dtemp2=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for dtemp2 could not be allocated\n");
 
  for (it=0;it<ny;it++) d[it]=0;

  //bfdesign(fpass,apass,fstop,astop,&npoles,&f3db);
  //fprintf(stderr,"npoles=%d,f3db=%f\n",npoles,f3db);

  for (ih=0;ih<nhcsp;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;

    for (it=0;it<nt;it++) dtemp[it]=0;
    for (it=0;it<nt;it++) dtemp2[it]=0;

    for (itau=0;itau<nt;itau++){
      if (fabs(m[itau])>1e-10){
	//fprintf(stderr,"m[%d]=%f\n",itau,m[itau]);
        moveout=hxh*velin2[itau];
	time=2*sqrt(t[itau]*t[itau]/4+moveout);
	itfloat=time/dt;
	itint=(int) floor(itfloat);
	a1=1-(itfloat-itint);
	a2=itfloat-itint;
	//j=ihxnt+itint;
	if ((itint+1) < nt ) {	
	  dtemp[itint]+=a1*m[itau];
	  dtemp[itint+1]+=a2*m[itau];
	}
      }      
    }
    //bflowpass (npoles,f3db,nt,dtemp,dtemp2);
    for (itau=0;itau<nt;itau++) d[ihxnt+itau]=dtemp[itau];
  }
  
  free1float(dtemp);
  free1float(dtemp2);


  return;
}


void wtcglsmig(float *t, float *vel, float *h,float *x,float *b,float *Qp, float tol, float step, float eps1, float eps2, int itercg)
  
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;
  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls
  extern int nt,nhcsp;
  int ny=nt*nhcsp;
  int nx=nt;
  //////////////////////////////////////////////////////////////////////  


  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1float(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1float(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1float(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(nx))==NULL)
     err("cannot allocate memory for gcv\n");   

  fprintf(stderr,"eps1=%f,eps2=%f,\n",eps1,eps2);
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));
  for (i=0;i<ny;i++) r[i]=b[i];
  /////////////////////////////////////////////////////////////
  ////// This call must be adapted in different wtcgsls versions
  migration(s,t,h,vel,r);
  //////////////////////////////////////////////////////////////
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Qp[i];
    q[i]=q1[i]/Qp[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){
    /////////////////////////////////////////////////////////////
    ////// This call must be adapted in different wtcgsls versions
    modelling(z,t,h,vel,Az);
    //////////////////////////////////////////////////////////////
    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < tol ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    //resold=resid;
    //////////////////////////////////////////////////////////////
    ////// This call must be adapted in different wtcgsls versions
    migration(s,t,h,vel,r);
    //////////////////////////////////////////////////////////////

    for(i=0;i<nx;i++){
      q1[i]=s[i]/Qp[i];
      q[i]=q1[i]/Qp[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    //fprintf(stderr,"j=%d,beta=%e\n",j,beta);
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(ny,r,r))/normb;
    fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }

  fprintf(stderr,"j=%d\n",j);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(Az);
  free1float(r);
  free1float(z1);
  free1float(z);
  free1float(x1);
  free1float(s);
  free1float(q1);
  free1float(q);        
  return;
}














