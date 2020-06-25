#include "radonfk.h"


void stoltz_wtcgls(float **data, float **model, float **Wd, float *h, int nh,  float *t, 
		   int nt,  float vel, inv_par inv)
{

  int  j;
  float **Wm; // Model weights
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;
  float sigmam;
  float sigmad=0;
  float *J;
  float **datap; // predicted data

  //////////////////////////////////////////////////////////////////////////
  // Set pointer to function for Conjugate Gradients
  void (*oper) (float **datain, float **dataout, int nt, int nh, float *t, 
		float *h, float vel, int adj);

  fprintf(stderr,"nt=%d\n",nt);  
  // Set Fourier transform parameters 

  oper=stoltzop2;

  Wm=ealloc2float(nt,nh);
  J=ealloc1float(inv.iter_end+1);
  datap=ealloc2float(nt,nh);
  zero_array(model,nh,nt);
  //  wdmask(Wd,nt,nh);  

  fprintf(stderr,"quantil1=%f, quantil2=%f\n",quantil1,quantil2);
  
  for (j=1;j<=inv.iter_end;j++){
    if (j==2) deviations(model[0],nh*nt,data[0],nh*nt,inv.norm,quantil1,
			 quantil2,&sigmam,&sigmad);
    fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
    weights(model[0],nt*nh,inv.norm,sigmam,Wm[0],j);    
    //wmmask(Wm,nt,nh);  
    J[j-1]=wpcgnr_mig(oper,model,data,nt,nh,t,h,vel,Wm,Wd,inv);
  }

  free1float(J);
  free2float(Wm);
  free2float(datap);

  return;
}

//////////////////////////////////////////////////////////////////////////////////////
/* 
     This function solves the system of equations 
     (WmT FH WdT M^{-1} Wd F Wm) m =  FH WdT Wd M^{-1} d
     oper is an operator implemented with the function
     where M is a  preconditioner acting on data space.
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);
     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )
     Wd allows to downweight bad data by using large values. 
     In general is a diagonal size(ndata) 
     M is the preconditioner. Diagonal size(nmodel) Large values
     penalize, small values focus the solution to a desired model value.
     M changes the null space, W does not. 
     Hence prior information about the model must be implemented by M. 
     
     Taken from Yousef Saad, pag 260: 
     Iterative methods for sparse linear systems
     W has bee added to the original algorihm 
     Daniel Trad - UBC March-2000
    
  */

#define plot 0
float wpcgnr_mig(void (*oper)  (float **data, float **model, int nt, int nh, float *t,
				float *h, float vel, int adj),
		 float **x, float **b,int nt, int nh, float *t, float *h,  
		 float vel, float **Wm, float **Wd, inv_par inv)
		 
 
{
  float normb, beta, alpha, alphanum, alphaden;
  int k,j=0,ih,it;
  float J; // Cost function and square residuals

  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.
  // rc Weighted residuals to compute cost function
  // rc Weighted solution to compute cost function

  float **r2,**r,**g,**s,**z,**w,*rho,*eta,rhold, **rc, **xc;
  float **dp; // test array for plots
  int nx=nt*nh;
  int ny=nx;  
  float rhobefore;
  float dt=t[1]-t[0];

  g=ealloc2float(nt,nh);
  s=ealloc2float(nt,nh);
  z=ealloc2float(nt,nh);
  r=ealloc2float(nt,nh);
  r2=ealloc2float(nt,nh);
  w=ealloc2float(nt,nh);
  rho=ealloc1float(inv.itercg+1);
  eta=ealloc1float(inv.itercg+1);
  xc=ealloc2float(nt,nh); 
  rc=ealloc2float(nt,nh);
  dp=ealloc2float(nt,nh);
  
  fprintf(stderr,"nt=%d\n",nt);
  fprintf(stderr,"inv.itercg=%d;inv.step=%f\n",inv.itercg,inv.step);
  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);
  inv.restart=1;
 
  if (0) ximageplotgather(Wm,nh,nt,dt,"Wm.su",0,"perc=100 legend=1");
  if (0)  xplotgather(b,nh,nt,dt,"datains.su",0,"bclip=1.20259 wclip=-1.175 xbox=600 legend=1&");
  
  if (inv.restart){
    zero_array(x,nh,nt);
    CequalAxB(r,Wd,b,nt,nh);
  }
  else{
    (*oper) (r,x,nt,nh,t,h,vel,0);
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++) 
	r[ih][it]=Wd[ih][it]*(b[ih][it]-r[ih][it]);
  }
  
  CequalAxB(r2,r,Wd,nt,nh);
  (*oper) (r2,g,nt,nh,t,h,vel,1);
  CequalAxB(z,g,Wm,nt,nh);
  
  normb=dot(nx,z[0],z[0]);
  rho[0]=1; 
     
  for (ih=0;ih<nh;ih++)
    for (it=0;it<nt;it++) 
      s[ih][it]=z[ih][it];

  rhold=inv.eps*2;
  k=0;
  while((rhold>inv.eps)&&(k<inv.itercg)){
    k++;
    
    (*oper) (w,s,nt,nh,t,h,vel,0);
       
    CequalAxB(w,w,Wd,nt,nh);
    alphanum=dot(nx,z[0],g[0]);
    alphaden=dot(ny,w[0],w[0]);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=inv.step;

    for (it=0;it<nt;it++){
      for (ih=0;ih<nh;ih++) x[ih][it]+=(alpha*s[ih][it]);

      for (ih=0;ih<nh;ih++){
	r[ih][it]-=(alpha*w[ih][it]);  
	r2[ih][it]=r[ih][it]*Wd[ih][it];    
      }

    }
    fprintf(stderr,"alpha=%e\n",alpha);
    
    (*oper) (r2,g,nt,nh,t,h,vel,1);
    
    CequalAxB(z,g,Wm,nt,nh);
    if (k>1) rhobefore=rho[k-1];
    else rhobefore=10^6;
    rho[k]=dot(nx,z[0],z[0])/normb;
    if ((rho[k]>rhobefore)&&(k>1)){
      fprintf(stderr,"rho[%d]=%e,===>rhobefore==%e\n",k,rho[k],rhobefore);
      break;
    }
    fprintf(stderr,"resm[%d]=%e,===TEST ===>res[%d]==%e\n",k,rho[k],k,dot(ny,r[0],r[0]));
    beta=dot(nx,z[0],g[0])/alphanum;
    CequalApluskxB(s,z,beta,s,nt,nh);    
    eta[k]=dot(nx,x[0],x[0]);

    if (plot) xplotgather(x,nh,nt,dt,"xiter.su",k,"perc=97 xbox=0 x1end=2.5 &");
    if (0)  xplotgather(r2,nh,nt,dt,"riter.su",k,"bclip=1.20259 wclip=-1.175 xbox=600 legend=1&");
    if (0){
      stoltzopinv2(dp,x,nt,nh,t,h,vel);
      xplotgather(dp,nh,nt,dt,"dp.su",k,"perc=97");
    }


  }
  
  if (0) save_vector(&rho[1],k,"rhofile");

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2
  
  /*****Cost Function *******/
  CequalAxB(xc,x,Wm,nt,nh);
  CequalAxB(rc,r,Wd,nt,nh);

  J=dot(nx,xc[0],xc[0])+dot(ny,rc[0],rc[0]);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);

  free2float(dp);
  free2float(xc);
  free2float(rc);
  free1float(eta);
  free1float(rho);
  free2float(w);
  free2float(r2);
  free2float(r);
  free2float(z);
  free2float(s);
  free2float(g);

  return(J);
}




void CequalAxB(float **C,float **A,float **B,int n1,int n2)
{
     int i1,i2;
     // Element wise multiplication of two matrices

     for (i2=0;i2<n2;i2++)
         for (i1=0;i1<n1;i1++)
             C[i2][i1]=A[i2][i1]*B[i2][i1];
          
     return;
}

void CequalApluskxB(float **C, float **A, float k, float **B, int n1,int n2)
{
     int i1,i2;
     // Element wise multiplication of two matrices

     for (i2=0;i2<n2;i2++)
         for (i1=0;i1<n1;i1++)
             C[i2][i1]=A[i2][i1]+k*B[i2][i1];
          
     return;
}

/* mask function designed to mask out energy mapping to unphysical areas
   of the model space */
void wmmask(float **Wm, int nt, int nq)
{

  int it, iq;
  // Define min trace to which signal is expected to map
  int limit=(int) (nq/3); 
  int slope;
  for (it=0;it<nt;it++){
    slope=(int) (30.0*(it/(1.0*nt)));
    for (iq=0;iq<(limit-slope);iq++){
      Wm[iq][it]*=0.1;
      Wm[nq-iq-1][it]*=0.1;
    }
  }
  return;
}


/* mask function designed to mask out energy mapping to unphysical areas
   of the model space */
void wdmask(float **Wd, int nt, int nh)
{

  int it, iq;
  // Define min trace to which signal is expected to map
  int limit=(int) (nh/2.1); 
  int slope;
  for (it=0;it<nt;it++){
    slope=(int) (160.0*(it/(1.0*nt)));
    if (limit<=slope) return;
    for (iq=0;iq<(limit-slope);iq++){
      Wd[iq][it]*=0.1;
      Wd[nh-iq-1][it]*=0.1;
    }
  }
  return;
}


void xplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[120];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 hbox=900 wbox=700 %s\n",s,s,num,s2);
  system(buf);
  return;
}

/* test for static variable used in fputtr */
void plot_after_stretch(float **d, int nh, int nt, int dt, char *s1, char *s2)
{
  xplotgather(d,nh,nt,dt,s1,0,s2);
  return;
}

void ximageplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[120];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 hbox=900 wbox=700 %s legend=1 \n",s,s,num,s2);
  system(buf);
  return;
}






