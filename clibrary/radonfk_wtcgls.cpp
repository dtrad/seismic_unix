#include "radonfk.h"
//#include <iostream.h>
#define plot 0
#define CLIP 1
#define VERBOSE 1
#define PLOTPAR "hbox=600 wbox=500 legend=1 &"

void stoltz_wtcgls(float **data, float **model, float **Wd, float *h, int nh,  float *t, 
		   int nt,  float vel, inv_par inv, float Wmthreshold, float **Wm)
{

  int  j;

  float quantil1=inv.eps1;
  float quantil2=inv.eps2;
  float sigmam;
  float sigmad=0;
  float *J;
  int verbose = VERBOSE ;
  //////////////////////////////////////////////////////////////////////////
  // Set pointer to function for Conjugate Gradients
  void (*oper) (float **datain, float **dataout, int nt, int nh, float *t, 
		float *h, float vel, int adj);

  fprintf(stderr,"nt=%d\n",nt);  
  oper=stoltzop2;
  J=ealloc1float(inv.iter_end+1);
  zero_array(model,nh,nt);
  //  wdmask(Wd,nt,nh);  

  fprintf(stderr,"quantil1=%f, quantil2=%f\n",quantil1,quantil2);
  
  for (j=1;j<=inv.iter_end;j++){
    if (j==2) deviations(model[0],nh*nt,data[0],nh*nt,inv.norm,quantil1,
			 quantil2,&sigmam,&sigmad);

    if (verbose){
	 fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
	 fprintf(stderr,"iteration number %d \n", j);
    }
    weights_test(model,nt,nh,inv.norm,sigmam,Wm,Wd,j,Wmthreshold);    
    //    weights(model[0],nt*nh,inv.norm,sigmam,Wm[0],j);    
    //wmmask(Wm,nt,nh);  
    if (j>1) xplotgather(Wm,nh,nt,t[1]-t[0],"WmBeforeLS.su",j,"perc=100 legend=1"); 
    J[j-1]=wpcgnr_mig(oper,model,data,nt,nh,t,h,vel,Wm,Wd,inv);
  }
  // If weights are used to design passmute zone it may be good
  // to calculate the weights once more after the last iteration
  // This makes the weights sharper though, so it may be worst.

  if (1) weights_test(model,nt,nh,inv.norm,sigmam,Wm,Wd,j,Wmthreshold);
  free1float(J);

  return;
}

void weights_test(float **m, int nt, int nh, int norm, float sigmam, float **Wm, float **Wd, int iter, float Wmthreshold)
  
{
     // This weighting function incorpores a new hyperparameter that seems
     // to overwrite the other ones, so it may be in the future the only
     // high res parameter required: Wmthreshold.
     // Wmthreshold defines the maximum possible model weight,
     // it is like a roof for the weights,
     // the lower it is, the lower the  resolution. 
     // It is like if in a mountain range we cut the peaks
     // at some height, and fill with the material all the valleys.
     // Also, the concept of areal weights is incorporated into this
     // weighting schemes.
     // In my experience, this is the most successful weighting scheme
     // so far for high res
     //
     // Daniel Trad - UBC, August 2002,

     int ih, it, ii, iih;
     float temp;
     int smooth=5;     // smoothing in the time direction
     int hsmooth = 3;  // smoothing in the offset direction
     int verbose = VERBOSE;
      
     if (iter==1){ 
       for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wm[ih][it]=1;
       return;
     }
     if (verbose) fprintf(stderr,"norm=%d,sigmam=%f,Wmthreshold=%f\n",norm,sigmam,Wmthreshold);
     if (norm==1) 
       for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wm[ih][it]=MAX(fabs(m[ih][it]),sigmam);
     else if(norm==0){
       for (ih=hsmooth;ih<nh-hsmooth;ih++){
	 for (it=smooth;it<nt-smooth;it++){
	   // Add all energy inside a window smooth x hsmooth
	   temp=0;
	   for (iih=-hsmooth;iih<hsmooth;iih++)
	     for (ii= -smooth;ii<smooth;ii++)
	       temp=temp+fabs(m[ih+iih][it+ii]);
	   // if the energy inside the window is larger than the roof
	   // keep it to the roof
	   Wm[ih][it]=MIN(Wmthreshold,(sigmam*sigmam+temp*temp));
	   // Next line is intended to up weight zero traces
	   // filling the near offset gap. It does not act if not zero traces
	   // inside the near offsets.
	   if ((Wd[ih][it]==0)&&(fabs(nh/2-ih)<5)) Wm[ih][it]=Wmthreshold;
	 }
       }
     }
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


float wpcgnr_mig(void (*oper)  (float **data, float **model, int nt, int nh, float *t,
				float *h, float vel, int adj),
		 float **x, float **b,int nt, int nh, float *t, float *h,  
		 float vel, float **Wm, float **Wd, inv_par inv)
		 
 
{
  float normb, beta, alpha, alphanum, alphaden;
  int k,ih,it;
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
  if (0)  xplotgather(b,nh,nt,dt,"datains.su",0,"clip=1 xbox=600 wbox=500 hbox=600 legend=1&");
  
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
    
    fprintf(stderr,"------------------------------iteration = %d\n",k);
    (*oper) (w,s,nt,nh,t,h,vel,0);
    if (plot)  xplotgather(w,nh,nt,dt,"witer.su",k,PLOTPAR);       
    CequalAxB(w,w,Wd,nt,nh);
    alphanum=dot(nx,z[0],g[0]);
    alphaden=dot(ny,w[0],w[0]);

    //if (alphaden < 0.) err("alphaden=%e\n",alphaden) break;
    
    fprintf(stderr,"alphanum=%e,alphaden=%e,k=%d\n", alphanum,alphaden,k);

    alpha=alphanum/alphaden;
    fprintf(stderr,"alpha=%e, alpha_scaled=%e\n",alpha,alpha*inv.step);    
    alpha*=inv.step;

    for (it=0;it<nt;it++){
      for (ih=0;ih<nh;ih++) x[ih][it]+=(alpha*s[ih][it]);

      for (ih=0;ih<nh;ih++){
	r[ih][it]-=(alpha*w[ih][it]);  
	r2[ih][it]=r[ih][it]*Wd[ih][it];    
      }

    }

    
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
    if (0) xplotgather(s,nh,nt,dt,"s.su",k,"xbox=0 hbox=600 wbox=500 legend=1 &");
    if (0) xplotgather(x,nh,nt,dt,"xiter.su",k,"xbox=0 &");
    if (0)  xplotgather(r2,nh,nt,dt,"riter.su",k,"wbox=500 xbox=600 hbox=600 legend=1&");
    if (0){
      stoltzopinv2(dp,x,nt,nh,t,h,vel);
      xplotgather(dp,nh,nt,dt,"dp.su",k,"perc=97");
    }


  }
  
  if (1) save_vector(&rho[1],k,"rhofile");

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
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  //sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 %s\n",s,s,num,s2);
  
  sprintf(buf,"suximage < %s title=%s%d legend=1 wbox=500 %s \n",s,s,num,s2);
  fprintf(stderr,buf);  
  system(buf);
  return;
}

void xplotdiff(float **d1, float **d2, int nh, int nt, float dt, char *s,  char *s2)
{
  char buf[180];
  int ih,it;
  float **d3=0;
  d3= ealloc2float(nt,nh);
  for (ih=0;ih<nh;ih++) for(it=0;it<nt;it++) d3[ih][it]=d1[ih][it]-d2[ih][it];

  save_gather(d3,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s legend=1 wbox=500 %s \n",s,s,s2);
  system(buf);
  
  free2float(d3);

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
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 hbox=900 wbox=700 %s legend=1 \n",s,s,num,s2);
  system(buf);
  return;
}






