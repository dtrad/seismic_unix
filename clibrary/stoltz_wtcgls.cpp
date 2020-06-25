#include "interpfk.h"
#include "segy.h"


void stoltz_wtcgls(float **data, float **model, float *h, int nh,  float *t, int nt,
		   float vel, inv_par inv, complex **F, complex **FLS, float *wavelet, int nw)
{

  register int it;
  int  j, ih;
  float **Wm; // Model weights
  float **Wd;// Data weights
  int filtout=0;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;
  float sigmam;
  float sigmad;
  float *J;
  int plot=1;
  float **datap; // predicted data
  
  //////////////////////////////////////////////////////////////////////////
  // Set pointer to function for Conjugate Gradients
  void (*oper) (float **datain, float **dataout, int nt, int nh, float *t, 
		float *h, float vel, complex **F, complex **FLS, float *wavelet, int nw, 
		int adj);

  
  // Set Fourier transform parameters 

  oper=stoltzop;

  Wm=ealloc2float(nt,nh);
  J=ealloc1float(inv.iter_end+1);
  Wd=ealloc2float(nt,nh); 
  datap=ealloc2float(nt,nh);
  zero_array(model,nh,nt);
  
  //if (taperflag==1) taper(data,nt,nh,ntaper,0);
  fprintf(stderr,"++++++++++++++++++1\n");
  // This is a filter for outliers or dominant bad data
  if (filtout){
    float qup=quest(0.99,nh*nt,data[0]);
    float qmean=quest(0.50,nh*nt,data[0]);
    for (ih=0;ih<nh;ih++) 
      for (it=0;it<nt;it++) 
	if (fabs(data[ih][it])>qup) Wd[ih][it]=qmean/qup;
	else Wd[ih][it]=1.0;
  }
  else for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1.0;

  //AequalB(d2,data,nt,nh);
  // Compute the adjoint 
  if (0){
    stoltzop(data,model,nt,nh,t,h,vel,F,FLS,wavelet,nw,1);
    if (plot==1){
      save_gather_test(model,nh,nt,0.004,"migrated1.su");
      system("suxwigb < migrated1.su perc=100 title=migrated1 & ");
    }	 
  }
  
  fprintf(stderr,"quantil1=%f, quantil2=%f\n",quantil1,quantil2);
  for (j=1;j<=inv.iter_end;j++){
       fprintf(stderr,"iter = %d\n",j);
       if (j==2) deviations(model[0],nh*nt,data[0],nh*nt,inv.norm,quantil1,
			    quantil2,&sigmam,&sigmad);
       fprintf(stderr,"sigmam=%f, sigmad=%f\n",sigmam,sigmad);
       weights(model[0],nt*nh,inv.norm,sigmam,Wm[0],j);    
       TRACE;
       J[j-1]=wpcgnr_mig(oper,model,data,nt,nh,t,h,vel,Wm,Wd,inv,F,FLS,wavelet,nw);
       TRACE;
       if (plot==1){
	    save_gather(model,nh,nt,0.004,"migrated2.su");
	    system("suxwigb < migrated2.su perc=100 title=migrated2 &");
       }
       if (plot==1){
	    stoltzop(model,datap,nt,nh,t,h,vel,F,FLS,wavelet,nw,0);
	    save_gather_test(datap,nh,nt,0.004,"unmigrated2.su");
	    system("suxwigb < unmigrated2.su perc=100 title=unmigrated2 & ");
       }		
  }
  
  //radon(model[0],data[0],L,nt,nh,nq,dt,d2,m2,nfft,fmax,0);
  free1float(J);
  free2float(Wm);
  free2float(Wd);
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


float wpcgnr_mig(void (*oper)  (float **datain, float **dataout, int nt, int nh, 
				float *t, float *h, float vel, complex **F, complex **FLS, 
				float *wavelet, int nw, int adj),
		 float **x, float **b,int nt, int nh, float *t, float *h, float vel, 
		 float **Wm, float **Wd, inv_par inv, complex **F, complex **FLS, 
		 float *wavelet, int nw)
 
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
  int nq=nh;
  int nx=nt*nq;
  int ny=nt*nh;  
  TRACE;


  g=ealloc2float(nt,nq);
  s=ealloc2float(nt,nq);
  z=ealloc2float(nt,nq);
  r=ealloc2float(nt,nh);
  r2=ealloc2float(nt,nh);
  w=ealloc2float(nt,nh);
  rho=ealloc1float(inv.itercg+1);
  eta=ealloc1float(inv.itercg+1);
  xc=ealloc2float(nt,nq); 
  rc=ealloc2float(nt,nh);
  

  fprintf(stderr,"inv.itercg=%d;inv.step=%f\n",inv.itercg,inv.step);
  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);
  inv.restart=1;
  TRACE;
  if (inv.restart){
    zero_array(x,nh,nt);
    CequalAxB(r,Wd,b,nt,nh);
  }
  else{
    (*oper) (r,x,nt,nh,t,h,vel,F,FLS,wavelet,nw,0);
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++) 
	r[ih][it]=Wd[ih][it]*(b[ih][it]-r[ih][it]);
  }
  
  CequalAxB(r2,r,Wd,nt,nh);
  (*oper) (r2,g,nt,nh,t,h,vel,F,FLS,wavelet,nw,1);
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
    (*oper) (s,w,nt,nh,t,h,vel,F,FLS,wavelet,nw,0);
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

    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++){
	x[ih][it]+=(alpha*s[ih][it]);
	r[ih][it]-=(alpha*w[ih][it]);  
	r2[ih][it]=r[ih][it]*Wd[ih][it];    
      }
    }
    fprintf(stderr,"alpha=%e\n",alpha);
    if (1){
      save_gather(s,nh,nt,0.004,"modelx.su");
      system("suxwigb < modelx.su perc=100 title=modelx ");
    }

    (*oper) (r2,g,nt,nh,t,h,vel,F,FLS,wavelet,nw,1);
    CequalAxB(z,g,Wm,nt,nh);
    
    rho[k]=dot(nx,z[0],z[0])/normb;
    
    fprintf(stderr,"resm[%d]=%e,===TEST ===>res[%d]==%e\n",k,rho[k],k,dot(ny,r[0],r[0]));
    beta=dot(nx,z[0],g[0])/alphanum;
    CequalApluskxB(s,z,beta,s,nt,nh);    
    eta[k]=dot(nx,x[0],x[0]);
  }


  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2
  
  /*****Cost Function *******/
  CequalAxB(xc,x,Wm,nt,nh);
  CequalAxB(rc,r,Wd,nt,nh);

  J=dot(nx,xc[0],xc[0])+dot(ny,rc[0],rc[0]);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);

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

void save_gather_test(float **d, int nh, int nt, float dt, char* name)
{
  segy tr;
  
  int  itr;
  FILE* fp;
  
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr.tracl=itr+1;
      tr.dt=(int) (dt*1e6);
      tr.ns=nt;
      tr.ntr=nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  //fflush(fp);
  fclose(fp);

  return;
  
}














