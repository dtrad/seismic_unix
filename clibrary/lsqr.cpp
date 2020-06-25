#include "su.h"
#include "clibrarytd.h"
#include <math.h>
void lsqr(float *t,float *qaxis,float *h,float *x,float *b,float tol,int reorth)
{
  /*
    % LSQR Solution of least squares problems by 
    % Lanczos bidiagonalization with/without reorthogonalization.
    % 
    % [x,rho,eta,U,B,V,nit] = wlsqr(A,W,b,k,tol,reorth)
    %
    % Input - A - The Matrix
    %         W - Weighting matrix
    %	  b - RHS vectors
    %	  k - Number of iterations
    %	  tol - Stopping criteria 0=GCV, else misfit
    % Output - x - Solution
    %	   rho - misfit
    %	   eta - model norm
    %          U,B,V - The matrixes from the Bidiagonalization
    %	   nit - Number of iterations till converge  
    % Reference: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
    % sparse linear equations and sparse least squares", ACM Trans.
    % Math. Software 8 (1982), 43-71.  */
  
  float pythag(float a, float b);
  float normb,nit,beta,alpha,alpha_old,beta_old,temp;
  float rrho,c1,s1,theta,rho_bar,phi,phi_bar;
  float delta, gamma_bar, rhs, z_bar,gamma,c2,s2,z,xnorm;
  int k,i,j,in,num,iter,iter2;
  extern int itercg,nx,ny,iter_end;
  extern float step, eps2;
  float **U,**B,**V,*W,*u,*v,*w,*q,*p,*r,*temp2,*eta,*rho,*gcv;
  int maxn;
  maxn=(nx > ny) ? nx : ny;

  if ((U=alloc2float(itercg,ny))==NULL)   //==> l(ny x itercg)
    err("cannot allocate memory for U\n");
  if ((V=alloc2float(itercg,nx))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for V\n");
  if ((B=alloc2float(2,itercg+1))==NULL)
     err("cannot allocate memory for B\n");
  if ((W=alloc1float(nx))==NULL) 
    err("cannot allocate memory for W\n");
  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for q\n");
  if ((v=alloc1float(nx))==NULL)
    err("cannot allocate memory for v\n");
  if ((u=alloc1float(ny))==NULL)
    err("cannot allocate memory for u\n");
  if ((w=alloc1float(nx))==NULL)
    err("cannot allocate memory for w\n");
  if ((r=alloc1float(nx))==NULL)
    err("cannot allocate memory for r\n");
  if ((p=alloc1float(ny))==NULL)
    err("cannot allocate memory for p\n");
  if ((temp2=alloc1float(maxn))==NULL)
    err("cannot allocate memory for temp2\n");
  if ((eta=alloc1float(itercg))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(itercg))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(itercg))==NULL)
    err("cannot allocate memory for gcv\n");
  
  for (i=0;i<nx;i++) x[i]=0.;
  if (reorth) fprintf(stderr,"Lanczos with reorthonalization\n");
  else fprintf(stderr,"Lanczos without reorthonalization\n");
  for (iter2=1;iter2<=iter_end;iter2++){ 
    for (i=0;i<nx;i++) 
      if (sqrt(fabs(x[i]))>1e-3) W[i]=1./sqrt(fabs(x[i]));
      else W[i]=eps2;
    if (iter2==1) {
      for(i=0;i<nx;i++) x[i]=0.;
      normb=sqrt(dot(ny,b,b));
    }
    
    
    nit=itercg;
    for(i=0;i< ny;i++) for(j=0;j<nit;j++) U[i][j]=0;
    for(i=0;i< nx;i++) for(j=0;j<nit;j++) V[i][j]=0;      
    c2 = -1; s2 = 0; xnorm = 0; z = 0;
    // Prepare for LSQR iteration.
    for (i=0;i<nx;i++) {v[i] = 0.; x[i] = 0.; beta = normb;} 
    if (beta==0) err("Right-hand side must be nonzero\n");
    for (i=0;i<ny;i++) u[i] = b[i]/beta;
    
    if (reorth)  for (i=0;i<ny;i++) U[i][0] = u[i];
    radonopi(q,t,h,qaxis,u,W);
    //rstack(t,qaxis,h,q,u,1);  //q=LH*u;  
    for (i=0;i<nx;i++) r[i] = q[i]/W[i] - beta*v[i];
    alpha = sqrt(dot(nx,r,r));
    for (i=0;i<nx;i++) v[i] = r[i]/alpha; 
    if (reorth)  for (i=0;i<nx;i++) V[i][0] = v[i]; 
    phi_bar = beta; rho_bar = alpha; 
    for (i=0;i<nx;i++) w[i] = v[i];     
    
    //Perform Lanczos bidiagonalization with/without reorthogonalization.
    for (iter=1;iter<nit;iter++){ 
      alpha_old = alpha; beta_old = beta;
      
      // Compute A*v - alpha*u.
      for (i=0;i<nx;i++) q[i] = v[i]/W[i];
      radonop(q,t,h,qaxis,p);
      //rstack(t,qaxis,h,q,p,0); //p=L*q;
      for (i=0;i<ny;i++) p[i]-=alpha*u[i];
      B[iter-1][1] = alpha;
      
      // Orthogonalize (if needed)
      if (reorth)
	for (temp=0,j=0;j<=iter-1;j++){
	  for (i=0;i<ny;i++) temp+=U[i][j]*p[i]; 
	  for (i=0;i<ny;i++) temp2[i]=temp*U[i][j]; 
	  for (i=0;i<ny;i++) p[i]=p[i]-temp2[i];
	}
      beta = sqrt(dot(ny,p,p));
      for (i=0;i<ny;i++) u[i] = p[i]/beta;
      B[iter-1][0] = beta;
      
      // Compute A'*u - beta*v.
      radonopi(q,t,h,qaxis,u,W);
      //rstack(t,qaxis,h,q,u,1);
      
      for (i=0;i<nx;i++) r[i] = q[i]/W[i] - beta*v[i];
      alpha = sqrt(dot(nx,r,r));
      for (i=0;i<nx;i++) v[i] = r[i]/alpha; 
      if (reorth)  for (i=0;i<nx;i++) V[i][0] = v[i]; 
      phi_bar = beta; rho_bar = alpha; 
      for (i=0;i<nx;i++) w[i] = v[i];     
      
      // Orthogonalize (if needed)
      if (reorth)   
	for (temp=0,j=0;j<=iter-1;j++){
	  for (i=0;i<nx;i++) temp+=V[i][j]*r[i]; 
	  for (i=0;i<nx;i++) temp2[i]=temp*V[i][j]; 
	  for (i=0;i<nx;i++) r[i]=r[i]-temp2[i];
	} 
      alpha = sqrt(dot(nx,r,r));
      for (i=0;i<nx;i++) v[i] = r[i]/alpha;     
      
      // Store U and V (if needed)
      if (reorth){
	for (i=0;i<ny;i++) U[i][iter] = u[i]; 
	for (i=0;i<nx;i++) V[i][iter] = v[i]; 
      }
      // Construct and apply orthogonal transformation.
      rrho = pythag(rho_bar,beta); 
      c1 = rho_bar/rrho;
      s1 = beta/rrho; 
      theta = s1*alpha; 
      rho_bar = -c1*alpha;
      phi = c1*phi_bar; 
      phi_bar = s1*phi_bar;
      
      // Compute solution norm and residual norm 
      
      delta = s2*rrho; 
      gamma_bar = -c2*rrho; 
      rhs = phi - delta*z;
      z_bar = rhs/gamma_bar; 
      eta[iter-1] = pythag(xnorm,z_bar);
      gamma = pythag(gamma_bar,theta);
      c2 = gamma_bar/gamma; 
      s2 = theta/gamma;
      z = rhs/gamma; 
      xnorm = pythag(xnorm,z);
      rho[iter-1] = fabs(phi_bar)/normb;
      fprintf(stderr,"%d iter, phi_bar=%e, itercg=%d\n",iter,phi_bar,itercg);  
      fprintf(stderr,"%d iter. Misfit = %f, Model Norm = %f\n",iter,rho[iter-1],eta[iter-1]);
      
      // Update the solution.
      for(i=0;i<nx;i++){
	x[i] = x[i] + step*(phi/rrho)*w[i]; 
	w[i] = v[i] - step*(theta/rrho)*w[i];
      }
      
      // Check for convergence
      if ((tol==0)&&(iter>1)){ // GCV criteria
	in = iter+1; //length(rho);
	for (temp=ny+1,i=0;i<iter;i++){
	  temp-=1;//([m:-1:m-in+1].^2);
	  gcv[i] = (rho[i]*rho[i])/(temp*temp);
	}
	fprintf(stderr,"gcv[%d]=%e,gcv[%d]=%e\n",iter-2,gcv[iter-2],
		iter-1,gcv[iter-1]);
	if (gcv[iter-2] < gcv[iter-1]){ 
	  fprintf(stderr,"GCV Criteria was reached in iteration %d\n",iter-1);
	  nit = iter-1;
	  //if (reorth){
	  //for (i=0;i<ny;i++) U[i][j] = u[i]; 
	  //for (i=0;i<nx;i++) V[i][j] = v[i]; 
	  //V = V(:,1:nit); 
	  //U = U(:,1:nit+1); 
	  //}
	  for (i=0;i<nx;i++) x[i] = x[i]/W[i];
	  return;
        } 
      }
      else{
	if (rho[iter-1]<tol){
	  fprintf(stderr,"Discrep was reached in iteration %d\n",iter-1);
	  nit = iter-1;
	  //if reorth, V = V(:,1:nit); U = U(:,1:nit+1);  end;
	  for (i=0;i<nx;i++) x[i] = x[i]/W[i];
	  return; 
	}
      }
      //if reorth, V = V(:,1:nit); U = U(:,1:nit+1);  end;
      for (i=0;i<nx;i++) x[i] = x[i]/W[i];
    }
    nit = iter-1;
    //return;
  }
  
  
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(temp2);
  free1float(p);
  free1float(r);
  free1float(w);
  free1float(u);
  free1float(v);
  free1float(q); 
  free1float(W);
  free2float(B);
  free2float(V); 
  free2float(U);       
  return;
}

float pythag(float a, float b)
{
  return(sqrt(a*a+b*b));
}












