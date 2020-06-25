#include "su.h"
#include "clibrary.h"
#include <math.h>
void lsqrc(complex **L, complex **LH,complex *x,complex *b,float tol,
	   float eps1,int reorth,complex **U,complex **B,complex **V,
	   complex *W,complex *u,complex *v,complex *w,complex *q,
	   complex *p, complex *r,complex *temp2, 
	   float *eta,float *rho,float *gcv)
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
  
  float pythag2(float a, float b);
  float normb,nit,beta,alpha,alpha_old,beta_old,tempr;
  float rrho,c1,s1,theta,rho_bar,phi,phi_bar;
  float delta, gamma_bar, rhs, z_bar,gamma,c2,s2,z,xnorm;
  int k,i,j,in,num,iter,iter2;
  extern int nq,nh,itercg,iter_end;
  extern float step;
  complex temp, czero;
  int nx=nq;
  int ny=nh;
  czero.r=czero.i=0;  
  //for (i=0;i<nx;i++) x[i]=czero;


  for (i=0;i<nx;i++) 
    if (sqrt(abs(x[i])) >1.e-3) W[i]=1./sqrt(abs(x[i]));
    else W[i]=eps1;
  
  for(i=0;i<nx;i++) x[i]=czero;
  normb=sqrt(rcdot(ny,b,b));
  
  
  
  nit=itercg;
  for(i=0;i< ny;i++) for(j=0;j<nit;j++) U[i][j]=czero;
  for(i=0;i< nx;i++) for(j=0;j<nit;j++) V[i][j]=czero;      
  c2 = -1; s2 = 0; xnorm = 0; z = 0;
  // Prepare for LSQR iteration.
  for (i=0;i<nx;i++) {v[i] = czero; x[i] = czero; beta = normb;} 
  if (beta==0) err("Right-hand side must be nonzero\n");
  for (i=0;i<ny;i++) u[i] = b[i]/beta;
    
  if (reorth)  for (i=0;i<ny;i++) U[i][0] = u[i];
  Atimesx(q,LH,u,nq,nh);   //q=LH*u;  
  for (i=0;i<nx;i++) r[i] = q[i]/W[i] - beta*v[i];
  alpha = sqrt(rcdot(nx,r,r));
  for (i=0;i<nx;i++) v[i] = r[i]/alpha; 
  if (reorth)  for (i=0;i<nx;i++) V[i][0] = v[i]; 
  phi_bar = beta; rho_bar = alpha; 
  for (i=0;i<nx;i++) w[i] = v[i];     
    
  //Perform Lanczos bidiagonalization with/without reorthogonalization.
  for (iter=1;iter<nit;iter++){ 
    alpha_old = alpha; beta_old = beta;
    
    // Compute A*v - alpha*u.
    for (i=0;i<nx;i++) q[i] = v[i]/W[i];
    Atimesx(p,L,q,nh,nq);    //p=L*q;
    for (i=0;i<ny;i++) p[i]-=alpha*u[i];
    B[iter-1][1] = alpha;
    
    // Orthogonalize (if needed)
    if (reorth)
      for (temp=czero,j=0;j<=iter-1;j++){
	for (i=0;i<ny;i++) temp+=conjg(U[i][j])*p[i]; 
	for (i=0;i<ny;i++) temp2[i]=temp*U[i][j]; 
	for (i=0;i<ny;i++) p[i]=p[i]-temp2[i];
      }
    beta = sqrt(rcdot(ny,p,p));
    for (i=0;i<ny;i++) u[i] = p[i]/beta;
    B[iter-1][0] = beta;
      
    // Compute A'*u - beta*v.
    Atimesx(q,LH,u,nq,nh);
    //displayA(u,nh);
    for (i=0;i<nx;i++) r[i] = q[i]/W[i] - beta*v[i];
    alpha = sqrt(rcdot(nx,r,r));
     
    for (i=0;i<nx;i++) v[i] = r[i]/alpha; 
    if (reorth)  for (i=0;i<nx;i++) V[i][0] = v[i]; 
    phi_bar = beta; rho_bar = alpha; 
    for (i=0;i<nx;i++) w[i] = v[i];     
      
    // Orthogonalize (if needed)
    if (reorth)   
      for (temp=czero,j=0;j<=iter-1;j++){
	for (i=0;i<nx;i++) temp+=conjg(V[i][j])*r[i];
	for (i=0;i<nx;i++) temp2[i]=temp*V[i][j]; 
	for (i=0;i<nx;i++) r[i]=r[i]-temp2[i];
      } 
    alpha = sqrt(rcdot(nx,r,r));
    for (i=0;i<nx;i++) v[i] = r[i]/alpha;     
    
    // Store U and V (if needed)
    if (reorth){
      for (i=0;i<ny;i++) U[i][iter] = u[i]; 
      for (i=0;i<nx;i++) V[i][iter] = v[i]; 
    }
    // Construct and apply orthogonal transformation.
    rrho = pythag2(rho_bar,beta); 
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
    eta[iter-1] = pythag2(xnorm,z_bar);
    gamma = pythag2(gamma_bar,theta);
    c2 = gamma_bar/gamma; 
    s2 = theta/gamma;
    z = rhs/gamma; 
    xnorm = pythag2(xnorm,z);
    rho[iter-1] = fabs(phi_bar)/normb;
    //fprintf(stderr,"%d iter, phi_bar=%e, itercg=%d\n",iter,phi_bar,itercg);  
    fprintf(stderr,"%d iter. Misfit = %f, Model Norm = %f\n",iter,rho[iter-1],eta[iter-1]);
   
    // Update the solution.
    for(i=0;i<nx;i++){
      x[i] = x[i] + step*(phi/rrho)*w[i]; 
      w[i] = v[i] - step*(theta/rrho)*w[i];
    }
    //fprintf(stderr,"phi=%f,rrho=%f,theta=%f,nx=%d\n",phi,rrho,theta,nx);

    // Check for convergence
    if ((tol==0)&&(iter>1)){ // GCV criteria
      in = iter+1; //length(rho);
      for (tempr=ny+1,i=0;i<iter;i++){
	tempr-=1;//([m:-1:m-in+1].^2);
	  gcv[i] = (rho[i]*rho[i])/(tempr*tempr);
      }
      //fprintf(stderr,"gcv[%d]=%e,gcv[%d]=%e\n",iter-2,gcv[iter-2],iter-1,gcv[iter-1]);
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
  return;
}

float pythag2(float a, float b)
{
  return(sqrt(a*a+b*b));
}












