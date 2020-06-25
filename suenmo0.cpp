/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suradtdi:  $Date: December 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
void enmo(float *d,float *t,float *p,float *v,float *m, int adj);
float testadjenmo(float *t, float *q, float *v);
float dot(int n, float *a, float *b);
void wtcglsenmo(float *t, float *p, float *v, float *x,float *b,float *Qp, float tol, float step, float eps1, float eps2, int itercg);

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SURADTDI -Inverse High Resolution Hyperbolic Radon transform          ", 
"	                                                          	",
" 	   								",
" suradtd < stdin > stdout offsetfile=  [no optional parameters]	",
" 									",
" Required parameters:							",
" offsetfile=	ascii file of new offset values. If there are no        ",
" 		changes use                                             ",
"               sugethw output=geom key=offset < sudata>offsetfile	",
" stdin: Radon model as it is produced by suhrrt2.                      ", 
"        It is a sufile with header.                                    ",  
"        This sufile must have the radon parameter value in             ",
"        key=f2                                                         ",
"                                                                       ",
" stdout: is a  sufile with the header and traces. The offset is        ",
"         given by offset file. Except offset, all other header words   ",
"         copied from original traces. As the field geometry do not     ",
"         exist for interpolated traces the header must corrected if    ",
"         interpolation is preformed. If original and final nh          ",
"         are the same, for example for multiple removal, the header    ",
"         is preserved.                                                 ",
" Optional parameters:              					",
" rtmethod=3      1-LRT 2-PRT 3-HRT                                     ",
"									",
" Example : # Inverse Radon Transform                                   ",
"  suradtdi offsetfile=sudata.off   < sudatarad > sudatarec             ",
"                                                                       ",
"                                                                       ",
"									",
NULL};

/* Credits:
 *      Daniel Trad
 *
 * Trace header fields accessed: ns, dt, key=f2
 *
 */
/**************** end self doc ***********************************/

segy tr;
int nt,nv,np;
float dt,dv,dp;

int main(int argc, char **argv)
{

	int i, j;
	float *d, *m, temp; 
        float *p, *t, *v, vmin, vmax, t0=0.;
        float *Wm; // Model Weight
	extern int nt,nh,np;
	extern float dt,dv,dp;
	int model;
        float eps,eps1,eps2;
        int iter_end, itercg;
        float step;
        int testadj; // !=0 test adjoint

	// smoothing
	int smooth; // =1 apply time smoothing
	int nl=3;  //  npoints left hand side
	int nr=3;  //  npoints left hand side
	int flags=2;  // 1 rectangular, 2 triangular
	//////////////////////////////////////////////
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Get parameters */
	if (!getparint("model", &model))  model =1;     
	if (!getparint("smooth", &smooth))  smooth =0;     
	if (!getparint("nv", &nv))  nv =50;
	if (!getparfloat("dv", &dv))  dv =50;
	if (!getparfloat("vmin", &vmin))  vmin =1000;
	if (!getparfloat("eps1", &eps1))  eps1 = 1;
	if (!getparfloat("eps2", &eps2))  eps2 = 1;
	if (!getparfloat("eps", &eps))  eps = 0;
	if (!getparint("iter_end", &iter_end))  iter_end = 1;
	if (!getparint("itercg", &itercg))  itercg = 10;
	if (!getparfloat("step", &step))  step = .9;
	if (!getparint("testadj", &testadj))  testadj = 0;
	/* Get info from first trace */
	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	if (!tr.ntr) err("ntr header field must be set"); 

	np =(int) tr.ntr;	
        nt = (int) tr.ns;
	dt = ((float) tr.dt)/1000000.0;

        if ((!nt)||(!nv)||(!np)) err("Error reading nt,nh or nq");

	
	// Allocate memory for data and model
  
	if ((d=ealloc1float(np*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
	
	if ((m=ealloc1float(nv*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");

	if ((Wm=ealloc1float(nv*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");	
	if ((p=ealloc1float(np))==NULL)
	  fprintf(stderr,"***Sorry, space for p could not be allocated\n");
  
	if ((v=ealloc1float(nv))==NULL)
	  fprintf(stderr,"***Sorry, space for v could not be allocated\n");

	if ((t=ealloc1float(nt))==NULL)
	  fprintf(stderr,"***Sorry, space for t could not be allocated\n");         
	// Because we want to use same struct array for data and model 
	// the maximun number of traces between them will be taken.

       	 
	j=0;
       	/* Loop over traces */
	do {
                register int i;
		p[j]= (float) tr.f2;

		for (i=0;i<nt;i++){
		        d[j*nt+i]=(float) tr.data[i];
		}
		j++;
                if (j > np) err("Number of traces > %d\n",np);     
	} while (gettr(&tr));

        np=j;
	fprintf(stderr,"np=%d\n",np);
        fprintf(stderr,"nt%d,dt %f\n ",nt,dt);

	for (j=0;j<nv;j++) v[j]=vmin+j*dv;
	for (i=0;i<nt;i++) t[i]=t0+i*dt;
	if (testadj) testadjenmo(t,p,v);
	enmo(m,t,p,v,d,1);

	for (j=1;j<=iter_end;j++){
	  for (i=0;i<nv*nt;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
	  wtcglsenmo(t,p,v,m,d,Wm,eps,step,eps1,eps2,itercg);  
	}
	

        enmo(m,t,p,v,d,0);  

        if (model){
	  if (smooth) smoothing(m,nt,nv,nl,nr,flags);
        
	  j=0;
	  do{
	    tr.f2=(float) v[j];
	    tr.ntr=nv;
	    
	    for (i=0;i<nt;i++)
	      tr.data[i]=m[j*nt+i];
	    puttr(&tr);
	    j++;
	  } while(j<nv);
	}
	else{
	  if (smooth) smoothing(d,nt,np,nl,nr,flags);
        
	  j=0;
	  do{
	    tr.f2=(float) p[j];
	    tr.ntr=np;
	    for (i=0;i<nt;i++)
	      tr.data[i]=d[j*nt+i];
	    puttr(&tr);
	    j++;
	  } while(j<np);
	}
	
        fprintf(stderr,"nt%d, nv %d \n ",nt,nv);

	free1float(Wm);
	free1float(t);
	free1float(v);
	free1float(p);
	free1float(m);
	free1float(d);            
	return EXIT_SUCCESS;
}

void enmo(float *m,float *t,float *p,float *v,float *d, int adj)
{
  register int it;
  int ip;
  int iv;
  float cost;
  float sint;
  float time;
  int itime;
  extern float dt;
  extern int nt;
  extern int nv;
  extern int np;
  float *ftime;
  float *dtemp;
  float *dint;

  if ((ftime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
 
  if (adj) for (it=0;it<(nv*nt);it++) m[it]=0;
  else for (it=0;it<(np*nt);it++) d[it]=0;

  for (ip=0;ip<np;ip++)
    for (iv=0;iv<nv;iv++){
      sint=p[ip]*v[iv];
      cost=sqrt(1-sint*sint);
      for (it=0;it<nt;it++){
	if (adj){
	  ftime[it]=t[it]/dt*cost;
	  dtemp[it]=d[ip*nt+it];
	}
	else{
	  if (cost>1e-5){ 
	    ftime[it]=t[it]/dt/cost;
	    dtemp[it]=m[iv*nt+it];
	  }
          else{
	    ftime[it]=0;
            dtemp[it]=0;
	  }	
	} 
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,ftime,dint);
      if (adj) for (it=0;it<nt;it++) m[iv*nt+it]+=dint[it];
      else for (it=0;it<nt;it++) d[ip*nt+it]+=dint[it];
      
    }
  if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
  else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
  
  free1float(dtemp);
  free1float(ftime);
  free1float(dint);  
  
  return;

}


float testadjenmo(float *t, float *p, float *v)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  extern int np;
  extern int nv;
  extern int nt;

  int nx=nv*nt;
  int ny=np*nt;


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
  
  enmo(mr2,t,p,v,dr1,1);
  enmo(mr1,t,p,v,dr2,0); 

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0){
    fprintf(stderr,"test=%f\n",dp1/dp2);
    return(dp1/dp2);
  }
  else{
    fprintf(stderr,"test=%f\n",dp2);
    return(dp2);
  }

  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}

void wtcglsenmo(float *t, float *p, float *v, float *x,float *b,float *Qp, float tol, float step, float eps1, float eps2, int itercg)
 
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;
  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls
  extern int nt,np,nv;
  int ny=nt*np;
  int nx=nt*nv;
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
  enmo(s,t,p,v,r,1);
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
    enmo(z,t,p,v,Az,0);
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
    enmo(s,t,p,v,r,1);
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







