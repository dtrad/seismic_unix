/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suradtdi:  $Date: December 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "clibrary.h"
#include "header.h"

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), float *t, float *h, float *q, float *vel,int nt, int nh, int nq);

float dot(int n, float *a, float *b);

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv);

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float *vel,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b, float *Wm, float *Wd, float *vel,float tol, float step, float eps1, float eps2, int itercg);

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUENMO - Velocity Analysis with sparseness for Slant Stack            ", 
"	                                                          	",
" 	   								",
" suenmo  < stdin > stdout                                      	",
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
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;
int method;
int main(int argc, char **argv)
{

	int it, j, ip, iv;
        register int i;
	float *d, *m, temp; 
        float *p, *t, *v, vmin, vmax, t0=0.;
        float *Wm; // Model Weights
	float *Wd; // residual Weights
	extern int nt,nh,np;
	extern float dt,dv,dp;
	int model;
        float eps,eps1,eps2;
        int iter_end, itercg;
        float step;
        int testadj; // !=0 test adjoint
	int taper;
	void (*radonellip) (float *m, float *d, unsigned short **index, 
			    int adj, int nt, int nh, int nq, int nsparse);


	/// Velocity Trend
        float *tvel;
	float *vel;
        float *velint;
	int ntvel;
	int nvel;
	int itv;
	
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
	if (!getparint("taper", &taper))  taper = 0;
        if (!getparint("method",&method)) method = 0; // sparse
	/* Introduce velocity trend to apply Hyp Vel Filtering */
	ntvel = countparval("tvel");
	if (ntvel==0) ntvel = 1;
	tvel = ealloc1float(ntvel);
	if (!getparfloat("tvel",tvel)) tvel[0] = 0.0;
	nvel = countparval("vel");
	if (nvel==0) nvel = 1;
	if (nvel!=ntvel) err("number of tmig and vmig must be equal");
	vel = ealloc1float(nvel);
	if (!getparfloat("vel",vel)) vel[0] = 2000.0;
	for (itv=1; itv<ntvel; ++itv)
		if (tvel[itv]<=tvel[itv-1])
			err("tvel must increase monotonically");
	////////////////////////////////////////////////////////

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

	if ((Wd=ealloc1float(np*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");
	
	if ((p=ealloc1float(np))==NULL)
	  fprintf(stderr,"***Sorry, space for p could not be allocated\n");
  
	if ((v=ealloc1float(nv))==NULL)
	  fprintf(stderr,"***Sorry, space for v could not be allocated\n");

	if ((t=ealloc1float(nt))==NULL)
	  fprintf(stderr,"***Sorry, space for t could not be allocated\n");    

	if ((velint=ealloc1float(nt))==NULL)
	  fprintf(stderr,"*Sorry, space for velint could not be allocated\n"); 


     
	// Because we want to use same struct array for data and model 
	// the maximun number of traces between them will be taken.

	headerfp = etmpfile();
	if (verbose) warn("using tmpfile() call");  

	/* Loop over traces */
	ip=0;
	do {
	  p[ip]= (float) tr.f2;
	  efwrite(&tr,HDRBYTES,1,headerfp);		
	  for (i=0;i<nt;i++){
	    d[ip*nt+i]=(float) tr.data[i];
	  }
	  ip++;
	  if (ip > np) err("Number of traces > %d\n",np);     
	} while (gettr(&tr));

	erewind(headerfp); 
        np=ip;

	if (verbose) fprintf(stderr,"np=%d,nt%d,dt %f\n",np,nt,dt);

	for (iv=0;iv<nv;iv++) v[iv]=vmin+iv*dv;
	for (i=0;i<nt;i++) t[i]=t0+i*dt;
	/* Create axis for velocities */
	intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);
	//for (it=0;it<nt;it++) fprintf(stderr,"velint[%d]=%f,t[%d]=%f\n",
	//			      it,velint[it],it,t[it]);


	if (testadj) testadjop(enmo1,t,p,v,velint,nt,np,nv);

	enmo1(m,t,p,v,d,velint,1,nt,np,nv);

	for (i=0;i<np*nt;i++) Wd[i]=1.;
        // Test: For data with gaps could be useful to taper the p=0 domain
	// where the gaps produce zeros data.
        if (taper) 
	  for (ip=(np/2-3);ip<np/2+3;ip++)
	    for (it=0;it<nt;it++)
	      Wd[ip*nt+it]=0.3;
	
	
	for (j=1;j<=iter_end;j++){
	  for (i=0;i<nv*nt;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
	  wtcgls(enmo1,nt,np,nv,t,p,v,m,d,Wm,Wd,velint,eps,step,eps1,eps2,itercg);  
	}
	
        enmo1(m,t,p,v,d,velint,0,nt,np,nv);  

        if (model){
	  if (smooth) smoothing(m,nt,nv,nl,nr,flags);
	  FILE *myfilep;
	  if((myfilep=fopen("model.su","w"))==NULL)
	    err("cannot open file=%s\n","model.su");        
	  iv=0;
	  do{
	    tr.f2=(float) v[iv];
	    tr.ntr=nv;
	    
	    for (i=0;i<nt;i++)
	      tr.data[i]=m[iv*nt+i];
	    fputtr(myfilep,&tr);
	    iv++;
	  } while(iv<nv);
          fclose(myfilep);
	}

	if (smooth) smoothing(d,nt,np,nl,nr,flags);
	ip=0;
	do{
	  efread(&tr,HDRBYTES,1,headerfp);
	  tr.f2=(float) p[ip];
	  tr.ntr=np;
	  for (i=0;i<nt;i++)
	    tr.data[i]=d[ip*nt+i];
	  puttr(&tr);
	  ip++;
	} while(ip<np);
	
	
        fprintf(stderr,"nt%d, nv %d \n ",nt,nv);

	free1float(velint);
        free1float(Wd);
	free1float(Wm);
	free1float(t);
	free1float(v);
	free1float(p);
	free1float(m);
	free1float(d);
	free1float(tvel);
	free1float(vel);
	
	efclose(headerfp);
           
	return EXIT_SUCCESS;
}

void enmo1_old(float *m,float *t,float *p,float *v,float *d,int adj, int nt, int np, int nv)
{
  register int it;
  int ip;
  int iv;
  float cost;
  float sint;
  float time;
  int itime;
  float dt=t[1]-t[0];
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
      
      if (fabs(sint)<0.9){
	cost=sqrt(1-(sint*sint));
	//fprintf(stderr,"cost=%f,iv=%d\n",cost,iv);
	//fprintf(stderr,"sint=%f,iv=%d\n",sint,iv);	
	for (it=0;it<nt;it++){
	  if (adj){
	    ftime[it]=t[it]/dt*cost;
	    dtemp[it]=d[ip*nt+it];
	  }
	  else{
	    ftime[it]=t[it]/dt/cost;
	    dtemp[it]=m[iv*nt+it];	
	  } 
	}
	ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,ftime,dint);
	if (adj) for (it=0;it<nt;it++) m[iv*nt+it]+=dint[it];
	else for (it=0;it<nt;it++) d[ip*nt+it]+=dint[it];	
      }
    }
  if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
  else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
  
  free1float(dtemp);
  free1float(ftime);
  free1float(dint);  
  
  return;

}

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv)
{
  register int it;
  int ip;
  int iv;
  float cost;
  float sint;
  float time;
  int itime;
  float dt=t[1]-t[0];
  float ftime;
  float dint;
  float a1;
  float a2;


  if (adj) for (it=0;it<(nv*nt);it++) m[it]=0;
  else for (it=0;it<(np*nt);it++) d[it]=0;

  for (it=0;it<nt;it++){
    for (ip=0;ip<np;ip++){
      for (iv=0;iv<nv;iv++){
	sint=p[ip]*(v[iv]+vel[it]);
      	if (fabs(sint)<0.95){
	  cost=sqrt(1-(sint*sint));
	  ftime=t[it]/dt*cost;
       	  if (adj){	
	    ints8r(nt,1.0,0,&d[ip*nt],0.0,0.0,1,&ftime,&dint);
	    m[iv*nt+it]+=dint;
	  }
	  else{
	    itime=(int) floor(ftime);
	    a2=ftime-itime;
	    a1=1-a2;
	    if (itime < nt)     d[ip*nt+itime]+=a1*m[iv*nt+it];
	    if ((itime+1) < nt) d[ip*nt+itime+1]+=a2*m[iv*nt+it];
	  } 
	}
      }
    }    
  }
  
  if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
  else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
 
  return;

}

int rad_ellip_sp(float *t,float *p,float *v, float *vel,int nt, int np, int nv)
{
  unsigned int it;
  int ip;
  int iv;
  float cost;
  float sint;
  float time;
  int itime;
  float dt=t[1]-t[0];
  float ftime;
  unsigned int j=0;
  for (it=0;it<nt;it++){
    for (ip=0;ip<np;ip++){
      for (iv=0;iv<nv;iv++){
	sint=p[ip]*(v[iv]+vel[it]);
      	if (fabs(sint)<0.95){
	  cost=sqrt(1-(sint*sint));
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (itime < nt) j++; 
	}
      }
    }    
  } 
  return(j);
}

void rad_ellip_sp(float *t,float *p,float *v, float *vel,int nt, int np, int nv, unsigned short **index, int nsparse)
{
  unsigned int it;
  int ip;
  int iv;
  float cost;
  float sint;
  float time;
  int itime;
  float dt=t[1]-t[0];
  float ftime;

  unsigned int j=0;
  for (it=0;it<nt;it++){
    for (ip=0;ip<np;ip++){
      for (iv=0;iv<nv;iv++){
	sint=p[ip]*(v[iv]+vel[it]);
      	if (fabs(sint)<0.95){
	  cost=sqrt(1-(sint*sint));
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (itime < nt){ 
	    index[0][j]=ip*nt+itime;
	    index[1][j]=iv*nt+it;
	    j++; 
	}
      }
    }    
  } 
  return;
}




void radonellip(float *m, float *d, unsigned short **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned long j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;

  if (!adj){
    for (j=0;j<ny;j++) d[j]=0;
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
  }
  else{
    for (j=0;j<nx;j++) m[j]=0;
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
  }
  return;
}









