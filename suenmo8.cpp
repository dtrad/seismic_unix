/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suradtdi:  $Date: December 1999  */

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "clibrary.h"
#include "header.h"
#include <time.h>
#include "inversion_par.h"

void irreg_velocity_grid(int nv, int nt, float pervmin, float dperv, float *t, float *vel,
			 float **vgrid, int centralq, int nreg);

int rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv,  unsigned short **nullspace, float *critangle);

void rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle);

void rad_ellip_sp_LI(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle);

void rad_ellip_sp_old(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, unsigned short **nullspace, float *critangle);

void radonellip(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

void radonellip_LI(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), float *t, float *h, float *q, float **vel,int nt, int nh, int nq);

float dot(int n, float *a, float *b);

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv);

void enmo1(float *m,float *t,float *p,float *v,float *d, float **vgrid,int adj, int nt, int np, int nv);

void rms2intvel(int n,float *t0, float *vs, float *v, float *h);

void myintlin (int nin, float xin[], float yin[], float yinl, float yinr, 
	int nout, float xout[], float yout[]);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

int tapercentre(float **data, int nh, int nt, int ntaper, int ntapertime, int centre);

void create_elliptical_taper(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, float **Wd, float *critangle);

void elliptical_stack_sinc(float *m,float *t,float *p,float *v,float *d, 
			   float **vel,int adj, int nt, int np, int nv);


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
inv_par inv;
segy tr;

char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
  FILE *myfilep;
  cwp_String modelfile=""; /* output sufile for the model */ 
  int j, ip, iv;
  register int it;
  float **d, **m, temp, **daux; 
  unsigned short **nullspace; 
  float *p, *t, *v, t0=0.;
  float **Wm; // Model Weights
  float **Wd; // residual Weights
  int nt,nh,np,nv;
  float dt,dv,dp;
  int model;
  int itercgfin;
  int testadj; // !=0 test adjoint

  int method; // =0 sparse; =1 operator
  int nsparse;
  unsigned int **index; // the head of the monster
  size_t sizeus=sizeof(unsigned int);

  // Taper
  int taper;
  int ntaper;
  int ntapertime;
  
  /// Velocity Trend
  float *tvel;
  float *vel;
  float *velint;
  int ntvel;
  int nvel;
  int itv;
  float **vgrid; // grid for velocity vgrid[iv][it]
  int nx;
  int ny;
  int centralq;
  float pervmin;
  int nreg; // defines how many traces besides the velocity trend have reg sampling

  // smoothing
  int smooth; // =1 apply time smoothing
  int nl=3;  //  npoints left hand side
  int nr=3;  //  npoints left hand side
  int flags=2;  // 1 rectangular, 2 triangular
  
  //////////////////////////////////////////////
  inv_par inv;
  inv.restart=1;

  //////////////////////////////////////////////////////////////////////////
  // Define pointers to functions to be used with operators

  void (*create_vgrid) (int nv, int nt, float pervmin, float dperv, float *t, float *vel,
			 float **vgrid, int centralq, int nreg);

  int (*count_sparse_elements) (float *t,float *p,float *v, float **vgrid,int nt, int np, 
				int nv, unsigned short **nullspace, float *critangle);

  void (*create_sparse_matrix) (float *t,float *p,float *v, float **vgrid,int nt, int np, 
				int nv, unsigned int **index, int nsparse, float *critangle);

  void (*elliptical_operator) (float *m, float *d, unsigned int **index, int adj, int nt, int nh,
			       int nq, int nsparse);

  void (*elliptical_operator_sinc) (float *m,float *t,float *p,float *v,float *d, 
				    float **vel,int adj, int nt, int np, int nv);


  ///////////////////////////////////////////////////////////////////////////
  // Define functions to use according to the method 
  create_vgrid=irreg_velocity_grid;
  count_sparse_elements=rad_ellip_sp;
  create_sparse_matrix=rad_ellip_sp_LI;
  elliptical_operator=radonellip_LI;
  elliptical_operator_sinc=elliptical_stack_sinc;

  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  
  /* Get parameters */
  if (!getparint("model", &model))  model =1;     
  if (!getparint("smooth", &smooth))  smooth =0;     


  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 0;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("itercgfin",&itercgfin)) itercgfin = 10 ;
  if (!getparfloat("step", &inv.step))  inv.step = .9;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("taper", &taper))  taper = 0;
  if (!getparint("method",&method)) method = 0; // sparse
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparint("norm", &inv.norm))  inv.norm =1; 
  if (!getparint("ntaper", &ntaper))  ntaper =15;
  if (!getparint("ntapertime", &ntapertime))  ntapertime = 200 ; 
  //////////////////// Vgrid ////////////////////
  if (!getparint("nv", &nv))  nv =50;
  if (!getparint("centralq",&centralq)) centralq=nv/2;
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dv", &dv))  dv =50;
  if (!getparint("nreg",&nreg)) nreg = 5;
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

  float *critang;
  float *critangle;
  float *hinterval;
  float *vinterval;

  critang=ealloc1float(nvel+10);
  critangle=ealloc1float(nt+1);
  vinterval=ealloc1float(nvel+1);
  hinterval=ealloc1float(nvel+1);
  rms2intvel(nvel,tvel,vel,vinterval,hinterval);
  critang[0]=0.99;

  for (itv=0;itv<nvel;itv++){ 
    critang[itv]=vel[itv]/vinterval[itv];
    fprintf(stderr,"vinterval[%d]=%f,hinterval[%d]=%f,critang[%d]=%f\n",
    	    itv,vinterval[itv],itv,hinterval[itv],itv,critang[itv]);
  }
  
  // Allocate memory for data and model
  d=ealloc2float(nt,np);
  daux=ealloc2float(nt,np);
  m=ealloc2float(nt,nv);
  Wm=ealloc2float(nt,nv);
  Wd=ealloc2float(nt,np);
  p=ealloc1float(np);
  v=ealloc1float(nv+1);
  t=ealloc1float(nt);
  velint=ealloc1float(nt);
  vgrid=ealloc2float(nt,nv+1);
 
  if ((nullspace=(unsigned short**) alloc2(nt,np,sizeof(unsigned short)))==NULL)
    fprintf(stderr,"Cannot allocate nullspace\n");
  
  // Because we want to use same struct array for data and model 
  // the maximun number of traces between them will be taken.
  
  headerfp = etmpfile();
  if (verbose) warn("using tmpfile() call");  
  
  /* Loop over traces */
  ip=0;
  do {
    p[ip]= (float) tr.f2;
    efwrite(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) d[ip],(const void *) tr.data,nt*FSIZE);
    ip++;
    if (ip > np) err("Number of traces > %d\n",np);     
  } while (gettr(&tr));
  
  erewind(headerfp); 
  np=ip;

  nx=nv*nt;
  ny=np*nt;
  
  if (verbose) fprintf(stderr,"np=%d,nt%d,dt %f\n",np,nt,dt);
  for (it=0;it<nt;it++) t[it]=t0+it*dt;
  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);
  create_vgrid(nv,nt,pervmin,dv,t,velint,vgrid,centralq,nreg);  
    /* 
       float dv2=dv*dv;
       for (iv=0;iv<nv;iv++)
       for (it=0;it<nt;it++){
       vgrid[iv][it]=velint[it]*velint[it]+(iv-((int) nv/2))*dv2;
       }      
    */
  //    Interpolate the sin of the critical angle
  for (it=0;it<nvel;it++)  fprintf(stderr,"tvel[%d]=%f,critang[%d]=%f\n",it,tvel[it],it,critang[it]);
  TRACE;
  if(1){
    it=0;
    for(itv=0;itv<nvel;itv++){
      while((t[it]<tvel[itv])&&(itv<nvel)){
	critangle[it]=critang[itv];
	it++;
      }
    }
  }  
  // Cancel critangle filter by making it equal to 1
  for (it=0;it<nt;it++) critangle[it]=0.99;
  TRACE;

  //intlin(nvel,tvel,critang,critang[0],critang[nvel-1],nt,t,critangle);
  free1float(critang);
  //  for (it=0;it<nt;it++)  fprintf(stderr,"critangle[%d]=%f\n",it,critangle[it]);  
  TRACE;
  free1float(vinterval);
  free1float(hinterval);
  TRACE;
  // Assign data weights
  for (ip=0;ip<np;ip++) for (it=0;it<nt;it++) Wd[ip][it]=1;
  // Test: For data with gaps could be useful to taper the p=0 domain
  // where the gaps produce zeros data.

  if (taper==1||taper==3){ 
    TRACE;
    int centre;
    float min=fabs(p[0]);
    for (ip=0;ip<np;ip++){
      if(fabs(p[ip])<min){
	min=fabs(p[ip]);
	centre=ip;
      }
    }
    tapercentre(Wd,np,nt,ntaper,ntapertime,centre);
    if (taper==1){
      save_gather(Wd,np,nt,dt,"Wd");
      system("suxwigb < Wd clip=1.5  title=\"plotgather\" &");
    }
  }
  
  float *wdvec;
  if ((wdvec=alloc1float(np+1))==NULL)
    err("cannot allocate wdvec\n");
  for (ip=0;ip<np;ip++) wdvec[ip]=Wd[ip][0];
  save_vector(&wdvec[0],np,"wd");
  free1float(wdvec);  
  int flagsparse=1;
  // Faster sparse matrix multiplication method==0

  // Count number of non zero elements for the operator index
  nsparse=count_sparse_elements(t,p,v,vgrid,nt,np,nv,nullspace,critangle);
  fprintf(stderr,"Nonzero elements for index=%d\n",nsparse);
  // Allocate the big monster
  if ((index=(unsigned int **) alloc2(nsparse+10,2,sizeof(unsigned int)))==NULL)
    err("Cannot allocate the monster index\n");
  create_sparse_matrix(t,p,v,vgrid,nt,np,nv,index,nsparse,critangle);
  
  if (0) if (testadj) testadjop(elliptical_operator,index,nt,np,nv,nsparse);
    
  if(taper==2){
    float **mtemp=ealloc2float(nt,nv);
    for (iv=0;iv<nv;iv++) memset((void *)mtemp[iv],(int) '\0',nt*FSIZE);
    for (it=0;it<nt;it++) mtemp[nv/2][it]=1;
    elliptical_operator(mtemp[0],Wd[0],index,0,nt,np,nv,nsparse);
    save_gather(Wd,np,nt,dt,"Wd");
    system("suxwigb < Wd clip=1.5  title=\"plotgather\" &");
    free2float(mtemp);
  }
  if(taper==3){
    create_elliptical_taper(t,p,v,vgrid,nt,np,nv,Wd,critangle);
    save_gather(Wd,np,nt,dt,"Wd");
    system("suxwigb < Wd  title=\"plotgather\" &");
  }

  save_gather(d,np,nt,dt,"d1");
  system("suxwigb < d1  title=\"plotdata\" &");
  
  if (1) testadjop(elliptical_operator_sinc,t,p,v,vgrid,nt,np,nv);

  if (0){ 
    elliptical_stack_sinc(m[0],t,p,v,d[0],vgrid,1,nt,np,nv);
    save_gather(m,nv,nt,dt,"m1");
    system("suxwigb < m1  title=\"plotadjoint\" &");
    elliptical_stack_sinc(m[0],t,p,v,d[0],vgrid,0,nt,np,nv);
    save_gather(d,np,nt,dt,"d2");
    system("suxwigb < d2  title=\"plotdata\" &");
  }

  if (0){
    elliptical_operator(m[0],d[0],index,1,nt,np,nv,nsparse);
    save_gather(m,nv,nt,dt,"m2");
    system("suxwigb < m2  title=\"plotadjoint2\" &");
    elliptical_operator(m[0],d[0],index,0,nt,np,nv,nsparse);
    save_gather(d,np,nt,dt,"d2");
    system("suxwigb < d2  title=\"plotdata\" &");
  }

  if (0){
    elliptical_operator(m[0],d[0],index,1,nt,np,nv,nsparse);
    for (j=1;j<=inv.iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
      if (j==inv.iter_end) inv.itercg=itercgfin;      
      inv.J[j]=wpcgnr(elliptical_operator,nt,np,nv,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv);
    }
    elliptical_operator(m[0],daux[0],index,0,nt,np,nv,nsparse);   
  }

  if (1){
    elliptical_stack_sinc(m[0],t,p,v,d[0],vgrid,1,nt,np,nv);
    for (j=1;j<=inv.iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
      if (j==inv.iter_end) inv.itercg=itercgfin;      
      inv.J[j]=wpcgnr(elliptical_operator_sinc,nt,np,nv,t,p,v,m[0],d[0],Wd[0],Wm[0],vgrid,inv);
    }
    elliptical_stack_sinc(m[0],t,p,v,daux[0],vgrid,0,nt,np,nv);
  }

  free2((void **) index);
 
  for (ip=0;ip<np;ip++) 
    for (it=0;it<nt;it++){
      if (nullspace[ip][it]) d[ip][it]=daux[ip][it];
      else d[ip][it]=0;
    }


  if (model){
    if (smooth) smoothing(m[0],nt,nv,nl,nr,flags);
    FILE *myfilep;
    if((myfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);        
    iv=0;
    do{
      //tr.f2=(float) SGN(vgrid[iv][0])*sqrt(fabs(vgrid[iv][0]));
      tr.f2=(float) vgrid[iv][0];
      tr.ntr=nv;
      
      for (it=0;it<nt;it++)
	tr.data[it]=m[iv][it];
      fputtr(myfilep,&tr);
      iv++;
    } while(iv<nv);
    fclose(myfilep);
  }
  
  if (smooth) smoothing(d[0],nt,np,nl,nr,flags);
  ip=0;
  do{
    efread(&tr,HDRBYTES,1,headerfp);
    tr.f2=(float) p[ip];
    tr.ntr=np;
    for (it=0;it<nt;it++)
      tr.data[it]=d[ip][it];
    puttr(&tr);
    ip++;
  } while(ip<np);
  
  if (1) save_vector(&inv.J[1],inv.iter_end,"costfunc");
  
  fprintf(stderr,"nt%d, nv %d \n ",nt,nv);

  free1float(critangle);
  free2float(vgrid);
  free1float(velint);
  free2float(Wd);
  free2float(Wm);
  free1float(t);
  free1float(v);
  free1float(p);
  free2float(m);
  free2((void **) nullspace);
  free2float(daux);
  free2float(d);
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
	    ints8r(nt,1.0,0,&d[ip],0.0,0.0,1,&ftime,&dint);
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
  
  if (adj) 
    for (iv=0;iv<nv;iv++) 
      for (it=0;it<nt*nv;it++) 
	m[iv*nt+it]/=sqrt(np*nv);
  else 
    for (ip=0;ip<np;ip++)
      for (it=0;it<nt*np;it++) 
	d[ip*nt+it]/=sqrt(nv*np);
 
  return;

}

void enmo1(float *m,float *t,float *p,float *v,float *d, float **vgrid,int adj, int nt, int np, int nv)
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
	sint=p[ip]*p[ip]*vgrid[iv][it];
      	if (fabs(sint)<0.99){
	  cost=sqrt(1-sint);
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
  
  //if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
  //else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
 
  return;

}

int rad_ellip_sp(float *t,float *p,float *v,float **vgrid,int nt, int np, int nv)
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
	sint=p[ip]*p[ip]*vgrid[iv][it];
      	if (fabs(sint)<0.99){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (itime < nt) j++; 
	}
      }
    }    
  }
  return(j);
}


int rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned short **nullspace, float *critangle)
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

  for (ip=0;ip<np;ip++) 
    for (it=0;it<nt;it++) 
      nullspace[ip][it]=0;

  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    for (iv=0;iv<nv;iv++){
      for (it=0;it<nt;it++){
	sint=p[ip]*p[ip]*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (it < nt){ 
	    nullspace[ip][itime]=1;
	    j++;
	  }
	}
      }
    } 
  }  
  return(j);
}

int rad_ellip_sp_LI(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned short **nullspace, float *critangle)
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

  for (ip=0;ip<np;ip++) 
    for (it=0;it<nt;it++) 
      nullspace[ip][it]=0;

  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    for (iv=0;iv<nv;iv++){
      for (it=0;it<nt;it++){
	sint=p[ip]*p[ip]*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (it < nt){ 
	    nullspace[ip][itime]=1;
	    j++;
	  }
	}
      }
    } 
  }  
  return(j);
}

void create_elliptical_taper(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, float **Wd, float *critangle)
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

  //for (ip=0;ip<np;ip++) memset((void *) Wd[ip],(int) '\0',nt*FSIZE);

  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    for (iv=nv/2;iv<nv/2+1;iv++){
      for (it=0;it<nt;it++){
	sint=p[ip]*p[ip]*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) (ftime+0.5);
	  Wd[ip][itime]*=sin(cost*PI/2);
	}
	else Wd[ip][itime]*=sin(cost*PI/2);
      }
    } 
  }  
  return;
}

void rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle)
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
  memset((void *)index[0],(int)'\0',nsparse*sizeof(unsigned int));
  memset((void *)index[1],(int)'\0',nsparse*sizeof(unsigned int));  
  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    for (iv=0;iv<nv;iv++){
      for (it=0;it<nt;it++){
	sint=p[ip]*p[ip]*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (it < nt){ 
	    index[0][j]=ip*nt+itime;
	    index[1][j]=iv*nt+it;
	    j++;
	  }
	}
      }
    } 
  }  
  return;
}

void rad_ellip_sp_LI(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle)
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
  memset((void *)index[0],(int)'\0',nsparse*sizeof(unsigned int));
  memset((void *)index[1],(int)'\0',nsparse*sizeof(unsigned int));  
  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    for (iv=0;iv<nv;iv++){
      for (it=0;it<nt;it++){
	sint=p[ip]*p[ip]*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=10*t[it]/dt*cost;
          itime=(int) (ftime+0.5);
	  if (it < nt){ 
	    index[0][j]=ip*10*nt+itime;
	    index[1][j]=iv*nt+it;
	    j++;
	  }
	}
      }
    } 
  }  
  return;
}



void radonellip(float *m, float *d, unsigned int **index, int adj, int nt, int np, int nv, int nsparse)
{
  unsigned long j;
  unsigned int ny=np*nt;
  unsigned int nx=nv*nt;
  unsigned int it;
  fprintf(stderr,"radonellip\n");

  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
  }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
  }
  if (0){ 
    if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
    else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
  }
  m[0]=0;
  d[0]=0;
  //if (adj) for (it=0;it<nx;it++) if (fabs(m[it])>100) fprintf(stderr,"m[%d]=%f\n",it,m[it]); 
  return;
}


void radonellip_LI(float *m, float *d, unsigned int **index, int adj, int nt, int np, int nv, int nsparse)
{
  unsigned long j;
  unsigned int ny=np*nt;
  unsigned int nx=nv*nt;
  unsigned int it;
  int id1;
  int im1;
  float a;
  float b;

  fprintf(stderr,"radonellip_LI\n");
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++){
      b=0.1*index[0][j];
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	d[id1]+=(1.0-a)*m[im1];
	d[id1+1]+=a*m[im1];
      }      
    }
  }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++){
      b=0.1*index[0][j];
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	m[im1]+=(1-a)*d[id1]+a*d[id1+1];
      }
    }
  }
  if (0){ 
    if (adj) for (it=0;it<nt*nv;it++) m[it]/=sqrt(np*nv);
    else for (it=0;it<nt*np;it++) d[it]/=sqrt(nv*np);
  }
  m[0]=0;
  d[0]=0;

  return;
}

void rad_ellip_sp_old(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, unsigned short **nullspace, float *critangle)
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
  float pxp;
  int itt;

  for (ip=0;ip<np;ip++) 
    for (it=0;it<nt;it++) 
      nullspace[ip][it]=0;

  memset((void *)index[0],(int)'\0',nsparse*sizeof(unsigned int));
  memset((void *)index[1],(int)'\0',nsparse*sizeof(unsigned int));  
  unsigned int j=0,jr=0;

  for (ip=0;ip<np;ip++){
    pxp=p[ip]*p[ip];
    for (iv=0;iv<nv;iv++){
      for (it=0;it<nt;it++){
	sint=pxp*vgrid[iv][it];
	if (fabs(sint)<critangle[it]){
	  cost=sqrt(1-sint);
	  ftime=t[it]/dt*cost;
          itime=(int) floor(ftime+0.5);
	  if (it < nt){
	    itt=ip*nv*nt+iv*nt+it;
	    index[0][itt]=ip*nt+itime;
	    index[1][itt]=iv*nt+it;
	    nullspace[ip][itime]=1; 
	  }
	}
      }
    } 
  }  
  return;
}

int tapercentre(float **data, int nh, int nt, int ntaper, int ntapertime, int centre)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  float *tapertime;

  int k;
  float s;
  int it, ih;
  if (ntaper<=0) return EXIT_SUCCESS;
  taper = ealloc1float(ntaper);
  tapertime = ealloc1float(nt);

  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    //taper[k] = s*s;
    taper[k]= s;//pow(s,4);
    //fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  for (k = 0; k < ntapertime; ++k) {
    s = sin(k*PI/(2*nt));
    //taper[k] = s*s;
    tapertime[k]= s;
    //fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  for(ih=0;ih<ntaper;ih++){
    for (it=0;it<ntapertime;it++) data[centre+ih][it]*=taper[ih]*tapertime[it];
    for (it=0;it<ntapertime;it++) data[centre-ih][it]*=taper[ih]*tapertime[it]; 
  }
  free1float(taper);
  free1float(tapertime);

  return EXIT_SUCCESS;
}

void irreg_velocity_grid(int nq, int nt, float pervmin, float dperv, float *t, float *vel,
			 float **vgrid, int centralq, int nreg)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here dperv is perturbation around the central velocity law
   dperv is a parameter used to generate the increasing space
   pervmin defines the minimum perturbation */
     

  float *perv;
  int it;
  int iq;
  int nqh=centralq;
  float vaux;
  float *incperv;   // vector for increasing velocity perturbation
  int inp;          // index for increasing velocity perturbation
  float temp;
  fprintf(stderr,"pervmin=%f,nreg=%d,dperv=%f\n",pervmin,nreg,dperv);

  perv=ealloc1float(nq+1);
  incperv=ealloc1float(nq+1);

  perv[nqh]=0;
  memset((void *) perv,(int) '\0',nq*FSIZE);

  for (iq=0;iq<nq;iq++) incperv[iq]=(iq-nqh)*pervmin*dperv;;
  
  for (inp=nqh,iq=nqh;iq<nq;iq++){
    if (abs(iq-nqh)<(nreg)) perv[iq+1]=perv[iq]+pervmin;
    else perv[iq+1] = perv[iq] + incperv[inp++] + pervmin;
  }
  for (inp=nqh,iq=nqh-1;iq>=0;iq--){
    if (abs(nqh-iq)<nreg*4) perv[iq]=perv[iq+1]-pervmin;
    else perv[iq]=perv[iq+1]+incperv[inp--]-pervmin;
  }
  for (iq=0;iq<nq;iq++){
    for (it=0;it<nt;it++){
      vaux=perv[iq];
      if ((temp=vel[it]*vel[it]+vaux) > 0) vgrid[iq][it]=temp;
      else  vgrid[iq][it]=0;
    }
  }

  free1float(perv); 
  free1float(incperv); 

  return;
}

void elliptical_stack_sinc(float *m,float *t,float *p,float *v,float *d, 
			   float **vgrid,int adj, int nt, int np, int nv)
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
  float *vec;
  float *t2;
  float dt2;
  int nt2;

  vec=ealloc1float(nt);
  t2=ealloc1float(nt);

  if (!adj) memset((void *) d,(int)'\0',nt*np*FSIZE); 
  else memset((void *) m,(int)'\0',nt*nv*FSIZE);  

  unsigned int j=0,jr=0;
  TRACE;

  for (iv=0;iv<nv;iv++){
    for (ip=0;ip<np;ip++){
      memset((void *) vec,(int)'\0',nt*FSIZE);  
      sint=p[ip]*p[ip]*vgrid[iv][0];
      cost=sqrt(1-sint);
      dt2=dt*cost;
      for (it=0;it < nt;it++) t2[it]=it*dt2;
      nt2=(int) (nt*cost);

      if (!adj){
	  ints8r(nt,dt2,0,&m[nt*iv],0.0,0.0,nt,t,vec);
	  for (it=0;it<nt;it++) d[ip*nt+it]+=vec[it];
      }
      else{
	ints8r(nt,dt,0,&d[nt*ip],0.0,0.0,nt,t2,vec);
	for (it=0;it<nt;it++) m[iv*nt+it]+=vec[it];
      } 
    }  
  }
  TRACE;
  free1float(vec);
  free1float(t2);
  TRACE;
  return;
}

void elliptical_stack_sinc_orig(float *m,float *t,float *p,float *v,float *d, 
			   float **vgrid,int adj, int nt, int np, int nv)
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
  float *vec;
  float *t2;
  float dt2;
  int nt2;

  vec=ealloc1float(nt);
  t2=ealloc1float(nt);

  if (!adj) memset((void *) d,(int)'\0',nt*np*FSIZE); 
  else memset((void *) m,(int)'\0',nt*nv*FSIZE);  

  unsigned int j=0,jr=0;
  TRACE;

  for (iv=0;iv<nv;iv++){
    for (ip=0;ip<np;ip++){
      memset((void *) vec,(int)'\0',nt*FSIZE);  
      sint=p[ip]*p[ip]*vgrid[iv][0];
      cost=sqrt(1-sint);

      if (!adj){
	if (cost<1 || 1){
	  dt2=dt*cost;
	  for (it=0;it < nt;it++) t2[it]=it*dt2;
	  nt2=(int) (nt*cost);
	  ints8r(nt,dt2,0,&m[nt*iv],0.0,0.0,nt,t,vec);
	  for (it=0;it<nt;it++) d[ip*nt+it]+=vec[it];
	}
      }
      else{
	//fprintf(stderr,"it=%d,ip=%d,iv=%d\n",it,ip,iv);
	dt2=dt*cost;
	for (it=0;it < nt;it++) t2[it]=it*dt2;
	nt2=(int) (nt*cost);
	ints8r(nt2,dt,0,&d[nt*ip],0.0,0.0,nt,t2,vec);
	for (it=0;it<nt;it++) m[iv*nt+it]+=vec[it];
      } 
    }  
  }
  TRACE;
  free1float(vec);
  free1float(t2);
  TRACE;
  return;
}






































