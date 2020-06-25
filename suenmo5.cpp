/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* suradtdi:  $Date: December 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "clibrary.h"
#include "header.h"
#include <time.h>
#include "inversion_par.h"

int rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv,  unsigned short **nullspace, float *critangle);

void rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle);

void rad_ellip_sp_old(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, unsigned short **nullspace, float *critangle);

void radonellip(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), float *t, float *h, float *q, float **vel,int nt, int nh, int nq);

float dot(int n, float *a, float *b);

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv);

void enmo1(float *m,float *t,float *p,float *v,float *d, float **vgrid,int adj, int nt, int np, int nv);

void rms2intvel(int n,float *t0, float *vs, float *v, float *h);

void myintlin (int nin, float xin[], float yin[], float yinl, float yinr, 
	int nout, float xout[], float yout[]);

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
  float *p, *t, *v, vmin, vmax, t0=0.;
  float **Wm; // Model Weights
  float **Wd; // residual Weights
  int nt,nh,np,nv;
  float dt,dv,dp;
  int model;
  int itercgfin;
  int testadj; // !=0 test adjoint
  int taper;
  int method; // =0 sparse; =1 operator
  int nsparse;
  unsigned int **index; // the head of the monster
  size_t sizeus=sizeof(unsigned int);

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

  // smoothing
  int smooth; // =1 apply time smoothing
  int nl=3;  //  npoints left hand side
  int nr=3;  //  npoints left hand side
  int flags=2;  // 1 rectangular, 2 triangular
  
  //////////////////////////////////////////////
  inv_par inv;
  inv.restart=1;
  // Define function variables for method 0 (sparse) and 1 (operator)
  void (*radonellip0) (float *m, float *d, unsigned int **index, 
		       int adj, int nt, int nh, int nq, int nsparse);
  void (*radonellip1) (float *m,float *t,float *p,float *v,float *d,
		       float **vgrid,int adj, int nt, int np, int nv);

  radonellip0=radonellip;
  radonellip1=enmo1;
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);
  
  /* Get parameters */
  if (!getparint("model", &model))  model =1;     
  if (!getparint("smooth", &smooth))  smooth =0;     
  if (!getparint("nv", &nv))  nv =50;
  if (!getparfloat("dv", &dv))  dv =50;
  if (!getparfloat("vmin", &vmin))  vmin =1000;
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

  critang=ealloc1float(nvel);
  critangle=ealloc1float(nt);
  vinterval=ealloc1float(nvel);
  hinterval=ealloc1float(nvel);
  rms2intvel(nvel,tvel,vel,vinterval,hinterval);
  critang[0]=0.99;
  for (itv=0;itv<nvel;itv++){ 
    critang[itv]=vel[itv]/vinterval[itv];
    fprintf(stderr,"vinterval[%d]=%f,hinterval[%d]=%f,critang[%d]=%f\n",
	    itv,vinterval[itv],itv,hinterval[itv],itv,critang[itv]);
  }
  
  // Allocate memory for data and model
  
  if ((d=ealloc2float(nt,np))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");

  if ((daux=ealloc2float(nt,np))==NULL)
    err("Cannot allocate daux\n");

  if ((nullspace=(unsigned short**) alloc2(nt,np,sizeof(unsigned short)))==NULL)
    fprintf(stderr,"Cannot allocate nullspace\n");
  
  if ((m=ealloc2float(nt,nv))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
  
  if ((Wm=ealloc2float(nt,nv))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");
  
  if ((Wd=ealloc2float(nt,np))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");
  
  if ((p=ealloc1float(np))==NULL)
    fprintf(stderr,"***Sorry, space for p could not be allocated\n");
  
  if ((v=ealloc1float(nv+1))==NULL)
    fprintf(stderr,"***Sorry, space for v could not be allocated\n");
  
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");    
  
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"*Sorry, space for velint could not be allocated\n"); 
  
  if ((vgrid=ealloc2float(nt,nv+1))==NULL)
    fprintf(stderr,"*Sorry, space for vgrid could not be allocated\n"); 

  
  // Because we want to use same struct array for data and model 
  // the maximun number of traces between them will be taken.
  
  headerfp = etmpfile();
  if (verbose) warn("using tmpfile() call");  
  
  /* Loop over traces */
  ip=0;
  do {
    p[ip]= (float) tr.f2;
    efwrite(&tr,HDRBYTES,1,headerfp);		
    for (it=0;it<nt;it++){
      d[ip][it]=(float) tr.data[it];
    }
    ip++;
    if (ip > np) err("Number of traces > %d\n",np);     
  } while (gettr(&tr));
  
  erewind(headerfp); 
  np=ip;

  nx=nv*nt;
  ny=np*nt;
  
  if (verbose) fprintf(stderr,"np=%d,nt%d,dt %f\n",np,nt,dt);
  int nvh=(int) (nv/2);
  v[nvh]=0;
  for (iv=nvh;iv<nv;iv++) v[iv+1]=v[iv]*(1+dv)+vmin;
  for (iv=nvh-1;iv>=0;iv--) v[iv]=v[iv+1]*(1+dv)-vmin;
  //  for (iv=0;iv<nv;iv++) fprintf(stderr,"v[%d]=%f,vmin=%f\n",iv,v[iv],vmin);
  
  for (it=0;it<nt;it++) t[it]=t0+it*dt;
  /* Create axis for velocities */
  intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);
  //for (it=0;it<nt;it++) fprintf(stderr,"velint[%d]=%f,t[%d]=%f\n",
  //			      it,velint[it],it,t[it]);
  //    Interpolate the sin of the critical angle
  for (it=0;it<nvel;it++)  fprintf(stderr,"tvel[%d]=%f,critang[%d]=%f\n",it,tvel[it],it,critang[it]);

  if(1){
    it=0;
    for(itv=0;itv<nvel;itv++){
      while((t[it]<tvel[itv])&&(itv<nvel)){
	critangle[it]=critang[itv];
	it++;
      }
    }
  }  
     
  //intlin(nvel,tvel,critang,critang[0],critang[nvel-1],nt,t,critangle);
  free1float(critang);
  for (it=0;it<nt;it++)  fprintf(stderr,"critangle[%d]=%f\n",it,critangle[it]);  

  free1float(vinterval);
  free1float(hinterval);

  float vaux;
  for (iv=0;iv<nv;iv++)
    for (it=0;it<nt;it++){
      vaux=v[iv]*(1.+2*t[it]);
      vgrid[iv][it]=velint[it]*velint[it]+vaux*vaux+2*velint[it]*vaux;
    }      
    
  
  // Assign data weights
  for (ip=0;ip<np;ip++) for (it=0;it<nt;it++) Wd[ip][it]=1;
  // Test: For data with gaps could be useful to taper the p=0 domain
  // where the gaps produce zeros data.

  if (taper){ 
    float par=1.5;
    for (ip=(np/2-10);ip<np/2+10;ip++)
      for (it=0;it<nt;it++){
	Wd[ip][it]=0.1*(t[it]+par)*(t[it]+par);
	//fprintf(stderr,"Wd[%d][%d]=%f\n",ip,it,Wd[ip*nt+it]);
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
  if (method==0){ 
    // Count number of non zero elements for the operator index
    if (flagsparse) nsparse=rad_ellip_sp(t,p,v,vgrid,nt,np,nv,nullspace,critangle);
    else nsparse=nt*nv*np;
    fprintf(stderr,"Nonzero elements for index=%d\n",nsparse);
    // Allocate the big monster
    if ((index=(unsigned int **) alloc2(nsparse+10,2,sizeof(unsigned int)))==NULL)
      err("Cannot allocate the monster index\n");
    if (flagsparse) rad_ellip_sp(t,p,v,vgrid,nt,np,nv,index,nsparse,critangle);
    else  rad_ellip_sp_old(t,p,v,vgrid,nt,np,nv,index,nsparse,nullspace,critangle);
    if (testadj) testadjop(radonellip0,index,nt,np,nv,nsparse);
    radonellip0(m[0],d[0],index,1,nt,np,nv,nsparse);
    for (iv=0;iv<nv;iv++) for (it=0;it<nt;it++) Wm[iv][it]=1;
    modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
    int itemp=inv.itercg;inv.itercg=5;
    wpcgnr(radonellip0,nt,np,nv,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv);
    inv.itercg=itemp;
    
    for (j=1;j<=inv.iter_end;j++){
      // norm==1 ==> L1 , ==0  Cauchy else l2
      modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
      if (j==inv.iter_end) inv.itercg=itercgfin;      
      inv.J[j]=wpcgnr(radonellip0,nt,np,nv,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv);
    }

    radonellip0(m[0],daux[0],index,0,nt,np,nv,nsparse);   
    free2((void **) index);
    //radonellip1(m[0],t,p,v,daux[0],vgrid,0,nt,np,nv);
    for (ip=0;ip<np;ip++) 
      for (it=0;it<nt;it++)
	if (nullspace[ip][it]) d[ip][it]=daux[ip][it];
  }
  else if (method==1){
    enmo1(m[0],t,p,v,d[0],vgrid,1,nt,np,nv);
    if (testadj) testadjop(radonellip1,t,p,v,vgrid,nt,np,nv);
    for (j=1;j<=inv.iter_end;j++){
 	modelweight(m[0],nx,inv.norm,inv.eps1,Wm[0]);
	inv.J[j]=wpcgnr(radonellip1,nt,np,nv,t,p,v,m[0],d[0],Wd[0],Wm[0],vgrid,inv);
	radonellip1(m[0],t,p,v,d[0],vgrid,0,nt,np,nv);  
    }	 
  } 
  if (model){
    if (smooth) smoothing(m[0],nt,nv,nl,nr,flags);
    FILE *myfilep;
    if((myfilep=fopen(modelfile,"w"))==NULL)
      err("cannot open file=%s\n",modelfile);        
    iv=0;
    do{
      tr.f2=(float) sqrt(vgrid[iv][0]);
      //tr.f2=(float) 1./(vgrid[iv][0]);
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



void radonellip(float *m, float *d, unsigned int **index, int adj, int nt, int np, int nv, int nsparse)
{
  unsigned long j;
  unsigned int ny=np*nt;
  unsigned int nx=nv*nt;
  unsigned int it;

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






























