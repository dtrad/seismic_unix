/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUHRRT:  $Date: March 1999  */
#define NNX 228
#define NT 2048
#include "su.h"
#include "segy.h"
#include "header.h"
#include "clibrarytd.h"

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
int nt,nh,nq,nx,ny,rtmethod,itercg;
float dt,dh,dq;
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *headerfp;		/* fp for header storage file		*/
int verbose;

int main(int argc, char **argv)
{
	FILE *offsetfilep;
	int i,it,ih,iq, nn,flag;
	float *d, *m, temp; 
        float *q, *t, *h, *tempos, qmin, qmax, t0=0.;
	extern int nt,nh,nq,rtmethod;
	extern float dt,dh,dq;
        int method;
	cwp_String offsetfile=""; /* file containing positions */
	/// Velocity Trend
	float *tvel;
	float *vel;
	float *velint;
	float **vgrid;
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
	if (!getparint("method", &method))  method =0; 
	if (!getparint("rtmethod", &rtmethod))  rtmethod =3; /*HRT default*/
	if (!getparint("smooth", &smooth))  smooth =0;     
	if (!getparint("verbose", &verbose))  verbose =0;
	
	if (getparstring("offsetfile",&offsetfile)){
	  if ((tempos=ealloc1float(500))==NULL)
	   fprintf(stderr,"***Space for temppos could not be allocated\n");
          
          flag=0;  
	  fprintf(stderr,"New offset given by %s\n",offsetfile);  
	  if((offsetfilep=fopen(offsetfile,"r"))==NULL)
		err("cannot open offset file=%s\n",offsetfile);
	  ih=0;
	  do{
	     nn=fscanf(offsetfilep,"%f",&temp); 
             tempos[ih]=temp;
             ih++;
	  }while(nn==1);
	  nh=ih-1;
          fprintf(stderr,"nh= %d\n",nh);
	  fclose(offsetfilep);
	}
        else {
	  flag=1;
          if (!getparint("nh", &nh))  
	    err("I need either offsetfile or original nh");
        }

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
	if (!tr.ntr){ 
             nq=NNX; 
             fprintf(stderr,"ntr header field is not set");
             fprintf(stderr,"nq set to %d\n",nq);
        } 
	else nq=(int) tr.ntr;
	
        nt     = (int) tr.ns;
	dt   = ((float) tr.dt)/1000000.0;
        if ((!nt)||(!nh)||(!nq)) err("Error reading nt,nh or nq");

	
	// Allocate memory for data and model
  
	if ((d=ealloc1float(nh*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for d could not be allocated\n");
	
	if ((m=ealloc1float(nq*nt))==NULL)
	  fprintf(stderr,"***Sorry, space for m could not be allocated\n");
	
	if ((q=ealloc1float(nq))==NULL)
	  fprintf(stderr,"***Sorry, space for q could not be allocated\n");
  
	if ((h=ealloc1float(nh))==NULL)
	  fprintf(stderr,"***Sorry, space for h could not be allocated\n");

	if ((t=ealloc1float(nt))==NULL)
	  fprintf(stderr,"***Sorry, space for t could not be allocated\n");    

	if ((velint=ealloc1float(nt))==NULL)
	  fprintf(stderr,"*Sorry, space for velint could not be allocated\n"); 

	if (method==8){
	  vgrid=ealloc2float(nt,nq);
	  for (iq=0;iq<nq;iq++) for(it=0;it<nt;it++)
	    vgrid[iq][it]=velint[it]*velint[it]+q[iq]*q[iq]+2*velint[it]*q[iq];
	} 


     
	headerfp = etmpfile();
	if (verbose) warn("using tmpfile() call");  
	
	iq=0;
       	/* Loop over traces */
	do {
       		register int i;
		efwrite(&tr,HDRBYTES,1,headerfp);    
		q[iq]= (float) tr.f2;
                if ((flag==1)&&(iq<nh)) h[iq]=(float) tr.offset; 					
		for (i=0;i<nt;i++){
		        m[i+nt*iq]=(float) tr.data[i];
		}
		iq++;     
	} while (gettr(&tr));
	erewind(headerfp);   
        nq=iq;
	fprintf(stderr,"nq=%d\n",nq);


        if (flag==0) {
	  for (ih=0;ih<nh;ih++) h[ih]=tempos[ih];
          free1float(tempos);
        }
	for (i=0;i<nt;i++) t[i]=t0+i*dt;
	/* Create axis for velocities */
	intlin(ntvel,tvel,vel,vel[0],vel[nvel-1],nt,t,velint);

	if (verbose) 
	  for (it=0;it<nt;it++) 
	    fprintf(stderr,"velint[%d]=%f,t[%d]=%f\n",it,velint[it],it,t[it]);

	if (method==8){
	  vgrid=ealloc2float(nt,nq);
	  for (iq=0;iq<nq;iq++) for(it=0;it<nt;it++)
	    vgrid[iq][it]=velint[it]*velint[it]+q[iq]*q[iq]+2*velint[it]*q[iq];
	} 

        nx=nq*nt;
        ny=nh*nt;
	//////////////////////////////////////////////////////////////////////
	
	void (*radon)(float *,float *,float *,float *,float *,int,int,int,int);
	void (*radon2)(float *,float *,float *,float *,float *,float *,int,int,int,int);    
	void (*radon3)(float *,float *,float *,float *,float *,float **,int,int,int,int);

	if (rtmethod==3){
	  radon=radonhyp;
	  radon2=radonhyp;
	  radon3=radonhyp;
	}
	else if (rtmethod==2) radon=radonpar;
	else if (rtmethod==1) radon=radonlin;


	if (method==7) radon2(m,t,h,q,d,velint,0,nt,nh,nq);
	else if (method==8) radon3(m,t,h,q,d,vgrid,0,nt,nh,nq);
	else radon(m,t,h,q,d,0,nt,nh,nq);

	if (smooth) smoothing(d,nt,nh,nl,nr,flags);
	warn("Inside suradtdi\n");

	ih=0;
	do{
	  if (ih<nq)  efread(&tr,HDRBYTES,1,headerfp);
	  tr.offset=(int) h[ih];
       	  tr.ntr=nh;
          tr.f2=0;
	  for (i=0;i<nt;i++)
	    tr.data[i]=d[i+ih*nt];
	  puttr(&tr);
	  ih++;
	 } while(ih<nh);
       
        if (verbose) fprintf(stderr,"nt%d, nh %d \n ",nt,nh);

	if (method==8) free2float(vgrid);
	free1float(velint);   
	free1float(t);
	free1float(h);
	free1float(q);
	free1float(m);
	free1float(d);
	free1float(vel);
	free1float(tvel);
	efclose(headerfp);
           
	return EXIT_SUCCESS;
}











