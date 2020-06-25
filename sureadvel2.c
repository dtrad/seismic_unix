/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* sureadvel: $Revision: 1.3 $ ; $Date: 2000/01/06 17:14:38 $	*/
/* Modified version to read promax velocity files */

#include "par.h"
#include "su.h"
#include "header.h"
#include "segy.h"
/*********************** self documentation ******************************/
char *sdoc[] = {
" 									",
" sureadvel2 - convert velocity sufile to par file format 		",
" 									",
" sureadvel2 < stdin > stdout vfile=	      				",
" 	vfile=v[nx][nz] to use with rayt2d	     			",
" Optional parameters:							",
" 	string1=\"par1=\"	first par string			",
" 	string2=\"par2=\"	second par string			",
"       factor=1                only takes every FACTOR traces          ",
" 	nz=300                  number of vertical cells	      	",
" 	nx=300                  number of horizontal cells	      	",
" 	dz=10                   vertical step size for vfile		",
"       step=20                 step*dt is the vertical time step       ",
" 									",
" This is a tool to convert velocities written in a sufile to parameter ",
" vectors in the form expected by getpar.  For example, if the input	",
" file is a sufile 							",
" then									",
"	sureadvel < input > output string1=tnmo string2=vnmo		",
" yields:								",
"	tnmo=t0,t1,...							",
"	vnmo=v0,v1,...							",
" 									",
NULL};
/**************** end self doc *******************************************/

/* Credits:
 *	CWP: Jack
        Modified Daniel Trad - UBC
 */


/* Caveat: A more general tool allowing n1 strings would be desirable. */
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
 

int
main(int argc, char **argv)
{
	int it, i1, i2, nt, step, ntr, ntr2, n1=0, n2 = 0;
	float x1, x2, dt, time;
	char *string1;
	char *string2;
	char *vfile;
	FILE *datafp1;
        FILE *datafp2;
	FILE *pvfile;           /* pointer to vfile                     */ 
	segy tr;
	int factor,itr; 
	float **vz;
	float *z;
	float dz;
	float dx;
	int nx,nz;
	unsigned ix,iz;
       	int jump;
	/* The following variables are required to calculate dx for vfile */
	float cdpspace;
	float cdpmin;
	float cdpmax;

	/* Hook up getpar */
	initargs(argc, argv);
	requestdoc(1);


	/* Get parameters and set up tmpfile */
	if (!getparstring("string1", &string1))	string1 = "tnmo";
	if (!getparstring("string2", &string2))	string2 = "vnmo";
	if (!getparstring("vfile",&vfile)) vfile= "vfile";
	if (!getparint("step",&step)) step=20;
	if (!getparint("factor",&factor)) factor=1; 
	if (!getparfloat("cdpspace",&cdpspace)) cdpspace=1;
	if (!getparfloat("dz",&dz)) dz=10.;
	if (!getparint("nz",&nz)) nz=300;
	if (!getparint("nx",&nx)) nx=300;

	/* We open four temp files:
	   tracefp and headerfp for the input file
	   datafp1 and datafp2 for time and velocities
	*/

	datafp1 = etmpfile();
	datafp2 = etmpfile();
	tracefp = etmpfile();
	headerfp = etmpfile();

	// Get info from first trace 
  	if (!gettr(&tr)) err("can't read first trace");
	if (!tr.dt) err("dt header field must be set");
	//if (!tr.cdp) err("cdp header field must be set");
	
	dt=tr.dt/1.e6;
	ntr=0;
        ntr2=0;
	fprintf(stderr,"dt=%f\n",dt);
	
       

	// Extract cdp from tr and save data for later pass over time and vel 
	printf("cdp=%d",(int) tr.cdp);
	cdpmin=tr.cdp;
	do{
	  if (fmod((float) ntr2,(float) factor)==0){ 
	    ntr++;
	    it=1;
	    nt=tr.ns;	  
	    fprintf(stderr,"cdp=%d,ntr=%d\n",tr.cdp,ntr);
	    if (ntr!=1) printf(",%d",(int) tr.cdp);
	    efwrite(&tr,HDRBYTES,1,headerfp);
	    efwrite(tr.data,FSIZE, nt, tracefp);
	  }
	  ntr2++;
	} while (gettr(&tr));    
	putchar('\n');
	cdpmax=tr.cdp;
	rewind(headerfp);
	rewind(tracefp);
	dx=(cdpmax-cdpmin)/(nx-1)*cdpspace;
	// First reading of datafile to construct velocity vz[x][z]
		
	vz=ealloc2float(nz,nx);
	z=ealloc1float(nz);
	if (ntr < nx){ 
	  nx=ntr;
	  warn("ntr < nx ==> nx set to the number of traces\n");
	}
	jump=(int) (floor(ntr/nx+0.5));

	for (z[0]=0,iz=1;iz<nz;iz++) z[iz]=z[iz-1]+dz;
	fprintf(stderr,"dz=%f;nz=%d,z[%d]=%f,jump=%d\n",dz,nz,nz-1,z[nz-1],jump);
	fprintf(stderr,"dx=%f;nx=%d\n",dx,nx);
	for (ix=0,itr=0;itr<ntr;itr++){
	  efread(tr.data,FSIZE,nt,tracefp);
	  efread(&tr,HDRBYTES,1,headerfp);
       	  if ((fmod(tr.tracl,jump)==0)&&(ix<nx)){
	    fprintf(stderr,"tr.tracl=%d\t",tr.tracl);
	    vz[ix][0]=tr.data[0];
	    for (it=0,iz=1;iz<nz;iz++){
	      it+=(2*dz/vz[ix][iz-1])/dt;
	      if (it<nt) vz[ix][iz]=tr.data[it];
	      else vz[ix][iz]=tr.data[nt-1];
	      fprintf(stderr,"vz[%d][%d]=%f\n",ix,iz,vz[ix][iz]);
	    }
	    ix++;
	    fprintf(stderr,"ix=%d,nx=%d\n",ix,nx);
	  }
	}
	// nx=ix-1;
	pvfile=efopen(vfile,"w");
	if (0) 
	  for (ix=0;ix<nx;ix++) 
	    for (iz=0;iz<nz;iz++) 
	      fprintf(stderr,"vz[ix][iz]=%f\n",vz[ix][iz]);
 
	efwrite(vz[0],FSIZE,nx*nz,pvfile);
	efclose(pvfile);


	rewind(headerfp);
	rewind(tracefp);
	
	for (itr=0;itr<ntr;itr++){
	  rewind(datafp1);
	  rewind(datafp2);
	  n1=0;
	  n2=0;
	  it=1;
	  efread(tr.data, FSIZE, nt, tracefp);
	  efread(&tr,HDRBYTES,1,headerfp);	  
          nt=tr.ns;
	  time=it*dt;
	  efwrite(&time, FSIZE, 1, datafp1);n1++;
	  efwrite(&tr.data[it], FSIZE, 1, datafp2);n2++;  
	  while (it < nt ){
	    while ((tr.data[it]==tr.data[it-1])&&(it<nt)){
	      it++;
	      // fprintf(stderr,"it=%d\n",it);
	    }
	    time=it*dt;
	    
	    efwrite(&time, FSIZE, 1, datafp1);n1++;
	    efwrite(&tr.data[it], FSIZE, 1, datafp2);n2++;
	    it+=step;
	  }

	  /* Rewind and get the time */
	  rewind(datafp1);
	  efread(&x1, FSIZE, 1, datafp1);
	  printf("%s=%g", string1, x1);
	  for (i1 = 2; i1 < n1; i1++) {
	    efread(&x1, FSIZE, 1, datafp1);
	    printf(",%g", x1);
	  }
	  putchar('\n');	
  
	  /* Rewind and get the x2's */
	  rewind(datafp2);
	  efread(&x2, FSIZE, 1, datafp2);
	  printf("%s=%g", string2, x2);
	  for (i2 = 2; i2 < n2; i2++) {
	    efread(&x2, FSIZE, 1, datafp2);
	    printf(",%g", x2);
	  }
	  putchar('\n');	  
	}
	
	fprintf(stderr,"%d cdps\n",ntr);
	fprintf(stderr,"Parameters required to use vfile in rayt2d\n");
	fprintf(stderr,"dx=%f,nx=%d\n",dx,nx);
 	fprintf(stderr,"dz=%f,nz=%d\n",dz,nz);

	free2float(vz);

	return EXIT_SUCCESS;
}
























