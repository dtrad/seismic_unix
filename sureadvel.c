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
" sureadvel - convert velocity sufile to par file format 		",
" 									",
" sureadvel <stdin >stdout 						",
" 									",
" Optional parameters:							",
" 	string1=\"par1=\"	first par string			",
" 	string2=\"par2=\"	second par string			",
" 									",
" This is a tool to convert velocities written in a sufile to parameter ",
" vectors in the form expected by getpar.  For example, if the input	",
" file is a sufile 							",
" then									",
"	sureadvel <input >output string1=tnmo string2=vnmo		",
" yields:								",
"	tnmo=t0,t1,...							",
"	vnmo=v0,v1,...							",
" 									",
NULL};
/**************** end self doc *******************************************/

/* Credits:
 *	CWP: Jack
        Modified Daniell Trad - UBC
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
	FILE *datafp1;
        FILE *datafp2;
	segy tr;
	int factor,itr; 
	
	/* Hook up getpar */
	initargs(argc, argv);
	requestdoc(1);


	/* Get parameters and set up tmpfile */
	if (!getparstring("string1", &string1))	string1 = "tnmo";
	if (!getparstring("string2", &string2))	string2 = "vnmo";
	if (!getparint("step",&step)) step=20;
	if (!getparint("factor",&factor)) factor=1; 
 
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
	if (!tr.cdp) err("cdp header field must be set");
	
	dt=tr.dt/1.e6;
	ntr=0;
        ntr2=0;
	fprintf(stderr,"dt=%f\n",dt);

	// Extract cdp from tr and save data for later pass over time and vel 
	printf("cdp=%d",(int) tr.cdp);
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

	return EXIT_SUCCESS;
}











