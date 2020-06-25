/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUREMOVE:  $Date: Julyh 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUSETNTR- set the sukeyword ntr to the right number of traces       ",
  "                                      				", 
  " susetntr< stdin > stdout [optional parameters]          		",
  "                                      				",
  " It can be used in the middle of a pipe because it uses temporary    ",
  " files rather than stdin or stdout.                            	",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
 **************** end self doc ***********************************/

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
segy tr;
static void closefiles(void);

int main(int argc, char **argv)
{
  int j;
  int ntr;
  int nt;			/* number of time samples */
  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/
  int verbose;		/* flag for echoing info		*/
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);

  /* Look for user-supplied tmpdir */
  if (!getparstring("tmpdir",&tmpdir) &&
      !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
  if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
    err("you can't write in %s (or it doesn't exist)", tmpdir);
  
  
  /* store traces and headers in tempfiles while getting a count */
  if (STREQ(tmpdir,"")) {
    tracefp = etmpfile();
    headerfp = etmpfile();
    if (verbose) warn("using tmpfile() call");
  } else { /* user-supplied tmpdir */
    char directory[BUFSIZ];
    strcpy(directory, tmpdir);
    strcpy(tracefile, temporary_filename(directory));
    strcpy(headerfile, temporary_filename(directory));
    /* Trap signals so can remove temp files */
    signal(SIGINT,  (void (*) (int)) closefiles);
    signal(SIGQUIT, (void (*) (int)) closefiles);
    signal(SIGHUP,  (void (*) (int)) closefiles);
    signal(SIGTERM, (void (*) (int)) closefiles);
    tracefp = efopen(tracefile, "w+");
    headerfp = efopen(headerfile, "w+");
    istmpdir=cwp_true;		
    if (verbose) warn("putting temporary files in %s", directory);
  }
  
  // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");
  nt = tr.ns;

  j = 0;
  do {
    ++j;
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data, FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);



  ntr=j;
  /* restore header fields and write output */
  for (j=0; j<ntr; ++j) {
    efread(&tr,HDRBYTES,1,headerfp);
    efread(tr.data,FSIZE,nt,tracefp);
    tr.ntr=ntr;
    puttr(&tr);
  }
  efclose(headerfp);
  if (istmpdir) eremove(headerfile);
  if (istmpdir) eremove(tracefile);

  return EXIT_SUCCESS;
}



/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}

















