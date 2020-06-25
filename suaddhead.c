/* Copyright (c) Colorado School of Mines, 1998.*/
/* All rights reserved.                       */

/* SUADDHEAD: $Revision: 1.15 $ ; $Date: 1997/07/28 22:36:46 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
" 									",
" SUADDHEAD - put headers on bare traces and set the tracl and ns fields",
" 									",
" suaddhead <stdin >stdout ns= ftn=0					",
" 									",
" Required parameter:							",
" 	ns=the number of samples per trace				",
" 									",
" Optional parameter:							",
" 	ftn=0		Fortran flag					",
" 			0 = data written unformatted from C		",
" 			1 = data written unformatted from Fortran	",
" 									",
" Trace header fields set: ns, tracl					",
" Use sushw/suchw to set other needed fields.				",
" 									",
" Caution: An incorrect ns field will munge subsequent processing.	",
" Note:    n1 and nt are acceptable aliases for ns.			",
" 									",
" Example:								",
" suaddhead ns=1024 <bare_traces | sushw key=dt a=4000 >segy_traces	",
" 									",
" This command line adds headers with ns=1024 samples.  The second part	",
" of the pipe sets the trace header field dt to 4 ms.			",
" 									",
NULL};

/* Credits:
 *	SEP: Einar Kjartansson
 *	CWP: Jack K. Cohen
 *
 * Trace header fields set: ns, tracl
 */

/**************** end self doc *******************************************/

segy tr;

int
main(int argc, char **argv)
{
	int ns;			/* number of samples			*/
	int ftn;		/* ftn=1 for Fortran			*/
	char junk[ISIZE];	/* to discard ftn junk  		*/
	cwp_Bool isreading=cwp_true;	/* true/false flag for while	*/


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get parameters */
	if (!getparint("n1", &ns)
	 && !getparint("nt", &ns)
	 && !getparint("ns", &ns))  err("must specify ns=");
	tr.ns = ns;
	if (!getparint("ftn", &ftn))	ftn = 0;
	if (ftn != 0 && ftn != 1)	err("ftn=%d must be 0 or 1", ftn);


	while (isreading==cwp_true) {
		static int tracl = 0;	/* one-based trace number */

		/* If Fortran data, read past the record size bytes */
		if (ftn) efread(junk, ISIZE, 1, stdin);

		/* Do read of data for the segy -- quit on <= 0 items */
		if (0 >= efread((char *) tr.data, FSIZE, ns, stdin)) break;
		tr.tracl = ++tracl;
		puttr(&tr);

		/* If Fortran data, read past the record size bytes */
		if (ftn) efread(junk, ISIZE, 1, stdin);
	}

	return EXIT_SUCCESS;
}
