/* Copyright (c) Colorado School of Mines, 1998.*/
/* All rights reserved.                       */

/* A2B: $Revision: 1.8 $ ; $Date: 1996/09/06 16:08:13 $	*/

#include "par.h"
#define ASCII_LINE 131072  /* 2^17 */

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" A2B - convert ascii floats to binary 				",
" 								",
" a2b <stdin >stdout outpar=/dev/tty 				",
" 								",
" Required parameters:						",
" 	none							",
" 								",
" Optional parameters:						",
" 	n1=2		floats per line in input file		",
" 								",
" 	outpar=/dev/tty	output parameter file, contains the	",
"			number of lines (n=)			",
" 								",
NULL};
/**************** end self doc ***********************************/

/* Credits:
 *	CWP: Jack K. Cohen, Dave Hale
 */


int
main(int argc, char **argv)
{
	char *outpar;		/* name of file holding output parfile	*/
	FILE *outparfp;		/* ... its file pointer			*/
	int n1;			/* number of floats per line		*/
	int n2 = 0;		/* number of lines in input file 	*/
	char buf[ASCII_LINE];	/* buffer for the ascii floats		*/
	float *x;		/* binary floats			*/


	/* Hook up getpar */
	initargs(argc, argv);
	requestdoc(1);


	/* Prevent floats from dumping on screen */
	switch(filestat(STDOUT)) {
	case BADFILETYPE:
		warn("stdout is illegal filetype");
		pagedoc();
	break;
	case TTY:
		warn("stdout can't be tty");
		pagedoc();
	break; 
	default: /* rest are OK */
	break;

	}


	/* Get parameters and do set up */
	if (!getparstring("outpar", &outpar))  outpar = "/dev/tty" ;
	outparfp = efopen(outpar, "w");

	if (!getparint("n1", &n1))  n1 = 2;
	x = ealloc1float(n1);


	/* Loop over data, converting to binary */
	while (gets(buf)) {
		char *p = buf;
		register int i1;

		for (i1 = 0; i1 < n1; ++i1) {
			while (*p == ' ' || *p == '\t')  ++p;
			if (1 != sscanf(p, "%f", x+i1))
				err("line #%d: scan failed on float #%d:\n%s",
					n2+1, i1+1, buf);
			while (*p != ' ' && *p != '\t' && *p != '\0')  ++p;
		}
		efwrite(x, FSIZE, n1, stdout);
		++n2;
	}


	/* Make par file */
	fprintf(outparfp, "n=%d\n", n2);


	return EXIT_SUCCESS;
}
