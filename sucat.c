/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVCAT: $Revision: 1.18 $ ; $Date: 2011/11/16 23:09:52 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SUCAT -  append one data set  				",
" sucat data1 data2 >stdout					",
" 								",
" Required parameters:						",
"        none							",
"                                                               ",
" 								",
NULL};

/* Credits:
 */
/**************** end self doc ***********************************/

segy intrace1;

int main(int argc, char **argv){

	FILE *fp1 = NULL;	/* file pointer for first file		*/
	FILE *fp2 = NULL;	/* file pointer for second file		*/
	int itr = 0;	/* number of trace being processed	*/
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1); /* two file args required */
        checkpars();

	/* Open two files given as arguments for reading */
	fp1 = efopen(argv[1], "r");
        if (argc > 2) fp2 = efopen(argv[2], "r");
        else{
            fprintf(stderr,"only one input given endl\n");
        }
	/* Loop over the traces */
	while (fgettr(fp1, &intrace1)){
		puttr(&intrace1); 
		++itr;
	}
	if (fp2)  
	while (fgettr(fp2, &intrace1)){
		puttr(&intrace1); 
		++itr;
	}
	  

	return(CWP_Exit());
}
