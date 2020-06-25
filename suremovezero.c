/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKILL: $Revision: 1.18 $ ; $Date: 2011/11/17 00:03:38 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 								",
  " SUCOUNTZERO - count zero traces					",
  " 								",
  " sucount <stdin >stdout [optional parameters]			",
  " 								",
  " Optional parameters:						",
  "	key=gx	header name to select traces to kill	        ",
  "	a=2	header value identifying min header to kill     ",
  "     b=      header value identifying max header to kill     ",
  "     rem=0   0 keep 1 remove zero traces								",
  " 	min= 		first trace to kill (one-based)		",
  " 	count=1		number of traces to kill 		",
  " 								",
  " Notes:							",
  "	If min= is set it overrides selecting traces by header.	",
  " 								",
  NULL};

/* Credits:
 *	CWP: Chris Liner, Jack K. Cohen
 *	header-based trace selection: Florian Bleibinhaus
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


segy tr;

int
main(int argc, char **argv)
{
  cwp_String key;		/* trace header			*/
  cwp_String type;	/* type for trace header	*/
  int index;		/* index of trace header	*/
  Value val;		/* trace header value		*/
  double dval,a,b;		/* trace header value		*/
  register int itr;	/* trace counter		*/
  int nt = 0;		/* number of time samples	*/
  int rem = 0;          /* keep (0) or remove (1) zero traces */
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1);


  /* Get parameters */
  if (!getparstring("key", &key)) key = "gx";
  if (!getpardouble("a", &a)) a = 0;
  if (!getpardouble("b", &b)) b = 0;
  if (!getparint("rem",&rem)) rem = 1;
  checkpars();

  /* Get type and index value */
  type  = hdtype(key);
  index = getindex(key);
  int countalive=0;
  int countdead=0;
  double sum=1;
  while (gettr(&tr)) {
    nt = tr.ns;
    gethval(&tr, index, &val);
    dval = vtod(type, val);
    if (( dval>=a )&&(dval<=b)){
      sum=0;
      for (int i=0;i<nt;i++) sum+=fabs(tr.data[i]);
      if (sum) countalive++;
      else     countdead++;
      

    }
    if (sum) puttr(&tr);
  }
  fprintf(stderr,"alive = %d , dead = %d \n", countalive, countdead);


  return(CWP_Exit());
}
