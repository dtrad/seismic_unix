#include "su.h"


void taper (int lxtaper, int lbtaper, 
	    int nx, int ix, int nt, float *trace, int inv)
/*****************************************************************************
Taper traces near left and right sides of trace window
******************************************************************************
Input:
lxtaper		length (in traces) of side taper
lbtaper		length (in samples) of bottom taper
nx		number of traces in window
ix		index of this trace (0 <= ix <= nx-1)
nt		number of time samples
trace		array[nt] containing trace

Output:
trace		array[nt] containing trace, tapered if within lxtaper of side
*****************************************************************************/
{
  int it;
  float xtaper;
  //inv=0;
  /* if near left side */
  if (ix<lxtaper) {
    xtaper = 0.54+0.46*cos(PI*(lxtaper-ix)/lxtaper);
    
    /* else if near right side */
  } else if (ix>=nx-lxtaper) {
    xtaper = 0.54+0.46*cos(PI*(lxtaper+ix+1-nx)/lxtaper);
    
    /* else x tapering is unnecessary */
  } else {
    xtaper = 1.0;
  }
  fprintf(stderr,"xtaper=%f\n",xtaper);
  
  /* if x tapering is necessary, apply it */
  if (xtaper!=1.0){
    if (!inv)
      for (it=0; it<nt; ++it)
	trace[it] *= xtaper;
    else
      if (xtaper>0.1)
	for (it=0; it<nt; ++it)
	  trace[it] /= (xtaper*1);
  }

  /* if requested, apply t tapering */
  for (it=MAX(0,nt-lbtaper); it<nt; ++it)
    trace[it] *= (0.54+0.46*cos(PI*(lbtaper+it+1-nt)/lbtaper));
}

