/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUREMOVE:  $Date: Julyh 2000  */

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUMUTEPAR- mute traces whose par is larger than given value         ", 
  " sumutepar < stdin > stdout [optional parameters]            	",
  "  cut=        threshold                                              ",
  " ntaper=      number of taper points                                 ",
  " left=1       =1 mutes less than (left part). =0 more than (right)   ", 
  "									",
  " Applies a mute starting at f2=cut with a taper of ntaper length     ",
  "                                   				        ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr;

int main(int argc, char **argv)
{
  int j;
  register int it;
  int nt;
  int nh;
  int ih;
  float sum;
  float *taper;
  int count;
  int ntaper;
  float dq=1;
  float cut;
  float dqtemp;
  int itaper;
  int left;
  float qtemp[2]={0,0};
  // Initialize 
  
  initargs(argc, argv);
  requestdoc(1);
  
  // Get parameters 
  //////////////////////////////////////////////////////////////////
  if (!getparint("ntaper",&ntaper)) ntaper=10;
  if (!getparfloat("cut",&cut)) cut=0;
  if (!getparint("left",&left)) left=1;
  
  taper=ealloc1float(ntaper);
  
  // Get info from first trace 
  if (!gettr(&tr)) err("can't get first trace");
  if (!tr.ns) err("ns header field must be set");
  nt = (int) tr.ns;
  dqtemp=tr.f2;
  /////////////////////////////////////////////////////////////////  
  for (ih=0;ih<ntaper;ih++){
    taper[ih]=1-sin((ih+1)*PI/(2*ntaper));
    fprintf(stderr,"taper[%d]=%f\n",ih,taper[ih]);
  }
  /********************* Get two first traces to find dq ************/
  for (j=0;j<2;j++){
    gettr(&tr);
    qtemp[j]=tr.f2;
  }
  dq=qtemp[1]-qtemp[0];
  rewind(stdin);
  gettr(&tr);
  /******************************************************************/

  count=0;
  j=0;
  do {   // Loop over traces
    j++;
    if (left){
      if (tr.f2<cut) { 
	itaper=(int) ((cut-tr.f2)/dq+0.5);
	//fprintf(stderr,"dq=%e,itaper=%d\n",dq,itaper);
	if (itaper < ntaper) for (it=0;it<nt;it++) tr.data[it]*taper[itaper];
	else for (it=0;it<nt;it++) tr.data[it]*=0;
      }
    }
    else{
      if (tr.f2>cut) { 
	itaper=(int) ((tr.f2-cut)/dq+0.5);
	//fprintf(stderr,"dq=%e,itaper=%d\n",dq,itaper);
	if (itaper < ntaper) for (it=0;it<nt;it++) tr.data[it]*taper[itaper];
	else for (it=0;it<nt;it++) tr.data[it]*=0;
      }
    }  
    puttr(&tr);
  }while (gettr(&tr));
  //fprintf(stderr,"number of good traces = %d\n", count);
  free1float(taper);
  return EXIT_SUCCESS;
}




















