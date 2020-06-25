/* Copyright (c) University of British Columbia, 2000.*/
/* All rights reserved.                       */

/* SUREMOVE:  $Date: Julyh 2000  */

#include "dan.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUTEST3D - simple tests with routines for 3 dimensions              ", 
  " sutests3d < stdin > stdout [optional parameters]          		",
  " 									",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
 **************** end self doc ***********************************/

segy tr;
void zero_array(float ***d, int n3, int n2, int n1);
void test3dto1d(float *d, int n1, int n2, int n3, int i1, int i2, int i3);
int main(int argc, char **argv)
{

  register int it;
  int ih=0, iv=0;
  int nt, nh, nv, ntr;
  float ***d;
  int i2,i3,n1;

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);

   
  // Get parameters 
  //////////////////////////////////////////////////////////////////
 // Get info from first trace 
  if (!getparint("nv", &nv))  nv = 2;
  if (!getparint("i3", &i3))  i3 = 0;

  if (!getparint("i2", &i2))  i2 = 10;
  if (!getparint("n1", &n1))  n1 = 20;

  if (!gettr(&tr)) err("can't get first trace");

  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  nt = (int) tr.ns;
  nh = (int) tr.ntr;



  d=ealloc3float(nt,nh,nv); 
  fprintf(stderr,"nt=%d,nh=%d,nv=%d\n",nt,nh,nv);
  /////////////////////////////////////////////////////////////////  
  for (iv=0;iv<nv;iv++){
    ih=0;
    if (iv < nv) erewind(stdin);
    while (gettr(&tr)){   // Loop over traces
      fprintf(stderr,"gathering data for ih=%d, iv=%d \n", ih, iv);
      for (it=0;it<nt;it++) d[iv][ih][it]=tr.data[it];
      ih++;
    }

    fprintf(stderr,"number of traces = %d\n", (ih)*(iv+1));
  }
  
  
  fprintf(stderr,"nt=%d,nh=%d,nv=%d\n",nt,nh,nv);

  ntr=nh*nv;
 


  int i1=25;
  

  for (i1=1;i1<n1;i1++){
    fprintf(stderr,"d[%d][%d][%d]=%f\n",i3,i2,i1,d[i3][i2][i1]);
    test3dto1d(*d[0],nt,nh,nv,i1,i2,i3);
  }
  //zero_array(d,nv,nh,nt); 

  for (iv=0;iv<nv;iv++)
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) tr.data[it]=d[iv][ih][it];
      tr.ntr=ntr;
      puttr(&tr);
      
    }


  free3float(d);
  
  return EXIT_SUCCESS;
}


void test3dto1d(float *d, int n1, int n2, int n3, int i1, int i2, int i3)
{
 
  fprintf(stderr,"d[%d][%d][%d]=%f\n",i3,i2,i1,d[i3 *n2*n1 + i2 * n1 + i1]);
}

















