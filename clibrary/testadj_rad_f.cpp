#include <math.h>
#include "su.h"
#include "clibrary.h"



float testadj_rad_f(complex **L, complex **LH)
{
  complex *dr1;
  complex *mr1;
  complex *dr2;
  complex *mr2;
  float dp1;
  float dp2;
  int it;
  extern int nh;
  extern int nq;
  extern int nt;
  


  if ((dr1=alloc1complex(nh))==NULL)
    fprintf(stderr,"***Sorry, space for dr1 could not be allocated\n");
  if ((mr1=alloc1complex(nq))==NULL)
    fprintf(stderr,"***Sorry, space for mr1 could not be allocated\n");
  if ((dr2=alloc1complex(nh))==NULL)
    fprintf(stderr,"***Sorry, space for dr2 could not be allocated\n");
  if ((mr2=alloc1complex(nq))==NULL)
    fprintf(stderr,"***Sorry, space for mr2 could not be allocated\n");
 
  for (it=0;it<nh;it++){
    dr1[it].r=frannor();
    dr1[it].i=frannor();
  }
  for (it=0;it<nq;it++){
    mr1[it].r=frannor();
    mr2[it].i=frannor();
  }

  Atimesx(dr2,L,mr1,nh,nq);
  Atimesx(mr2,LH,dr1,nq,nh);

  dp1=rcdot(nh,dr1,dr2);
  dp2=rcdot(nq,mr1,mr2);

  if (dp2!=0){
    fprintf(stderr,"testadj=%f\n",dp1/dp2);
    return(dp1/dp2);
  }else
    return(0.);
  
  free1complex(mr2);
  free1complex(dr2);
  free1complex(mr1);
  free1complex(dr1);

}



