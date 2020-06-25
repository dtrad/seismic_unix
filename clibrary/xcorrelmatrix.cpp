#include "su.h"
#include "pwdf.h"

void xcorrelmatrix(float **A, float **B, float **C, int n1, int n2)
{
  /* Correlation of two matrices in the 2D frequency domain */
  float **AA, **BB, **CC; 
  complex **AAkf, **BBkf, **CCkf;
  int i1,i2;
  int nk, nw;

  
  int n1h=(int) n1/2;
  int n2h=(int) n2/2;

  /* Pad with zeros around. In matlab would be
     datain1z=[ zeros(n1h,2*n2); zeros(n1,n2h) datain1  zeros(n1,n2h); zeros(n1h,2*n2)];
  */

  for (i2=0;i2<2*n2;i2++){     // zeros(n1h,2*n2)
    for (i1=0;i1<n1h;i1++){               
      BB[i2][i1]=0;
      CC[i2][i1]=0;
    }
  }

  for (i2=0;i2<n2h;i2++){     // zeros(n1,n2h)
    for (i1=0;i1<n1;i1++){
      BB[i2][i1+n1h]=0;
      CC[i2][i1+n1h]=0;
    }
  }     
    
  for (i2=0;i2<n2;i2++){          //datain1
    for (i1=0;i1<n1;i1++){              
      BB[i2+n2h][i1+n1h]=B[i2][i1];
      CC[i2+nh2][i1+n1h]=C[i2][i1];
    }
  }     
    
  for (i2=0;i2<n2h;i2++){     // zeros(n1,n2h)
    for (i1=0;i1<n1;i1++){
      BB[i2][i1+n1h]=0;
      CC[i2][i1+n1h]=0;
    }
  }     
  
  for (i2=0;i2<2*n2;i2++){     // zeros(n1h,2*n2)
    for (i1=0;i1<n1h;i1++){               
      BB[i2][i1+n1]=0;
      CC[i2][i1+n1]=0;
    }
  }

  fksize(n1,n2,&nw,&nk);
  fft2_xt2kf(BB, BBkf,n1,n2,nw,nk);
  fft2_xt2kf(CC, CCkf,n1,n2,nw,nk);
  Atimes_conjB_elem(AAkf,BBkf,,CCkf,n1,n2);

  
  


  free2float(AA);
  free2float(BB);
  free2float(CC);
  free2complex(AAkf);
  free2complex(BBkf);
  free2complex(CCkf);


}
