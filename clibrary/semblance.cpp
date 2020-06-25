#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void semblance( float *m, float *t, float *h, float *q, float *d, 
		int nt, int nh, int nq, float dt) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,*den,*suma;
  int nx, ny;
  nx=nt*nq;
  ny=nt*nh;
  suma=(float*) malloc((size_t) (nx+1)*sizeof(float));
  den=(float*) malloc((size_t) (nx+1)*sizeof(float));
  if (!den) printf("allocation failure in vector, semblance ");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq;iq++){
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	time=sqrt(pow(t[itau],2)+pow(h[ih],2)*q[iq]);
	it=time/dt;
	j=ih*nt+(int) floor(it);
	if ((it!=nt)&&(j<ny)&&(k<nx)){ 
	  m[k]=m[k]+d[j];
	  den[k]=den[k]+d[j]*d[j]; /**sincin(it-floor(it));*/
	}
      }
      
    }
  }
  int nwind=10;
  for (i=0;i<nx;i++) suma[i]=0;
  for (k=nwind;k<nx-nwind;k++)
    for (i=k-nwind;i<k+nwind;i++) suma[i]+=(m[i]*m[i]);
  for (k=0;k<nx;k++)  m[k]=suma[k];
  
  for (i=0;i<nx;i++) suma[i]=0;
  for (k=nwind;k<nx-nwind;k++) 
    for (i=k-nwind;i<k+nwind;i++) suma[i]+=(den[i]);
  for (k=0;k<nx;k++)  den[k]=suma[k];
  
  for (k=0;k<nx;k++) if (den[k]>1e-10) m[k]=(m[k])/(nh*den[k]); else m[k]==0;
  free(den);
  free(suma);
  return;
}






