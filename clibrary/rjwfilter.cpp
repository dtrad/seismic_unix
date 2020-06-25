#include "su.h"
#include "segy.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
#include "Complex.h"
/* Root j omega filter for Kirchoff 2D
The first routine is for a single trace
The second is for a two dimensional array 
 */
void rjwfilter(float *d,int nt, float dt)
{
  int it,  iw;          /* counters for time, and freq          */
  register float *rt;	/* real trace				*/
  register complex *ct; /* complex transformed trace		*/
  int nfft;		/* transform length			*/
  int nf;		/* number of frequencies		*/
  float dw;		/* output sample interval 2PIdf 	*/ 
  complex j;            /* Imaginary unit                       */
  j.r=0;j.i=1;
  const double  pi=acos(-1.);
  float onfft;          /* 1/nfft                               */
  float	const2 = sqrt(2.0); 
  float amp;
  float temp;
 
  nfft = npfaro(nt, LOOKFAC * nt);
  if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
  nf = nfft/2 + 1;
  dw = 2.0*PI/(nfft*dt);
  onfft = 1.0/nfft;

  if ((rt = ealloc1float(nfft))==NULL)
    err("cannot allocate memory for rt\n");
  if ((ct = ealloc1complex(nf))==NULL)
    err("cannot allocate memory for ct\n"); 

  fprintf(stderr,"In rjwfilter nfft=%d, nf=%d nt=%d, dt=%f, dw=%f, onfft=%f MIN=%d \n",nfft,nf,nt,dt,dw,onfft,MIN(nfft,nt));   
  
  for (it=0;it<MIN(nt,nfft);it++)
    rt[it]=d[it];
  if (nt<nfft) for (it=nt;it<nfft;it++)  rt[it]=0;    
  pfarc(1, nfft, rt, ct);
  /*  phase shifts PI/4 or sqrt(d/dt) */
  /* sqrt(i)=(1+i)/sqrt(2) */    
  for (iw=0;iw<nf;iw++){
    amp=(sqrt(dw*iw)/nfft);
    temp = (ct[iw].r-ct[iw].i)*amp*const2;
    ct[iw].i = (ct[iw].r+ct[iw].i)*amp*const2;
    ct[iw].r = temp; 
  }
  
  pfacr(-1, nfft, ct, rt);
  for (it = 0; it < MIN(nfft,nt); it++) 
    d[it] = rt[it]*onfft;
  if (nfft < nt) for (it = nfft; it < nt; it++) d[it] = 0;      
  
  
  free1complex(ct);
  free1float(rt);
      
  return;
}
  

void rjwfilter(float **d,int nt,int nh, float dt)
{
  int it, ih, iw;       /* counters for time, offset and freq   */
  register float *rt;	/* real trace				*/
  register complex *ct;/* complex transformed trace		*/
  int nfft;		/* transform length			*/
  int nf;		/* number of frequencies		*/
  float dw;		/* output sample interval in Hz		*/ 
  complex j;            /* Imaginary unit                       */
  j.r=0;j.i=1;
  const double  pi=acos(-1.);
  float onfft;          /* 1/nfft                               */
  float	const2 = sqrt(2.0); 
  float amp;
  float temp;
  
  nfft = npfaro(nt, LOOKFAC * nt);
  if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
  nf = nfft/2 + 1;
  dw = 2.0*PI/(nfft*dt);
  onfft = 1.0/nfft;

  if ((rt = ealloc1float(nfft))==NULL)
    err("cannot allocate memory for rt\n");
  if ((ct = ealloc1complex(nf))==NULL)
    err("cannot allocate memory for ct\n"); 

  fprintf(stderr,"In rjwfilter nfft=%d, nf=%d nh=%d nt=%d, dt=%f, dw=%f, onfft=%f MIN=%d \n",nfft,nf,nh,nt,dt,dw,onfft,MIN(nfft,nt));   
  
  for (ih=0;ih<nh;ih++){
    for (it=0;it<MIN(nt,nfft);it++)
      rt[it]=d[it][ih];
    if (nt<nfft) for (it=nt;it<nfft;it++)  rt[it]=0;    
    pfarc(1, nfft, rt, ct);
    for (iw=0;iw<nf;iw++){
      /* phase shifts PI/4 	*/
      amp=(sqrt(dw*iw)/nfft);
      temp = (ct[iw].r-ct[iw].i)*amp*const2;
      ct[iw].i = (ct[iw].r+ct[iw].i)*amp*const2;
      ct[iw].r = temp; 
    }
    pfacr(-1, nfft, ct, rt);
    for (it = 0; it < MIN(nfft,nt); it++) 
      d[it][ih] = rt[it];
      if (nfft < nt) for (it = nfft; it < nt; it++) d[it][ih] = 0;      
  }  
  
  free1complex(ct);
  free1float(rt);
      
  return;
}
  















