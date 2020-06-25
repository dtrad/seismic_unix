void radonopi(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau,ktime,itime;
  unsigned int j;
  float *tau,*mint,*dtemp,time,hxh,pxhxh,tt;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((tau=alloc1float(nh*nt))==NULL)
    err("cannot allocate memory for tau \n");
  if ((mint=alloc1float(nt))==NULL)
    err("cannot allocate memory for mint \n");
  if ((mtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (iq=0;iq<nq;iq++){
    pxhxh=hxh*q[iq];
    iqxnt=iq*nt;
    for (ktime=0,ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;    
      for (itime=0;itime<nt;itime++)
	if ((tt=t[itime]*t[itime]-pxhxh) >= 0.0){
          ktime+=1;
	  tau[ktime]=sqrt(tt);
	}
    }
    intlin(ktime,tau,d,d[0],d[ktime-1],nt,t,mint);
    for (itau=0;itau<nt;itau++){
      k=iqxnt+itau;        
      m[k]=mint[itau];            
    }
  }
  free1float(dtemp);
  free1float(tau);
  free1float(mint);
  return;
}
