#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
float costfunc(complex *d, complex *dc, complex *u,int nh,int nq,
float sigma,int norma)
{
      int i;
      float  r2, power,cost;

      r2=0;
      for(i=0;i<nh;i++)
         r2=r2+abs((d[i]-dc[i]))*abs((d[i]-dc[i]));
      
      power=0;
      if (norma==10) 
          for(i=0;i<nq;i++)
              power=power+log(1+real(conjg(u[i])*u[i])/sigma);            
      else if(norma==1) 
          for(i=0;i<nq;i++)
              power=power+abs(u[i]);
          
      
      cost=r2+power;

      //     write(0,*) 'Jtot,Jdata,Jmodel',costfunc,r2,power
      
      
      return(cost);
}
      

      




