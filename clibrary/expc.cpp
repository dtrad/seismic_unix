#include "Dcomplex.h"
#include "Complex.h"
complex expc(complex x)
{
complex z,i;
i.r=0;i.i=1;
z=exp(x.r)*(cos(x.i)+i*sin(x.i));
return(z);
}
dcomplex expc(dcomplex x)
{
dcomplex z,i;
i.r=0;i.i=1;
z=exp(x.r)*(cos(x.i)+i*sin(x.i));
return(z);
}











