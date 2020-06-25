#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int t2t0(float t,float h,float v, float dt);

int main(void)
{
  float dt=0.004;
  float t,h,v;
  int it0;
  puts("t,h,v?\n");
  scanf("%f %f %f",&t,&h,&v);
  printf("t=%f,h=%f,v=%f\n",t,h,v);
  it0=t2t0(t,h,v,dt);
  printf("it=%f,it0=%d,t=%f\n",t/dt+0.5,it0,it0*dt);
  return 0;
}
