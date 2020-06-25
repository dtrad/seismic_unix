#include <math.h>
void rms2intvel(int n,float *t0, float *vs, float *v, float *h)
{
  register int i;		/* counter				*/
  float t1, t2;		/* temporaries for one-way times	*/
  float v1, v2;		/* temporaries for stacking v's		*/
  float dt;		/* temporary for t0/2 difference	*/
  float dvh;		/* temporary for v*h difference		*/

  h[0] = 0.5 * vs[0] * t0[0];
  v[0] = vs[0];
  for (i = 1; i < n; i++) {
    t2 = 0.5 * t0[i]; t1 = 0.5 * t0[i-1];
    v2 = vs[i]; v1 = vs[i-1];
    dt = t2 - t1;
    dvh = v2*v2*t2 - v1*v1*t1;
    h[i] = sqrt(dvh * dt);
    v[i] = sqrt(dvh / dt);
  }


}
