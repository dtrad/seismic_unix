#include "su.h"

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t)
{
	static int indx=0;
	int it;
	float a1,a2;

	/* if before first cdp, constant extrapolate */
	if (cdpt<=cdp[0]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[0][it];
			oa1t[it] = oa1[0][it];
			oa2t[it] = oa2[0][it];
		      };
	
	/* else if beyond last cdp, constant extrapolate */
	} else if (cdpt>=cdp[ncdp-1]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[ncdp-1][it];
			oa1t[it] = oa1[ncdp-1][it];
			oa2t[it] = oa2[ncdp-1][it];
		      };
	
	/* else, linearly interpolate */
	} else {
		xindex(ncdp,cdp,cdpt,&indx);
		a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
		a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
		for (it=0; it<nt; ++it) {
			ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
			oa1t[it] = a1*oa1[indx][it]+a2*oa1[indx+1][it];
			oa2t[it] = a1*oa2[indx][it]+a2*oa2[indx+1][it];
		      };
	}
}

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt)
{
  static int indx=0;
  int it;
  float a1,a2;
  
  /* if before first cdp, constant extrapolate */
  if (cdpt<=cdp[0]) for (it=0; it<nt; ++it) ovvt[it] = ovv[0][it];
  /* else if beyond last cdp, constant extrapolate */
  else if (cdpt>=cdp[ncdp-1]) for (it=0; it<nt; ++it) ovvt[it] = ovv[ncdp-1][it];
  /* else, linearly interpolate */
  else {
    xindex(ncdp,cdp,cdpt,&indx);
    a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
    a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
    for (it=0; it<nt; ++it) ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
  }
}











