#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/cwp.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"

//complex cexp(complex);

void matrix_2(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod)
{
  //       Transformation matrix.
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         pmin - minimum parameter seek by the transform
  //         pmax - maximum parameter seek by the transform
  //            w - the normalized freq. at which the transform is evaluated
  //   near_offet - near offset trace in meters.
  //

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  L=F.WU has size np x nh such that m=L.u
  //	  LH=FH.WV has size nh x np such that u=LH.m
  //
  //     Note: I have eliminated dh from L because I use a formulation 
  //     that allows cancelation 
  //     of the dh terms. This allows us to use Cholesky instead of LU. 
  //     If Levinson is used (gauss_gauss_To.f) Wu is multiplied by F to get L.
  //		Daniel Trad- 22-02-99
 
        register int i;
	register int j;  
        complex  arg;
        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	

        for (j=0;j<nq;j++)
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];

	      l[i][j]=dq*cwp_cexp1(-1*arg);
	      lh[j][i]=dh[i]*cwp_cexp1(arg);
	} 	  	
        return;
}






