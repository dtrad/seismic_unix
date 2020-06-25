#include "su.h"
#include "clibrary.h"

complex exp(complex);

void matrix_4(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, float *g)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        register int i;
	register int j;  
        complex  arg;
        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;
        for (i=0;i<nh;i++) Aperture+=dh[i];        
        
        for (j=0;j<nq;j++){
	  R[j].r=0;
	  R[j].i=0;
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
	      else if (rtmethod==3)         
		arg.i=w*g[i]*q[j];
    	      l[i][j]=exp(-arg);
	      lh[j][i]=(1./Aperture)*dh[i]*exp(arg);
              R[j]+=lh[0][i]*l[i][j]; //Top row of LL=LH*L
	  }
	}
  	        
        return;
}

void matrix_3(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        register int i;
	register int j;  
        complex  arg;
        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;
        for (i=0;i<nh;i++) Aperture+=dh[i];        

        for (j=0;j<nq;j++){
	  R[j].r=0;
	  R[j].i=0;
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
    	      l[i][j]=exp(-arg);
	      lh[j][i]=(1./Aperture)*dh[i]*exp(arg);
              R[j]+=lh[0][i]*l[i][j]; //Top row of LL=LH*L
	      //fprintf(stderr,"lh[%d][%d]=(%f,%f)\n",j,i,lh[j][i].r,lh[j][i].i);
	  }
	}
  	//TRACE;        
        return;
}


void matrix_3(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod,float *Wd)
{
//This version add Wd (inv(covariance matrix) = Wd^T Wd) to the R computation

        register int i;
	register int j;  
        complex  arg;
        float Aperture=0.0;
  
	for (i=0;i<nh;i++) Aperture+=dh[i];
        //fprintf(stderr,"Aperture=%f\n", Aperture);
        //Aperture=nh;
       
        for (j=0;j<nq;j++){
	  R[j].r=0;
	  R[j].i=0;
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
    	      l[i][j]=Wd[i]*exp(-arg);
	      lh[j][i]=Wd[i]/nh*exp(arg);
              R[j]+=lh[0][i]*l[i][j]; //Top row of LL=LH*L
	      //R[j]+=lh[0][i]*l[i][j];
	  }
	}
  	        
        return;
}




void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod)
  // This version does not compute R
{

        register int i;
	register int j;  
        complex  arg;
        float Aperture=0.0;
        for (i=0;i<nh;i++) Aperture+=dh[i]; 
        
        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
    	      l[i][j]=exp(-arg);
	      lh[j][i]=(1./Aperture)*dh[i]*exp(arg);
	  }
	}
  	        
        return;
}

void matrix_3(complex **l,float *pos,float *q,int nh,int nq,float w, int rtmethod)
  // This version does not compute R nor lh. Is call for hrrti.cpp for d=Lm
{

        register int i;
	register int j;  
        complex  arg; 
        
        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
    	      l[i][j]=exp(-arg);
	  }
	}
  	        
        return;
}
void matrix_3(complex **l,float *pos,float *q,int nh,int nq,float w, int rtmethod, float *g)
  // This version does not compute R nor lh. Is call for hrrti.cpp for d=Lm
{

        register int i;
	register int j;  
        complex  arg; 
        
        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
	      else if (rtmethod==3)
		arg.i=w*g[i]*q[j];
    	      l[i][j]=exp(-arg);
	  }
	}
  	        
        return;
}

void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, int flag)
  // This version does not compute R
{

        register int i;
	register int j;  
        complex  arg;


        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
	      l[i][j]=exp(-arg)/nq;
	      lh[j][i]=exp(arg)/nh;
	      //l[i][j]=exp(-arg)/sqrt(nq*nh);
	      //lh[j][i]=exp(arg)/sqrt(nq*nh);
	  }
	}
  	        
        return;
}


void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float *Wd, float dq, int rtmethod, int flag)
  // This version does not compute R and uses the factorization of data 
  // Covariance matrix Wd (such that inv(Cd)=Wd^T*Wd
{

        register int i;
	register int j;  
        complex  arg;
        //float Aperture=0.0;
        //for (i=0;i<nh;i++) Aperture+=dh[i]; 
        //float scale=1/sqrt(fabs(Aperture)*nq);
	

        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
    	      //l[i][j]=scale*dh2[i]*exp(-arg);
	      //lh[j][i]=scale*dh2[i]*exp(arg);
	      l[i][j]=Wd[i]*exp(-arg)/sqrt(nq*nh);
	      lh[j][i]=Wd[i]*exp(arg)/sqrt(nh*nq);

	  }
	}
  	        
        return;
}

void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float *Wd, float dq, int rtmethod, float *g)
  // This version does not compute R and uses the factorization of data 
  // Covariance matrix Wd (such that inv(Cd)=Wd^T*Wd
  // The difference with previous version is that here a
  // pseudo hyperbolic Forward and Adjoint transforms are
  // computed by using a refernce depth.
  // Ref. Foster and Mosther, Geophysics Vol57 No3 (1992)
{

        register int i;
	register int j;  
        complex  arg;
        //float Aperture=0.0;
        //for (i=0;i<nh;i++) Aperture+=dh[i]; 
        //float scale=1/sqrt(fabs(Aperture)*nq);
	

          

        for (j=0;j<nq;j++){
	  for (i=0;i<nh;i++){
              
              arg.r=0;
              if (rtmethod==1)
		arg.i=w*pos[i]*q[j];
              else if (rtmethod==2)
		arg.i=w*pos[i]*pos[i]*q[j];
	      else if (rtmethod==3)         
		arg.i=w*g[i]*q[j];
	      
    	      //l[i][j]=scale*dh2[i]*exp(-arg);
	      //lh[j][i]=scale*dh2[i]*exp(arg);
	      l[i][j]=Wd[i]*exp(-arg)/sqrt(nq*nh);
	      lh[j][i]=Wd[i]*exp(arg)/sqrt(nh*nq);

	  }
	}
	
        return;
}

void matrix_5(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, float *g)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99
 
        register int i;
	register int j;  
        complex  arg,arg2;
	complex exparg,exparg2;
	float dq0=q[1]-q[0];
	float q0=q[0];
        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;
        for (i=0;i<nh;i++) Aperture+=dh[i];        
	for (j=1;j<nq;j++) R[j].r=R[j].i=0;
	arg.r=0;
	arg2.r=0;
        for (i=0;i<nh;i++){
	  if (rtmethod==1){
	    arg.i=w*pos[i]*dq0;
	    arg2.i=w*pos[i]*q0;
	  }
	  else if (rtmethod==2){
	    arg.i=w*pos[i]*pos[i]*dq0;
	    arg2.i=w*pos[i]*pos[i]*q0;
	  }
	  else if (rtmethod==3){
 	    arg.i=w*g[i]*dq0;        
	    arg2.i=w*g[i]*q0;
	  }
	  exparg=exp(arg);
	  exparg2=exp(-arg2);
	  l[i][0]=exparg2;
	  lh[0][i]=exparg2;
	  R[0]=lh[0][i]*l[i][0]; 	  
	  for (j=1;j<nq;j++){
	    l[i][j]=l[i][j-1]*exparg;;
	    lh[j][i]=(1./Aperture)*dh[i]*lh[j-1][i]*exparg;
	    R[j]+=lh[0][i]*l[i][j]; //Top row of LL=LH*L
	  }
	}
  	        
        return;
}



void matrix_5(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod)
{

  //       Transformation matrices. L, LH and R= top row of LH*L
  //       This matrix relates the cmp gather and the velocity
  //       gather in the f-x space.

  //       Input parameters:
  //
  //       np   - number of parameters= number of traces of the velocity gather

  //         nh   - number of traces of the CMP
  //         q - radon parameter 
  //         pos  - offset
  //            w - the normalized freq. at which the transform is evaluated
  //         dh  - delta offset.
  //        rtmethod      1  LRT 2 PRT
  //        

  //       Out parameter:
  //
  //       Notes:
  //
  //       The parameter p in the velocity gather is the slowness
  //	  LH=FH.WU has size np x nh such that m=LH.u
  //	  L=F.WV has size nh x np such that u=L.m
  //      R= top row of LH*L 
  //		Daniel Trad- 22-02-99


        register int i;
	register int j;  
        complex  arg,arg2;
	complex exparg,exparg2;
	float dq0=q[1]-q[0];
	float q0=q[0];
        //      In main: 
	//      l=alloc2complex(nq,nh);  // ==> L(nh x nq)
	//	lh=alloc2complex(nh,nq); // ==> LH(nq x nh)	
        
        //dh[0]=dh[1];
        float Aperture=0.0;
        for (i=0;i<nh;i++) Aperture+=dh[i];        
	for (j=1;j<nq;j++) R[j].r=R[j].i=0;
	arg.r=0;
	arg2.r=0;
        for (i=0;i<nh;i++){
	  if (rtmethod==1){
	    arg.i=w*pos[i]*dq0;
	    arg2.i=w*pos[i]*q0;
	  }
	  else if (rtmethod==2){
	    arg.i=w*pos[i]*pos[i]*dq0;
	    arg2.i=w*pos[i]*pos[i]*q0;
	  }

	  exparg=exp(arg);
	  exparg2=exp(-arg2);
	  l[i][0]=exparg2;
	  lh[0][i]=exparg2;
	  R[0]=lh[0][i]*l[i][0]; 	  
	  for (j=1;j<nq;j++){
	    l[i][j]=l[i][j-1]*exparg;;
	    lh[j][i]=(1./Aperture)*dh[i]*lh[j-1][i]*exparg;
	    R[j]+=lh[0][i]*l[i][j]; //Top row of LL=LH*L
	  }
	}
  	        
        return;
}









