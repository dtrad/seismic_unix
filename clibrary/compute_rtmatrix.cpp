void compute_rtmtarix(complex *F, int maxfreq, float df, float *h, int nh, 
		      float *q, int nq, float *t, 
		      int nt, int fmax, float parmute)
{
  /* Variables with capital letters corresponds to matrices */
    

  int iq,ih,freq;
  float w=0,wa;
  float *g=0;  
  int rtmethod=2;
  float depth=2000; // Not used
  complex **L;    // RT oprator
  complex **RT=0; // LS RT operator
  float *mute=0; // mute operator 
  complex *rtm=0; // LS RT operator after mute 
  complex **RTM=0;
  float *stackop=0;

  g=ealloc1float(nh);
  L=ealloc2complex(nq,nh); 
  RT=ealloc2complex(nh,nk);
  mute=ealloc1float(nh);
  stackop=ealloc1float(nh);
  RTM=ealloc2complex(nq,nq);
  rtm=ealloc1complex(nq);

  radon_moveout(h,g,nh,rtmethod,depth); 
  mute_vector(mute,q,nq,parmute);

  stack_vector(stackop,nh);

  for (freq=1;freq<maxfreq;freq++){
     w=2*PI*freq*df;
     radon_matrix(L,g,q,nh,nq,w);
     LSmatrix(L,RT,nh,nq,eps); 
     xtimesA(rtm,mute,RT,nq,nh);
     Atimesx(RTM,L,rtm,nh,nq);
     xtimesy(finalop,stackop,RTM,nh);
     for (ih=0;ih<nh;ih++) F[freq][ih]=finalop[ih];
  }

  free1complex(rtm);
  free1complex(rtms);
  free2complex(FLS);
  free2complex(L);
  free1float(g);
  
  return;
}

void LSmatrix(complex **F, complex **FLST, int nh, int nk, float eps)
{
  int ik;
  complex **FF=0;
  complex **FT=0;

  FT=ealloc2complex(nh,nk);
  FF=ealloc2complex(nk,nk);


  AHtimesB(FF,F,F,nh,nk,nh,nk);
  // regularization
  for (ik=0;ik<nk;ik++) FF[ik][ik]+=eps;
  transpose(F,FT,nh,nk);

  inverse_matrix_multiply(nk,FF,nh,nk,FT,FLS);  


  free2complex(FT);
  free2complex(FF);


  return;
}

void mute_vector(float *mute,float *q, int nq, float parmute)
{
  for (iq=0;iq<nq;iq++)
    if (q[iq]<=parmute) mute[iq]=0;
    else mute[iq]=1;
  
  return;
}


void stack_vector(float *stackop, int nh)
{
  for (ih=0;ih<nh;ih++) stackop[ih]=1./nh;
  return;
} 

void AHtimesB(complex **C, complex **A, complex **B, int nr1, int nc1, int nr2, int nc2)
{
  // C= A^H x B (H means Hermitian) 
  int i,j,k, ni=0;
  int nr=nc1;
  int nc=nc2;
  if (nr1==nr2) ni=nr1;
  else err("nr1 must be equal to nr2");
  


  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      C[i][j].r=0;
      C[i][j].i=0;
      for (k=0;k<ni;k++)
	C[i][j]=C[i][j]+conjg(A[k][i])*B[k][j];
    }
  } 
  return;
}

void AtimesBH(complex **C,complex **A,complex **B,int nr1, int nc1, int nr2, int nc2)
{
  // C= A x B^H (H means Hermitian) 
  int i,j,k, ni=0;
  int nr=nr1;
  int nc=nr2;
  if (nc1==nc2) ni=nc1;
  else err("nr1 must be equal to nr2");
  

  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      C[i][j].r=0;
      C[i][j].i=0;
      for (k=0;k<ni;k++)
	C[i][j]=C[i][j]+A[i][k]*conjg(B[j][k]);
    }
  } 
  return;
}

void transpose(complex **F, complex **FT, int nr, int nc)
{
  int i,j;

  for (i=0;i<nr;i++)
    for (j=0;j<nc;j++)
      FT[j][i]=conjg(F[i][j]);
  return;
}


