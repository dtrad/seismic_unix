c-----------------------------------------------------------------------------
 
      SUBROUTINE FFT(LX,X,ISIGNI)
 

c     Claerbout's FORK.FOR 

c     Input parameters:
c
c       X(LX)        - complex vector to be tranformed or inverse-transformed
c       ISIGNI.eq.-1 - means forward transform
c       ISIGNI.eq. 1 - means inverse transform
c
c     Output parameters:
c
c     X(LX)   - complex vector containing the discrete Fourier tranform
c               or the data depending on ISIGNI
c
c     Notes:
c
c       LX is an FFT number e.g., LX=J**2 where J is an integer.
c
c       The original f X(lx) is destroyed and replaced
c       by the Fourier transform or the Inverse Fourier Transform
c       depending on ISIGNI.




      

      COMPLEX   * 16 CARG,CW,CTEMP,X(lx)
      REAL      * 8  SIGNI,sc
      signi=dfloat(isigni)
      J=1
      DO 30 I=1,LX
      IF(I.GT.J) GO TO 10
      CTEMP=X(J)
      X(J)=X(I)
      X(I)=CTEMP
10    M=LX/2
20    IF(J.LE.M) GO TO 30
      J=J-M
      M=M/2
      IF(M.GE.1) GO TO 20
30    J=J+M
      L=1
40    ISTEP=2*L
      DO 50 M=1,L
      CARG=(0.d0,1.d0)*(3.14159265d0*SIGNI*dfloat(M-1))/dfloat(L)
      CW=CDEXP(CARG)
      DO 50 I=M,LX,ISTEP
      CTEMP=CW*X(I+L)
      X(I+L)=X(I)-CTEMP
50    X(I)=X(I)+CTEMP
      L=ISTEP
      IF(L.LT.LX) GOTO 40
      if(isigni.eq.-1) goto 60
      sc=1.d0/dsqrt(dfloat(lx)) 
      DO 70 I=1,LX
70    X(I)=X(I)*sc 
60    RETURN
      END

