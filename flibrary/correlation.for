      SUBROUTINE CORRELATION (N,LAG,MODE,X,Y,R0,R)
C
C   This program computes either the unbiased or biased complex correlation
C   estimates between complex data sample arrays X and Y.  If X=Y, then the
C   autocorrelation is computed.
C
C   Input Parameters:
C
C     N    - Number of data samples in arrays X and Y (integer)
C     LAG  - Number of correlation lags to compute [ lags from  0  to  LAG
C            are computed and stored in R0 and R(1) to R(LAG) ] (integer)
C     MODE - Set to 0 for unbiased correlation estimates; otherwise, biased
C            correlation estimates are computed (integer)
C     X    - Array of complex data samples X(1) through X(N)
C     Y    - Array of complex data samples Y(1) through Y(N)
C
C   Output Parameters:
C
C     R0   - Complex correlation estimate for lag 0
C     R    - Array of complex correlation estimates for lags 1 to LAG
C
C   Notes:
C
C
C     External arrays X,Y must be dimensioned .GE. N and array R must be
C     dimensioned .GE. LAG in the calling program.
C
      COMPLEX X(1),Y(1),R(1),R0,SUM
      DO 30 K=0,LAG
        NK=N-K
        SUM=(0.,0.)
        DO 10 J=1,NK
10        SUM=SUM+X(J+K)*CONJG(Y(J))
        IF (K .NE. 0)  GO TO 20
        R0=SUM/FLOAT(N)
        GO TO 30
20      IF (MODE .EQ. 0)  R(K)=SUM/FLOAT(N-K)    ! Eq.  (5.9)
        IF (MODE .NE. 0)  R(K)=SUM/FLOAT(N)      ! Eqs. (5.13),(5.19)
30      CONTINUE
      RETURN
      END


c
