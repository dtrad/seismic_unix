      SUBROUTINE Toeplitz(M,T0,T,Z,X,ISTAT)
C      call toeplitz(n2-n1,r0,r,Z,xx,ISTAT)
C   Solves the set of complex linear simultaneous equations
C                           TX = Z
C   by a variation of the Levinson algorithm.  T is a complex  M+1
C   by  M+1  Hermitian Toeplitz matrix, Z is the known right-hand-
C   side complex column vector of  M+1  elements,  and  X  is the
C   solution vector of  M+1  complex elements.
C
C   Input Parameters:
C
C     M  - Order of matrix T (integer)
C     T0 - Scalar corresponding to real matrix element t(0)
C          (This element must be real due to Hermitian symmetry.)
C     T  - Array of  M  complex matrix elements t(1),...,t(M)
C          from the left column of the Toeplitz matrix
C     Z  - Array of M+1 complex elements of the right-hand-side
C          vector.   Program element Z(k+1) corresponds to text
C          element z(k), for k=0 to k=M
C
C   Output Parameters:
C
C     X  - Array of  M+1  complex elements of solution vector.
C          Program element X(k+1) corresponds to text element
C          x(k), for k=0 to k=M
C     ISTAT - Integer status indicator at time of exit
C             0 for normal exit
C             1 if P=0. (singular matrix)
C
C   Notes:
C
C   External array T must be dimensioned .GE. M and arrays X,Z must
C   be dimensioned .GE. M+1 in the calling program.  Internal array
C   A must be dimensioned .GE. M .
C
      COMPLEX T(1),X(1),Z(1),A(228)
      COMPLEX TEMP,SAVE,ALPHA,BETA
      REAL P,T0
      P=T0
      ISTAT=1
      IF (P .EQ. 0.)  RETURN
C   Handle  M=0  as a special case
      X(1)=Z(1)/T0                            ! Eq. (3.175)
      IF (M .LE. 0)  RETURN
C
C   Main recursion
C
      K=0
100   K=K+1
      SAVE=T(K)
      BETA=X(1)*T(K)
      IF (K .EQ. 1)  GO TO 20
      DO 10 J=1,K-1
        SAVE=SAVE+A(J)*T(K-J)                 ! Eq. (3.136)
10      BETA=BETA+X(J+1)*T(K-J)               ! Eq. (3.173)
20    TEMP=-SAVE/P
      P=P*(1.-REAL(TEMP)**2-AIMAG(TEMP)**2)   ! Eq. (3.158)
      IF (P .LE. 0.)  RETURN
30    A(K)=TEMP                               ! Eq. (3.139)
      ALPHA=(Z(K+1)-BETA)/P                   ! Eq. (3.174)
      IF (K .EQ. 1)  GO TO 50
      KHALF=K/2
      DO 40 J=1,KHALF
        KJ=K-J
        SAVE=A(J)
        A(J)=SAVE+TEMP*CONJG(A(KJ))           ! Eqs. (3.147),(3.157)
        IF (J .EQ. KJ)  GO TO 40
        A(KJ)=A(KJ)+TEMP*CONJG(SAVE)          ! Eqs. (3.147),(3.157)
40      CONTINUE
50    X(K+1)=ALPHA
      DO 60 J=1,K
60      X(J)=X(J)+ALPHA*CONJG(A(K-J+1))       ! Eq. (3.171)
      IF (K .LT. M)  GO TO 100
      ISTAT=0
      RETURN
      END
