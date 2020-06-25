c-----------------------------------------------------------------------------
      
 
      SUBROUTINE HERM (M,T0,T,Z,X,ISTAT)
 
C   Solves the set of complex * 16 linear simultaneous equations
C                           TX = Z
C   by a variation of the Levinson algorithm.  T is a complex * 16 M+1
C   by  M+1  Hermitian Toeplitz matrix, Z is the known right-hand-
C   side complex * 16column vector of  M+1  elements,  and  X  is the
C   solution vector of  M+1  complex * 16elements.
C
C   Input Parameters:
C
C     M  - Order of matrix T (integer)
C     T0 - Scalar corresponding to real matrix element t(0)
C          (This element must be real due to Hermitian symmetry.)
C     T  - Array of  M  complex * 16matrix elements t(1),...,t(M)
C          from the left column of the Toeplitz matrix
C     Z  - Array of M+1 complex * 16elements of the right-hand-side
C          vector.   Program element Z(k+1) corresponds to text
C          element z(k), for k=0 to k=M
C
C   Output Parameters:
C
C     X  - Array of  M+1  complex * 16elements of solution vector.
C          Program element X(k+1) corresponds to text element
C          x(k), for k=0 to k=M
C     ISTAT - Integer status indicator at time of exit
C             0 for normal exit
C             1 if P=0. (dsingular matrix)
C
C   Notes:
C
C   External array T must be dimensioned .GE. M and arrays X,Z must
C   be dimensioned .GE. M+1 in the calling program.  Internal array
C   A must be dimensioned .GE. M .
C


      parameter  (nnt=228)
      COMPLEX   * 16    T(nnt),X(nnt),Z(nnt),A(nnt)
      COMPLEX   * 16    TEMP,SAVE,ALPHA,BETA
      REAL      * 8     P,T0
      P=T0
      ISTAT=1
      IF (P .EQ. 0.d0) RETURN
	
C   Handle  M=0  as a special case
      X(1)=Z(1)/dcmplx(T0,0.d0)
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
        SAVE=SAVE+A(J)*T(K-J)
10      BETA=BETA+X(J+1)*T(K-J)
20    TEMP=-SAVE/P
c     P=P*(1.- DREAL(TEMP)**2-DIMAG(TEMP)**2)
      p=p*(1.d0-dreal(temp*dconjg(temp)))
      IF (P .LE. 0.d0)  RETURN
30    A(K)=TEMP
      ALPHA=(Z(K+1)-BETA)/dcmplx(P,0.d0)
      IF (K .EQ. 1)  GO TO 50
      KHALF=K/2
      DO 40 J=1,KHALF
        KJ=K-J
        SAVE=A(J)
        A(J)=SAVE+TEMP*DCONJG(A(KJ))
        IF (J .EQ. KJ)  GO TO 40
        A(KJ)=A(KJ)+TEMP*DCONJG(SAVE)
40      CONTINUE
50    X(K+1)=ALPHA
      DO 60 J=1,K
60      X(J)=X(J)+ALPHA*DCONJG(A(K-J+1))
      IF (K .LT. M)  GO TO 100
      ISTAT=0
      RETURN
      END

    
c-----------------------------------------------------------------------------

