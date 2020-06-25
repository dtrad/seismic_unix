
      SUBROUTINE CHOLESKY(M,EPS,A,B,ISTAT)
C
C This program solves a Hermitian symmetric set of complex linear
C simultaneous equations using the Cholesky decomposition method.
C The solution replaces the original contents of array B. Contents
C of array A are destroyed after this routine is called.
C
C                    AX = B
C
C   Input Parameters:
C        
C      M   - Order of the matrix (number of linear equations)
C      EPS - Epsilon (quantity for testing loss of significance;
C            depends on machine precision; suggest 1.E-15)
C      A   - Array of complex matrix elements stored columnwise
C            (i.e., A(1,1) is stored as A(1), A(1,2) as A(2),
C            A(2,2) as A(3), etc.  Only the top triangular part
C            of the A matrix is stored since the other half is
C            obtained by Hermitian symmetry)
C      B   - Array of complex elements of right-hand-side vector
C
C
C   Output Parameters:
C
C      B   - Complex solution X vector stored in place of B vector
C      ISTAT - Integer status indicator at time of exit
C              0 for normal exit
C              -1 if matrix is singular
C              +K if there is loss of numerical significance or if
C                 a nonpositive-definite matrix detected at pivot K
C
C   Notes:
C
C     External array A must be dimensioned .GE. M(M+1)/2 and array B
C     must be dimensioned .GE. M in the calling program.
C
c     ********************************************************
c     Remember to change the dimension of A is you change  nnx        
c     dim(A)=nnx*(nnx+1)/2
c     ********************************************************
	
      parameter    (nnx=228) 
      COMPLEX *16  A(nnx*(nnx+1)),B(nnx),SUM
      real    * 8  tiny,eps,dpiv 
C
C   Factor into triangular and diagonal form   !  Eq. (3.76)
C
      ISTAT=0
      KPIV=0
      DO 100 K=1,M
        KPIV=KPIV+K
        IND=KPIV
        LEND=K-1
        TINY=dABS(EPS*dREAL(A(KPIV)))
        DO 100 I=K,M
          SUM=(0.d0,0.d0)
          IF (LEND .EQ. 0)  GO TO 40
          LPIV=KPIV
          DO 30 L=1,LEND
            LPIV=LPIV+L-K-1
30          SUM=SUM+dREAL(A(LPIV))*A(IND-L)*dCONJG(A(KPIV-L))
40        SUM=A(IND)-SUM
          IF (I .NE. K)  GO TO 80
C
C   Test for negative pivot element and loss of significance
C
          IF (dREAL(SUM) .GT. TINY)  GO TO 90
          IF (dREAL(SUM) .GT. 0.d0)  GO TO 70
          ISTAT=-1
          RETURN
70        IF (ISTAT .GT. 0)  GO TO 90
          ISTAT=K
90        A(KPIV)=dCMPLX(dREAL(SUM),0.d0)
          DPIV=1.d0/dREAL(SUM)
          GO TO 100
80        A(IND)=SUM*DPIV
100       IND=IND+I
C
C   Back solution for intermediate column vector solution  ! Eq. (3.74)
C
      KPIV=1
      DO 200 K=2,M
        KPIV=KPIV+K
        SUM=B(K)
        DO 210 J=1,K-1
210       SUM=SUM-B(K-J)*dCONJG(A(KPIV-J))
200     B(K)=SUM
C
C   Back solution for final column vector solution    !  Eq. (3.75)
c
      KPIV=(M*(M+1))/2
      B(M)=B(M)/dREAL(A(KPIV))
      DO 300 K=M,2,-1
        KPIV=KPIV-K
        IND=KPIV
        SUM=B(K-1)/dREAL(A(KPIV))
        DO 310 J=K,M
          IND=IND+(J-1)
310       SUM=SUM-B(J)*A(IND)
300     B(K-1)=SUM
      RETURN
      END