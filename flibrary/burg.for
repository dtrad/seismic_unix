      SUBROUTINE BURG (N,IP,X,P,A,ISTAT)
C
C   Program to estimate the complex autoregressive parameters by
C   the Burg algorithm.
C
C   Input Parameters:
C
C     N  - Number of data samples
C     IP - Order of autoregressive process
C     X  - Array of complex data samples X(1) through X(N)
C
C   Output Parameters:
C
C     P  - Real variable representing driving noise variance
C     A  - Array of complex autoregressive parameters A(1) to A(IP)
C     ISTAT - Integer status indicator at time of exit
C             0 for normal exit
C             1 for numerical ill-conditioning (P < 0)
C
C   Notes:
C
C     External array X must be dimensioned .GE. N and array A
C     must be dimensioned .GE. IP in the calling program.
C     Internal arrays EF and EB must be dimensioned .GE. N .
C
      parameter  (nnt=2048)

      COMPLEX X(nnt),A(nnt),EF(nnt),EB(nnt),NUM,SAVE1,SAVE2
      REAL P,DEN,TEMP
C
C   Initialization
C
      ISTAT=0
      P=0.
      DO 10 J=1,N
10      P=P+REAL(X(J))**2+AIMAG(X(J))**2   ! Eq. (8.12)
      DEN=P*2.
      P=P/N
      IF (IP .EQ. 0)  RETURN
      DO 20 J=1,N
        EF(J)=X(J)
20      EB(J)=X(J)                         ! Eq. (8.11)
      TEMP=1.
      K=0
C
C   Main recursion
C
100   K=K+1
      NUM=(0.,0.)
      DO 30 J=K+1,N
30      NUM=NUM+EF(J)*CONJG(EB(J-1))
      DEN=TEMP*DEN-REAL(EF(K))**2-AIMAG(EF(K))**2
     *    -REAL(EB(N))**2-AIMAG(EB(N))**2  ! Eq. (8.10)
      SAVE1=-2.*NUM/DEN                    ! Eq. (8.14)
      TEMP=1.-REAL(SAVE1)**2-AIMAG(SAVE1)**2
      P=P*TEMP                             ! Eq. (8.4)
      IF (TEMP .GT. 0.)  GO TO 40
      ISTAT=1
      RETURN
40    A(K)=SAVE1
      IF (K .EQ. 1)  GO TO 60
      KHALF=K/2
      DO 50 J=1,KHALF
        KJ=K-J
        SAVE2=A(J)
        A(J)=SAVE2+SAVE1*CONJG(A(KJ))      ! Eq. (8.2)
        IF (J .EQ. KJ)  GO TO 50
        A(KJ)=A(KJ)+SAVE1*CONJG(SAVE2)     ! Eq. (8.2)
50      CONTINUE
60    IF (K .EQ. IP)  RETURN
      DO 70 J=N,K+1,-1
        SAVE2=EF(J)
        EF(J)=SAVE2+SAVE1*EB(J-1)          ! Eq. (8.7)
70      EB(J)=EB(J-1)+CONJG(SAVE1)*SAVE2
      GO TO 100
      END
 


