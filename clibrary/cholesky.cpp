
void cholesky(complex M, float EPS,complex A,complex B, int ISTAT)
  /*
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
  */	
  //    complex  A(nnx*(nnx+1)),B(nnx),sum
float    tiny, dpiv 
int kpiv, k, i, ind, lend
complex sum, 
      //   Factor into triangular and diagonal form   !  Eq. (3.76)

istat=0;
kpiv=0
for (k=0;k<m;k++)
    kpiv=kpiv+k;
    ind=kpiv;
    lend=k-1
    tiny=abs(eps*real(A(kpiv)));
    for (i=k;k<m;i++){
        sum.r=0; sum.i=0;
        if (lend!= 0){
          lpiv=kpiv;
          for (l=0;l<lend;l++){
              lpiv=lpiv+l-k-1
30            sum=sum+real(A(lpiv))*A(ind-L)*conjg(A(kpiv-L))
	      }
	    }
40        sum=A(ind)-sum
          if (i==k){
C
C   Test for negative pivot element and loss of significance
C
          IF (dREAL(sum) .GT. TINY)  GO TO 90
          IF (dREAL(sum) .GT. 0.d0)  GO TO 70
          ISTAT=-1
          RETURN
70        IF (ISTAT .GT. 0)  GO TO 90
          ISTAT=K
90        A(kpiv)=dCMPLX(dREAL(sum),0.d0)
          dpiv=1.d0/dREAL(sum)
          GO TO 100
80        A(ind)=sum*dpiv
          }
100       ind=ind+I
      }
}
C
C   Back solution for intermediate column vector solution  ! Eq. (3.74)
C
      kpiv=1
      for 200 K=2,M
        kpiv=kpiv+K
        sum=B(K)
        for 210 J=1,K-1
210       sum=sum-B(K-J)*dCONJG(A(kpiv-J))
200     B(K)=sum
C
C   Back solution for final column vector solution    !  Eq. (3.75)
c
      kpiv=(M*(M+1))/2
      B(M)=B(M)/dREAL(A(kpiv))
      for 300 K=M,2,-1
        kpiv=kpiv-K
        ind=kpiv
        sum=B(K-1)/dREAL(A(kpiv))
        for 310 J=K,M
          ind=ind+(J-1)
310       sum=sum-b[j]*A(ind)
300     B(K-1)=sum
      return
	}