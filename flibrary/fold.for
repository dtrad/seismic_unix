

!C FOLD IS CONVOLUTION
      SUBROUTINE FOLD (LA,A,LB,B,LC,C)
      DIMENSION A(2),B(2),C(2)
      LC=LA+LB-1
      !CALL ZERO(LC,C)
      DO 10 I=1,LA
      DO 10 J=1,LB
      K=I+J-1
10    C(K)=A(I)*B(J)+C(K)
      RETURN
      END
