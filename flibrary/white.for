      SUBROUTINE WHITE(NOISE,SDEV,MEAN)
!C   ****************************************
      REAL MEAN,NOISE
	integer*4 iseed
	iseed=1

      A=0.
      DO 1 J=1,12
      Z=RAN(iseed)
1     A=A+Z
      NOISE=(A-6.)*SDEV+MEAN
      RETURN
      END
