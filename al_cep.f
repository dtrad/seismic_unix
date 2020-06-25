c   ************************************************************
c  al_cep computes a smoothed spectrum using the cepstrum

	subroutine al_cep(x,nx,nft,lcut,water,amp,ceps)

c   inputs are

c   x(nx) is the input amp. spec.
c   nft is the usual for fft computation
c   lcut is the cut-off in the cepstrum

c   outputs are

c   amp(nft) is the output smoothed spectrum
c   ceps(nft) is the real cepstrum

c   ************************************************************

	real x(nx),amp(nft),c(1024),ceps(1024)
	complex*8 cx(1024),work(1024)

c   at this stage do the homorphic transformation
c   find the powerspectrum


c   find the log of the spectrum and the cepstrum


	do  i=1,nx
	c(i)=abs(x(i))
	enddo
	do i=nx+1,nft/2+1
	c(i)=water
	enddo
	do i=nft/2+2,nft
	c(i)=c(nft-i+2)
	enddo

	do i=1,nft
	cx(i)=log(c(i))
	enddo

	call fft(cx,nft,1,work)

	do 5 i=1,nft
	c(i)=cx(i)/float(nft)
5	ceps(i)=c(i)

c   ceps(i) is now the cepstrum

	do 11 i=lcut,nft-(lcut-2)
11	c(i)=0.0
	do 7 i=1,nft
7	cx(i)=c(i)
	call fft(cx,nft,-1,work)
	do 8 i=1,nft
	c(i)=cx(i)
	c(i)=exp(c(i))
8	amp(i)=c(i)
c8	 amp(i)=sqrt(amp(i))


	return
	end


      SUBROUTINE FFT(DATA,N,ISIGN,WORK)
C FFT IS FOURG DESCRIBED IN UBC FOURT
C CONSIDER DATA AS COMPLEX*8 INPUT AND OUTPUT
C ISIGN=-1 FOR FORWARD TRANSFORM
C ISIGN=+1 FOR INVERSE TRANSFORM (AND DIVIDE DATA BY N)
      DIMENSION DATA(1),WORK(1),IFACT(32)
      TWOPI=ABS(ATAN2(0.0,-1.0))*2.0
      IF(ISIGN.LT.0)TWOPI=-TWOPI
      IF=0
      NPART=N
      DO 50 ID=1,N,2
      IDIV=ID
      IF(ID-1)10,10,20
   10 IDIV=2
   20 IQUOT=NPART/IDIV
      IF(NPART-IDIV*IQUOT)40,30,40
   30 IF=IF+1
      IFACT(IF)=IDIV
      NPART=IQUOT
      GOTO20
   40 IF(IQUOT-IDIV)60,60,50
   50 CONTINUE
   60 IF(NPART-1)80,80,70
   70 IF=IF+1
      IFACT(IF)=NPART
   80 NFACT=IF
      IP0=2
      IP3=IP0*N
      IWORK=1
      I3REV=1
      DO 110 I3=1,IP3,IP0
      WORK(IWORK)=DATA(I3REV)
      WORK(IWORK+1)=DATA(I3REV+1)
      IP2=IP3
      DO 100 IF=1,NFACT
      IP1=IP2/IFACT(IF)
      I3REV=I3REV+IP1
      IF(I3REV-IP2)110,110,90
   90 I3REV=I3REV-IP2
  100 IP2=IP1
  110 IWORK=IWORK+IP0
      IWORK=1
      DO 120 I3=1,IP3,IP0
      DATA(I3)=WORK(IWORK)
      DATA(I3+1)=WORK(IWORK+1)
  120 IWORK=IWORK+IP0
      IF=0
      IP1=IP0
  130 IF(IP1-IP3)140,240,240
  140 IF=IF+1
      IFCUR=IFACT(IF)
      IP2=IP1*IFCUR
      THETA=TWOPI/FLOAT(IFCUR)
      SINTH=SIN(THETA/2.0)
      ROOTR=-2.0*SINTH*SINTH
      ROOTI=SIN(THETA)
      THETA=T
      THETA=TWOPI/FLOAT(IP2/IP0)
      SINTH=SIN(THETA/2.0)
      WSTPR=-2.0*SINTH*SINTH
      WSTPI=SIN(THETA)
      WMINR=1.
      WMINI=0.
      DO 230 I1=1,IP1,IP0
      IF(IFCUR-2)150,150,170
  150 DO 160 I3=I1,IP3,IP2
      J0=I3
      J1=I3+IP1
      TEMPR=WMINR*DATA(J1)-WMINI*DATA(J1+1)
      TEMPI=WMINR*DATA(J1+1)+WMINI*DATA(J1)
      DATA(J1)=DATA(J0)-TEMPR
      DATA(J1+1)=DATA(J0+1)-TEMPI
      DATA(J0)=DATA(J0)+TEMPR
  160 DATA(J0+1)=DATA(J0+1)+TEMPI
      GOTO220
  170 IWMAX=IP0*IFCUR
      DO 210 I3=I1,IP3,IP2
      I2MAX=I3+IP2-IP1
      WR=WMINR
      WI=WMINI
      DO 200 IWORK=1,IWMAX,IP0
      I2=I2MAX
      SUMR=DATA(I2)
      SUMI=DATA(I2+1)
  180 I2=I2-IP1
      TEMPR=SUMR
      SUMR=WR*SUMR-WI*SUMI+DATA(I2)
      SUMI=WR*SUMI+WI*TEMPR+DATA(I2+1)
      IF(I2-I3)190,190,180
  190 WORK(IWORK)=SUMR
      WORK(IWORK+1)=SUMI
      TEMPR=WR
      WR=WR*ROOTR-WI*ROOTI+WR
  200 WI=TEMPR*ROOTI+WI*ROOTR+WI
      IWORK=1
      DO 210 I2=I3,I2MAX,IP1
      DATA(I2)=WORK(IWORK)
      DATA(I2+1)=WORK(IWORK+1)
  210 IWORK=IWORK+IP0
  220 TEMPR=WMINR
      WMINR=WMINR*WSTPR-WMINI*WSTPI+WMINR
  230 WMINI=TEMPR*WSTPI+WMINI*WSTPR+WMINI
      IP1=IP2
      GOTO130
  240 RETURN
      END





