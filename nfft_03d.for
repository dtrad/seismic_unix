

	program nfft_03d

c	nfft_03d is like nfft_03c but tries aCauchy window
c	04/01/01
c	tries an exponential
c	also a Cauchy

c	nfft does fast DFT
c	Uniform, Gaussian and Poisson sampling

c	harmonic data as input
c	2 filters are used

c	the parameters that I am trying so far are:
c	dt for i/p =.5
c	nt=100, m=100
c	fr1=.1, fr2=.6
c	for Poisson (2 and 2) (1,.5 and .1 appear to work just as well)
c	for Gaussian filter (.2 and 10 )
c	for Cauchy filter, at Poisson (.1 and 11, very tough), the
c	parameters are var=.4, lwind=100 give great results
c	also	       var=.2, lwind=100 give great results
c	for Hanning filter  (?)
c	decon length 101
c	resampling dt=.5

c	for resampling at .1 change G filter to (.05 and 4-10)
c	take care with the decon. The best appears to be 301 (max 501)
c	still  works well

c	!!! the NFFT and DFT results are virtually identical.

c	!!! there is something kaput with the Hanning

c	!!! the very interesting thing is the way the resampled
c	    time series compares with the original. Different look
c	    same spectrum



C   THIS PROGRAM COMPUTES THE F.T. USING  DFT
C   FULL FORWARD AND INVERSE TRANSFORMS ARE COMPUTED
c	the input is 2 frequencies, one not aliased, the other yes
c	a good choice is n=100, m=400, fr1=.25. fr2=.6
c	the Nyquist being at .5

c	dft_4 allows for chatter in the sampling (different chatter
c	to that in dft_5p and the dft_w programs - less chatter


c	and more irregular sampling)

	parameter(lp=10000)

c      IMPLICIT REAL*8(A-H,O-Z)
	REAL*4 QI(lp),BI(lp),bid(lp),t(lp),tn(lp)
	real*4 ps_q(lp),ps_quc(lp),ps_qu(lp),ps_qn(lp),ps_qnf(lp)
	real*4 ps_fft(lp)
	real*4 tsample(lp),tsample_n(lp),qi_n(lp),qi_u(lp)
	real*4 yplt(lp,10),w(1000),wr(lp)
	real lwind
	complex ai,DFT(lp),sum,cw(lp),work(lp),cwr(lp)
	integer isample(lp),isample_n(lp)

C
C   INPUT THE PARAMETERS

	write(*,*)'desired DT for I/P signal (0.5) qi, dt  =?'
	write(*,*)
	read(*,*)dti


	ai=cmplx(0.0,1.0)
	pi=acos(-1.)
	PI2=pi*2.

	write(*,*)
	write(*,*)'# of pnts in trace, N =?; # of pnts in DFT, M =?'
	write(*,*)
	read(*,*)N,M


c ****** Nr is the number of points in the irregular spacing

c ****** i/p the 2 frequencies
	write(*,*)
	write(*,*)'fr1 =?; fr2 =?'
	write(*,*)
	read(*,*)fr1,fr2


c ****** choose the type of pdf for the chatter
c	if ichat=0 it is uniform
c	if ichat=1 it is Gaussian
c	if ichat=2 it is Poisson
	write(*,*)'choose pdf, 0=uniform, 1=Gaussian, 2=Poisson'
	write(*,*)'ichat =?'
	write(*,*)
	read(*,*)ichat


		if(ichat.eq.1)then
	write(*,*)
	write(*,*)'sdev =?; seed =?'
	write(*,*)
	read(*,*)sdev,seed
	start=ran(seed)
		endif

		if(ichat.eq.0)then
	write(*,*)
	write(*,*)'max sdt =?; iseed =?'
	write(*,*)
	read(*,*)sdev,iseed
		endif

		if(ichat.eq.2)then
	write(*,*)
	write(*,*)'sampr =? (~1/Dt); iseed =?'
	write(*,*)
	read(*,*)sdev,iseed
		endif

c ****** choose the type of band-limiting filter
c	if iband=0 it is exponential
c	if iband=1 it is Gaussian
c	if iband=2 it is Cauchy
	write(*,*)'choose filter, 0=Exponential, 1=Gaussian, 2=Cauchy'
	write(*,*)'iband =?'
	write(*,*)
	read(*,*)iband
		do i=1,n
		w(i)=0.
		enddo

	if(iband.eq.1)then
	write(*,*)'variance of W, var =? ; cut-off length, lwind =?'
	write(*,*)
	read(*,*)var,lwind
	do i=1,n
	w(i)=exp(-(i-1)**2/var)
	enddo
	write(*,*)(w(i),i=1,lwind)
	pause '  '
	endif

c	nexth

	  if(iband.eq.0)then
	  write(*,*)'cut-off length, lwind =?, var =?'
	  write(*,*)
	  read(*,*)lwind,var
	  DO  I=1,lwind
c	   w(I)=0.5*(1.+COS(PI*FLOAT(I-1)/float(lwind)))
	  w(I)=exp(-abs(FLOAT(I-1)/var))
	  enddo
	  write(*,*)(w(i),i=1,lwind)
	  pause '  '
	  endif

	if(iband.eq.2)then
	write(*,*)'variance of W, var =? ; cut-off length, lwind =?'
	write(*,*)
	read(*,*)var,lwind
	do i=1,lwind
	w(i)=1./(1.+(i-1)**2/var)
	enddo
	write(*,*)(w(i),i=1,lwind)
	pause '  '
	endif


C   QI(I) ARE THE TIME SERIES WITH UNIFORM SAMPLING
C   QI_N(I) ARE THE TIME SERIES WITH NON-UNIFORM SAMPLING
C   QI_U(I) ARE THE TIME SERIES WITH UNIFORM RE-SAMPLING
c

c ****** the DFT is defined as
c	 DFT(k)=Sum(x(n)*exp(-i2*pi*dt*dk*nk)

C   INPUT THE MODEL PARAMETERS


c	create some  harmonics

	ndi=n/dti
	if(sdev.eq.0)then
	do i=1,ndi
	tn(i)=i-1
	enddo
	endif

	do i=1,Ndi
	t(i)=(i-1)
	tn(i)=i-1
	if(sdev.ne.0.and.ichat.eq.1)then
	call white(sdt,sdev,0.)
	tn(i)=tn(i)+sdt
	if(i.eq.1)tn(i)=0.
	endif
	 if(sdev.ne.0.and.ichat.eq.0)then
	   sdt=arand(iseed)
	   sdt=sdt*sdev*2.-sdev
	   tn(i)=tn(i)+sdt
	   if(i.eq.1)tn(i)=0.
	 endif
	  if(sdev.ne.0.and.ichat.eq.2)then
	      xx=arand(iseed)
	      sdt=-log(1.-xx)/sdev
	      tn(i)=tn(i)+sdt
	      if(i.eq.1)tn(i)=0.
	  endif

	qi(i)=sin(2*pi*t(i)*dti*fr1+.65)
	qi(i)=qi(i)+sin(2*pi*t(i)*dti*fr2+0.0)
	qi_n(i)=sin(2*pi*tn(i)*fr1+.65)
	qi_n(i)=qi_n(i)+sin(2*pi*tn(i)*fr2+0.0)

	enddo

c ****** now create the sampling arrays for plotting

			do i=1,n
			isample(i)=20.*t(i)+0.5
			isample_n(i)=20.*tn(i)+0.5
			enddo

		nn=n*10
		do i=1,nn
		tsample(i)=0.
		tsample_n(i)=0.
		enddo

			do i=1,n
			tsample(isample(i)+1)=1.
			tsample_n(isample_n(i)+1)=1.
			enddo
			tsample_n(1)=1.


	loop=0
c	 if(loop.eq.0)go to 2000	COMMENT 1
c	 if(loop.eq.0)go to 2000
1000	continue
	loop=loop+1


	write(*,*)'desired DT =?'
	write(*,*)
	read(*,*)dt

c ****** determine the number of points, Nr, for the regular
c	 interpolation

	Nr=N/Dt+0.5
	write(*,*)'total # of regular points  nr=',nr
	pause '  '
	mr=nr

	jjj=0
	k1=1
c	 k2=4*lwind
	k2=2/dt*lwind
c ****** the convolution loop
c ****** this band-limits and samples at equal intervals of Dt
		do i=1,nr
		t(i)=dt*(i-1)
		sumr=0.
	if(i.ge.2.and.abs(tn(k)-t(i)).gt.lwind)then
c	 k1=i/2-lwind
c	 k2=i/2+2*lwind
	k1=i*dt-lwind
	k2=i*dt+2*lwind
	endif

c	nexth

		do 101 k=k1,k2
		if(abs(tn(k)-t(i)).gt.lwind)go to 101
	      if(iband.eq.1)then
		sumr=sumr+qi_n(k)*exp(-(t(i)-tn(k))**2/var)
	      endif
	      if(iband.eq.0)then
c		 w(I)=0.5*(1.+COS(PI*((t(i)-tn(k)))/float(lwind)))
		w(I)=exp(-abs((t(i)-tn(k))/var))
		sumr=sumr+qi_n(k)*w(i)
	      endif
	      if(iband.eq.2)then
c		 w(I)=0.5*(1.+COS(PI*((t(i)-tn(k)))/float(lwind)))
		w(I)=1./(1.+(t(i)-tn(k))**2/var)
		sumr=sumr+qi_n(k)*w(i)
	      endif
101		continue
		qi_u(i)=sumr

		enddo

	call remdc(qi,ndi,qn)
	call remdc(qi_u,nr,qn)
	call remdc(qi_n,nr,qn)

	call norm(qi,ndi,qn)
	call norm(qi_u,nr,qn)
	call norm(qi_n,nr,qn)

c ****** for comparison the PS of q

		do i=1,mr
		dft(i)=0.
		dft(i)=qi(i)
		enddo
		call fft(dft,mr,-1,work)
		do i=1,mr
		ps_q(i)=dft(i)*conjg(dft(i))
		enddo

		call norm(ps_q,mr,pn)

c ****** now do the FFT and the deconvolution.
	Dk=1./(FLOAT(Nr)*Dt)

c ****** 1st the FFT of the Nr samples, qu

	do i=1,mr
	dft(i)=0.
	dft(i)=qi_u(i)
	enddo
	call fft(dft,mr,-1,work)
		do i=1,mr
		ps_quc(i)=dft(i)*conjg(dft(i))
		enddo
		call norm(ps_quc,mr,pn)


c ****** 2nd the FFT of the weighting filter, w(i)

		if(iband.eq.1)then
		do i=1,lwind
		w(i)=exp(-((i-1)*dt)**2/var)
		enddo
		do i=lwind+1,nr
		w(i)=0
		enddo
		endif

		    if(iband.eq.0)then
		    do i=1,lwind
c		     w(I)=0.5*(1.+COS(PI*FLOAT(I-1)*dt/float(lwind)))
		w(I)=exp(-abs((i-1)*dt)/var)
		    enddo
		    do i=lwind+1,nr
		    w(i)=0
		    enddo
		    endif

		if(iband.eq.2)then
		do i=1,lwind
		w(i)=1./(1.+((i-1)*dt)**2/var)
		enddo
		do i=lwind+1,nr
		w(i)=0
		enddo
		endif


c ****** prepare for FFT
		 do  i=nr/2+2,nr
		 w(i)=w(nr-i+2)
		 enddo


	do i=1,mr
	cwr(i)=0.
	cwr(i)=w(i)
	enddo
	call fft(cwr,mr,-1,work)
	do k=1,mr
	wr(k)=cwr(k)
	enddo
	write(*,*)'FFT of w(i)'
	write(*,100)(wr(i),i=1,mr/2+1)
100	format(10f7.3)


c ****** 3rd the deconvolution in the f domain
	write(*,*)
	write(*,*)'the # of FFT points for w(i) =?'
	write(*,*)
	read(*,*)mr0

			do i=1,mr0
c			 if(wr(i).lt.0.2)wr(i)=wr(i)+0.2
			DFT(i)=DFT(i)/Wr(i)
			enddo

c   do the flipping around the Nyquist
	do  i=mr/2+2,mr
	dft(i)=conjg(dft(mr-i+2))
	enddo
		do i=1,mr
		ps_qu(i)=dft(i)*conjg(dft(i))
		enddo
		call norm(ps_qu,mr,pn)

c ****** now reconstruct in the time domain
c	 again, use the FFT


	call fft(dft,mr,1,work)
	do k=1,mr
	bi(k)=dft(k)/mr
	enddo

	iper=10
       call TAPR1D(bi,Nr,IPER)

	call norm(bi,nr,qn)


c ****** now compare with the DFT results of tn

      DO  k=1,Mr
	sum=0
c	 do j=1,Nr
	do j=1,N
	sum=sum+qi_n(j)*cexp(-ai*pi2*tn(j)*dk*(k-1))
	enddo
      DFT(k)=sum
	enddo

		do i=1,mr
		ps_qn(i)=dft(i)*conjg(dft(i))
		enddo
		call norm(ps_qn,mr,pn)


c ****** taper to band-limit the transform
	iper=10
	call TAPR1C(dft,mr/2+1,IPER)

c   do the flipping around the Nyquist
	do 20 i=mr/2+2,mr
	dft(i)=conjg(dft(mr-i+2))
20	continue


c ****** now reconstruct in the time domain
c	 again, use the FFT

	call fft(dft,mr,1,work)
	do k=1,mr
	bid(k)=dft(k)/mr
	enddo

	call norm(bid,mr,qn)

c ****** check after the FFt of the resampled qn

	do k=1,mr
	dft(k)=bid(k)
	enddo
	call fft(dft,mr,-1,work)
		do i=1,mr
		ps_qnf(i)=dft(i)*conjg(dft(i))
		enddo
		call norm(ps_qnf,mr,pn)


c ****** check after the FFt of the  qn
	do k=1,mr
	dft(k)=qi_n(k)
	enddo
	call fft(dft,mr,-1,work)
		do i=1,mr
		ps_fft(i)=dft(i)*conjg(dft(i))
		enddo
		call norm(ps_fft,mr,pn)
c				next

c ****** PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

2000	continue
	if(loop.eq.0)then
	do  ii=1,nn
	yplt(ii,1)=tsample_n(ii)
	yplt(ii,2)=tsample(ii)
	enddo
	ns=2
	mp=nn+1
	ysep=.9
	mode=3
	delt=1.
	dbn=10
	xsiz=7.
	ylt=5.
	ysiz=(ylt-0.5)/float(ns+1)
c191	 format(a80)

c ****** now do the plotting setup

	call device(12,11,18,0,'oplotp ')
	call ploton

	call oldplt(0.)
	call newpen(14)

	t1=0.0
	call ftrace(yplt,lp,mp,ns,delt,mode,dbn,t1,Xsiz,Ysiz,ysep,
     1	' ,')

	call endplt

	call oldplt(0.)
	call worig(0.,0.)
	go to 1000
	endif

	do  ii=1,nr
	yplt(ii,1)=qi(ii)
	yplt(ii,2)=qi_u(ii)
	yplt(ii,3)=bi(ii)
	yplt(ii,4)=bid(ii)
	yplt(ii,5)=qi_n(ii)
c	 yplt(ii,4)=bi_1(ii)
	enddo
	ns=5
c	 mp=n+1
	mp=nr
		write(*,*)'# of points for plotting, mp =?'
		write(*,*)
		read(*,*)mp
	ysep=.9
	mode=1


	delt=1.
	dbn=10
	xsiz=7.
	ylt=5.
	ysiz=(ylt-0.5)/float(ns+1)
c191	 format(a80)

c ****** now do the plotting setup

	call device(12,11,18,0,'oplotp ')
	call ploton

	call oldplt(0.)
	call newpen(14)

	t1=0.0
	call ftrace(yplt,lp,mp,ns,delt,mode,dbn,t1,Xsiz,Ysiz,ysep,
     1	' ,')

	call endplt

	call oldplt(0.)
	call worig(0.,0.)


	do  ii=1,mr
	yplt(ii,1)=ps_q(ii)
	yplt(ii,2)=ps_quc(ii)
c	 yplt(ii,3)=ps_qn(ii)
	yplt(ii,3)=ps_qu(ii)
	yplt(ii,4)=ps_qnf(ii)
	yplt(ii,5)=ps_fft(ii)
	enddo
	ns=5
c	 mp=n+1
	mp=nr
	ysep=.9
	mode=1


	delt=1.
	dbn=10
	xsiz=7.
	ylt=5.
	ysiz=(ylt-0.5)/float(ns+1)
c191	 format(a80)

c ****** now do the plotting setup

	call device(12,11,18,0,'oplotp ')
	call ploton

	call oldplt(0.)
	call newpen(14)

	t1=0.0
	call ftrace(yplt,lp,mp,ns,delt,mode,dbn,t1,Xsiz,Ysiz,ysep,
     1	' ,')

	call endplt

	call oldplt(0.)
	call worig(0.,0.)


c ****** PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP


	stop
	end


C   ****************************************
      SUBROUTINE WHITE(NOISE,SDEV,MEAN)
C   ****************************************
      REAL MEAN,NOISE
      A=0.
      DO 1 J=1,12
      Z=RAN(0.)
1     A=A+Z
      NOISE=(A-6.)*SDEV+MEAN
      RETURN
      END


      REAL FUNCTION ARAND (IR)
******************************************************
*						     *
*     Generate uniform (0,1) pseudo-random numbers   *
*	  based on linear congruential method	     *
*						     *
*     At the first entry give any positive integer   *
*     less than 1664501 for IR			     *
*     But it must not be altered in the calling      *
*     program after the first entry		     *
*						     *
******************************************************
      IMPLICIT INTEGER*4 (I-N)
*
      PARAMETER (ML = 1664501)
      PARAMETER (IA = 1229, IC = 351750)
      PARAMETER (ANORM = 1.0 / ML)
*
      IR = MOD (IA * IR + IC, ML)
      ARAND = IR * ANORM
*
      RETURN
      END


      SUBROUTINE TAPR1C(XIN,NTT,IPER)
C   TAPR1c APPLIES A BELL OF IPER PERCENT TO COMPLEX DATA
C   XIN IS THE INPUT LENGTH NTT
C   XIN IS ALSO THE OUTPUT
      complex XIN(2)
      PI=3.141593
      RIX1=FLOAT(NTT*IPER)/100.
      IX1=RIX1
      IX2=NTT-IX1
      IF(IPER.EQ.0) GO TO 99
      DO 3 I=IX2,NTT
      xx=xin(i)
3     xin(I)=xx*0.5*(1.+COS(PI*FLOAT(IX2-I)
     # /FLOAT(IX1)))
99    RETURN
      END

      SUBROUTINE TAPR1D(XIN,NTT,IPER)
C   TAPR1D APPLIES A BELL OF IPER PERCENT TO DATA
C   XIN IS THE INPUT LENGTH NTT
C   YOUT IS THE OUTPUT
      DIMENSION XIN(NTT)
      PI=3.141593
      RIX1=FLOAT(NTT*IPER)/100.
      IX1=RIX1
      IX2=NTT-IX1
      IF(IPER.EQ.0) GO TO 99
      IXP1=IX1+1
      DO 2 I=1,IXP1
2     xin(I)=xin(I)*0.5*(1.+COS(PI*FLOAT(IXP1-I)
     # /FLOAT(IX1)))
      DO 3 I=IX2,NTT
3     xin(I)=xin(I)*0.5*(1.+COS(PI*FLOAT(IX2-I)
     # /FLOAT(IX1)))
99    RETURN
      END
