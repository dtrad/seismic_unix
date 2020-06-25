c       *********************************************************
c
c       PROGRAM: hr_ss_2.for
c
c	  SLANT STAKCS
c
c       Jun/1994
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c
c       E-mail: sacchi@geop.ubc.ca
c
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047
c
c	  Modified for non- evenly spaced data
c	  by Daniel Trad - December 98		
c       E-mail: dtrad@geop.ubc.ca
c       *********************************************************
c
c       GENERAL DESCRIPTION:
c
c       Given the seimic data  computes a slant stack. The S.Stack 
c       is then maped back to the data space. Undesired components
c       can be masked in the velocity gather before the reconstruction 
c       of the data. Doing so, we can isolate and keep only
c       part of the information contain in the original data.
c       
c
c       Input parameters:
c
c         The following parameters are read from unit 10 called
c         'hr_ss.in'
c
c         filein          - file that contains (unit 20):
c						d(nt,nh) - CMP gather
c         fileout_1       - file where the velocity gather is written (unit 30).
c         fileout_2       - file where the reconstructed CMP is written (unit 40).
c         filein          - file that contains offset axis:
c         esp1            - prewhitening in % (1-5% is ok)         
c         pmin            - minimum 1/velocity to invert
c         pmax            - maximum 1/velocity to invert
c         nq              - number of traces of the velocity gather
c         nh				- number of traces
c         nt				- number of samples per trace
c         dt				- sampling interval in sec
c         t0				- time where the window starts in sec
c                                   must start in t0=0sec.
c         nw              - number of time-velocity windows that we desire
c                           to mask in the velocity gather (mask = filter out). 
c                           If nw=0 the velocity gather is not filtered.
c                           The following 4 parameters define each window:
c         tw1(nw) tw2(nw) - time range where the velocity is masked.
c         pw1(nw) pw2(nw) - range of slopes to mask.    
c
c        Output parameters:
c
c          m(nt,nq)       - Velocity gather. This is written in unit 30 (fileout_1)
c          d_rec(nt,nq)   - reconstructed CMP from the velocity gather.
c                           This is written in unit 40 (fileout_2)
c
c                   
c        Notes:
c      
c        The output is organized as follows:
c        m(it,iq), it=1,nt, iq=1,nq
c                    iq=1 is the first trace of the velocity
c                    gather corresponding to q parameter qmin. qmin=1./vmax^2.
c                    iq=nq is the last trace of the velocity
c                    gather corresponding to q parameter qmax. qmax=1./vmin^2.
c                    where vmin and vmax are the minimum and maximum stacking
c                    velocity seek by the procedure.
c                    it=1 corresponds to 0sec. and it=nt to the last sample.
c                 
  

        parameter (nnt=512,nnx=128)

        real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax, pos(nnx)
        real * 8 dt,t0,dq,pmin,pmax 
        real * 8 tw1(10),tw2(10),pw1(10),pw2(10),cn(nnx)
         
        character * 20 filein,fileout_1,fileout_2,fileout_3

        open(10,file='hr_ss_2.in',status='unknown')

c------ read data
 
        read(10,*) filein      ! data
        read(10,*) fileout_1   ! velocity gather is written here.
        read(10,*) fileout_2   ! data computed from the velocity gather.
        read(10,*) fileout_3   ! offset file 
        read(10,*) eps1,eps2
        read(10,*) pmin,pmax
        read(10,*) nq,nh,nt
	  read(10,*) dt,t0
        read(10,*) nw 
        if (nw.ne.0) read(10,*) (tw1(i),tw2(i),pw1(i),pw2(i),i=1,nw)


        open(20,file=filein,form='binary',status='unknown')
        open(30,file=fileout_1,form='binary',status='unknown')
        open(40,file=fileout_2,form='binary',status='unknown')
        open(50,file=fileout_3,status='unknown')      
	
        do ih=1,nh
	    do it=1,nt
			read(20) d(it,ih)
		enddo
		read(50,*) pos(ih)
		cn(ih)=1/eps2  ! assign a value to the inverse of covariance
	  enddo

c------ maximum and minimum q parameter to invert

        qmin=pmin
        qmax=pmax
        dq=(qmax-qmin)/(nq-1) ! q interval

c-----  compute the slant stack

c      *************************************************

        call S_STACK_offset(d,nh,nt,pos,dt,
     #              m,nq,qmin,qmax,eps1,cn)

c      *************************************************
       

        write(30) nq,nt
        write(30) dq,dt
        write(30) qmin,t0
        do 12 iq=1,nq
12      write(30) (m(it,iq),it=1,nt)

 
c-----   THE UNDESIRE PARTS OF THE VELOCITY GATHER ARE
c-----   FILTERED HERE. THEN IT IS TRANSFORM BACK TO
c-----   THE OFFSET SPACE. 

        if (nw.eq.0) goto  1

        do iw=1,nw
        it1=idint(tw1(iw)/dt+1)
        it2=idint(tw2(iw)/dt+1)
        iq1= idint((pw1(iw)-qmin)/dq)
        iq2= idint((pw2(iw)-qmin)/dq)
        do it=it1,it2 
        do iq=iq1,iq2
        m(it,iq)=0.d0
        enddo
        enddo
        enddo


c-----  the velocity gather is used to recostruct the CMP gather
c-----  you can change the parametrs of the recostruction. 
c-----  i.e., you can use a different near_offset and 
c-----  and a different nunber of traces. This is used 
c-----  to extrapolate near and far offset traces.


c       ********************************

1       call sum_slope(d_rec,nh,nt,m,nq,
     #     qmin,qmax,dt,pos)

c       ********************************

c------ write the reconstructed CMP gather. d_rec(nt,nh)

        write(40) nh,nt
        write(40) 0.d0,dt
        write(40) 0.d0,0.d0
        do 14 ih=1,nh
14      write(40) (d_rec(it,ih),it=1,nt)


31      format(1f15.5)
32      format(1f15.5)
        stop
        end 
 

c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------

        subroutine  S_STACK_offset(d,nh,nt,pos,dt,
     #              m,nq,qmin,qmax,eps1,cn)
c
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047

c      Given a CMP gather compute the Velocity Stack
c
c      Input parameters:
c
c        d(nt,nh)   - CMP gather
c        nt         - number of samples of each trace
c        nh         - number of traces
c        nq         - number of velocity traces to estimate
c        qmin       - minimum q parameter  (sec^2/m^2)
c        qmax       - maximum q parameter  (sec^2/m^2)
c        eps1       - damping in % (1.-5% is o.k)
c        iter_end   - maximum number of iterations (5 is o.k.)
c
c      Output paramters:
c 
c        m(nt,nq)   - VELOCITY GATHER
c 
c      Notes:
c  
c        The array m(nt,nq) is  organized as follows:
c        m(it,iq), it=1,nt, iq=1,nq
c                    iq=1 is the first trace of the velocity
c                    gather corresponding to q parameter qmin=1./vmax^2.
c                    iq=nq is the last trace of the velocity
c                    gather corresponding to a q parameter qmax=1./vmin^2.
c                    Where vmin and vmax are the min and max velocities
c                    that you expect to see.
c

 
        parameter    (nnx=128, nnt=512)

        complex * 16 l(nnx,nnx), lh(nnx,nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx), pos(nnx), cn(nnx)
        real * 8     tmax,dt,dq,qmin,qmax,eps1
        real * 8     eps2,f_low, f_high, fn, w ,pi

        pi=4.d0*datan(1.d0)
        
        Tmax=dt*(nt-1)
        dq=(qmax-qmin)/dfloat(nq-1) !sampling interval of the variable q

c------ compute the FFT number

        do 30 j=1,20
        nfft=2**j
        if(nfft.ge.nt) goto 40
30      continue
40      continue
        n2=nfft/2+1

c------ Transform data to f-x 

        call FX_go (nh,nt,nfft,d,aux2,-1)

c------ Define the band in which the inversion is carried out

        fn=1./(2.d0*dt)
        f_low=0.
        f_high=fn
        kl=f_low*dfloat(nfft)*dt+1
        kh=f_high*dfloat(nfft)*dt+1
        kl=kl+2
        kh=kh-2

c------ Make the inversion for each frequency

        do 1000 if=kl,kh

        
c------ dc contains the data at frequcny if

        do ih=1,nh
        dc(ih)=aux2(if,ih)  
        enddo

        w=2.d0*pi*dfloat(if-1)/(dfloat(nfft)*dt) ! freq in rad/sec

c------ compute the linear operator L, that linkd the data with the velocity gather

	  call matrix_off(L,lh,qmin,qmax,pos,nq,nh,w)
        
c------ solve dc=L.mc using Cauchy-Gauss model

        call GAUSS_GAUSS_OFF(L,lh,dc,mc,nh,nq,eps1,cn)

        do iq=1,nq
        aux2(if,iq)=mc(iq)
        enddo

1000    continue

c------ impose symmetry to the frequency axis

        call FX_symmetry (aux2,nq,nt,nfft,kl,kh)

c------ transform the freq-q space to tau-q space

        call FX_go (nq,nt,nfft,m,aux2, 1)   


       
        return 
        end 


c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------

        subroutine FX_go (nx,nt,nfft,stx,sfx,index)


c       Transform the data to the FX domain
c       and viceversa depending on index
c
c       Input parameters:

c         stx(nt,nx)      -  if index.eq.-1  time-offset  data
c         stf(nfft,nx)    -  if index.eq. 1  freq-offset data
c         nfft            -  lenght of the transform

c       Output parameters:

c         stf(nfft,nx)    - if index.eq.-1  freq-offset
c         stx(nt,nx)      - if index.eq. 1  time-offset
c         nfft            - lenght of the transform
c
c       Notes:
c
c        The input/output changes according to index:
c
c          index = -1 TX ----> FX
c          index =  1 FX ----> TX



        parameter   (nnt=512, nnx=128)
        complex * 16    sfx(nnt,nnx), aux(nnt)
        real    * 8     stx(nnt,nnx)

        if(index.eq.-1)  then

        do 140 ix=1,nx
        do 130 it=1,nt
130     aux(it)=dcmplx(stx(it,ix),0.d0)
        do 120 it=nt+1,nfft
120     aux(it)=dcmplx(0.d0, 0.d0)

        call fft(nfft,aux,-1)

        do 110 if=1,nfft
110     sfx(if,ix)=aux(if)
140     continue
         
        endif

        if(index.eq.1) then
        do 240 ix=1,nx
        do 220 if=1,nfft
220     aux(if)=sfx(if,ix)

        call fft(nfft,aux, 1)

        do 230 it=1,nfft
230     stx(it,ix)=dreal(aux(it))
240     continue
        
        endif

        return
        end


c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------


        subroutine FX_symmetry (sfx,nx,nt,nfft,iflow,ifhigh)


c       Imposes symmetry to the FX domain
c
c       Input parameters:
c
c         sfx(nfft/2+1,nx) - signal in the FX. Only positive freqs.
c
c
c       Output parameters:

c         sfx(nfft,nx)   - signal in the FX with conjugate symmetric
c                          freq. axis.



        parameter        (nnt=512, nnx=128)
        complex   * 16   sfx(nnt,nnx)


c-------makes the f-h or f-q symmetric

        do 200 ix=1,nx
        do 180 if=1,iflow
180     sfx(if,ix)=dcmplx(0.d0, 0.d0)
        do 190 if=ifhigh+1,nfft/2+1
190     sfx(if,ix)=dcmplx(0.d0, 0.d0)
200     continue

        do 210 ix=1,nx
        do 210 if=nfft/2+2,nfft
210     sfx(if,ix) = dconjg(sfx(nfft-if+2,ix))

        

        return
        end

c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
 
      SUBROUTINE FFT(LX,X,ISIGNI)
 

c     Claerbout's FORK.FOR 

c     Input parameters:
c
c       X(LX)        - complex vector to be tranformed or inverse-transformed
c       ISIGNI.eq.-1 - means forward transform
c       ISIGNI.eq. 1 - means inverse transform
c
c     Output parameters:
c
c     X(LX)   - complex vector containing the discrete Fourier tranform
c               or the data depending on ISIGNI
c
c     Notes:
c
c       LX is an FFT number e.g., LX=J**2 where J is an integer.
c
c       The original f X(lx) is destroyed and replaced
c       by the Fourier transform or the Inverse Fourier Transform
c       depending on ISIGNI.




      parameter (nnt=512)

      COMPLEX   * 16 CARG,CW,CTEMP,X(nnt)
      REAL      * 8  SIGNI,sc
      signi=dfloat(isigni)
      J=1
      DO 30 I=1,LX
      IF(I.GT.J) GO TO 10
      CTEMP=X(J)
      X(J)=X(I)
      X(I)=CTEMP
10    M=LX/2
20    IF(J.LE.M) GO TO 30
      J=J-M
      M=M/2
      IF(M.GE.1) GO TO 20
30    J=J+M
      L=1
40    ISTEP=2*L
      DO 50 M=1,L
      CARG=(0.d0,1.d0)*(3.14159265d0*SIGNI*dfloat(M-1))/dfloat(L)
      CW=CDEXP(CARG)
      DO 50 I=M,LX,ISTEP
      CTEMP=CW*X(I+L)
      X(I+L)=X(I)-CTEMP
50    X(I)=X(I)+CTEMP
      L=ISTEP
      IF(L.LT.LX) GOTO 40
      if(isigni.eq.-1) goto 60
      sc=1.d0/dsqrt(dfloat(lx)) 
      DO 70 I=1,LX
70    X(I)=X(I)*sc 
60    RETURN
      END


c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
        
        subroutine sum_slope(d,nh,nt,v,nq,
     #     qmin,qmax,dt,pos)

c       Compute the data  by summation over slopes
c
c       Input parameters:
c      
c           v(nt,np)  -  velocity gather
c           dt        -  time interval in sec.
c           pos       -  offset axis in meters
c           np        -  number of traces of the velocity gather 
c           qmin      -  minimum q parameter of the velocity gather(s^2/m^2)
c           qmax      -  maximum q parameter of the velocity gather(s^2/m^2) 
c
c       Output parameters:
c
c         d(nt,nh)    - cdp gather 

        parameter (nnt=512, nnx=128)

        real * 8 d(nnt,nnx),qmin,qmax,q,dt,v(nnt,nnx) 
        real * 8 h,tau,t,aux,pos(nnx),dq

        dq=(qmax-qmin)/(nq-1)
        do 100 ih=1,nh 
        h=pos(ih)
        do 100 it=1,nt
        d(it,ih)=0.d0
        t=(it-1)*dt
        do 200 iq=2,nq
        q=qmin+(qmax-qmin)*dfloat(iq-1)/dfloat(nq-1)
        tau=t-h*q
        if(tau.gt.0.d0) then  
        tau=tau/dt
        if(tau.ge.1d0) then
        i1=idint(tau-tw1/dt)
        i2=i1+1
        aux=v(i1,iq)+
     #        (v(i2,iq)-v(i1,iq))*(tau-dfloat(i1))/dfloat(i2-i1)
        else
        aux=0.d0
        endif
        d(it-1,ih)=d(it-1,ih)+aux
        endif 
200     continue
100     continue

        return
        end



c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------

      SUBROUTINE CHOLESKY (M,EPS,A,B,ISTAT)
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
	
      parameter    (nnx=128) 
      COMPLEX *16  A(8256),B(nnx),SUM
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



c------------------------------------------------

      subroutine GAUSS_GAUSS_OFF(l,lh,d,model,nd,np,eps1,cn)
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047
c	  Modified to include not evenly spaced offset
c       Daniel Trad-- 29-12-98	
c       E-mail: dtrad@geop.ubc.ca

c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Gauss-Gauss regularization.

c     Input Parameters:
c
c         L(np,nd) - Linear operator (complex)
c         LH(nd,np) - Linear operator (complex)

c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).


c     Output Parameters:
c
c         m(np) - vector of unkonw parameters
c
c     Notes:
c
c     The damping eps1 represents a percentage of the diagonal
c     of the matrix L.L^H.
c
c     remenber that since we are solving an underdetermine problem
c     np>nd.
c
c     This version uses CHOLESKY method to retrieve the solution. 

      parameter    (nnx=128) 

      complex *16  a(nnx*(nnx+1)),model(nnx),L(nnx,nnx),LH(nnx,nnx) 
      complex *16  d(nnx),aux(nnx),lambda(nnx),sum 
      real    * 8  power(nnx),cn(nd),eps1

c---- Initial model is the adjoint

  
      do  i=1,nd
      aux(i)=d(i)
      enddo

c-----compute the matrix to invert
c-----eps1 is the damping 1% is ok

      n=0
      do 200 j=1,nd 
      do 200 i=1,j  
      sum=dcmplx(0.d0,0.d0)     
      do 300 k=1,np
300   sum=sum+(l(k,j))*lh(i,k)
      n=n+1
      if(i.eq.j) A(n)=sum*(1.d0+eps1/100.d0)
      if(i.ne.j) A(n)=sum
200   continue

      call  CHOLESKY (nd,1.d-15,A,aux,ISTAT)

      do  i=1,nd
      lambda(i)=aux(i)
      enddo

      do i=1,np
		sum=dcmplx(0.d0,0.d0) 
		do  k=1,nd
			sum=sum+l(i,k)*lambda(k)*Cn(k)
		enddo 
		model(i)=sum
      enddo
     
      return
      end






c-----------------------------------------------------------------------------

      subroutine matrix_off(L,lh,pmin,pmax,pos,np,nh,w)

c       Transformation matrix.
c       This matrix relates the cmp gather and the velocity
c       gather in the f-x space.

c       Input parameters:
c
c         np   - number of parameters= number of traces of the velocity gather

c         nh   - number of traces of the CMP
c         pmin - minimum parameter seek by the transform
c         pmax - maximum parameter seek by the transform
c            w - the normalized freq. at which the transform is evaluated
c   near_offet - near offset trace in meters.
c

c       Out parameter:
c
c       Notes:
c
c       The parameter p in the velocity gather is the slowness
c	  L=F.WU has size np x nh such that m=L.u
c	  LH=FH.WV has size nh x np such that u=LH.m

        parameter (nnx=128)

        complex * 16 l(nnx,nnx), lh(nnx,nnx)
        real * 8     dp,dh,arg,w,pmin,pmax,p,h
        real * 8     pos(nnx)

        dp=(pmax-pmin)/dfloat(np-1)
        do 100 j=1,np
        p=pmin+(pmax-pmin)*dfloat(j-1)/dfloat(np-1)
        do 100 i=1,nh
        h=pos(i)
	
        if (i.eq.1) then
		dh=pos(i+1)-pos(i)
	  else
		dh=pos(i)-pos(i-1)
	  endif
	   
        arg=w*h*p
	  l(j,i)=dh*cdexp(dcmplx(0.d0,arg))
100     lh(i,j)=dp*cdexp(dcmplx(0.d0,-arg))
    	  	
        return
        end

c---------------------------------------------------
