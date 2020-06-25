c       *********************************************************
c
c       PROGRAM: hr_vs.for
c
c       HIGH RESOLUTION SLANT STAKCS
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
c
c       *********************************************************
c
c       GENERAL DESCRIPTION:
c
c  -                           nh - number of traces
c                              nt - number of samples per trace
c                              dh - trace interval in meters
c                              dt - sampling interval in sec
c                     near_offset - position of the first trace in meters
c                              t0 - time where the window starts in sec
c                                   must start in t0=0sec.
c                        d(nt,nh) - CMP gather
c         fileout_1       - file where the velocity gather is written (unit 30).
c         fileout_2       - file where the reconstructed CMP is written (unit 40).
c         esp1            - prewhitening in % (1-5% is ok)         
c         pmin            - minimum 1/velocity to invert
c         pmax            - maximum 1/velocity to invert
c         nq              - number of traces of the velocity gather
c         iter_end        - maximum number of iterations (5 is 0.k.)
c         nw              - number of time-velocity windows that we desire
c                           to mask in the velocity gather (mask = filter out). 
c                           If nw=0 the velocity gather is not filtered.
c                           The following 4 parameters define each window:
c         tw1(nw) tw2(nw) - time range where the velocity is masked.
c         pw1(nw) pw2(nw) - range of slopes to mask.    

c        Output parameters:

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
c         The traces of the CMP (data) are considered equally spaced
c         the minumum  receiver-source distance is near_offset (meters)
c         The interval between traces of the CMP is dh. 
c  
c         Everything is done by subr. 'V_STACK' the rest of the main is
c         for INPUT/OUTPUT and to compute the reconstructed CMP, this 
c         is done by subr. 'sum_velocity' 

        parameter (nnt=512,nnx=128)

        real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax
        real * 8 near_offset,dt,dh,t0,dq,pmin,pmax 
        real * 8 tw1(10),tw2(10),pw1(10),pw2(10)
         
        character * 20 filein,fileout_1,fileout_2

        open(10,file='hr_ss.in',status='unknown')

c------ read data
 
        read(10,*) filein      ! data
        read(10,*) fileout_1   ! velocity gather is written here.
        read(10,*) fileout_2   ! data computed from the velocity gather.
        read(10,*) eps1,eps2
        read(10,*) pmin,pmax
        read(10,*) nq
        read(10,*) iter_end 
        read(10,*) nw 
        read(10,*) (tw1(i),tw2(i),pw1(i),pw2(i),i=1,nw)

        open(20,file=filein,status='unknown')
        open(30,file=fileout_1,form='binary',status='unknown')
        open(40,file=fileout_2,form='binary',status='unknown')
      
	  read(20,*) nh,nt
	  read(20,*) dh,dt
        read(20,*) near_offset,t0
	
        do 10 it=1,nt
10      read(20,*) (d(it,ih),ih=1,nh)

c------ maximum and minimum q parameter to invert

        qmin=pmin
        qmax=pmax
        dq=(qmax-qmin)/(nq-1) ! q interval

c-----  compute the slant stack

c      *************************************************

        call S_STACK(d,nh,nt,near_offset,dh,dt,
     #              m,nq,qmin,qmax,eps1,eps2,iter_end)

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
     #     qmin,qmax,dt,dh,near_offset)

c       ********************************

c------ write the reconstructed CMP gather. d_rec(nt,nh)

        write(40) nh,nt
        write(40) dh,dt
        write(40) near_offset,0.d0
        do 14 ih=1,nh
14      write(40) (d_rec(it,ih),it=1,nt)


31      format(1f15.5)
32      format(1f15.5)
        stop
        end 
 

c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------

        subroutine  S_STACK(d,nh,nt,near_offset,dh,dt,
     #              m,nq,qmin,qmax,eps1,eps2,iter_end)
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
c        near_offset- position of the near offset trace in meters.
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

        complex * 16 l(nnx,nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx)
        real * 8     tmax,dt,dq,dh,near_offset,qmin,qmax,eps1
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

        call matrix  (L,qmin,qmax,nq,nh,dh,w,near_offset)

c------ solve dc=L.mc using Cauchy-Gauss model

        call GAUSS_GAUSS(L,dc,mc,nh,nq,eps1)

c       call CAUCHY_GAUSS_1(L,dc,mc,nh,nq,
c    #           eps1,iter_end)

c       call CAUCHY_GAUSS_2(l,dc,mc,nh,nq,eps1,iter_end,eps2)

c------ aux2(if,iq) contains the freq-q space 
 
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
 
        subroutine t_2(nh,nt,dt,nts,dts,d,index)

c       Index.eq. 1 then stretch the time axis
c       Index.eq.-1 then unstretch the time axis
       

c       when INDEX=1 this is the IN/OUT

c       Input parameters:
c
c         d(nt,nh)    - data before tranformation
c               ds    - original sampling rate (sec)
c
c        Output parameters
c
c         d(nts,nh)  - data after tranformation
c         dts        - sampling rate in the stretched axis (sec**2)
c
c        Note:
c
c        When index=-1 the IN/OUT are interchanged

  
        parameter (nnt=512, nnx=128)
 
        real * 8 d(nnt,nnx),aux(nnt),dt,dts,t

        if(index.eq. 1) n=nts
        if(index.eq.-1) n=nt
        
        do 100 ih=1,nh
        do 200 i=1,n
        if(index.eq. 1) t=(dsqrt(dfloat(i)*dts))/dt
        if(index.eq.-1) t=     ((dfloat(i)*dt)**2)/dts
        if(t.ge.1.d0) then
        i1=idint(t)
        i2=i1+1
        aux(i)=d(i1,ih)+
     #            (d(i2,ih)-d(i1,ih))*(t-dfloat(i1))/dfloat(i2-i1)
        else
        aux(i)=0.d0
        endif
200     continue
        do 300 i=1,n
300     d(i,ih)=aux(i)
100     continue
                
        return
        end

    
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
        
        subroutine sum_slope(d,nh,nt,v,nq,
     #     qmin,qmax,dt,dh,near_offset)

c       Compute the data  by summation over slopes
c
c       Input parameters:
c      
c           v(nt,np)  -  velocity gather
c           dt        -  time interval in sec.
c           dh        -  space interval in meters
c           np        -  number of traces of the velocity gather 
c           qmin      -  minimum q parameter of the velocity gather(s^2/m^2)
c           qmax      -  maximum q parameter of the velocity gather(s^2/m^2) 
c         near_offset -  near offset in meters
c
c       Output parameters:
c
c         d(nt,nh)    - cdp gather 

        parameter (nnt=512, nnx=128)

        real * 8 d(nnt,nnx),qmin,qmax,q,dt,dh,v(nnt,nnx) 
        real * 8 h,tau,t,aux,near_offset,dq

        dq=(qmax-qmin)/(nq-1)
        do 100 ih=1,nh 
        h=dfloat(ih-1)*dh+near_offset
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

      subroutine CAUCHY_GAUSS_2(l,d,model,nd,np,eps1,iter_end,eps3)
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047


c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Cauchy-Gauss regularization.

c     Input Parameters:
c
c         L(nd,np) - Linear operator (complex)
c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).
c         iter_end - maximum number of iteration (5 is ok)

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
c     This version uses CHOLESKY method iteratively
c     retrieve teh solution. 

      parameter    (nnx=128) 

      complex *16  a(8256),model(nnx),d(nnx),L(nnx,nnx),sum 
      complex *16  aux(nnx),lambda(nnx) 
      real    * 8  power(nnx),p(nnx),pm
      real    * 8  eps1,eps2,eps3,s,snew,sold,th

c---- Initial model is the adjoint

      th=0.000001d0

      do 10 i=1,np
      model(i)=dcmplx(0.d0,0.d0)
      do 10 j=1,nd
      model(i)=model(i)+dconjg(l(j,i))*d(j)
10    continue

      s=0.0
      pm=0.0
      do i=1,np
       power(i)=dreal(model(i)*dconjg(model(i)))
        pm=pm+power(i)
        enddo
         pm=pm/np
            
         if(eps3.eq.0.d0) eps2=pm
         if(eps3.ne.0.d0) eps2=eps3*pm
         do i=1,np
         p(i)=1.d0+power(i)/eps2
         s=s+dlog(1.d0+power(i)/eps2)
      enddo
       sold=s

      do 1000 iter=1,iter_end

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
300   sum=sum+p(k)*dconjg(l(j,k))*l(i,k)
      n=n+1
      if(i.eq.j) a(n)=sum*(1.d0+eps1/100.d0)
      if(i.ne.j) a(n)=sum
200   continue

      call  CHOLESKY (nd,1.d-15,A,aux,ISTAT)

      do  i=1,nd
      lambda(i)=aux(i)
      enddo

      do i=1,np
      sum=dcmplx(0.d0,0.d0) 
       do  k=1,nd
        sum=sum+lambda(k)*dconjg(l(k,i))
         enddo 
        model(i)=p(i)*sum
       enddo

      pm=0.d0
      do i=1,np
       power(i)=dreal(model(i)*dconjg(model(i)))
       pm=pm+power(i)
      enddo 
      if(eps3.eq.0.d0) then 
      pm=pm/np
      eps2=pm
      endif

c----- If the relative variation of the cost function is le<th
c----- returns 

      s=0.d0 
      do i=1,np
       p(i)=1.d0+power(i)/eps2
        s=s+dlog(1.d0+power(i)/eps2)
         enddo 
          snew=s
           if(dabs(snew-sold)/dabs(snew).lt.th)  return 
            sold=s

1000    continue

        return
        end






c-----------------------------------------------------------------------------
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


      parameter  (nnt=128)
      COMPLEX   * 16    T(nnt),X(nnt),Z(nnt),A(nnt)
      COMPLEX   * 16    TEMP,SAVE,ALPHA,BETA
      REAL      * 8     P,T0
      P=T0
      ISTAT=1
      IF (P .EQ. 0.d0)  RETURN
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
c-----------------------------------------------------------------------------
        

      subroutine CAUCHY_GAUSS_1(l,d,m,nd,np,
     #           eps1,iter_end)

c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267   
c       Tel (h): (604)-224-0949   
c       Fax:     (604)-822-6047


c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Cauchy-Gauss regularization.
c-----USES LEVINSON

     
c     Input Parameters:
c
c         L(nd,np) - Linear operator (complex)
c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).
c         iter_end - maximum number of iteration (5 is ok)

c     Output Parameters:    
c
c         m(np) - vector of unkonw parameters
c
c      
c     Notes:
c         
c     The damping eps1 represents a percentage of the diagonal 
c     of the Toeplitz matrix which is solved in each iteration.
c         
c     remenber that since we are solving an underdetermine problem 
c     np>nd.
c

      parameter    (nnt=128,nnf=512)

      complex *16  t(nnt),m(nnt),d(nnt),L(nnt,nnf),sum 
      complex *16  aux(nnt),lambda(nnt) 
      real    * 8  eps1,eps,power(nnf),p(nnf)
      real    * 8  eps2,s,snew,sold,th,t0,norm    

      th=0.001
   
      eps=1.d-8

c---- Initial model is computed using the adjoint

      do 10 i=1,np
      m(i)=dcmplx(0.d0,0.d0)
      do 10 j=1,nd
      m(i)=m(i)+dconjg(l(j,i))*d(j)
10    continue

      do i=1,np
      power(i)=(dreal(m(i)*dconjg(m(i))))
      enddo

      iter_max=iter_end
      norm=0.d0
      do  i=1,np
      norm=norm+power(i)
      enddo 
      eps2=norm/np

c---- compute nonlinear weights and Cauchy norm S.

      s=0.d0
      do  i=1,np
      p(i)=1.d0+power(i)/eps2
      s=s+dlog(1.d0+power(i)/eps2)
      enddo
      sold=s

c---- start to minimize the Cauchy-Gauss cost function

      do 1000 iter=1,iter_max

       write(*,*) iter

      do  i=1,nd
      aux(i)=d(i)
      enddo

c---- compute the Toeplitz matrix
c---- the first lag is loaded in t(1), the second in t(2),... 
c---- to is the diagonal term of the Toeplitz form (the zero lag).
 
      do 200 i=1,nd
      sum=dcmplx(0.d0,0.d0)
      do 300 k=1,np
300   sum=sum+p(k)*dconjg(l(1,k))*l(i,k)
      if(i.eq.1.) then
      t0=sum*(1.d0+eps1/100.)   ! damping
      endif
      if(i.ne.1) t(i-1)=sum
200   continue

c---- solve the Teopliz form using Levinson

      ntot=nd-1
      call HERM (ntot,T0,T,aux,lambda,ISTAT)

c---- update the Lagrange multipliers

      do 400 i=1,np
      sum=dcmplx(0.d0,0.d0) 
      do 450 k=1,nd
450   sum=sum+lambda(k)*dconjg(l(k,i))
400   m(i)=p(i)*sum

c---- estimate new weighs and the Cauhcy cretrion, S.

      do i=1,np
      power(i)=(dreal(m(i)*dconjg(m(i))))
      enddo

      norm=0.d0
      do  i=1,np
      norm=norm+power(i)
      enddo 
      eps2=norm/np

c---- return if the Cauchy criterion does not change.

      s=0.d0
      do  i=1,np
      p(i)=1.d0+power(i)/eps2
      s=s+dlog(1.d0+power(i)/eps2)
      enddo
      snew=s
      if(dabs((snew-sold)/snew).lt.th) return
      sold=snew 


1000  continue

      return
      end

c-----------------------------------------------------------------
c-----------------------------------------------------------------
c-----------------------------------------------------------------

        subroutine matrix(L,pmin,pmax,np,nh,dh,w,near_offset)

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
c       The parameter p in the velocity gather is the square solowness
c       p= 1/(velocity)^2


   

        parameter (nnx=128)

        complex * 16 l(nnx,nnx)
        real * 8     dp,arg,dh,w,pmin,pmax,p,h,near_offset,sn

        dp=(pmax-pmin)/dfloat(np-1)
        sn=dsqrt(dp*dh)

        do 100 j=1,np
        p=pmin+(pmax-pmin)*dfloat(j-1)/dfloat(np-1)
        do 100 i=1,nh
        h=dfloat(i-1)*dh+near_offset
        arg=w*h*p
100     l(i,j)=cdexp(dcmplx(0.d0,-arg))*sn

        return
        end





c---------new new new 
      subroutine GAUSS_GAUSS(l,d,m,nd,np,eps1)

c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267   
c       Tel (h): (604)-224-0949   
c       Fax:     (604)-822-6047


c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Gauss-Gauss regularization (damped least squares).
c-----USES LEVINSON

     
c     Input Parameters:
c
c         L(nd,np) - Linear operator (complex)
c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).

c     Output Parameters:    
c
c         m(np) - vector of unkonw parameters
c
c     Notes:
c         
c     The damping eps1 represents a percentage of the diagonal 
c     of the Toeplitz matrix which is solved in each iteration.
c         
c     remenber that since we are solving an underdetermine problem 
c     np>nd.
c

      parameter    (nnt=128,nnf=512)

      complex *16  t(nnt),m(nnt),d(nnt),L(nnt,nnf),sum 
      complex *16  aux(nnt)
      real    * 8  eps1,eps,t0

      eps=1.d-8

c---- compute the Toeplitz matrix
c---- the first lag is loaded in t(1), the second in t(2),... 
c---- to is the diagonal term of the Toeplitz form (the zero lag).
 
      do 200 i=1,np
      sum=dcmplx(0.d0,0.d0)
      do 300 k=1,nd
300   sum=sum+dconjg(l(k,1))*l(k,I)
      if(i.eq.1.) then
      t0=sum*(1.d0+eps1/100.)   ! damping
      endif
      if(i.ne.1) t(i-1)=sum
200   continue

      do i=1,np
      aux(i)=0.
      do j=1,nd
      aux(i)=aux(i)+conjg(l(j,i))*d(j)
      enddo
      enddo

c---- solve the Teopliz form using Levinson

      ntot=np-1
      call HERM (ntot,T0,T,aux,m,ISTAT)


      return
      end









