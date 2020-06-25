c  Gap filling program for Shots and CMP Gathers 

c  This is the gap filling program that uses AR
c  spectral modeling to retrieve a continuous gap.

c  This program read data and parameters from
c  a file called ``gap.in''. 

c  parameters read from gap.in

c  filein  : a binary file containing the cdp or shot
c  fileout : a binary file containing the cdp or shot
c  nt,nx   : number of samples and number of channels of the cmp
c  n1,n2   : from where to where the gap extends
c  ip      : lenght of the peo filter 

c the ouput goes into another binary called 
c fileout. If the output file is fileone='myfile'
c you can use SU to see it
c
c xwigb n1=800<myfile perc=99&
c where n1=nt=800 is the number of samples 


c TO PROCESS REAL DATA: 
c     I think that the best is to read segy files
c     with gaps fill in with zeros traces and
c     know from where to where the gap goes (n1,n2).
c     In this way your header is fix and you don't
c     have to create it for the new
c     traces in the gap.



 
c (c) M.D.SACCHI, Bayesian Solutions Ltd.
c                 Edmonton, AB, Canada
c                 (403)-438-4336
c                 sacchi@phys.ualberta.ca
 


      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)
      character  * 20 filein,fileout 

      open(10,file='gap.in')
      read(10,'(a20)') filein
      read(10,'(a20)') fileout
      read(10,*) nt,nx
      read(10,*) dt,dx
      read(10,*) n1,n2
      read(10,*) ip

      call read_bin(filein,nt,nx,data)

c add a gap here

      do ix=n1,n2
      do it=1,nt
       data(it,ix)=0
        enddo
         enddo

c  Write the data with gap

       call gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)
       call write_bin(fileout,nt,nx,data)

       stop
       end

c     -------------------------------------------

      subroutine gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)

      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx),mu
      complex    data_fx(nnt,nnx)
      complex    x(nnx),a(nnx)
      complex    x2(nnx),a2(nnx)


      nf=2048
      if(nf.lt.nt) pause 'Check NF'

 
      call FX_go (nx,nt,nf,data,data_fx,-1)

c Do it for all freqs in the band

      iflow = 5
      ifhigh = nf/2+1-5

      do if=iflow,ifhigh
        do ix=1,nx
         x(ix)=data_fx(if,ix)
           enddo
        call BURG (N1-1,IP,X,P,A,ISTAT)
       do ix=1,nx-n2
        x2(ix)=data_fx(if,ix+n2)
         enddo
          call BURG (Nx-n2,IP,X2,P,A2,ISTAT)
         do i=1,ip
          a(i)=0.5*(a(i)+a2(i))
           enddo
       call gapf(x,nx,n1,n2,ip,a)   
           do ix=1,nx
            data_fx(if,ix)=x(ix)
             enddo
          enddo

      call  FX_symmetry (data_fx,nx,nt,nf,iflow,ifhigh)
      call  FX_go (nx,nt,nf,data,data_fx,1)

      return
      end

c     -------------------------------------------------------
c
c     This is the good one. 
c     Works ok for real and complex data. 
c     We do not consider any type of symetry,
c     the conitinuity of the gap is required.

      subroutine gapf(x,n,n1,n2,ip,a)   

c   1D GAP filling.

      parameter  (nnt=2048,nnx=128)
      complex    x(nnt),a(nnx),r(nnx),rr
      complex    g(nnx),xx(nnt),b(nnt),c(nnx),z(nnx),R0C
      real       r0

      
c   Form the peo from the one-step-ahead Pred. filter
      
      g(1)=cmplx(1.,0.0)
      do i=1,ip
      g(i+1)=a(i)
      enddo

c   Autocorrelation of the peo      
      call correlation(ip+1,ip+1,1,g,g,R0C,R)
      R0 = real(R0C)

c   Compute righ-side term 
      perc=0.
      R0=R0*(1.+perc/100.)
      do j=n1,n2
      b(j+1-n1)=cmplx(0.,0.)
      do m=-ip+n1,n1-1
      l=m-j
      if(l.eq.0) rr=cmplx(r0,0.)
      if(l.gt.0) rr=r(m-j)
      if(l.lt.0) rr=conjg(r(j-m))
      b(j+1-n1)=b(j+1-n1)+x(m)*conjg(rr)
      enddo
      enddo

      do j=n1,n2
      c(j+1-n1)=cmplx(0.,0.)
      do m=n2+1,ip+n2
      l=m-j
      if(l.eq.0) rr=cmplx(r0,0.)
      if(l.gt.0) rr=r(m-j)
      if(l.lt.0) rr=conjg(r(j-m))
      c(j+1-n1)=c(j+1-n1)+x(m)*conjg(rr)
      enddo
      enddo

      do i=1,n2-n1+1
      z(i)=-b(i)-c(i)
      enddo

c   Solve for the gap

      call toeplitz(n2-n1,r0,r,Z,xx,ISTAT)

c   Fill the gap

      do i=n1,n2
      x(i)=xx(i-n1+1)
      enddo
      return
      end 



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
 




      SUBROUTINE CORRELATION (N,LAG,MODE,X,Y,R0,R)
C
C   This program computes either the unbiased or biased complex correlation
C   estimates between complex data sample arrays X and Y.  If X=Y, then the
C   autocorrelation is computed.
C
C   Input Parameters:
C
C     N    - Number of data samples in arrays X and Y (integer)
C     LAG  - Number of correlation lags to compute [ lags from  0  to  LAG
C            are computed and stored in R0 and R(1) to R(LAG) ] (integer)
C     MODE - Set to 0 for unbiased correlation estimates; otherwise, biased
C            correlation estimates are computed (integer)
C     X    - Array of complex data samples X(1) through X(N)
C     Y    - Array of complex data samples Y(1) through Y(N)
C
C   Output Parameters:
C
C     R0   - Complex correlation estimate for lag 0
C     R    - Array of complex correlation estimates for lags 1 to LAG
C
C   Notes:
C
C
C     External arrays X,Y must be dimensioned .GE. N and array R must be
C     dimensioned .GE. LAG in the calling program.
C
      COMPLEX X(1),Y(1),R(1),R0,SUM
      DO 30 K=0,LAG
        NK=N-K
        SUM=(0.,0.)
        DO 10 J=1,NK
10        SUM=SUM+X(J+K)*CONJG(Y(J))
        IF (K .NE. 0)  GO TO 20
        R0=SUM/FLOAT(N)
        GO TO 30
20      IF (MODE .EQ. 0)  R(K)=SUM/FLOAT(N-K)    ! Eq.  (5.9)
        IF (MODE .NE. 0)  R(K)=SUM/FLOAT(N)      ! Eqs. (5.13),(5.19)
30      CONTINUE
      RETURN
      END


c
      SUBROUTINE Toeplitz(M,T0,T,Z,X,ISTAT)
C
C   Solves the set of complex linear simultaneous equations
C                           TX = Z
C   by a variation of the Levinson algorithm.  T is a complex  M+1
C   by  M+1  Hermitian Toeplitz matrix, Z is the known right-hand-
C   side complex column vector of  M+1  elements,  and  X  is the
C   solution vector of  M+1  complex elements.
C
C   Input Parameters:
C
C     M  - Order of matrix T (integer)
C     T0 - Scalar corresponding to real matrix element t(0)
C          (This element must be real due to Hermitian symmetry.)
C     T  - Array of  M  complex matrix elements t(1),...,t(M)
C          from the left column of the Toeplitz matrix
C     Z  - Array of M+1 complex elements of the right-hand-side
C          vector.   Program element Z(k+1) corresponds to text
C          element z(k), for k=0 to k=M
C
C   Output Parameters:
C
C     X  - Array of  M+1  complex elements of solution vector.
C          Program element X(k+1) corresponds to text element
C          x(k), for k=0 to k=M
C     ISTAT - Integer status indicator at time of exit
C             0 for normal exit
C             1 if P=0. (singular matrix)
C
C   Notes:
C
C   External array T must be dimensioned .GE. M and arrays X,Z must
C   be dimensioned .GE. M+1 in the calling program.  Internal array
C   A must be dimensioned .GE. M .
C
      COMPLEX T(1),X(1),Z(1),A(100)
      COMPLEX TEMP,SAVE,ALPHA,BETA
      REAL P,T0
      P=T0
      ISTAT=1
      IF (P .EQ. 0.)  RETURN
C   Handle  M=0  as a special case
      X(1)=Z(1)/T0                            ! Eq. (3.175)
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
        SAVE=SAVE+A(J)*T(K-J)                 ! Eq. (3.136)
10      BETA=BETA+X(J+1)*T(K-J)               ! Eq. (3.173)
20    TEMP=-SAVE/P
      P=P*(1.-REAL(TEMP)**2-AIMAG(TEMP)**2)   ! Eq. (3.158)
      IF (P .LE. 0.)  RETURN
30    A(K)=TEMP                               ! Eq. (3.139)
      ALPHA=(Z(K+1)-BETA)/P                   ! Eq. (3.174)
      IF (K .EQ. 1)  GO TO 50
      KHALF=K/2
      DO 40 J=1,KHALF
        KJ=K-J
        SAVE=A(J)
        A(J)=SAVE+TEMP*CONJG(A(KJ))           ! Eqs. (3.147),(3.157)
        IF (J .EQ. KJ)  GO TO 40
        A(KJ)=A(KJ)+TEMP*CONJG(SAVE)          ! Eqs. (3.147),(3.157)
40      CONTINUE
50    X(K+1)=ALPHA
      DO 60 J=1,K
60      X(J)=X(J)+ALPHA*CONJG(A(K-J+1))       ! Eq. (3.171)
      IF (K .LT. M)  GO TO 100
      ISTAT=0
      RETURN
      END





      subroutine fork(lx,x,signi)
C
      complex carg,cexp,cw,ctemp,x(2048)
      j=1
      sc=(1./lx)
      do 30 i=1,lx
      if(i.gt.j) go to 10
      ctemp=x(j)*sc
      x(j)=x(i)*sc
      x(i)=ctemp
10    m=lx/2
20    if(j.le.m) go to 30
      j=j-m
      m=m/2
      if(m.ge.1) go to 20
30    j=j+m
      l=1
40    istep=2*l
      do 50 m=1,l
      carg=(0.,1.)*(3.14159265*signi*(m-1))/l
      cw=cexp(carg)
      do 50 i=m,lx,istep
      ctemp=cw*x(i+l)
      x(i+l)=x(i)-ctemp
50    x(i)=x(i)+ctemp
      l=istep
      if(l.lt.lx) go to 40
      if(signi.eq.1.) goto 60
      do 70 i=1,lx
70    x(i)=x(i)*lx
60    return
      end


      subroutine write_bin(filename,nt,nx,x)
      parameter  (nnt=2048,nnx=128)
      real       x(nnt,nnx)
      character * 20  filename
      open(unit=11,file=filename,access='direct',recl=4*nt)
      do ix=1,nx
       write(11,rec=ix,err=2000)(x(it,ix),it=1,nt)
        enddo
         close(11)
2000  continue
 
        return
         end

c     -----------------------------------------------------------------

      subroutine read_bin(filename,nt,nx,x)
 
      parameter  (nnt=2048,nnx=128)
      real       x(nnt,nnx)
      character * 20 filename
 
      open(unit=12,file=filename,access='direct',recl=4*nt)
 
       do ix=1,nx
        read(12,rec=ix,err=2000) (x(it,ix),it=1,nt)
          enddo
 
2000  continue
         close(12)
 
      return
      end

c     -----------------------------------------------------------------
 



        subroutine FX_go (nx,nt,nfft,stx,sfx,index)

c------ index=-1 TX ---> FX
c------ index= 1 FX----> TX

        parameter  (nnt=2048,nnx=128)
        complex         sfx(nnt,nnx), aux(nnt)
        real            stx(nnt,nnx)
        if(index.eq.-1)  then
        do 140 ix=1,nx
        do 130 it=1,nt
130     aux(it)=cmplx(stx(it,ix),0.0)
        do 120 it=nt+1,nfft
120     aux(it)=cmplx(0.0, 0.0)
        call fork(nfft,aux,-1.)
        do 110 if=1,nfft
110     sfx(if,ix)=aux(if)
140     continue
        endif
        if(index.eq.1) then
        do 240 ix=1,nx
        do 220 if=1,nfft
 
220     aux(if)=sfx(if,ix)
        call fork(nfft,aux, 1.)
        do 230 it=1,nfft
230     stx(it,ix)=real(aux(it))
240     continue
        endif
        return
        end
 
 
c       -------------------------------------------
 
        subroutine FX_symmetry (sfx,nx,nt,nfft,iflow,ifhigh)

        parameter       (nnt=2048,nnx=128)
        complex          sfx(nnt,nnx)
        do 200 ix=1,nx
        do 180 if=1,iflow
180     sfx(if,ix)=cmplx(0.0, 0.0)
        do 190 if=ifhigh+1,nfft/2+1
190     sfx(if,ix)=cmplx(0.0, 0.0)
200     continue
        do 210 ix=1,nx
        do 210 if=nfft/2+2,nfft
210     sfx(if,ix) = conjg(sfx(nfft-if+2,ix))
        return
        end
 




c-----------------------------------------------------------------



        subroutine nmoc(pos,data,nt,nh,dt,t0,dh,
     #          near_offset,tv,v,nv,index)
       
c------ index .eq. 1  ==> data  --->data_c
c------ index .eq.-1  ==> data_c--->data

        parameter (nnt=2048,nnx=128)

        real  data(nnt,nnx)
        real  data_c(nnt,nnx),v(10),tv(10)
        real  t0,t,dt,dh,tau,h,near_offset,vi(nnt)
        real  pos(nnx)

        if(index.ne.1.and.index.ne.-1) then
        write(*,*)'Check index in subr. nmoc'
        return
        endif



        it0=ifix(t0/dt)
        do iv=1,nv-1
        it1=ifix(tv(iv)/dt)+1
        it2=ifix(tv(iv+1)/dt)+1
        do i=it1,it2+1
        vi(i-it0)=v(iv)+(v(iv+1)-v(iv))*(i-it1)/(it2-it1)
        enddo
        enddo


        do it=1,nt
        do ih=1,nh
        data_c(it,ih)=0.0
        enddo
        enddo

        if(index.eq.1) then

        do it=2,nt
        tau=(it-1)*dt+t0
        do  ih=1,nh
        h=pos(ih)
        t=sqrt(tau**2+(h**2)/vi(it)**2)-t0
        t=t/dt
        i1=ifix(t)
        i2=i1+1
        data_c(it-1,ih)=data(i1,ih)+
     #     (data(i2,ih)-data(i1,ih))*(t-float(i1))/float(i2-i1)
        enddo
        enddo

        endif

        if(index.eq.-1) then

        do it=2,nt
        tau=(it-1)*dt+t0
        do  ih=1,nh
        h=pos(ih)
        t=tau**2-(h**2)/vi(it)**2
        if(t.ge.0.d0) then
        t=sqrt(t)-t0
        t=t/dt
        i1=ifix(t)
        i2=i1+1
        data_c(it-1,ih)=data(i1,ih)+
     #     (data(i2,ih)-data(i1,ih))*(t-float(i1))/float(i2-i1)
        endif
        enddo
        enddo

        endif

         do it=1,nt
         do ih=1,nh
         data(it,ih)=data_c(it,ih)
         enddo
         enddo

        return
        end




 
 
      function power(data,nt,nx,n1,n2)
      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)


           p=0.
               do it=1,nt
                 do ix=n1,n2
                  p=p+data(it,ix)**2
                   enddo
                    enddo
               power=p/(nt*(n2-n1+1))

               return
               end

                  

      subroutine smooth(data,nt,nx,n1,n2)
      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)
      real       aux(nnt,nnx)


        do it=1,nt
         do ix=n1,n2
          aux(it,ix)=data(it,ix)+data(it,ix-1)+data(it,ix+1)
           enddo
           enddo

        do it=1,nt
         do ix=n1,n2
          data(it,ix)=aux(it,ix)/3.
           enddo
           enddo

               return
               end



      subroutine LMSC(d,nd,mu,w,nw)

      parameter  (nnt=2048,nnx=128)

      complex   w(nnx),dh(nnt),d(nnt)
      real      mu

      do i=2,nw                     !Initail filter
       w(i)=0.
        enddo
         w(1)=1.


      do k=nw+1,nd

      dh(k)=0.
       do m=1,nw                    !forward prediction stage
         dh(k)=dh(k)+w(m)*d(k-m)
           enddo

           do m=1,nw
            w(m)=w(m)+2.*mu*(d(k)-dh(k))*conjg(d(k-m))
             enddo

             enddo

 
           do i=1,nw               ! filter in PEO form
            w(i)=-w(i)
             enddo 

           return
          end
       

      
