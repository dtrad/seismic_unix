c
c Fast Radon Parabolic Transform
c
c -------------------------------------------------
c This code computes the high resolution parabolic
c Radon transform using cg+fft
c -------------------------------------------------
c
c M.D.Sacchi, Dep. of Physics, UofA.
c 
c -------------------------------------------------

         parameter      (n1x=1024)

         real           m(n1x,n1x)
         real           d(n1x,n1x)
         real           w(n1x)
         integer        iqw(10)

         character      * 30 fileout1
         character      * 30 fileout2
         character      * 30 fileout3
         character      * 30 fileout4

         common/data/iqw

         open (11,file='input')

         read(11,*) nh
         read(11,*) io
         read(11,*) tol
         read(11,111) fileout1
         read(11,111) fileout2
         read(11,111) fileout3
         read(11,111) fileout4
111      format(a30)

         dh = 5.
         hmin = 0.
         hmax = hmin+(nh-1)*dh

         nt = 512 
         dt = 4./1000.
         tmin = 0.

         nq = nh
         qmin = 0.
         qmax =250.*dt
         dq = (qmax-qmin)/(nq-1)

         write(*,*) qmin,dq
              
         do j=1,nt
          do i=1,nq
           m(i,j) = 0.
          enddo
         enddo


         call ricker(w,nw,15.,dt)

         iq = 20
         itau=80
         do i=1,nw
         m(iq,itau+i+nw/2) = m(iq,itau+i+nw/2)+w(i)
         enddo

         iq = 20
         itau=110
         do i=1,nw
         m(iq,itau+i+nw/2) = m(iq,itau+i+nw/2)+w(i)
         enddo

         iq = 30
         itau=110
         do i=1,nw
         m(iq,itau+i+nw/2) = m(iq,itau+i+nw/2)+w(i)
         enddo

         iq = 20
         itau=125
         do i=1,nw
         m(iq,itau+i+nw/2) = m(iq,itau+i+nw/2)-w(i)
         enddo

         call write_bin(n1x,n1x,fileout1,nq,nt,m)

         f1 = 2.
         f2 = 90. 

         call RT(d,nh,nt,dh,dt,hmin,tmin,m,nq,dq,qmin,1,f1,f2)
         CALL WRITe_bin(n1x,n1x,fileout2,nh,nt,d)

         t1 = secnds(0.0) 
         call RT(d,nh,nt,dh,dt,hmin,tmin,m,nq,dq,qmin,-1,f1,f2)
         t2 = secnds(0.0) 
         write(*,*) 'Time of adjoint operator:', t2-t1
         CALL WRITe_bin(n1x,n1x,fileout3,nq,nt,m)

         t1 = secnds(0.0)
         call RI(nh,nt,dh,dt,hmin,tmin,m,nq,dq,qmin,f1,f2,io,tol)
         t2 = secnds(0.0) 
         write(*,*) 'Time of inversion :', t2-t1
         call write_bin(n1x,n1x,fileout4,nq,nt,m)

         stop
         end


c ---------------------------------------------------------------------------


         subroutine RT(d,nh,nt,dh,dt,hmin,tmin,m,nq,dq,qmin,mode,
     #                  f1,f2)

         parameter      (n1x=1024)

         real           m(n1x,n1x)
         real           d(n1x,n1x)
         complex        mc(n1x,n1x)
         complex        dc(n1x,n1x)
         complex        aux(n1x),arg
         complex        wf(n1x),wi(n1x)

         call  prefft (nt, 0,nexp,wf)
         call  prefft (nt, 1,nexp,wi)

         pi = 4.*atan(1.)

         hmax = hmin+(nh-1)*dh
         if1 = f1*nt*dt+1
         if2 = f2*nt*dt+1

         if(if1.lt.1)      pause 'check f1'
         if(if2.gt.nt/2+1) pause 'check f2'
         if(if1.ge.if2)    pause 'check f1 and f2'

          do j=1,if1-1
           do i=1,nq
            mc(i,j)=cmplx(0.,0.)
           enddo
           do i=1,nh
            dc(i,j)=cmplx(0.,0.)
           enddo
          enddo

          do j=if2+1,nt/2+1
           do i=1,nq
            mc(i,j)=cmplx(0.,0.)
           enddo
           do i=1,nh
            dc(i,j)=cmplx(0.,0.)
           enddo
          enddo

          if(mode.eq.1) then 

          do i = 1,nq 
           do j = 1,nt
            aux(j) = cmplx(m(i,j),0.)
           enddo
            call fft (nt, 0,nexp,wf,aux)
           do j = if1,if2
            mc(i,j) = aux(j)
           enddo
          enddo

         do j = if1,if2
            w=2.*pi*(j-1)/(nt*dt)
          do k = 1,nh
            h = hmin+(k-1)*dh
            dc(k,j) = cmplx(0.,0.)
           do i = 1,nq
            q = qmin+(i-1)*dq
            arg = cmplx(0.,-w*q*(h/hmax)**2)
            dc(k,j) = dc(k,j) + mc(i,j)*cexp(arg)
           enddo 
          enddo 
         enddo 

         do j = 1,nt/2
          do k = 1,nh
           dc(k,nt-j+1) = conjg(dc(k,j+1))
          enddo
         enddo

         do i = 1,nh 
          do j = 1,nt
           aux(j) = dc(i,j)
          enddo
           call fft (nt, 1,nexp,wi,aux)
          do j = 1,nt
           d(i,j) = real(aux(j))
          enddo
         enddo

         else

         do i = 1,nh
          do j = 1,nt
           aux(j) = cmplx(d(i,j),0.)
          enddo
           call fft (nt, 0,nexp,wf,aux)
          do j = if1,if2
           dc(i,j) = aux(j)
          enddo
         enddo

       do j = if1,if2
           w=2.*pi*(j-1)/(nt*dt)
        do k = 1,nq
           q = qmin+(k-1)*dq
           mc(k,j) = cmplx(0.,0.)
         do i = 1,nh
           h = hmin+(i-1)*dh
           arg = cmplx(0., w*q*(h/hmax)**2)
          mc(k,j) = mc(k,j) + dc(i,j)*cexp(arg)
         enddo 
        enddo 
       enddo 

       do j = 1,nt/2
        do k = 1,nq
         mc(k,nt-j+1) = conjg(mc(k,j+1))
        enddo
       enddo

       do i = 1,nq 
        do j = 1,nt
         aux(j) = mc(i,j)
        enddo
         call fft (nt, 1,nexp,wi,aux)
        do j = 1,nt
         m(i,j) = real(aux(j))
        enddo
       enddo

       endif


      return
      end

c     -------------------------------------

      subroutine  circ_mult(r,n,f,g,nf,d)

c     Perform the multiplication of a Toeplitz form with a vector.
c     The procedure is done using a circulant matrix derived from
c     the first row of the Toeplitz form. The Toeplitz form, in fact,
c     is not needed. Everything you need is the firt row of the matrix.

c     INPUT  r(n): first row of the autocorrelation matrix (Teoplitz)
c            f(n): a filter
c              nf:   number of freq. samples to perform the circular conv.

c     OUTPUT g(n): the product of the toeplitz form with f

c     USE:   fork: the old fft by Claerbout.

c    NOTE   that T(1) must be real in order to  have an Hermitian 
c           form T

      parameter (n1x=1024)

      complex r(n1x),f(n1x),g(n1x),d(n1x)
      complex rc(n1x),fc(n1x),gc(n1x)
      complex zero,wf(n1x),wi(n1x),aux
           

       zero=cmplx(0.,0.)

      goto 111
      do i=1,n
       g(i) = zero
        do j=1,n
          if(i.ge.j) aux=f(j)*r(iabs(i-j+1))
          if(i.lt.j) aux=f(j)*conjg(r(iabs(j-i+1)))
         g(i)=g(i)+aux
        enddo
       enddo

       do i=1,n
        g(i) = g(i)+d(i)*f(i)
       enddo

        RETURN

111      continue

       if(nf.lt.2*n) pause 'Check nf in subr. CIRC_MULT'

       call prefft (nf, 0,nexp,wf)
       call prefft (nf, 1,nexp,wi)

       do i=n+1,nf         ! pad with zeros the autocorrelation
        r(i)=zero
       enddo

       do i=1,n-1          ! make things periodic to use the DFT
        r(nf-i+1)=conjg(r(i+1))
       enddo
   
       do i=1,nf             
        rc(i)=r(i)  
       enddo
               
       do i=1,n
        fc(i)=f(i)
       enddo

       do i=n+1,nf
        fc(i)=zero
       enddo

       call fft (nf, 0,nexp,wf,fc)
       call fft (nf, 0,nexp,wf,rc)

       do i=1,nf    
        gc(i)=rc(i)*fc(i)
       enddo

       call fft (nf, 1,nexp,wi,gc)

       do i=1,n
        g(i)=gc(i)            ! truncation - undo the zero padding
       enddo
             
       do i=1,n               ! add the diagonal term
        g(i)=g(i)+d(i)*f(i)
       enddo
              
        return
        end

c     -------------------------------------------

      function dot(n,x,y)

c    This is used by the cg program
c    Compute the inner product
c    dot=(x,y)

       parameter (nnx=1024)

       complex   x(nnx), y(nnx)
       complex   val ,dot

       val=cmplx(0.,0)
       do i=1,n 
        val = val + conjg(x(i))*y(i)
       enddo
       dot=val

       return
       end

c     -------------------------------------------
c     CG - SOLVERS  
c     -------------------------------------------

      subroutine cg(n,T,diag,x,b,nf,niter,rms,iter,tol)

c     Solution of a Teoplitz system using CG and Circulants.
c     Multiplications are done via FFT. The sub. circ_mult does
c     the multiplication (circular) of the Toeplitz form with
c     a vector. The CG uses one circular mult. per iteration. For
c     large problems these may lead to algorithm that is more
c     rapid than Levinson recursion.


c     INPUT:     n: maximum autocorrelation lag. 
c                T(n): first row of the Toeplitz matrix
c                x(n): initial solution to the system of eq. Tx=b    
c                      if you doubt set x=0.
c                b(n): righ hand term. Crosscorrelation term in Wienner filters
c                nf  :number of FFT samples for the circular multiplication.
c                niter : maximum number of iterations. niter.leq.nx
c                
c    OUTPUT:    x(n): is the filter or vector or unknowns 
c               rms : error at the final iteration
c               iter: final iteration number. At this point the rms error
c                     is below a preassigned tolerance

c    USE:       circ_mult: the subr. to do the circular  multiplication
c               dot      : a fucntion to compute inner products.


c     There are some usefull parameters internally defined. i.e. tol
c     this is tol for the rms, if the rms.le.tol, the cg stops an
c     returns the solution at the iteration iter.
c     eps is a regularization term for one of the divisions in the alogorithm


      parameter (n1x=1024)
      complex  T(n1x),x(n1x),b(n1x)
      complex  r(n1x),g(n1x),d(n1x),diag(n1x)
      real     rms,beta,r1,r3
      complex  r2,alpha,dot


      if(aimag(T(1)).ne.0.0) pause 'T(1) must be real'

      call  circ_mult(T,n,x,r,nf,diag)   ! initial sol. is used 
                                         ! to get an initial error
       do i=1,n
        r(i)=b(i)-r(i)
         enddo
          r1=real(dot(n,r,r))
       beta=0.

       do iter=1,niter                ! start the cg

         do i=1,n
          d(i)=r(i)+beta*d(i)
           enddo

           call circ_mult(T,n,d,g,nf,diag)

          r2=dot(n,d,g)

           alpha=r1/r2
             do i=1,n
              x(i)=x(i)+alpha*d(i)
               r(i)=r(i)-alpha*g(i)
                enddo

           r3=real(dot(n,r,r))
            beta=r3/r1
             r1=r3

             rms=sqrt(r3/n)

              if(rms.le.tol) return 

             enddo
             return
              end 


c  -----------------------------------------------------------

      subroutine write_bin(nxmax,nzmax,filename,nx,nz,model)

c Write out a binary, read it with xwigb or ximage (from su)

      real
     :     model(nxmax,nzmax)

      character
     :                filename*30

      NSGI = 1   ! use 4 for Sun/Linux
      open(unit=13,file=filename,access='direct',recl=NSGI*nz)

       do ix=1,nx
        write(13,rec=ix,err=2000)(model(ix,iz),iz=1,nz)
         enddo
        close(13)

2000  continue
      return
      end

      subroutine RI(nh,nt,dh,dt,hmin,tmin,m,nq,dq,qmin,
     #                  f1,f2,io,tol)

       parameter      (n1x=1024)

       real           m(n1x,n1x)
       complex        aux(n1x),arg1,arg2,a(n1x),T(n1x),b(n1x)
       complex        mc(n1x,n1x),diag(n1x)
       complex        wf(n1x),wi(n1x)

c  Prepare sin/cos tables for the fft

       call prefft (nt, 0,nexp,wf)
       call prefft (nt, 1,nexp,wi)

       pi = 4.*atan(1.)

       hmax = hmin+(nh-1)*dh

       if1 = f1*nt*dt+1
       if2 = f2*nt*dt+1

       if(if1.lt.1)      pause 'check f1'
       if(if2.gt.nt/2+1) pause 'check f2'
       if(if1.ge.if2)    pause 'check f1 and f2'

       do j=1,if1-1
        do i=1,nq
         mc(i,j)=cmplx(0.,0.)
        enddo
       enddo

       do j=if2+1,nt/2+1
        do i=1,nq
         mc(i,j)=cmplx(0.,0.)
        enddo
       enddo

c  Map data into fx

       do i = 1,nq 
        do j = 1,nt
         aux(j) = cmplx(m(i,j),0.)
        enddo
         call fft (nt, 0,nexp,wf,aux)
        do j = if1,if2
          mc(i,j) = aux(j)
        enddo
       enddo

c  Start to loop on the freq.

       do j = if1,if2

        w=2.*pi*(j-1)/(nt*dt)

          do i = 1,nq
           b(i) = mc(i,j)
          enddo

c  Define 1 row of the Toeplitz Matrix

          q1 = qmin

         do iq2 = 1,nq

          T(iq2) = cmplx(0.,0.)
           q2 = qmin+(iq2-1)*dq

          do k = 1,nh

            h = hmin+(k-1)*dh
            arg1 = cmplx(0.,-w*q1*(h/hmax)**2)
            arg2 = cmplx(0., w*q2*(h/hmax)**2)
            T(iq2) = T(iq2)+cexp(arg1)*cexp(arg2) 

          enddo
         enddo

c  Solve (T+D)a=b 

        do i=1,nq
         diag(i)=100.0
        enddo

        i1 = (120*dt-qmin)/dq+1.1
        i2 = (140*dt-qmin)/dq+1.1
        diag(20)=0.0001
        diag(30)=0.0001
        epsi = 10.01

        nf = 2*nq
ccccc   ni = nq
        ni = ifix(nq/5)
ccccc

c  Solve using a method selected by io

         if(io.eq.1) call ch(nq,T,diag,a,b,istat)
         if(io.eq.2) call cg(nq,T,diag,a,b,nf,ni,rms,iter,tol)
         if(io.eq.3) call to(nq,T,epsi,a,b,istat)


        do i = 1,nq
         mc(i,j) = a(i)
        enddo

        enddo   ! end freq loop

c   Map data back to tx
c     First impose symmetries 

        do j = 1,nt/2
         do k = 1,nq
          mc(k,nt-j+1) = conjg(mc(k,j+1))
         enddo
        enddo

        do i = 1,nq 
         do j = 1,nt
          aux(j) = mc(i,j)
         enddo
          call fft (nt, 1,nexp,wi,aux)
         do j = 1,nt
          m(i,j) = real(aux(j))
         enddo
        enddo

        return
        end

        subroutine ricker(w,nw,f,dt)
        real  w(1024),pi,alfa,beta,f,dt

        pi=4.0*atan(1.d0)       
        nw = 2.5/f/dt
        nc=nw/2+1
        do 100 i=1,nw
        alfa=(i-nc)*pi*f*dt
        beta=alfa**2
        w(i)=(1.0-2.0*beta)*exp(-beta)
100     continue
        return
        end

C
      subroutine prefft (n,mode,nexp,w)
C
C   Input Parameters:
C
C     N     - Number of data samples to be processed (integer-must be a
C             power of two)
C     MODE  - Set to 0 for discrete-time Fourier series (Eq. 2.C.1) or
C             1 for inverse (Eq. 2.C.2)
C
C   Output Parameters:
C
C     NEXP  - Indicates power-of-2 exponent such that N=2**NEXP .
C             Will be set to -1 to indicate error condition if N
C             is not a power of 2 (this integer used by sub. FFT)
C     W     - Complex exponential array
C
C   Notes:
C
C     External array W must be dimensioned .GE. N by calling program.
C
      complex w(1),c1,c2
      nexp=1
5     nt=2**nexp
      if (nt .ge. n)  go to 10
      nexp=nexp+1
      go to 5
10    if (nt .eq. n)  go to 15
      nexp=-1                 ! error:  n is not a power of 2
      return
15    s=8.*atan(1.)/float(nt)
      c1=cmplx(cos(s),-sin(s))
      if (mode .ne. 0)  c1=conjg(c1)
      c2=(1.,0.)
      do 20 k=1,nt
        w(k)=c2
20      c2=c2*c1
      return
      end


      subroutine fft (n,mode,nexp,w,x)
C
C   Input Parameters:
C
C     N,MODE,NEXP,W - See parameter list for subroutine PREFFT
C     X             - Array of N complex data samples, X(1) to X(N)
C
C   Output Parameters:
C
C     X - N complex transform values replace original data samples
C         indexed from k=1 to k=N, representing the frequencies
C         (k-1)/NT hertz
C
C   Notes:
C
C     External array X must be dimensioned .GE. N by calling program.
C
      parameter(n1x=1024)
      complex x(n1x),w(n1x),c1,c2
      mm=1
      ll=n
      do 70 k=1,nexp
        nn=ll/2
        jj=mm+1
        do 40 i=1,n,ll
          kk=i+nn
          c1=x(i)+x(kk)
          x(kk)=x(i)-x(kk)
40        x(i)=c1
        if (nn .eq. 1) go to 70
        do 60 j=2,nn
          c2=w(jj)
          do 50 i=j,n,ll
            kk=i+nn
            c1=x(i)+x(kk)
            x(kk)=(x(i)-x(kk))*c2
50          x(i)=c1
60        jj=jj+mm
        ll=nn
        mm=mm*2
70      continue
      nv2=n/2
      nm1=n-1
      j=1
      do 90 i=1,nm1
        if (i .ge. j)  go to 80
        c1=x(j)
        x(j)=x(i)
        x(i)=c1
80      k=nv2
85      if (k .ge. j)  go to 90
        j=j-k
        k=k/2
        go to 85
90      j=j+k
      if (mode .eq. 0)  s=1.
      if (mode .ne. 0)  s=1./float(n)
      do 100 i=1,n
100     x(i)=x(i)*s
      return
      end



      subroutine ch(m,T,diag,x,y,istat)

c     solve (T+diag)x=b
c     subroutine cholesky(m,eps,a,b,istat)
c
 
C This program solves a Hermitian symmetric set of complex linear
C simultaneous equations using the Cholesky decomposition method.
C The solution replaces the original contents of array B. Contents
C of array A are destroyed after this routine is called.
 
C                    AX = B
 
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
      parameter (n1x=1024)
      complex a(n1x*n1x),b(n1x),sum
      complex T(n1x),diag(n1x),x(n1x),y(n1x)
C
C   Factor into triangular and diagonal form   !  Eq. (3.76)
C

c**
      eps=1.e-15
      do i=1,m
      b(i) = y(i)
      enddo
c**
c Pack the Toeplitz form in a Upper diag matrix, this
c is what Cholesky needs.
c**
      k=1
      do i=1,m
      do j=1,i
      if(i.ne.j) A(k)=conjg(T(i-j+1))
      if(i.eq.j) A(k)=conjg(T(i-j+1)) + conjg(diag(i))
      k=k+1
      enddo
      enddo
c**

      istat=0
      kpiv=0
      do 100 k=1,m
        kpiv=kpiv+k
        ind=kpiv
        lend=k-1
        tiny=abs(eps*real(a(kpiv)))
        do 100 i=k,m
          sum=(0.,0.)
          if (lend .eq. 0)  go to 40
          lpiv=kpiv
          do 30 l=1,lend
            lpiv=lpiv+l-k-1
30          sum=sum+real(a(lpiv))*a(ind-l)*conjg(a(kpiv-l))
40        sum=a(ind)-sum
          if (i .ne. k)  go to 80
C
C   Test for negative pivot element and loss of significance
C
          if (real(sum) .gt. tiny)  go to 90
          if (real(sum) .gt. 0.)  go to 70
          istat=-1
          return
70        if (istat .gt. 0)  go to 90
          istat=k
90        a(kpiv)=cmplx(real(sum),0.)
          dpiv=1./real(sum)
          go to 100
80        a(ind)=sum*dpiv
100       ind=ind+i
C
C   Back solution for intermediate column vector solution  ! Eq. (3.74)
C
      kpiv=1
      do 200 k=2,m
        kpiv=kpiv+k
        sum=b(k)
        do 210 j=1,k-1
210       sum=sum-b(k-j)*conjg(a(kpiv-j))
200     b(k)=sum
C
C   Back solution for final column vector solution    !  Eq. (3.75)
C
      kpiv=(m*(m+1))/2
      b(m)=b(m)/real(a(kpiv))
      do 300 k=m,2,-1
        kpiv=kpiv-k
        ind=kpiv
        sum=b(k-1)/real(a(kpiv))
        do 310 j=k,m
          ind=ind+(j-1)
310       sum=sum-b(j)*a(ind)
300     b(k-1)=sum

      do i=1,m
      x(i) = b(i)
      enddo

      return
      end


      subroutine to (n,T,epsi,x,z,istat)
C
c   MODIFIED, check t 
C   Solves the set of complex linear simultaneous equations
C             suitably named, .ins file to 
%% control the generation of files from your source file; this file 
%% should contain your own preambles for the files it generates, not 
%% those in the standard .ins files. 
%% 
%% The names of the source files used are shown above. 
%% 
%% 
%% 
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Dou,...,t(M)
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
      parameter (n1x=1024)
      complex t(n1x),x(n1x),z(n1x),a(n1x)
      complex toep(n1x)
      complex temp,save,alpha,beta
      real p,t0,epsi
 
c**
      M=N-1
      t0=real(T(1))+epsi
      do i=1,m
      toep(i) = (T(i+1))
      enddo
c**

      p=t0
      istat=1
      if (p .eq. 0.)  return
C   Handle  M=0  as a special case
      x(1)=z(1)/t0                            ! eq. (3.175)
      if (m .le. 0)  return
C
C   Main recursion
C
      k=0
100   k=k+1
      save=toep(k)
      beta=x(1)*toep(k)
      if (k .eq. 1)  go to 20
      do 10 j=1,k-1
        save=save+a(j)*toep(k-j)                 ! eq. (3.136)
10      beta=beta+x(j+1)*toep(k-j)               ! eq. (3.173)
20    temp=-save/p
      p=p*(1.-real(temp)**2-aimag(temp)**2)   ! eq. (3.158)
      if (p .le. 0.)  return
30    a(k)=temp                               ! eq. (3.139)
      alpha=(z(k+1)-beta)/p                   ! eq. (3.174)
      if (k .eq. 1)  go to 50
      khalf=k/2
      do 40 j=1,khalf
        kj=k-j
        save=a(j)
        a(j)=save+temp*conjg(a(kj))           ! eqs. (3.147),(3.157)
        if (j .eq. kj)  go to 40
        a(kj)=a(kj)+temp*conjg(save)          ! eqs. (3.147),(3.157)
40      continue
50    x(k+1)=alpha
      do 60 j=1,k
60      x(j)=x(j)+alpha*conjg(a(k-j+1))       ! eq. (3.171)
      if (k .lt. m)  go to 100
      istat=0
      return
      end
