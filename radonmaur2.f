        real   dd(512,2048),mm(512,2048)

        character * 20 filename

        real    offset(512),w(2048)
        complex Wf(2048),Wi(2048)

        if1 = 25
        if2 = 105

        nt = 1024
        nf = 1024
        nx = 64
        nq = 64

        dt = 4/1000.
        q0 = 0*dt
        dq = dt*30

        pi = 4.*atan(1.)

        do i=1,nf
        w(i) = (i-1)*2.*pi/(nf-1)/dt
        enddo

        do i=1,nx
        offset(i) = 10. + (i-1)*5.
        enddo

        ofmax = offset(nx)

        do i=1,nx
        offset(i) = (offset(i)/ofmax)**2
        enddo

        do j=1,nt
        do i=1,nq
        mm(i,j) = 0.
        enddo
        enddo
        mm(4,111) =-1. 
        mm(9,218) = 1. 
        mm(8,318) =-2. 
        mm(2,418) = 1. 
        mm(5,618) =-1. 

        call  prefft (nf,0,nexp,wf)
        call  prefft (nf,1,nexp,wi)

c        call timer(it1) 
        call radon(mm,dd,nx,nq,nt,nf,
     #            if1,if2,w,offset,q0,dq,Wf,Wi,nexp,'f')
c        call timer(it2)
c        write(*,*) it2-it1
c        call timer(it1) 

        filename = 'outd1'
        call read_write('w',filename,dd,nx,nt)
        call read_write('w','adj1',mm,nq,nt)
        do j=1,nt
        do i=1,nq
        mm(i,j) = 0.
        enddo
        enddo

        call radon_inv(mm,dd,nx,nq,nt,nf,
     #              if1,if2,w,offset,q0,dq,Wf,Wi,nexp)

        filename = 'out'
        call read_write('w',filename,mm,nq,nt)

        call radon(mm,dd,nx,nq,nt,nf,
     #            if1,if2,w,offset,q0,dq,Wf,Wi,nexp,'f')

        filename = 'outd2'
        call read_write('w',filename,dd,nx,nt)
        stop
        end


        subroutine radon(mm,dd,nx,nq,nt,nf,
     #                 if1,if2,w,offset,q0,dq,Wf,Wi,nexp,flag)

        real    dd(512,2048),mm(512,2048)
        real    offset(512),w(2048)

        complex D(512,2048),M(512,2048)
        complex Wf(2048),Wi(2048)
        
        character * 1 flag

c flag='f' means forward 
        if(flag.eq.'f') then 
         call  tx2fx (mm,M,nq,nt,nf,if1,if2,wf,nexp,'f')
         call  fx2fq (M,D,nx,nq,nf,if1,if2,w,
     #                   q0,dq,offset,'f')
         call  tx2fx (dd,D,nx,nt,nf,if1,if2,wi,nexp,'i')
        endif

            
c flag='f' means adjoint 
        if(flag.eq.'a') then 
         call  tx2fx (dd,D,nx,nt,nf,if1,if2,wf,nexp,'f')
         call  fx2fq (M,D,nx,nq,nf,if1,if2,w,
     #                   q0,dq,offset,'a')
         call  tx2fx (mm,M,nq,nt,nf,if1,if2,wi,nexp,'i')
        endif

        return
        end

c      ------------------------------------------------------------

        subroutine  tx2fx(dd,D,nx,nt,nf,if1,if2,w,nexp,conj)

        parameter   (nxmax=512,ntmax=2048)

        character * 1 conj

        real*4          
     :              dd(nxmax,ntmax)
        complex
     :              D(nxmax,ntmax),
     :              aux(ntmax), w(ntmax)

        if(conj.eq.'f') then 

        do i = 1,nx 

         do j = 1,nt
          aux(j) = cmplx(dd(i,j),0.)
           enddo

         do j = nt+1,nf
          aux(j) = cmplx(0.,0.)
           enddo

        call fft (nf,0,nexp,w,aux)

         do j = 1,if1-1
          D(i,j) = cmplx(0.,0.)
           enddo

         do j = if1,if2
          D(i,j) = aux(j)
           enddo

         do j = if2+1,nf/2+1 
          D(i,j) = cmplx(0.,0.)
           enddo

         do j = nf/2+2,nf
          D(i,j) = conjg(D(i,nf-j+2))
           enddo

        enddo

        else

        do i = 1,nx 
        
        do j = 1,if1-1 
         aux(j) = cmplx(0.,0.)
          enddo
         
         do j = if1,if2 
          aux(j) = D(i,j)
           enddo

         do j = if2+1,nf/2+1 
          aux(j) = cmplx(0.,0.)
           enddo

         do j = nf/2+2,nf
          aux(j) = conjg(aux(nf-j+2))
           enddo

         call fft (nf,1,nexp,w,aux)

         do j = 1,nt
          dd(i,j) = real(aux(j))
           enddo

        enddo

        endif

        return
        end


c      ------------------------------------------------------------


      subroutine prefft (n,mode,nexp,w)
C
C   Input Parameters:
C
C     N     - Number of data samples to be processed (integer-must be a
C             power of two)
C     MODE  - Set to 0 for discrete-time Fourier series  or
C             1 for inverse
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
      complex w(2048),c1,c2
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

c      ------------------------------------------------------------

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
      parameter(n1x=2048)
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

c  -----------------------------------------------------------

        subroutine fx2fq(M,D,nx,nq,nf,if1,if2,w,
     #                   q0,dq,offset,conj)

        complex  M(512,2048),D(512,2048)
        complex  C0, DELTA_C0
        real     offset(512),w(2048)

        character * 1 conj

        if(conj.eq.'a') then

        do j = if1,if2
        do i = 1,nq
        M(i,j) = cmplx(0.,0.)
        enddo
        enddo

        do j = if1,if2
        do k = 1,nx
        C0 = cexp(cmplx(0.,w(j)*offset(k)*(q0-dq)))
        DELTA_C0 = cexp(cmplx(0.,w(j)*offset(k)*dq))
        do i = 1,nq
        C0 = C0*DELTA_C0
        M(i,j) = M(i,j) + D(k,j)*C0
        enddo
        enddo
        enddo

        endif

        if(conj.eq.'f') then

        do j = if1,if2
        do i = 1,nx
        D(i,j) = cmplx(0.,0.)
        enddo
        enddo
    
c       do j = if1,if2
c       do i = 1,nq
c       do k = 1,nx
c       C0 = cexp(cmplx(0.,-w(j)*offset(k)*(q0+dq*(i-1))))
c       D(k,j) = D(k,j) + M(i,j)*C0
c       enddo
c       enddo
c       enddo

        do j = if1,if2
        do k = 1,nx
        C0 = cexp(cmplx(0.,-w(j)*offset(k)*(q0-dq)))
        DELTA_C0 = cexp(cmplx(0.,-w(j)*offset(k)*dq))
        do i = 1,nq
        C0 = C0*DELTA_C0
        D(k,j) = D(k,j) + M(i,j)*C0
        enddo
        enddo
        enddo

        endif

        return
        end


      subroutine read_write(io,filename,d,nx,nt)
c     read  is io='r'
c     write is io='w'
      real d(512,2048)
      character * 20 filename
      character * 1, io
      if (io.eq.'r') then
      open(unit=10,file=filename,access='direct',recl=4*nt)
      do ix=1,nx
      read(10,rec=ix,err=2000) (d(ix,it),it=1,nt)
      enddo
      close(10)
      endif
      if (io.eq.'w') then
      open(unit=10,file=filename,access='direct',recl=4*nt)
      do ix=1,nx
      write(10,rec=ix,err=2000) (d(ix,it),it=1,nt)
      enddo
      close(10)
      endif
2000  continue
      return
      end

        subroutine radon_inv(x,y,nh,np,nt,nf,
     #                 if1,if2,w,offset,q0,dq,Wf,Wi,nexp)

     
        real    y(512,2048),x(512,2048)
        real    ss(512,2048),s(512,2048)
        real    r(512,2048),g(512,2048)
        real    offset(512),w(2048)

        complex Wf(2048),Wi(2048)
        


      itmax=50
 
      do it=1,nt
      do ip=1,np
      x(ip,ih)=0.
      enddo
      enddo
 
 
      do it=1,nt
      do ih=1,nh
      r(ih,it)=y(ih,it)
      enddo
      enddo
 
      call radon(g,r,nh,np,nt,nf,
     #     if1,if2,w,offset,q0,dq,Wf,Wi,nexp,'a')
 
      do it=1,nt
      do ip=1,np
      s(ip,it)=g(ip,it)
      enddo
      enddo
           gammam=0.
           do it=1,nt
           do ip=1,np
           gammam = gammam + g(ip,it)*g(ip,it)
           enddo
           enddo
 
      do iter=1,itmax
 
      call radon(s,ss,nh,np,nt,nf,
     #          if1,if2,w,offset,q0,dq,Wf,Wi,nexp,'f')
 
           do it=1,nt
           do ih=1,nh
           den = den + ss(ih,it)*ss(ih,it)
           enddo
           enddo
 
           alpha=gammam/den
 
           do it=1,nt
           do ip=1,np
           x(ip,it)=x(ip,it)+alpha*s(ip,it)
           enddo
           enddo
 
 
           do it=1,nt
           do ih=1,nh
           r(ih,it)=r(ih,it)-alpha*ss(ih,it)
           enddo
           enddo
 
      call radon(g,r,nh,np,nt,nf,
     #          if1,if2,w,offset,q0,dq,Wf,Wi,nexp,'a')
 
           gamma=0.
           do it=1,nt
           do ip=1,np
           gamma = gamma + g(ip,it)*g(ip,it)
           enddo
           enddo
           beta=gamma/gammam
           gammam=gamma
 
           do it=1,nt
           do ip=1,np
           s(ip,it)=g(ip,it)+beta*s(ip,it)
           enddo
           enddo
 
           e=0.
           do it=1,nt
           do ih=1,nh
           e=e+r(ih,it)*r(ih,it)
           enddo
           enddo
           rms=sqrt(e/(nh*nt))
 
           write(*,*) 'CGLS rms:',rms
           enddo
           return
           end
 

