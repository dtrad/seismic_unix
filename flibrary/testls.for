      program testls
      implicit none
      integer n,i
      parameter (n=3)
      real*8 A(n,n), x(n), b(n), x0(n)
      data A /1.61000,0.79800,0.85000,0.79800,1.11890,0.45000,
     &0.85000,0.45000,1.05000/
      data b  /4.0224,3.6661,2.3700/
      data x0 /1.20000,2.30000,0.30000/
      data x /0, 0, 0/
      call ls(A,x,b,n)
      do i=1,n
       print*,x0(i),x(i),A(i,1),A(i,2),A(i,3),b(i)
      enddo
      end

      subroutine ls(A,x,b,n)
c     input 
c     x initial guess
c     b rhs
c     n dimension
c     Define function Atimesx(x,b,n) ==> b=A_(nxn)*x
      implicit none
      integer kmax, k, i, n, nmax
      parameter (nmax=1000)
      real*8      A(n,n),x(n), b(n), r(nmax), p(nmax), w(nmax)
      real*8      alpha, beta, rho, rho2, dot, tol, bb, den

 
      kmax=n
      tol=1.e-13

      call Atimesx(r,A,x,n,n)

      do i=1,n
        r(i)=b(i)-r(i)
      enddo

      rho=dot(n,r,r)
      bb=sqrt(dot(n,b,b))

      rho=dot(n,r,r)

      do while((sqrt(rho).gt.tol*bb).and.(k.lt.kmax))
           k=k+1
           if (k.eq.1) then
              do i=1,n
                 p(i)=r(i)
              enddo
           else 
              beta=rho/rho2
              do i=1,n
                 p(i)=r(i)+beta*p(i)
              enddo                                         
           endif
           call Atimesx(w,A,p,n,n)
           den=dot(n,p,w)  
           alpha=rho/den

           do i=1,n
           x(i)=x(i)+alpha*p(i)      
           enddo

           do i=1,n
           r(i)=r(i)-alpha*w(i)
           enddo

           rho2=rho
           rho=dot(n,r,r)
        enddo
        return
        end

      function dot(n,x,y)

c     Compute the inner product
c     dot=(x,y)

         real*8  x(n), y(n)
         real*8  dot
         real*8  val

         val=0.0d0
         do i=1,n 
         val = val + x(i)*y(i)
         enddo

         dot=val

         return
         end

      subroutine Atimesx(b,A,x,nr,nc)
      implicit none
      integer i,j, nr, nc
      real*8 A(nr,nc), x(nc), b(nr)

      do i=1,nr
         b(i)=(0.d0,0.d0)
         do j=1,nc
            b(i)=b(i)+A(i,j)*x(j)
         enddo
      enddo
      return
      end











