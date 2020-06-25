      subroutine ls(x,b,n)
c     input 
c     x initial guess
c     b rhs
c     n dimension
c     Define function Atimesx(x,b,n) ==> b=A_(nxn)*x
      implicit none
      integer kmax, k, i, n
      real*8      x(n), b(n), r(n), p(n), w(n)
      real*8      alpha, beta, rho, rho2, dot, tol, bb, den

 
      kmax=n/5 
      tol=1.e-7

      call Atimesx(x,r,n)

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
           call Atimesx(p,w,n)
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













