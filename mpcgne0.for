      program mpcgne0
c     Test for mpcgne.for. 
c     A matrix A and b are read from input.txt 
c     and solve for x such that
c
c     b=A*x
c
c     MI is the preconditioning matrix with size m x m ( model size squared)
c
c     If the weight is low the
c     model is penalized. MI corresponds to the inverse model Covariance.
c     Daniel Trad - May 2000
      implicit none
      integer nmax, mmax, itercgmax, n, m
      parameter (nmax=1000,mmax=1000, itercgmax=100)
      parameter (n=5,m=6)
      
      real*8 A(nmax,mmax), MI(mmax,mmax), x(m), b(n)
      integer i,j,itercg,restart
      real*8 tol,step
      
      tol=1e-5
      restart=0
      step=1
      itercg=m

      open(10,file='input.txt',status='unknown')

      do i=1,n
          read(10,*) (A(i,j), j=1,m)
      enddo

      do i=1,n
         do j=1,m
           write(*,*) A(i,j)
         enddo
      enddo

      do i=1,m
         do j=1,m
            MI(i,j)=0
         enddo
      enddo

      do i=1,m 
         MI(i,i)=1
      enddo

      MI(6,6)=0

      do i=1,n
          read(10,*) b(i)
          write(*,*) b(i)
      enddo

      call mpcgne(A,n,m,x,b,MI,tol,step,itercg,restart)
      
       do i=1,m
          write(*,*) x(i)
      enddo     
    
      close(10)

      end

      subroutine mpcgne(A,n,m,x,b,MI,tol,step,itercg,restart)
c     MCGNE preconditioned conjugate Gradient Normal error CGNE
c     Reference: Yousef Saad: Iterative method for sparse linear systems
c     Given the problem b= A x
c     It solves  A MI A^T u = b with u= MI A^T u 
c     which corresponds to the middle preconditioned CGNE
c     A A^T u = b with x=A^T u
c     This is useful to modify the nullspace of the transformation
c     as the solution is free from the nullspace of MI A but not A
c     MI is the inverse preconditioner or inverse model covariance
c     Daniel Trad - May 2000 
      implicit none
      integer nmax, mmax, itercgmax
      parameter (nmax=1000,mmax=1000, itercgmax=100)
      integer i,k,n,m,  itercg, restart
      real*8 A(nmax,mmax), MI(mmax,mmax), x(m), b(n)
      real*8 u(nmax),g(mmax),s(nmax)
      real*8 r(nmax),w(nmax),gtemp(mmax),xtemp(mmax)
      real*8 rho(itercgmax),eta(itercgmax),rhold,tol
      real*8 alphanum, alphaden, alpha, dot
      real*8 beta, step,  Jcost
      real*8 normb
c     Temp variables
c     r  Conjugate Gradient (Residual)
c     g  Gradient
c     z  Preconditioned gradient
c     s  Search direction
c     w  Conjugate search direction
c     M  precondtioner on data space
c     Wd model and data weights respectively.

      integer nx,ny
      nx=m
      ny=n
      normb=dot(ny,b,b)

      do i=1,ny 
         u(i)=0. 
         s(i)=b(i)
         r(i)=b(i) 
      enddo
  
      k=0 
      rhold=tol*2 
      do while ((rhold.gt.tol).AND.(k.lt.itercg))
         k=k+1 
         call Atimesx(s,A,g,1,n,m)
         call Atimesx(gtemp,MI,g,0,m,m)
         call Atimesx(w,A,gtemp,0,n,m) 
    
         alphanum=dot(ny,r,r) 
         alphaden=dot(ny,w,s) 
         alpha=alphanum/alphaden 
    
c         if (alphaden .lt. 0.) 
c         write(*,*) 'alphaden=',alphaden 

         if (alphaden .lt. 1e-15 ) then 
            write(*,*) 'alphanum=',alphanum,'alphaden=',alphaden,
     #           'alpha=',alpha
         endif
         
         do i=1,ny 
            u(i)=u(i)+alpha*s(i)
            r(i)=r(i)-alpha*w(i)
         enddo  
         
         rho(k)=dot(ny,r,r)
         rho(k)=rho(k)/normb 

         beta=dot(ny,r,r)
         beta=beta/alphanum
 
         write(*,*) 'rho(',k,')=',rho(k)
         do i=1,ny 
            s(i)=r(i)+beta*s(i) 
         enddo
   
         eta(k)=dot(ny,u,u) 
      
      enddo


c     Because the number of iterations is acting as regularization 
c     parameter the cost function is || AM^{-1}m-b||^2
      
      Jcost=dot(ny,r,r)   
      call Atimesx(u,A,xtemp,1,n,m) 
      call Atimesx(x,MI,xtemp,0,m,m) 

      write(*,*) 'iter=',k,'Jcost=',Jcost       

      return
      end

      subroutine Atimesx(b,A,x,adj,n,m)
c     It calculates b = A*x if adjoint == 0 
c     or x = A^T b if adjoint == 1
c     References Claerbout: Processing vs Inversion (WWW), chapter 5
c     Daniel Trad - May 2000
      implicit none
      integer nmax, mmax
      parameter (nmax=1000, mmax=1000)
      integer i,j,adj,n,m
      real*8 A(nmax,mmax), x(*), b(*)
    

      if (adj.eq.0) then
         do i=1,n 
            b(i)=0
         enddo
      else
         do i=1,m 
            x(i)=0
         enddo
      endif
      
      
      do i=1,n
         do j=1,m
            if (adj.eq.0) then
               b(i)=b(i)+A(i,j)*x(j)
            else 
               x(j)=x(j)+A(i,j)*b(i)
            endif   
         enddo    
      enddo

      return
      end

      function dot(n,x,y)
c     Compute the inner product
c     dot=(x,y) for real*8 x,y
      implicit none
      integer n,i
      real*8 x(n), y(n), dot, val
      
      val=0.0d0

      do i=1,n 
         val = val +  x(i) * y(i)
      enddo

      dot=val
      
      return
      end













