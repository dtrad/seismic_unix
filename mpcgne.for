      program mpcgne0
 
      implicit none
      integer nmax, mmax, itercgmax, n, m
      parameter (nmax=100,mmax=100, itercgmax=100)
      parameter (n=5,m=6)

      real*8 A(nmax,mmax), MI(mmax,mmax), x(m), b(n)
      integer i,j,itercg,restart,nn,mm
      real*8 tol,step
      
      tol=1e-5
      restart=0
      step=1
      itercg=m
c     print,* 'pppppppp'     
      open(10,file='input.txt',form='formatted',status='unknown')
c     open(10,file='./pppp')
c      read(10,*) nn,mm
c      write(*,*) nn,mm
      do i=1,n
         write(*,*) i
        read(10,*) (A(i,j), j=1,m)
        write(*,*) (A(i,j), j=1,m)   
      enddo

      do i=1,n
         do j=1,m
            MI(i,j)=0
         enddo
         MI(i,i)=1
      enddo

      do i=1,n
           read(10,*) b(j)  
           write(*,*) b(j)
      enddo

      call mpcgne(A,n,m,x,b,MI,tol,step,itercg,restart)
      
    
      close(10)
      stop
      end

      subroutine mpcgne(A,n,m,x,b,MI,tol,step,itercg,restart)
      implicit none
      integer nmax, mmax, itercgmax
      parameter (nmax=100,mmax=100, itercgmax=100)
      integer i,j,k,adj,n,m,  itercg, restart
      real*8 A(nmax,mmax), MI(mmax,mmax), x(m), b(n)
      real*8 u(nmax),g(mmax),s(nmax)
      real*8 r(nmax),w(nmax)
      real*8 rho(itercgmax),eta(itercgmax),rhold,tol
      real*8 alphanum, alphaden, alpha, dot
      real*8 betanum, betaden, beta, step,  Jcost
      real*8 normb
c     Temp pointers
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
         call Atimesx(g,MI,g,1,m,m)
         call Atimesx(w,A,g,0,n,m) 
    
         alphanum=dot(ny,r,r) 
         alphaden=dot(ny,w,s) 
         alpha=alphanum/alphaden 
    
         if (alphaden .lt. 0.) write(*,*) 'alphaden=',alphaden 

         if (alphaden .lt. 1e+10 ) then 
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
 
         write(*,*) 'rho(',k,')=',rho(k),'beta=',beta     
         do i=1,ny 
            s(i)=r(i)+beta*s(i) 
         enddo
   
         eta(k)=dot(ny,u,u) 
      
      enddo


c     Because the number of iterations is acting as regularization 
c     parameter the cost function is || AM^{-1}m-b||^2
      
      Jcost=dot(ny,r,r)   
      call Atimesx(u,A,x,1,n,m) 
      call Atimesx(x,MI,x,1,m,m) 

      write(*,*) 'iter=',k,'Jcost=',Jcost       

      return
      end

      subroutine Atimesx(b,A,x,adj,n,m)
      implicit none
      integer nmax, mmax
      parameter (nmax=1000, mmax=1000)
      integer i,j,adj,n,m
      real*8 A(nmax,mmax), x(m), b(n)


      if (adj.eq.0) then
         do i=1,n 
            b(i)=0
         enddo
      else
         do i=1,m 
            x(i)=0
         enddo
      endif
      
      
      do i=0,n
         do j=0,m
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




