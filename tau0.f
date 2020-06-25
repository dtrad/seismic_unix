

c  Tau-v 
c  M.D.Sacchi, Dept. of Physics, UofA...

c  Inversion of Velocity stacks via conjugate gradients.
  
c  This is a tau-p with hyperbolic travel-time curve
	
	include "c:\mauricio\flibrary\read_bin_d2s.for"
	parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)
	real      x(nnx),y(nnx),h(n1x),t(n1t),p(n1x)
     	character * 30 filein

c   Read parameters form file 11

      open (11,file='tau_pp.in')      ! file with parameters 

      read(11,1) filein               ! binary file with seismic data
      read(11,*) nt                   ! number of samples per trace
      read(11,*) dt                   ! sampling interval
      read(11,*) t0                   ! time in sec of the first sample
      read(11,*) nh                   ! number of traces to read 
      read(11,*) dh                   ! offset interval 
      read(11,*) h_near               ! offset of the first trace
      read(11,*) np                   ! number of velocities to estimate
      read(11,*) p0                   ! minimum velocity m/s
      read(11,*) pmax                 ! Maximum velocity
1     format(a20)

c   Offset

      do i=1,nh
      h(i) = h_near +(i-1)*dh
      enddo

c   Velocity axis 

      dp = (pmax-p0)/(np-1)
      do i=1,np
      p(i) = p0+(i-1)*dp
      enddo
       
c   Tiem axis 

      do i=1,nt
      t(i)=t0+(i-1)*dt
      enddo

c  Read binary data 

c     call read_bin (filein,nt,nh,y)
	call read_bin_d2s(filein,nt,nh,y)	

c  Write data into inp; only the first nh trace   

      call write_bin('inp',nt,nh,y)
      call write_ascii('inp.dat',nt,nh,y) 

c  Compute a velocity panel via the adjoint

      call slant(dt,t,nt,p,np,h,nh,x,y,1)
      call write_bin('adj',nt,np,x)


c  Invert the velocity panel using conjugate gradients 

      call ls(dt,t,nt,p,np,h,nh,x,y)
      call write_bin('inv',nt,np,x)


c  Compute the predicted data 

      call slant(dt,t,nt,p,np,h,nh,x,y,0)
      call write_bin('pred',nt,nh,y)

      stop
      end


c     --------------------------------------------


      subroutine slant(dt,t,nt,p,np,h,nh,m,d,conj)

c  Compute velocity panels  when conj = 0
c  Compute CMP gathers when      conj = 1 

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

      real        d(nnx)
      real        m(nnx),h(n1x)
      real        t(n1t),p(n1x)

      integer conj

      if(conj.eq.1) call clean(m,np*nt)
      if(conj.eq.0) call clean(d,nh*nt)
      
       do ih = 1,nh

       do ip=1,np

       do itau=1,nt

       k=(ip-1)*nt+itau

       time=sqrt (t(itau)**2+(h(ih)/p(ip))**2)
       it=1+0.5+time/dt
       j=(ih-1)*nt+it
       if(it.le.nt.and.it.ge.1) then
        if(conj.eq.1) m(k)=m(k)+d(j)
        if(conj.eq.0) d(j)=d(j)+m(k)
       endif

       enddo
       enddo

      enddo


      return
      end

c     ----------------------------------------------

      function dot(n,x,y)

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

c     Compute the inner product
c     dot=(x,y)

         real  x(nnx), y(nnx)
         real  dot
         real  val

         val=0.0d0
         do i=1,n 
         val = val + x(i)*y(i)
         enddo

         dot=val

         return
         end


c     -------------------------------------------
      subroutine read_bin(filename,nt,nx,x)

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

      real   x(nnx)
      real   d(n1x,n1t)
      character * 20 filename

      open(unit=10,file=filename,access='direct',recl=4*nt)

            do ix=1,nx
             read(10,rec=ix,err=2000) (d(ix,it),it=1,nt)
              enddo

         do ix=1,nx
          do it=1,nt
           k=(ix-1)*nt+it
            x(k)=d(ix,it)
             enddo
              enddo

2000  continue
         close(10)

      return
      end 


c     -------------------------------------------
      subroutine write_bin(filename,nt,nx,x)

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

      real      x(nnx)
      real      d(n1x,n1t)
      character * 20 filename


      open(unit=12,file=filename,access='direct',recl=4*nt)


         do ix=1,nx
          do it=1,nt
           k=(ix-1)*nt+it
            d(ix,it)=real(x(k))
             enddo
              enddo


            do ix=1,nx
             write(12,rec=ix,err=2000) (d(ix,it),it=1,nt)
              enddo

2000  continue
       close(12)

      return
      end 
c-----------------------------------
      subroutine write_ascii(filename,nt,nx,x)

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

      real      x(nnx)
      real      d(n1x,n1t)
      character * 20 filename


      open(unit=15,file=filename)


         do ix=1,nx
          do it=1,nt
           k=(ix-1)*nt+it
            d(ix,it)=real(x(k))
             enddo
              enddo

	   
            do it=1,nt
             write(15,'(40f12.4)') (d(ix,it),ix=1,nx)
              enddo

2000  continue
       close(15)

      return
      end 


c     -------------------------------------------
      subroutine clean(x,nx)

      parameter (n1x=1600,n1t=2048,nnx=n1t*n1x)

      real        x(nnx)

      do i=1,nx
      x(i) = 0.0
      enddo
   
      return
      end


c     -------------------------------------------
      subroutine ls(dt,t,nt,p,np,h,nh,x,y)
      parameter (n1x=1000,n1t=2048,nnx=n1t*n1x)

      real      tol
      real      x(nnx), y(nnx), r(nnx)
      real      g(nnx), s(nnx), ss(nnx)
      real      alpha, beta, gamma, gammam, dot
      real      rms,e,den

      real      t(n1t),p(n1x),h(n1x)

      nx=nt*np
      ny=nt*nh 

      itmax=6
      tol=1.e-7

      do i=1,nx
      x(i)=0.
      enddo


      do i=1,ny
      r(i)=y(i)
      enddo

      call slant(dt,t,nt,p,np,h,nh,g,r,1)

      do i=1,nx
      s(i)=g(i)
      enddo

      gammam=dot(nx,g,g)

      do iter=1,itmax

      call slant(dt,t,nt,p,np,h,nh,s,ss,0)

           den=dot(ny,ss,ss)
           alpha=gammam/den

           do i=1,nx
           x(i)=x(i)+alpha*s(i)
           enddo

           do i=1,ny
           r(i)=r(i)-alpha*ss(i)
           enddo

      call slant(dt,t,nt,p,np,h,nh,g,r,1)

           gamma=dot(nx,g,g)
           beta=gamma/gammam
           gammam=gamma
           do i=1,nx
           s(i)=g(i)+beta*s(i)
           enddo

           e=dot(ny,r,r)

           rms=sqrt(e/ny)

           write(*,*) 'CGLS rms:',rms

           enddo

           return
           end










