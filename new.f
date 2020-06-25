c     PRIMAL DUAL LOG-LOG BARRIER LP algorithm
c     =========================================

c     This version solves a LS problem at each barrier
c     iteration using a CG solver   

      parameter   (nnx=30000,nbig=120000)
      real        x(nnx),c(nnx),b(nnx)
      real        gamma,delta,A(nbig)
      real        bp(nnx),wave(21)
      integer     irow(nbig),icol(nbig)

c   set inital variables
c   primal-dual LP log-barrier method

      nh=12
      dh=60.
      nv=12
      dv=10.
      v0=2000.
      t0=0.
      dt=0.002
      nt=200
      h0=10


C  GET MATRIX IN SPARSE FORMAT ...

      call NMO(nh,dh,h0,nt,dt,t0,nv,dv,v0,A,icol,irow,na)

 
      nw = 21
      do i=1,nw
      wave(i) = exp(-6.*(float(i-10)/11.)**2)    
      wave(i) = 1.
      enddo

      np=nv*nt
      nd=nh*nt

      write(*,*) 'np:',np
      write(*,*) 'nd:',nd
      write(*,*)' # of elements:', na

      do i=1,np
      x(i)=0.
      enddo

c   Make a velocity model

      do it=1,nw
      k=(nv/2-1)*nt+nt/3+it
      x(k)=wave(it)
      k=(nv/2+4-1)*nt+nt/3+5+it
      x(k)=wave(it)
      k=(nv/2-4-1)*nt+nt/3-9+it
      x(k)=wave(it)
      enddo

      call matmul0(0,na,A,np,x,nd,b,irow,icol)

       idum=-39
       do i=1,nd
       b(i)=b(i)+gasdev(idum)*0.0
       enddo

      delta=1.e0
      gamma=1.e-7


       open(unit=4,file='data',access='direct',recl=4*nt)
        do ih=1,nh
         write(4,rec=ih,err=2000) (b((ih-1)*nt+it),it=1,nt)
          enddo

     
       open(unit=7,file='in',access='direct',recl=4*nt)
        do iv=1,nv
         write(7,rec=iv,err=2000) (x((iv-1)*nt+it),it=1,nt)
          enddo

      do i=1,np
      x(i)=0.
      c(i)=1.
      enddo
 

      call PD_lp_bm(np,nd,na,c,A,x,b,irow,icol,
     #              iter,delta,gamma,PI,DI,DG)

      write(*,*)' # of iterations of PD_lp: ', iter
      write(*,*)' Primal Infeasibility :', PI
      write(*,*)' Dual Infeasibility   :', DI
      write(*,*)' Dualilty Gap         :', DG

c   Predicted data 

      call  matmul0 (0,na,A,np,x,nd,bp,irow,icol)

      error=sqrt(dot(nd,b,bp)/nd)
      write(*,*) ' MSE:', error

        open(unit=6,file='out',access='direct',recl=4*nt)
        do iv=1,nv
         write(6,rec=iv,err=2000) (x((iv-1)*nt+it),it=1,nt)
          enddo
 

2000  stop
      end


c     -----------------------------------------------------

      subroutine PD_lp_bm(np,nd,na,c,A,x,b,irow,icol,
     #                    iter,delta,gamma,PI,DI,DG)

c     Minimize c'x +gamma^2x'x+p'p
c     subject Ax+delta.p=b, x.ge.0

      parameter   (nnx=30000,iter_max=4,nbig=1200000,
     #             PI_t=1.e-1,DI_t=1.e-1,DG_t=1.e-1)

      real        a(nbig),x(nnx),c(nnx),b(nnx)
      real        D(nnx),v(nnx),t(nnx),r(nnx),z(nnx),y(nnx)
      real        rhs(nnx)
      real        mu,gamma,delta,rhod,rhop,x1,z1
      real        xp(nnx),yp(nnx)
      real        dx(nnx),dy(nnx),dz(nnx)
      integer     irow(nbig),icol(nbig)

c     set inital variables

      do i=1,np
      x(i)=1.
      z(i)=1.
      enddo
      do i=1,nd
      y(i)=1.
      enddo

       mu=5.e-2
       iter=0

       do k=1,iter_max

          iter=iter+1

c--1
           call  matmul0(1,na,A,np,xp,nd,y,irow,icol)

          do i=1,np
          t(i)=c(i)+gamma*gamma*x(i)-z(i)-xp(i)
          v(i)=mu-x(i)*z(i)
          d(i)=1./(z(i)/x(i)+gamma*gamma)
          enddo

c--2
       call  matmul0(0,na,A,np,x,nd,yp,irow,icol)

          do i=1,nd
          r(i)=b(i)-yp(i)-delta*delta*y(i)
          enddo

c     right side term of the system

          do i=1,np
          xp(i)=d(i)*(v(i)/x(i)-t(i))
          enddo

c--3
       call  matmul0 (0,na,A,np,xp,nd,yp,irow,icol)

          do i=1,nd
          rhs(i)=r(i)-yp(i)
          enddo

c     solve ADA'dy=rhs

       call CG_solver(A,D,na,np,nd,irow,icol,dy,rhs,delta,itercg)

c    set dx and dz 

c--4
       call  matmul0(1,na,A,np,xp,nd,dy,irow,icol)

          write(*,*) itercg

          do i=1,np
          dx(i)=d(i)*xp(i)+d(i)*(v(i)/x(i)-t(i))
          dz(i)=v(i)/x(i)-z(i)*dx(i)/x(i)
          enddo

c    primal and dual step sizes 
       
           x1=x(1)+dx(1)
           z1=z(1)+dz(1)
           kx=1
           kz=1

          do i=2,np
           if(z1.ge.(z(i)+dz(i))) then
               z1=z(i)+dz(i)
               kz=i
             endif
           if(x1.ge.(x(i)+dx(i))) then 
                x1=x(i)+dx(i)
                kx=i
             endif 
           enddo

               if(x1.ge.0.) step_p=1.
               if(x1.lt.0.) step_p=-x(kx)/dx(kx)

               if(z1.ge.0.) step_d=1.
               if(z1.lt.0.) step_d=-z(kz)/dz(kz)

             rhop=.99*step_p
             rhod=.99*step_d

             do i=1,np
              x(i)=x(i)+rhop*dx(i)
              z(i)=z(i)+rhod*dz(i)
             enddo

             do i=1,nd
              y(i)=y(i)+rhod*dy(i)
             enddo
               
              mu=(1.-amin1(rhop,rhod,0.99))*mu
              write(*,*) '*******'
              write(*,*) rhop,rhod,amin1(rhop,rhod,0.99)
c        primal-dual conditions are tested here

               PI=dot(nd,r,r)/(dot(np,x,x)+1.)
               DI=dot(np,t,t)/(dot(nd,y,y)+1.)
               DG=dot(np,z,x)/(dot(np,x,x)+1.)

               write(*,*)' # of iterations of PD_lp: ', iter
               write(*,*)' Primal Infeasibility :', PI
               write(*,*)' Dual Infeasibility   :', DI
               write(*,*)' Dualilty Gap         :', DG              

          if(PI.lt.PI_t.and.DI.lt.DI_t.and.DG.lt.DG_t) goto 11

          enddo

   
11        return
          end



c        -------------------------------------
 
         function dot( n, x, y)
         parameter (nnx=30000)
         real  x(nnx), y(nnx)
         real  val ,dot
 
         val=0.0
         do i=1,n
         val = val + x(i)*y(i)
         enddo
         dot=val
 
         return
         end
 
 
c     ----------------------------------------

          subroutine matmul0(conj,nb,bb,nx,x,ny,y,irow,icol)
 
c    Sparse matrix multiplication....
 
          parameter (nnx=30000, nny=30000, nbig=120000)
          integer conj,irow(nbig),icol(nbig)
          real bb(nbig), x(nnx), y(nnx)


          if(conj.eq.0) then

           do iy=1,ny
             y(iy)=0.0
               enddo
            do k=1,nb
             y(irow(k))=y(irow(k))+bb(k)*x(icol(k))
              enddo

          else

          do ix=1,nx
           x(ix)=0.0
            enddo
          do k=1,nb
           x(icol(k))=x(icol(k))+bb(k)*y(irow(k))
            enddo

          endif

          return
          end
 
 
c     ----------------------------------------------------------

      subroutine  CG_solver(A,D,na,np,nd,irow,icol,x,b,delta,iter)

c     Solve the system
c                        ADA'x=y whenre A is sparse and D is
c                        diagonal      

      parameter   (nnx=30000,nbig=120000,iter_max=200,tol=1.e-3)
      real        a(nbig),x(nnx),D(nnx),b(nnx)
      real        r(nnx),p(nnx),w(nnx),xp(nnx)
      integer     irow(nbig),icol(nbig)

         do i=1,nd
         x(i)=0.
         r(i)=b(i)
         enddo

         rho=dot(nd,r,r)

         do iter=1,iter_max

           if(iter.eq.1) then
             do i=1,nd
             p(i)=r(i)
           enddo

           else
             beta=rho/rho_old
             do i=1,nd
             p(i)=r(i)+beta*p(i)
             enddo
           endif
         
c        xp=A'p
            call  matmul0 (1,na,A,np,xp,nd,p,irow,icol)

c        D.xp
            do i=1,np
            xp(i)=D(i)*xp(i)
            enddo

            call  matmul0 (0,na,A,np,xp,nd,w,irow,icol)
               
c       Regularization

            do i=1,nd
            w(i)=w(i)+p(i)*delta**2
            enddo

             alpha=rho/dot(nd,p,w)
             
               do i=1,nd
               x(i)=x(i)+alpha*p(i)
               r(i)=r(i)-alpha*w(i)
               enddo
               rho_old=rho
               rho=dot(nd,r,r)
                 if(sqrt(rho/nd).le.tol) return

                 enddo  ! main loop

                   return
                   end


             


c     ----------------------------------------------------------------


      subroutine NMO(nh,dh,h0,nt,dt,t0,nv,dv,v0,A,icol,irow,ndata)

c   The NMO operator is load in the matrix A,icol,irow,with ndata elements.

      parameter (nnx=30000, nny=30000, nbig=120000, tol=1.e-5)

      real       dt,t0,v0,dv,h0,dh,t,tau,v,h
      real       a(nbig)
      integer    icol(nbig),irow(nbig)

      ny=nt*nh
      nx=nt*nv
      do itau=1,nt
      tau=t0+(itau-1)*dt
      do ih=1,nh
      h=h0+(ih-1)*dh
      do iv=1,nv
      v=v0+(iv-1)*dv
      t=(sqrt(tau*tau+h*h/(v*v))-t0)/dt
      it=1+0.5+t
      l=l+1
      j =(ih-1)*nt+it
      k =(iv-1)*nt+itau
      irow(l)=j
      icol(l)=k
      a(l)=1.
      enddo   
      enddo
      enddo   
      ndata=l
      return
      end



      function gasdev(idum)
      integer idum
      real gasdev
      integer iset
      real fac,gset,rsq,v1,v2,ran1
      save iset,gset
      data iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end
      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
      real ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     *ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      return
      end






