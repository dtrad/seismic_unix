      subroutine graddesc(L,LH,d,u,nx,np,eps1,iter_end,dh,dq,freq,Jtot)
      implicit none   
      integer i,j,k,kk, iter_end, iter, nnx, nnp, np, nx, freq,laststep
      parameter (nnx=228, nnp= 228) 
 
      complex*16 u(nnp),d(nnx),L(nnx,nnp),LH(nnp,nnx),g(nnp), DD(nnp)
      complex*16 gtemp(nnp),gtemp2(nnp),gtemp3(nnp), gtemp4(nnp)
      complex*16 alfanum, alfaden, alfa, dtemp(nnp), cdot, dc(nnx) 
      real*8 eps1, eps, sigman, dh(nnx), dq, value
      real*8 costfunc, Jtot(*), noise, power, resid, power2(nnp)
      real*8 pmax, pmin, eps2, pmaxold, residold, Jtotlast, Jtotprev
      integer norma
    
c     In this routine d=Lu and u=LHd to keep Thorson convention
c     Based in the Thorson'thesis  algorithm
c     Look that the rest of program is opposite
      eps=1d-5
      sigman=1d-5
      norma=1  ! Huber
      norma=10 ! Cauchy
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Unconstrained gradient step
c-----     LHL=LH*L 
c      write(0,*) eps, iter_end, nx, np, dh(2), freq

c-----Apply weights to L and LH

      do j=1,nx
         do i=1,np
           LH(i,j)=LH(i,j)*dh(j)
         enddo
      enddo 
  

c---------------------------------
c----       dtemp=LH*d
      call cAtimesx(gtemp2,LH,d,np,nx,nnp,nnx)  
      call cxequaly(u,gtemp2,np)

c----       Minimum for u=(eps,eps)
      do i=1,np
            power2(i)=dreal(dconjg(u(i))*u(i))
            if(power2(i).lt.eps) then
            u(i)=dcmplx(eps,eps)
            power2(i)=dreal(dconjg(u(i))*u(i))
            endif   
      enddo
c----       Sigma=min(u(i)**2)

      pmin=power2(1)
      do i=1,np
              if (power2(i).lt.pmin) pmin=power2(i)
      enddo
      pmin=pmin*1d-3
      if (pmin.gt.1.d-6) then 
         eps2=pmin
      else 
         eps2=1d-6
      endif
      eps2=pmin
        
      write(0,*) eps2, freq
c-------------------------------------------------
c------     gtemp=LH*L*u
      call cAtimesx(dtemp,L,u,nx,np,nnx,nnp)
      call cAtimesx(gtemp,LH,dtemp,np,nx,nnp,nnx)      

c-----      g=LH*L*u-LH*d 
      call cxminusy(g,gtemp,gtemp2,np)


c------    alfa=g'*g/(g'L'Lg)
      alfanum=cdot(g,g,np)
      call cAtimesx(dtemp,L,g,nx,np,nnx,nnp)
      call cAtimesx(gtemp3,LH,dtemp,np,nx,nnp,nnx)
      alfaden=(cdot(g,gtemp3,np))
      alfa=alfanum/(alfaden+dcmplx(eps,0.d0))
     
      do i=1,np
           u(i)=u(i)-alfa*g(i)
      enddo

c----Estimate Qp
      iter=1
      Jtotlast=0
      Jtotprev=10
      do while ((iter.le.iter_end).and.(Jtotlast.le.Jtotprev))
         power=0.d0
         do i=1,np
              power2(i)=dreal(dconjg(u(i))*u(i))
              power=power+power2(i)
         enddo
         pmax=power2(1)
         do i=1,np
              if (power2(i).gt.pmax) pmax=power2(i)
         enddo
         pmin=power2(1)
         do i=1,np
              if (power2(i).lt.pmin) pmin=power2(i)
         enddo
         if (pmin.gt.1.d-6) eps2=pmin
         if (pmin.lt.1.d-6) eps2=1d-6
c        write(0,*) 'eps2',eps2
c        power=power/np
c        eps2=power
         do i=1,np
             if (norma.eq.10) then
                DD(i)=dcmplx(1.d0/(eps2+power2(i)),0.d0)+eps
             elseif (norma.eq.1) then   
                DD(i)=dcmplx(1.d0/eps2/dsqrt(power2(i)),0.d0)
             endif   
c            write(0,*) DD(i)
c            DD(i)=(1.d0,0.d0) ! Test for Gauss
         enddo
         

c---- Gradient steps


c---- NP loop

         resid=10.d0
         residold=100.d0
         k=np+1
         power=1.d0

         laststep=np-3
         if (iter.gt.4) laststep=np-30
         if (iter.gt.8) laststep=np-50
         if (iter.eq.iter_end) laststep=1

         do while ((k.ne.1).and.(resid.gt.1e-10).and.(k.gt.laststep))
             k=k-1
            call cAtimesx(dtemp,L,u,nx,np,nnx,nnp)
 
            call cAtimesx(gtemp,LH,dtemp,np,nx,nnp,nnx)

            call cxtimesy(gtemp4,DD,u,np)

c           do kk=1,np
c            write(0,*) gtemp(kk),gtemp4(kk)
c           enddo

            call cxplusy(gtemp,gtemp,gtemp4,np)

c-------(L'L+D)u-L'd

            call cxminusy(g,gtemp,gtemp2,np) !new gradient

c-----   alfa=g'*g/(g'(L'L+D)g)

            alfanum=cdot(g,g,np)
            call cAtimesx(dtemp,L,g,nx,np,nnx,nnp)
            call cAtimesx(gtemp3,LH,dtemp,np,nx,nnp,nnx)
            call cxtimesy(gtemp4,DD,g,np)
            call cxplusy(gtemp3,gtemp3,gtemp4,np)
            alfaden=(cdot(g,gtemp3,np))
            alfa=alfanum/(alfaden+dcmplx(eps,0.d0))

            residold=resid
            resid=0.d0
            do i=1,np
              u(i)=u(i)-alfa*g(i)
              resid=resid+dconjg(alfa*g(i))*(alfa*g(i))
            enddo
c           resid=dreal(alfa)*dreal(cdot(g,g,np))
            
            power=dreal(cdot(u,u,np))
            
c           write(0,*) resid,power
         enddo                  !loop for k

         call cAtimesx(dc,L,u,nx,np,nnx,nnp)

         Jtot(iter)=costfunc(d,dc,u,nx,np,eps2,norma)
         write(0,*) 'Jiter', Jtot(iter), iter,k, freq
         Jtotlast=Jtot(iter)
         if (iter.eq.1) then
                 Jtotprev=Jtot(iter)+10
         else
                 Jtotprev=Jtot(iter-1)
         endif 
 
         iter=iter+1
      enddo  !loop for niter 

      return
      end
