      subroutine conjgrad(L,LH,d,u,nx,np,eps1,iter_end)
      implicit none   
      integer i,j,k, iter_end, iter, nnx, nnp, np, nx
      parameter (nnx=228, nnp= 228) 
 
      complex*16 u(nnp),d(nnx),L(nnx,nnp),LH(nnp,nnx),g(nnp)
      complex*16 LHL(nnp,nnp),gtemp(nnp),gtemp2(nnp)
      complex*16 alfanum, alfaden, alfa 
      real*8 eps1, sigman, DD(nnp)

c     In this routine d=Lu and u=LHd to keep Thorson convention
c     Based in the Thorson'thesis  algorithm
c     Look that the rest of program is opposite

      sigman=eps1

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Unconstrained gradient step
c-----     LHL=LH*L 
 
      do i=1,np
         do j=1,np
            LHL(i,j)=(0.d0,0.d0)
            do k=1,nx
               LHL(i,j)=LH(i,k)*L(k,j)
            enddo
         enddo
      enddo
c-----      g=LH*L*u-LH*d
      do i=1,np
         do j=1,np
           gtemp(i)=LHL(i,j)*u(j)
         enddo
      enddo
      

      do i=1,np
         do j=1,nx
           gtemp2(i)=LH(i,j)*d(j)
         enddo
      enddo
  
      do i=1,np
           g(i)=gtemp(i)-gtemp2(i)
      enddo

c------    alfa=g'*g/(g'L'Lg)
           
      do i=1,np
           alfanum=dconjg(g(i))*g(i)
      enddo


      do i=1,np
           alfaden=dconjg(g(i))*gtemp(i)
      enddo

      alfa=alfanum/alfaden
      
      do i=1,np
           u(i)=u(i)-alfa*g(i)
      enddo

c----Estimate Qp

      do i=1,np
           DD(i)=sigman/zabs(u(i))
      enddo

c---- Gradient steps


      do i=1,np
         do j=1,np
            LHL(i,j)=LHL(i,j)+DD(i)
         enddo
      enddo   

c---- NP loop

      do iter=1,iter_end
      do k=np,1,-1

c-------(L'L+D)u-L'd

        do i=1,np
           do j=1,np
             gtemp(i)=LHL(i,j)*u(j)  ! (L'L+D)u (it changes for every k)
           enddo
        enddo
      

        do i=1,np
             g(i)=gtemp(i)-gtemp2(i) !g=(L'L+D)u-L'd
        enddo       


c------    alfa=g'*g/(g'L'Lg)
           
        do i=1,np
           alfanum=dconjg(g(i))*g(i)
        enddo


        do i=1,np
           alfaden=dconjg(g(i))*gtemp(i)
        enddo

        alfa=alfanum/alfaden
      
        do i=1,np
           u(i)=u(i)-alfa*g(i)
        enddo

       enddo !loop for k
      enddo  !loop for niter 
      return
      end
