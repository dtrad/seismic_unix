      function costfunc(d,dc,u,nx,np,sigma,norma)

      implicit none
      integer  nx, np,i, norma
      complex*16 u(*), d(*), dc(*)
      real*8 sigma,r2, costfunc, power

      r2=0
      do i=1,nx
         r2=r2+zabs((d(i)-dc(i)))**2
      enddo
      power=0.d0
      if (norma.eq.10) then
          do i=1,np
              power=power+dlog(1+dreal(dconjg(u(i))*u(i))/sigma)
          enddo
            
      elseif(norma.eq.1) then
          do i=1,np
              power=power+dcmplx(zabs(u(i)),0.d0)
          enddo
      endif
      costfunc=r2+power

c     write(0,*) 'Jtot,Jdata,Jmodel',costfunc,r2,power
      
      
      return
      end

      

      
