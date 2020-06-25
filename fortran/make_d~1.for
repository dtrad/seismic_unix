c This is a program to generate hyperbolas
c to test the gap filling technique.
c
c f77 -O2 xx.f -o xx.exe and then run xx.exe 
c this program generates a binary file that can 
c be read with SU

c to display it do :   xwigb n1=800<data&
c then run the gap filling algorithm.


c (c) M.D.SACCHI, Bayesian Solutions Ltd.
c                 Edmonton, AB, Canada
c                 (403)-438-4336
c                 sacchi@phys.ualberta.ca


      real       d(200,800)
      real       w(800)
      integer    ix(10),iz(10)



      dt=0.002
      dx=10.
      dz=10.

      nt=800
      nr=100

      v = 2800.
      f = 25.

      call ricker(f,dt,w,nw)
        
      is = nr/2 
      nscat = 8

      ix(1) = 40
      iz(1) = 20

      ix(2) = 50
      iz(2) = 30

      ix(3) = 70
      iz(3) = 60

      ix(4) = 40
      iz(4) = 43

      ix(5) = 46
      iz(5) = 40

      ix(6) = 41
      iz(6) = 48

      ix(7) = 41
      iz(7) = 78

      ix(8) = 51
      iz(8) = 88


      do i=1,nt
      do ir = 1,nr
      d(ir,it) = 0.
      enddo
      enddo

      do ir = 1,nr

      do k = 1,nscat
      T = sqrt(( (ir-ix(k))*dx )**2 + (iz(k)*dz)**2 )/v+
     :    sqrt(( (is-ix(k))*dx )**2 + (iz(k)*dz)**2 )/v
      it = ifix(T/dt)+1

      do iw = 1,nw
      if(it+iw.le.nt) d(ir,it+iw) = d(ir,it+iw) + w(iw)
      enddo

      enddo
      enddo

      call write_bin('data',nt,nr,d)
     
      stop
      end


      subroutine ricker(f,dt,w,nw)

      parameter (nnw=800)
      real      w(nnw)

c   generate a ricker source wavelet of
c   central freq. f
c
c   input:  f: central f in hz
c          dt: sampling in sec
c
c   output w(nw): the wavelet
c             nw: is the length which is
c                 computed in the program


      pi=4.*atan(1.)
      nw=4./(f*dt)+1
      do i=1,nw
      alpha=(i-nw/2)*pi*f*dt
      beta=alpha**2
      w(i)=(1.-2.*beta)*exp(-beta)
      enddo
      return
      end





      subroutine write_bin(filename,nt,nx,x)

      real      x(200,800)
      character * 20 filename

      open(unit=11,file=filename,access='direct',recl=4*nt)

           do ix=1,nx
            write(11,rec=ix,err=2000) (x(ix,it),it=1,nt)
             enddo

2000  continue
      close(11)
      return
      end



