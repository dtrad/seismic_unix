      subroutine cxminusy(z,x,y,n)
      implicit none
      integer i, n
      complex*16 z(n), x(n), y(n)
     
      do i=1,n
            z(i)=x(i)-y(i)
      enddo

      return
      end
