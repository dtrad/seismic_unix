      subroutine cxequaly(x,y,n)
      implicit none
      integer i, n
      complex*16  x(n), y(n)
     
      do i=1,n
            x(i)=y(i)
      enddo

      return
      end
