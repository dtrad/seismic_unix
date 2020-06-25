      subroutine Atimesx(b,A,x,nr,nc)
      implicit none
      integer i,j, nr, nc
      complex*16 A(nr,nc), x(nc), b(nr)

      do i=1,nr
         b(i)=(0.d0,0.d0)
         do j=1,nc
            b(i)=b(i)+A(i,j)*x(j)
         enddo
      enddo
      return
      end
