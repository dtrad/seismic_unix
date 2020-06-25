      subroutine Atimesx(A,x,b,nr,nc)
      implicit none
      complex*16 A(nr,nc), x(nc), b(nr)
      integer i,j
      do i=1,nr
         b(i)=(0.d0,0.d0)
         do j=1,nc
            b(i)=b(i)+A(i,j)*x(j)
         enddo
      enddo
      return
      end
