      subroutine cAtimesx(b,A,x,nr,nc,nnr,nnc)
      implicit none
      integer i,j, nr, nc, nnr, nnc, nnx
      parameter(nnx=228)
      complex*16 A(nnx,nnx), x(nnx), b(nnx)

      do i=1,nr
         b(i)=dcmplx(0.d0,0.d0)
         do j=1,nc
           b(i)=b(i)+A(i,j)*x(j)
         enddo
      enddo
      return
      end
