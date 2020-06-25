      program testlib
      implicit none
      integer i,n,m,j, ii, jj
      parameter(n=5, m=3)
      complex*16 x(n), y(n), z(n), cdot, dot, A(n,m)
      real val
      do i=1,n
         val=real(i)
         x(i)=dcmplx(val,val)
      enddo
     
      do i=1,n
         ii=i**2
         jj=i-1
         y(i)=dcmplx(ii,jj)
      enddo
      do i=1,n
        do j=1,m
          A(i,j)=dcmplx(i,j)
         enddo
       enddo

c      call cAtimesx(z,A,x,n,m)
      call cxminusy(z,x,y,n)

c      print*,dot
      do i=1,n
         print*,i,z(i)
      enddo
      stop
      end
