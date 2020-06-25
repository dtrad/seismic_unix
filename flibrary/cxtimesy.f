      subroutine cxtimesy(z,x,y,n)

c     Compute the inner product
c     dot=(x,y) for complex x,y
         implicit none
         integer n,i
         complex*16  x(n), y(n), z(n)
                 
         do i=1,n 
            z(i) =  x(i)*y(i)
         enddo

         return
         end


