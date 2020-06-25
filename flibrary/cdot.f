      function cdot(x,y,n)

c     Compute the inner product
c     dot=(x,y) for complex x,y
         implicit none
         integer n,i
         complex*16  x(n), y(n), cdot, val
        
         val=(0.0d0,0.d0)
         do i=1,n 
            val = val +dconjg( x(i))*y(i)
         enddo

         cdot=val

         return
         end


