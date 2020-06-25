      subroutine smooth(data,nt,nx,n1,n2)
      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)
      real       aux(nnt,nnx)


        do it=1,nt
         do ix=n1,n2
          aux(it,ix)=data(it,ix)+data(it,ix-1)+data(it,ix+1)
           enddo
           enddo

        do it=1,nt
         do ix=n1,n2
          data(it,ix)=aux(it,ix)/3.
           enddo
           enddo

               return
               end
