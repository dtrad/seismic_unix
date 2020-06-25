      function power(data,nt,nx,n1,n2)
      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)


           p=0.
               do it=1,nt
                 do ix=n1,n2
                  p=p+data(it,ix)**2
                   enddo
                    enddo
               power=p/(nt*(n2-n1+1))

               return
               end