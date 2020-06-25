c-------------------------------------------------------------
       subroutine norma(x,nt,nh)

c------Normalization of y
       
       parameter (nnt=2048,nnx=128)

       real * 8 x(nnt,nnx)
       real * 8 xmax

       xmax=x(1,1)
       do it=1,nt
       do ih=1,nh
       if(dabs(x(it,ih)).gt.xmax) xmax=dabs(x(it,ih))
       enddo
       enddo
       
       do it=1,nt
       do ih=1,nh
       x(it,ih)=x(it,ih)/xmax
       enddo
       enddo
       return
       end
