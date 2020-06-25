      subroutine LMSC(d,nd,mu,w,nw)

      parameter  (nnt=2048,nnx=128)

      complex   w(nnx),dh(nnt),d(nnt)
      real      mu

      do i=2,nw                     !Initail filter
       w(i)=0.
        enddo
         w(1)=1.


      do k=nw+1,nd

      dh(k)=0.
       do m=1,nw                    !forward prediction stage
         dh(k)=dh(k)+w(m)*d(k-m)
           enddo

           do m=1,nw
            w(m)=w(m)+2.*mu*(d(k)-dh(k))*conjg(d(k-m))
             enddo

             enddo

 
           do i=1,nw               ! filter in PEO form
            w(i)=-w(i)
             enddo 

           return
          end
       

      
