
      subroutine forward(l,nd,np,d,m)

c     Matrix multiplication L.m=d
c     L,m and d are complex
c
c     Input parameters:
c
c       L(nd,np)  - linear operator
c       m(np)     -  complex vector of np elements
c
c     Output  parameters:
c
c       d(nd)     - complex vector of nd elements


      parameter    (nnx=128)
      complex * 16 m(nnx),d(nnx),L(nnx,nnx)

      do j=1,nd
      d(j)=dcmplx(0.d0,0.d0)
      do i=1,np
      d(j)=d(j)+l(j,i)*m(i)
      enddo
      enddo

      return
      end


