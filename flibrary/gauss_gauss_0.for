
      subroutine GAUSS_GAUSS_0(l,lh,d,model,nd,np,eps1)
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047
c	  Modified to include not evenly spaced offset
c       Daniel Trad-- 29-12-98	
c       E-mail: dtrad@geop.ubc.ca

c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Gauss-Gauss regularization.

c     Input Parameters:
c
c         L(np,nd) - Linear operator (complex)
c         LH(nd,np) - Linear operator (complex)

c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).


c     Output Parameters:
c
c         m(np) - vector of unkonw parameters
c
c     Notes:
c
c     The damping eps1 represents a percentage of the diagonal
c     of the matrix L.L^H.
c
c     remenber that since we are solving an underdetermine problem
c     np>nd.
c
c     This version uses CHOLESKY method to retrieve the solution. 

      parameter    (nnx=228) 

      complex *16  a(nnx*(nnx+1)),model(nnx),L(nnx,nnx),LH(nnx,nnx) 
      complex *16  d(nnx),aux(nnx),lambda(nnx),sum
      real    * 8  power(nnx),eps1
	nnxx=nnx
c---- Initial model is the adjoint

  
      do  i=1,nd
      aux(i)=d(i)
      enddo

c-----compute the matrix to invert
c-----eps1 is the damping 1% is ok

      n=0
      do 200 j=1,nd 
      do 200 i=1,j  
      sum=dcmplx(0.d0,0.d0)     
      do 300 k=1,np
300   sum=sum+lh(i,k)*l(k,j)
      n=n+1
      if(i.eq.j) A(n)=sum*(1.d0+eps1/100.d0)
      if(i.ne.j) A(n)=sum
200   continue

      call  CHOLESKY2(nd,1.d-15,A,aux,ISTAT,nnxx)

      do  i=1,nd
      lambda(i)=aux(i)
      enddo

      do i=1,np
		sum=dcmplx(0.d0,0.d0) 
		do  k=1,nd
			sum=sum+l(i,k)*lambda(k)
		enddo 
		model(i)=sum
      enddo
     
      return
      end
