      subroutine gauss_gauss_t0(l,lh,d,model,nd,np,eps1,dh)
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047
c	  Modified Daniel Trad, Department of earth and Ocean Sciences, UBC.	
c       January, 30. 1999.

c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Cauchy-Gauss regularization with Levinson recursion.

c     Input Parameters:
c
c         L(nd,np) - Linear operator (complex)
c         d(nd)    - data (complex)
c         eps1     - pre-whitening in % (1% is ok).
c         iter_end - maximum number of iteration (5 is ok)

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
c     This version uses Levinson recursion method iteratively
c     retrieve the solution. 
C     Because the offset interval is not constant the problem to solve is 
c        m=inv(L.LH+Qp).L d
c	where L and LH are defined as 
c	   m=Ld and d=LH m	
c      This is done in two steps
c     1- Solve (L.LH+Qp)x=Ld
c     Then m=x
c	Remember that Qp is diagonal and Qp=sigma^2 1/(1+(|v|/sigma)^2

      parameter    (nnx=228) 

      complex *16  a(nnx), model(nnx),d(nnx),L(nnx,nnx),lh(nnx,nnx)
      complex *16  aux(nnx),lambda(nnx),sum, ll(nnx,nnx)  
      real    * 8  power(nnx),p(nnx),pm, a0
      real    * 8  eps1,eps2,s,snew,sold,th, dh(nnx)

c---- Right hand side is the adjoint
	
      do j=1,nd
		do i=1,np
			ll(i,j)=dh(j)*l(i,j)
		enddo
	enddo
	 
      do i=1,np
		model(i)=dcmplx(0.d0,0.d0)
		do j=1,nd
			model(i)=model(i)+ll(i,j)*d(j)
		enddo
		aux(i)=model(i)
	enddo



c-----compute the matrix to invert
c-----eps1 is the damping 1% is ok
      
       do 200 j=1,np 
         sum=dcmplx(0.d0,0.d0)
	   do 300 k=1,nd		!wrong 300  --> sum=sum+l(1,k)*lh(k,j) 
300			sum=sum+ll(j,k)*lh(k,1)        !This is (LH.inv(Qp).L) 
	   if(j.eq.1) a0=dreal(sum*(1.d0+eps1/100.d0))  ! Add Cn term
	   if(j.ne.1) a(j-1)=sum
200	 continue

       ntot=np-1

       call HERM (ntot,a0,a,aux,model,ISTAT)
       if (ISTAT.NE.0) print*,"Problems in herm. ISTAT.NE.0" 


        return
        end

c-----------------------------------------------------------------
