      subroutine CAUCHY_GAUSS_0(l,lh,d,model,nd,np,eps1,iter_end)
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047


c     Solve an underdetermined  linear inverse problem, L.m=d,
c     using the Cauchy-Gauss regularization.

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
c     This version uses CHOLESKY method iteratively
c     retrieve teh solution. 
C     Because np> nd the problem to solve is 
c        m=inv(Qp)L(Cn+LH.inv(Qp).L).d
c      This is done in two steps
c     1- Solve x=(Cn+LH.inv(Qp).L).d
c     Then m=inv(Qp).Lx
c	Remember that Qp is diagonal and inv(Qp)=(1/sigma^2)(1+(|v|/sigma)^2

      parameter    (nnx=228) 

      complex*16 a(nnx*(nnx+1)),model(nnx),d(nnx),L(nnx,nnx),lh(nnx,nnx) 
      complex*16 aux(nnx),lambda(nnx),sum
      real*8 power(nnx),p(nnx),pm
      real*8 eps1,eps2,s,snew,sold,th
	integer nnxx
	nnxx=nnx

c---- Initial model is the adjoint

      th=0.00001d0

c      do 10 i=1,np
c		model(i)=dcmplx(0.d0,0.d0)
c		do 10 j=1,nd
c			model(i)=model(i)+l(i,j)*d(j)
c10    continue

      s=0.0
      pm=0.0
      do i=1,np
		power(i)=dreal(model(i)*dconjg(model(i)))
		pm=pm+power(i)
      enddo
      
	s=power(1)
      pm=pm/np
      eps2=pm
      
	do i=1,np
         if(power(i).gt.s) s=power(i)
      enddo
      
	smin=power(1)
      
	do i=1,np
         if(power(i).le.smin) smin=power(i)
      enddo
         
      if(smin.gt.1.d-6) eps2=smin
      if(smin.le.1.d-6) eps2=1.d-6
c      if(s.gt.1.d-6) eps2=s
c      if(s.le.1.d-6) eps2=1.d-6
c     write(*,*)' max: ',s, ' min: ',smin, ' eps2 : ', eps2
       	
      do i=1,np
         p(i)=1.d0+power(i)/eps2
         s=s+dlog(1.d0+power(i)/eps2)
      enddo
      sold=s
     
      do 1000 iter=1,iter_end

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
300   sum=sum+p(k)*lh(i,k)*l(k,j)        !This is (LH.inv(Qp).L)
      n=n+1
      if(i.eq.j) a(n)=sum*(1.d0+eps1/100.d0)  ! Add Cn term
      if(i.ne.j) a(n)=sum
200   continue

      call  CHOLESKY2(nd,1.d-15,A,aux,ISTAT,nnxx) !Solution to (LH.inv(Qp).L).x=d



      do  i=1,nd
      lambda(i)=aux(i)
      enddo

      do i=1,np
      sum=dcmplx(0.d0,0.d0) 
       do  k=1,nd
        sum=sum+l(i,k)*lambda(k)            !m=inv(Qp).L.x
         enddo 
        model(i)=p(i)*sum
       enddo

      pm=0.d0
      do i=1,np
       power(i)=dreal(model(i)*dconjg(model(i)))
       pm=pm+power(i)
      enddo 
      pm=pm/np
      eps2=pm

c----- If the relative variation of the cost function is le<th
c----- returns 

      s=0.d0 
      do i=1,np
       p(i)=1.d0+power(i)/eps2
        s=s+dlog(1.d0+power(i)/eps2)
         enddo 
          snew=s
           if(dabs(snew-sold)/dabs(snew).lt.th)  return 
            sold=s

1000    continue

        return
        end

c-----------------------------------------------------------------
