      subroutine CAUCHY_GAUSS_T0(l,lh,d,model,nd,np,eps1,iter_end)
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

      parameter    (nnx=128) 

      complex *16  a(nnx), model(nnx),d(nnx),L(nnx,nnx),lh(nnx,nnx) 
      complex *16  aux(nnx),aux2(nnx),lambda(nnx),sum
      real    * 8  power(nnx),p(nnx),pm, a0
      real    * 8  eps1,eps2,s,snew,sold,th

c---- Initial model is the adjoint

      th=0.000001d0
	
      do i=1,np
		model(i)=dcmplx(0.d0,0.d0)
		do j=1,nd
			model(i)=model(i)+l(i,j)*d(j)
		enddo
		aux(i)=model(i)
	enddo

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

c     write(*,*)' max: ',s, ' min: ',smin, ' eps2 : ', eps2
      	
      do i=1,np
         p(i)=(1.d0+power(i)/eps2)
         s=s+dlog(1.d0+power(i)/eps2)
      enddo
      sold=s
     
      do 1000 iter=1,iter_end

c-----compute the matrix to invert
c-----eps1 is the damping 1% is ok
	 if (iter.gt.1) print*,iter
      
      do 200 j=1,np 
         sum=dcmplx(0.d0,0.d0)
	   if (iter.lt.4) then
			aux2(j)=aux(j)
	   else
			aux2(j)=aux(j) * 1/p(j)     
	   endif
    		!model(j)=dcmplx(0.d0,0.d0)
	   do 300 k=1,nd
c wrong 300  --> sum=sum+l(1,k)*lh(k,j) 
300			sum=sum+l(j,k)*lh(k,1)        !This is (LH.inv(Qp).L)
	   if(j.eq.1) a0=dreal(sum*(1.d0+eps1/100.d0))  ! Add Cn term
	   if(j.ne.1) a(j-1)=sum
200   continue

      ntot=np-1

      call HERM (ntot,a0,a,aux2,model,ISTAT)
      if (ISTAT.NE.0) print*,"Problems in herm. ISTAT.NE.0" 
c     call  CHOLESKY (nd,1.d-15,A,aux,ISTAT)  ! Solution to (LH.inv(Qp).L).x=d

c      do  i=1,nd
c      lambda(i)=aux(i)
c      enddo

c      do i=1,np
c      sum=dcmplx(0.d0,0.d0) 
c       do  k=1,nd
c        sum=sum+l(i,k)*lambda(k)            !m=inv(Qp).L.x
c         enddo 
c        model(i)=p(i)*sum
c       enddo


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
c           if(dabs(snew-sold)/dabs(snew).lt.th)  return 
            sold=s

1000    continue

        return
        end

c-----------------------------------------------------------------
