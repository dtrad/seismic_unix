c     -------------------------------------------------------
c
c     This is the good one. 
c     Works ok for real and complex data. 
c     We do not consider any type of symetry,
c     the conitinuity of the gap is required.

      subroutine gapf(x,n,n1,n2,ip,a,perc)   

c   1D GAP filling.

      parameter  (nnt=2048,nnx=228)
      complex    x(nnt),a(nnx),r(nnx),rr
      complex    g(nnx),xx(nnt),b(nnt),c(nnx),z(nnx),R0C
      real       r0,perc
      
c   Form the peo from the one-step-ahead Pred. filter
      
      g(1)=cmplx(1.,0.0)
      do i=1,ip
		g(i+1)=a(i)
      enddo

c   Autocorrelation of the peo      
      call correlation(ip+1,ip+1,1,g,g,R0C,R)
      R0 = real(R0C)

c   Compute righ-side term 
              write(1,*)' perc=',perc
      R0=R0*(1+perc/100.)
      do j=n1,n2
	      b(j+1-n1)=cmplx(0.,0.)
	      do m=-ip+n1,n1-1
			l=m-j
			if(l.eq.0) rr=cmplx(r0,0.)
			if(l.gt.0) rr=r(m-j)
			if(l.lt.0) rr=conjg(r(j-m))
			b(j+1-n1)=b(j+1-n1)+x(m)*conjg(rr)
		  enddo
      enddo

      do j=n1,n2
		c(j+1-n1)=cmplx(0.,0.)
		do m=n2+1,ip+n2
			l=m-j
			if(l.eq.0) rr=cmplx(r0,0.)
			if(l.gt.0) rr=r(m-j)
			if(l.lt.0) rr=conjg(r(j-m))
			c(j+1-n1)=c(j+1-n1)+x(m)*conjg(rr)
		enddo
      enddo

      do i=1,n2-n1+1
		z(i)=-b(i)-c(i)
      enddo

c   Solve for the gap

      call toeplitz(n2-n1,r0,r,Z,xx,ISTAT)
      if (ISTAT.eq.1) write(1,*) 'ISTAT in Toeplitz', ISTAT
c   Fill the gap

      do i=n1,n2
      x(i)=xx(i-n1+1)
      enddo
      return
      end 
