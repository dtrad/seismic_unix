c  Gap filling program for Shots and CMP Gathers 

c  This is the gap filling program that uses AR
c  spectral modeling to retrieve a continuous gap.

c  This program read data and parameters from
c  a file called ``gap.in''. 

c  parameters read from gap.in

c  filein  : a binary file containing the cdp or shot
c  fileout : a binary file containing the cdp or shot
c  nt,nx   : number of samples and number of channels of the cmp
c  n1,n2   : from where to where the gap extends
c  ip      : lenght of the peo filter 

c the ouput goes into another binary called 
c fileout. If the output file is fileone='myfile'
c you can use SU to see it
c
c xwigb n1=800<myfile perc=99&
c where n1=nt=800 is the number of samples 


c TO PROCESS REAL DATA: 
c     I think that the best is to read segy files
c     with gaps fill in with zeros traces and
c     know from where to where the gap goes (n1,n2).
c     In this way your header is fix and you don't
c     have to create it for the new
c     traces in the gap.



 
c (c) M.D.SACCHI, Bayesian Solutions Ltd.
c                 Edmonton, AB, Canada
c                 (403)-438-4336
c                 sacchi@phys.ualberta.ca

	  
	  include "c:\mauricio\flibrary\sfx_go.for"
	  include "c:\mauricio\flibrary\sfx_symmetry.for"
	  include "c:\mauricio\flibrary\fork.for"
	  include "c:\mauricio\flibrary\toeplitz.for" 
	  include "c:\mauricio\flibrary\burg.for"
	  include "c:\mauricio\flibrary\correlation.for"
	  include "c:\mauricio\flibrary\read_bin.for"
	  include "c:\mauricio\flibrary\write_bin.for"

      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx)
      character  * 20 filein,filein2,fileout 

      open(10,file='gap.in')
      read(10,'(a20)') filein
	read(10,'(a20)') filein2
      read(10,'(a20)') fileout
      read(10,*) nt,nx
      read(10,*) dt,dx
      read(10,*) n1,n2
      read(10,*) ip
	


	call read_bin(filein,nt,nx,data) 

c add a gap here

      do ix=n1,n2
		do it=1,nt
			data(it,ix)=0
		enddo
      enddo
	call write_bin(filein2,nt,nx,data) 
c      Write the data with gap
c	 subroutine gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)
       call gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)
	
      call write_bin(fileout,nt,nx,data) 
      
       stop
       end

c     -------------------------------------------

      subroutine gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)
c     call to
c			SUBROUTINE BURG (N,IP,X,P,A,ISTAT) 


      parameter  (nnt=2048,nnx=128)
      real       data(nnt,nnx),mu
      complex    data_fx(nnt,nnx)
      complex    x(nnx),a(nnx)
      complex    x2(nnx),a2(nnx)


      nf=2048
      if(nf.lt.nt) pause 'Check NF'

 
      call sFX_go (nx,nt,nf,data,data_fx,-1)

c Do it for all freqs in the band

      iflow = 5
      ifhigh = nf/2+1-5

      do if=iflow,ifhigh
        do ix=1,nx
           x(ix)=data_fx(if,ix)
        enddo

        call BURG (N1-1,IP,X,P,A,ISTAT)
        do ix=1,nx-n2
           x2(ix)=data_fx(if,ix+n2)
        enddo
        call BURG (Nx-n2,IP,X2,P,A2,ISTAT)
        do i=1,ip
          a(i)=0.5*(a(i)+a2(i))
        enddo
        call gapf(x,nx,n1,n2,ip,a)   
        do ix=1,nx
          data_fx(if,ix)=x(ix)
        enddo
      enddo

      call  sFX_symmetry (data_fx,nx,nt,nf,iflow,ifhigh)
      call  sFX_go (nx,nt,nf,data,data_fx,1)

      return
      end

c     -------------------------------------------------------
c
c     This is the good one. 
c     Works ok for real and complex data. 
c     We do not consider any type of symetry,
c     the conitinuity of the gap is required.

      subroutine gapf(x,n,n1,n2,ip,a)   

c   1D GAP filling.

      parameter  (nnt=2048,nnx=128)
      complex    x(nnt),a(nnx),r(nnx),rr
      complex    g(nnx),xx(nnt),b(nnt),c(nnx),z(nnx),R0C
      real       r0
      
c   Form the peo from the one-step-ahead Pred. filter
      
      g(1)=cmplx(1.,0.0)
      do i=1,ip
		g(i+1)=a(i)
      enddo

c   Autocorrelation of the peo      
      call correlation(ip+1,ip+1,1,g,g,R0C,R)
      R0 = real(R0C)

c   Compute righ-side term 
      perc=0.
      R0=R0*(1.+perc/100.)
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

c   Fill the gap

      do i=n1,n2
      x(i)=xx(i-n1+1)
      enddo
      return
      end 
