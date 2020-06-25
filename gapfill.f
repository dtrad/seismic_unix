	subroutine gapfill(datac,nt,nx,dt,dx,n1,n2,ip,perc)

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

	  

      parameter  (nnt=2048,nnx=228)
      real       data(nnt,nnx),datac(*),perc
	open(1,file='temp')

      write(1,*) nt, nx, dt, dx	, ip, n1, n2, perc

      do ix=1,nx
	 do it=1,nt    
            data(it,ix)=datac((it-1)+(ix-1)*nt)
	enddo
      enddo	

      call gapf_fx(data,nt,nx,dt,dx,n1,n2,ip,perc)
      
      do ix=1,nx
	 do it=1,nt
c	    write(1,*) data(it,ix)    
            datac((it-1)+(ix-1)*nt)=data(it,ix)
	enddo
      enddo		
      close(1)
      
       return
       end


