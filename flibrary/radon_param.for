      subroutine radon_param(fmax,xmin,xmax,dx,qmin,qmaxt,qmax,dq,nq
     #  ,rtmethod)

c Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin, 
c it computes dq and the  maximum theoretical qmaxt.
c Hence,  NMO must be adjusted such that q < qmaxt.
c qmax is the value actually used as detremined by nq and dq	
c The dx can be the average or maximum, it is under research.
c rtmethod=1 PRT
c rtmethod=2 LRT 
c Daniel Trad- UBC- 16-2-99
   
	implicit none
	real*8 fmax, xmin , xmax, dx, qmin, qmax, dq, qmaxt
	integer nq, rtmethod
	if (rtmethod.eq.1) then
	  dq= 1/(fmax*(xmax-xmin)**2)
	  dq=0.8*dq
	  qmax=qmin+dq*(nq-1)
	  qmaxt = 1/(2*fmax*(xmax-xmin)*dx)
	elseif (rtmethod.eq.2) then
	  dq= 1/(fmax*(xmax-xmin))
	  dq=0.8*dq
	  qmax=qmin+dq*(nq-1)   
	  qmaxt = 1/(fmax*dx)

	endif  
	return
	end











