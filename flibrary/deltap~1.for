	FUNCTION deltap_max (fmax,xmin,xmax)
	implicit none	
      real*8 fmax, xmax, xmin, deltap_max
	deltap_max = 1/(fmax*(xmax-xmin)**2)
	return
	end 