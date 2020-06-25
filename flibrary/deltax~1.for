	FUNCTION deltax_max (fmax,pmin,pmax,xmin,xmax)
	implicit none	
      real*8 fmax, xmax, xmin, deltax_max, pmin, pmax
	deltax_max = 1/(fmax*(xmax-xmin)*(pmax-pmin))
	return
	end 