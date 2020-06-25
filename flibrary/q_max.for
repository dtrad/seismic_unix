	FUNCTION q_max (fmax,xmin,xmax,dx)
	implicit none	
      real*8 fmax, xmax, xmin, dx, q_max
	q_max = 1/(fmax*(xmax-xmin)*dx)
	return
	end 