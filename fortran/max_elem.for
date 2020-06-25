	subroutine interval(x,lx,mx,ax)
	implicit none
	real*8 x(lx), mx, ax, dx
	integer i
		
      mx=0
	ax=0	
	do i=2,lx
		dx=x(i)-x(i-1)
		if (dx.gt.mx) mx=dx
		ax=ax+dx
	enddo
	ax=ax/lx
	return
	end	
