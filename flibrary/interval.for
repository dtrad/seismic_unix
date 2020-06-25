	subroutine interval(x,lx,mx,ax)
	implicit none
	
	integer i, lx
	real*8 x(lx), mx, ax, dx		
      mx=0
	ax=0	
	do i=2,lx
		dx=x(i)-x(i-1)
		if (dx.gt.mx) mx=dx
		ax=ax+dx
	enddo
	ax=ax/(lx-1)
	return
	end	
