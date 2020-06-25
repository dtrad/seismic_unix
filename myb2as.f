	program mya2b
	real a
	integer nt
	character*20 filein, fileout
	print*, 'filein'
	read(*,*) filein
	print*, 'fileout'
	read(*,*) fileout

	print*, 'nt'
	read(*,*) nt
	open(2,file=filein)
	open(1,file=fileout,form='unformatted')
	n=0
	do j=1,2000
	do i=1,nt
	read(1,end=10) a
	write(2,*) a
	enddo
	n=n+1
	print*,n
	enddo
10	close(1)
	close(2)
	end
