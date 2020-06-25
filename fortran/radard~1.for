	real*8 yplt(1024, 202)
	open(2,file='\daniel\radon\Dan_t.dat') 
	open(3,file='\daniel\radon\radar2zp.bin',form='binary') 
	
	read(2,*)ns,nt
709	format(2i5)

	write(*,*)'# of traces       =',ns
	write(*,*)'# of points/trace =',nt

	do 1111 kk=1,ns
	   read(2,*)(yplt(i,kk),i=1,nt)
1111	continue
	do 1112 kk=1,ns
	do 1113 j=nt+1,1024
	   yplt(j,kk)=0
1113	continue
1112	continue
	do k=1,100	
	do j=1,1024
	      write(3) yplt(j,k)
	enddo
	enddo	
	close(2)
      close(3)	
	end
