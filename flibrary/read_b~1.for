      subroutine read_bin_d2s(filename,nt,nx,x)

      parameter (n1x=1600,n1t=1000,nnx=n1t*n1x)

      real   x(nnx)
      real *8  d(n1x,n1t)
      character * 20 filename

      open(unit=10,file=filename,form='binary')

      do ix=1,nx
		do it=1,nt
             read(10,err=2000) d(ix,it)
          enddo
	enddo

      do ix=1,nx
          do it=1,nt
            k=(ix-1)*nt+it
            x(k)=sngl(d(ix,it))
          enddo
      enddo

2000  continue
      close(10)

      return
      end 
