      subroutine read_bin(filename,nt,nx,x)
	parameter (nnt=2048,nnx=128)
      real x(nnt,nnx)
      character *20 filename
 
      open(unit=12,file=filename,access='direct',recl=nt)
 
       do ix=1,nx
        read(12,rec=ix,err=2000) (x(it,ix),it=1,nt)

          enddo
	goto 2001 
2000  print*,'error reading file'
2001	continue
      close(12)
 
      return
      end