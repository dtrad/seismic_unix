      subroutine write_bin(filename,nt,nx,x)
      parameter (nnt=2048,nnx=128)
      real x(nnt,nnx)
      character * 20  filename
      open(unit=11,file=filename,access='direct',recl=nt)
      do ix=1,nx
       write(11,rec=ix,err=2000)(x(it,ix),it=1,nt)
        enddo
         close(11)
	goto 2001 
2000  print*, "error writing file"
2001	continue
 
        return
         end