
       real * 8 a
       character * 20 filein,fileout
	  print*,'filein'	
        read(*,*) filein
	  print*,'fileout'
        read(*,*) fileout

        open(20,file=filein,status='unknown')
        open(30,file=fileout,form='binary',status='unknown')

	  do while (eof(20).NE.(.true.))		
		  read(20,*,end=100) a
		  write(30) a
	  enddo	
100	  continue
	end