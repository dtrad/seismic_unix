      
C     driver for routine lubksb
	include "c:\daniel\radon\flibrary\ludcmp_c.for"	
      include "c:\daniel\radon\flibrary\lubksb_c.for"

	PROGRAM xlubksb
	INTEGER NP
      PARAMETER (NP=20)
      complex*16 p,a(NP,NP),b(NP,NP),c(NP,NP),x(NP)
      INTEGER j,k,l,m,n,indx(NP)
      CHARACTER txt*3
      open(7,file='c:\daniel\press\other\MATRX1.DAT',status='old')
      read(7,*)
10    read(7,*)
      read(7,*) n,m
      read(7,*)
      read(7,*) ((a(k,l), l=1,n), k=1,n)
      read(7,*)
      read(7,*) ((b(k,l), k=1,n), l=1,m)
	a(1,1)=a(1,1)+cmplx(0,1)
	a(2,1)=a(1,1)+cmplx(0,4)
	b(1,1)=b(1,1)+cmplx(0,4)
C     save matrix a for later testing
      do 12 l=1,n
        do 11 k=1,n
          c(k,l)=a(k,l)
11      continue
12    continue
C     do LU decomposition
      call ludcmp_c(c,n,NP,indx,p)
C     solve equations for each right-hand vector
      do 16 k=1,m
        do 13 l=1,n
          x(l)=b(l,k)
13      continue
        call lubksb_c(c,n,NP,indx,x)
C     test results with original matrix
        write(*,*) 'Right-hand side vector:'
        write(*,'(1x,6f12.6)') (b(l,k), l=1,n)
        write(*,*) 'Result of matrix applied to sol''n vector'
        do 15 l=1,n
          b(l,k)=0.0
          do 14 j=1,n
            b(l,k)=b(l,k)+a(l,j)*x(j)
14        continue
15      continue
        write(*,'(1x,6f12.6)') (b(l,k), l=1,n)
        write(*,*) '***********************************'
16    continue
      write(*,*) 'Press RETURN for next problem:'
      read(*,*)
      read(7,'(a3)') txt
      if (txt.ne.'END') goto 10
      close(7)
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .
