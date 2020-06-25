c===============================================================
        subroutine mute (x,Lm,Ln,m1,n1,m2,n2,m3,n3,m4,n4,y,IMUTE)

c  Zeroes a rectangular region of the matrix x(Lm,Ln),
c  taking into account that it has 2D Amplitude Spectrum simetries.
c  The region is defined by 4 points (counterclockwise) with
c  the following conditions:
c  Point 1 must be the left-most point (smallest column)
c  Point 2 must be the upper-most point (smallest row)
c  Point 3 must be the right-most point (largest column)
c  Point 4 must be the bottom-most point (largest row)
c
c  Matrix y(Lm,Ln) is the output.
c  IMUTE.eq.2: reflect the muting area with respect to the vertical axis.

	parameter       (mmax=256,nmax=256)
	real*8          x(mmax,nmax),y(mmax,nmax),
     #                  a1,a2,a3,a4,b1,b2,b3,b4,e
	data            e/0.001d0/

	nmini=min(n1,n2,n3,n4)
	nmaxi=max(n1,n2,n3,n4)
	mmini=min(m1,m2,m3,m4)
	mmaxi=max(m1,m2,m3,m4)
	if (n1.gt.nmini) then
	   write (*,*) 'Mute warning: 1 is not left-most point'
	   write (*,*) 'n1,n2,n3,n4: ',n1,n2,n3,n4
c	   return
	endif
	if (m2.gt.mmini) then
	   write (*,*) 'Mute warning: 2 is not upper-most point'
	   write (*,*) 'm1,m2,m3,m4: ',m1,m2,m3,m4
c	   return
	endif
	if (n3.lt.nmaxi) then
	   write (*,*) 'Mute warning: 3 is not right-most point'
	   write (*,*) 'n1,n2,n3,n4: ',n1,n2,n3,n4
c	   return
	endif
	if (m4.lt.mmaxi) then
	   write (*,*) 'Mute warning: 4 is not bottom-most point'
	   write (*,*) 'm1,m2,m3,m4: ',m1,m2,m3,m4
c	   return
	endif

	a1=dble(m2-m1)/dble(n2-n1+e)
	b1=dble(n2-n1)/dble(m2-m1-e)
	a2=dble(m3-m2)/dble(n3-n2+e)
	b2=dble(n3-n2)/dble(m3-m2+e)
	a3=dble(m4-m3)/dble(n4-n3-e)
	b3=dble(n4-n3)/dble(m4-m3+e)
	a4=dble(m1-m4)/dble(n1-n4-e)
	b4=dble(n1-n4)/dble(m1-m4-e)

	do n=1,Ln
	   m11=m1+idnint(a1*(n-n1))
	   m21=m2+idnint(a2*(n-n2))
	   m31=m3+idnint(a3*(n-n3))
	   m41=m4+idnint(a4*(n-n4))
	   do m=1,Lm/2+1
	      y(m,n)=x(m,n)
	      if (m.ge.m11.and.m.ge.m21.and.m.le.m31.and.m.le.m41) then
	         y(m,n)=0.d0
	         y(Lm-m+1,Ln-n+1)=0.d0
		 if (IMUTE.eq.2) then
		    y(Lm-m+1,n)=0.d0
		    y(m,Ln-n+1)=0.d0
		 endif
	      endif
	   enddo
	enddo
	return
	end
