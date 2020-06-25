      SUBROUTINE ludcmp_c(a,n,np,indx,d)
	implicit none
      INTEGER n,np,indx(n),NMAX
      complex*16 d,a(np,np),TINY,dumc,sum
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      real*8 aamax,vv(NMAX),dumr
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dumr=vv(i)*cdabs(sum)
          if (dumr.ge.aamax) then
            imax=i
            aamax=dumr
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dumc=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dumc
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dumc=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dumc
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .
