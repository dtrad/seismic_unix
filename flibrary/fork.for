      subroutine fork(lx,x,signi)
C
      complex carg,cexp,cw,ctemp,x(2048)
      j=1
      sc=(1./lx)
      do 30 i=1,lx
      if(i.gt.j) go to 10
      ctemp=x(j)*sc
      x(j)=x(i)*sc
      x(i)=ctemp
10    m=lx/2
20    if(j.le.m) go to 30
      j=j-m
      m=m/2
      if(m.ge.1) go to 20
30    j=j+m
      l=1
40    istep=2*l
      do 50 m=1,l
      carg=(0.,1.)*(3.14159265*signi*(m-1))/l
      cw=cexp(carg)
      do 50 i=m,lx,istep
      ctemp=cw*x(i+l)
      x(i+l)=x(i)-ctemp
50    x(i)=x(i)+ctemp
      l=istep
      if(l.lt.lx) go to 40
      if(signi.eq.1.) goto 60
      do 70 i=1,lx
70    x(i)=x(i)*lx
60    return
      end

