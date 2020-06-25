 
        subroutine t_2(nh,nt,dt,nts,dts,d,index)

c       Index.eq. 1 then stretch the time axis
c       Index.eq.-1 then unstretch the time axis
       

c       when INDEX=1 this is the IN/OUT

c       Input parameters:
c
c         d(nt,nh)    - data before tranformation
c               ds    - original sampling rate (sec)
c
c        Output parameters
c
c         d(nts,nh)  - data after tranformation
c         dts        - sampling rate in the stretched axis (sec**2)
c
c        Note:
c
c        When index=-1 the IN/OUT are interchanged

  
        parameter (nnt=2048, nnx=228)
 
        real * 8 d(nnt,nnx),aux(nnt),dt,dts,t

        if(index.eq. 1) n=nts
        if(index.eq.-1) n=nt
        
        do 100 ih=1,nh
        do 200 i=1,n
        if(index.eq. 1) t=(dsqrt(dfloat(i)*dts))/dt
        if(index.eq.-1) t=     ((dfloat(i)*dt)**2)/dts
        if(t.ge.1.d0) then
        i1=idint(t)
        i2=i1+1
        aux(i)=d(i1,ih)+
     #            (d(i2,ih)-d(i1,ih))*(t-dfloat(i1))/dfloat(i2-i1)
        else
        aux(i)=0.d0
        endif
200     continue
        do 300 i=1,n
300     d(i,ih)=aux(i)
100     continue
                
        return
        end
