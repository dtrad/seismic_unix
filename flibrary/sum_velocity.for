        
        subroutine sum_velocity(d,nh,nt,v,nq,
     #     qmin,qmax,dt,dh,near_offset)

c       Compute the CMP gather by summation over the velocity space
c
c       Input parameters:
c      
c           v(nt,np)  -  velocity gather
c           dt        -  time interval in sec.
c           dh        -  space interval in meters
c           np        -  number of traces of the velocity gather 
c           qmin      -  minimum q parameter of the velocity gather(s^2/m^2)
c           qmax      -  maximum q parameter of the velocity gather(s^2/m^2) 
c         near_offset -  near offset in meters
c
c       Output parameters:
c
c         d(nt,nh)    - cdp gather 

        parameter (nnt=2048, nnx=228)

        real * 8 d(nnt,nnx),qmin,qmax,q,dt,dh,v(nnt,nnx) 
        real * 8 h,tau,t,aux,near_offset,dq

        dq=(qmax-qmin)/(nq-1)
        do 100 ih=1,nh 
        h=dfloat(ih-1)*dh+near_offset
        do 100 it=1,nt
        d(it,ih)=0.d0
        t=(it-1)*dt
        do 200 iq=2,nq
        q=qmin+(qmax-qmin)*dfloat(iq-1)/dfloat(nq-1)  !dq=deltaq/nq
        tau=t**2-(h**2)*q								!q=qmin+dq.index
        if(tau.gt.0.d0) then  
        tau=dsqrt(tau)/dt
        if(tau.ge.1d0) then
        i1=idint(tau-tw1/dt)
        i2=i1+1
        aux=v(i1,iq)+
     #        (v(i2,iq)-v(i1,iq))*(tau-dfloat(i1))/dfloat(i2-i1)
        else
        aux=0.d0
        endif
        d(it-1,ih)=d(it-1,ih)+aux
        endif 
200     continue
100     continue

        return
        end
