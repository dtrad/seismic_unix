        
        subroutine sum_slope(d,nh,nt,pos,v,nq,qmin,qmax,dt,ipower)

c       Compute the data  by summation over slopes
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

        parameter (nnt=2048, nnx=128)

        real * 8 d(nnt,nnx),qmin,qmax,q,dt,v(nnt,nnx) 
        real * 8 h,tau,t,aux,dq,pos(nnx)
	  integer ipower

        dq=(qmax-qmin)/(nq-1)
        do 100 ih=1,nh 
        h=pos(ih)
        do 100 it=1,nt
        d(it,ih)=0.d0
        t=(it-1)*dt
        do 200 iq=1,nq
        q=qmin+(qmax-qmin)*dfloat(iq-1)/dfloat(nq-1)
        tau=t-(h**ipower)*q
        tau=tau/dt
        if(tau.gt.1.d0) then
        i1=idint(tau)
        i2=i1+1
        aux=v(i1,iq)+
     #        (v(i2,iq)-v(i1,iq))*(tau-dfloat(i1))/dfloat(i2-i1)
        else
        aux=0.d0
        endif
        d(it,ih)=d(it,ih)+aux
200     continue
100     continue

        return
        end


