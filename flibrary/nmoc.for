        subroutine nmoc(data,nt,nh,pos,dt,t0,tv,v,nv,data_c,index)
        
c------ index .eq. 1  ==> data  --->data_c
c------ index .eq.-1  ==> data_c--->data

        parameter (nnt=2048,nnx=228)

        real * 8 data(nnt,nnx),data_c(nnt,nnx),v(10),tv(10)
        real * 8 t0,t,dt,tau,h,vi(nnt),pos(nnx)

        write(*,*) index
        if(index.ne.1.and.index.ne.-1) then
        write(*,*)'Check index in subr. nmoc'
        return
        endif
 
        it0=idint(t0/dt)
        do iv=1,nv-1
        it1=idint(tv(iv)/dt)+1 
        it2=idint(tv(iv+1)/dt)+1
        do i=it1,it2+1
        vi(i-it0)=v(iv)+(v(iv+1)-v(iv))*(i-it1)/(it2-it1)
        enddo
        enddo

        if(index.eq.1) then

        do it=1,nt
        do ih=1,nh
        data_c(it,ih)=0.0
        enddo
        enddo

c------ Apply a NMOC with the intermediate velocity

        do it=2,nt
        tau=(it-1)*dt+t0 
        do  ih=1,nh
        h=pos(ih)
        t=dsqrt(tau**2+(h**2)/vi(it)**2)-t0 
        t=t/dt
        i1=idint(t)
        i2=i1+1
        data_c(it-1,ih)=data(i1,ih)+
     #     (data(i2,ih)-data(i1,ih))*(t-dfloat(i1))/dfloat(i2-i1)
        enddo
        enddo

        endif

c------- restablish the NMO

        if(index.eq.-1) then

        do it=1,nt
        do ih=1,nh
        data(it,ih)=0.0
        enddo
        enddo

        do it=2,nt
        tau=(it-1)*dt+t0 
        do  ih=1,nh
        h=pos(ih)
        t=tau**2-(h**2)/vi(it)**2
        if(t.ge.0.d0) then
        t=dsqrt(t)-t0
        t=t/dt
        i1=idint(t)
        i2=i1+1
        data(it-1,ih)=data_c(i1,ih)+
     #     (data_c(i2,ih)-data_c(i1,ih))*(t-dfloat(i1))/dfloat(i2-i1)
        endif
        enddo
        enddo
       
        endif 
       
        return
        end

c-------------------------------------------------------------