        subroutine  ss_stack_forw(pos,d,nh,nt,dt,m,nq,qmin,qmax)
c
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c       E-mail: sacchi@geop.ubc.ca
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047

c      Given a CMP gather compute the Velocity Stack
c
c      Input parameters:
c
c        d(nt,nh)   - CMP gather
c        nt         - number of samples of each trace
c        nh         - number of traces
c        near_offset- position of the near offset trace in meters.
c        nq         - number of velocity traces to estimate
c        qmin       - minimum q parameter  (sec^2/m^2)
c        qmax       - maximum q parameter  (sec^2/m^2)
c      Output paramters:
c 
c        m(nt,nq)   - VELOCITY GATHER
c 
c      Notes:
c  
c        The array m(nt,nq) is  organized as follows:
c        m(it,iq), it=1,nt, iq=1,nq
c                    iq=1 is the first trace of the velocity
c                    gather corresponding to q parameter qmin=1./vmax^2.
c                    iq=nq is the last trace of the velocity
c                    gather corresponding to a q parameter qmax=1./vmin^2.
c                    Where vmin and vmax are the min and max velocities
c                    that you expect to see.
c
 
        parameter    (nnx=228, nnt=2048)

        complex * 16 l(nnx,nnx),lh(nnx,nnx),pos(nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx)
        real * 8     tmax,dt,dq,dh,near_offset,qmin,qmax
        real * 8     f_low, f_high, fn, w ,pi

        pi=4.d0*datan(1.d0)
        
        Tmax=dt*(nt-1)
        dq=(qmax-qmin)/dfloat(nq-1) !sampling interval of the variable q

c------ compute the FFT number

        do 30 j=1,20
        nfft=2**j
        if(nfft.ge.nt) goto 40
30      continue
40      continue
        n2=nfft/2+1

        call FX_go (nq,nt,nfft,m,aux2,-1)

        fn=1./(2.d0*dt)
        f_low=0.
        f_high=fn
        kl=f_low*dfloat(nfft)*dt+1
        kh=f_high*dfloat(nfft)*dt+1
        kl=kl+2
        kh=kh-2

        do 1000 if=kl,kh

        do iq=1,nq
        mc(iq)=aux2(if,iq)  
        enddo

        w=2.d0*pi*dfloat(if-1)/(dfloat(nfft)*dt) ! freq in rad/sec
c       subroutine matrix_2(L,lh,pmin,pmax,pos,np,nh,w)
        call matrix_1(L,lh,qmin,qmax,pos,nq,nh,w)
      

        do ih=1,nh
        dc(ih)=dcmplx(0.d0,0.d0)
        do iq=1,nq
        dc(ih)=dc(ih)+LH(ih,iq)*mc(iq)
        enddo
        enddo
        
        do ih=1,nh
        aux2(if,ih)=dc(ih)
        enddo

1000    continue

c------ impose symmetry to the frequency axis

        call FX_symmetry (aux2,nh,nt,nfft,kl,kh)

c------ transform the freq-q space to tau-q space

        call FX_go (nh,nt,nfft,d,aux2, 1)   

        return 
        end 