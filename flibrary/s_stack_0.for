
       subroutine S_STACK_0(d,nh,nt,pos,dt,m,nq,qmin,qmax,eps1,
     #    niter,option,dq)
	
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
c        nq         - number of velocity traces to estimate
c        qmin       - minimum q parameter  (sec^2/m^2)
c        qmax       - maximum q parameter  (sec^2/m^2)
c        eps1       - damping in % (1.-5% is o.k)
c        iter_end   - maximum number of iterations (5 is o.k.)
c
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

        complex * 16 l(nnx,nnx), lh(nnx,nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx), pos(nnx), cn(nnx)
        real * 8     tmax,dt,dq,qmin,qmax,eps1
        real * 8     eps2,f_low, f_high, fn, w ,pi
	  integer option
		
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

c------ Transform data to f-x 

        call FX_go (nh,nt,nfft,d,aux2,-1)

c------ Define the band in which the inversion is carried out

        fn=1./(2.d0*dt)
        f_low=0.
        f_high=fn
        kl=f_low*dfloat(nfft)*dt+1
        kh=f_high*dfloat(nfft)*dt+1
        kl=kl+2
        kh=kh-2

c------ Make the inversion for each frequency

        do 1000 if=kl,kh

c        print*, 'freq index', if
c------ dc contains the data at frequcny if

        do ih=1,nh
        dc(ih)=aux2(if,ih)  
        enddo

        w=2.d0*pi*dfloat(if-1)/(dfloat(nfft)*dt) ! freq in rad/sec

c------ compute the linear operator L, that linkd the data with the velocity gather

	  call matrix_1(L,lh,qmin,qmax,pos,nq,nh,w)
        
c------ solve dc=L.mc using Cauchy-Gauss model
	  if (option.eq.1) then
        call GAUSS_GAUSS_0(L,lh,dc,mc,nh,nq,eps1)
	  elseif (option.eq.2) then
	  call CAUCHY_GAUSS_0(l,lh,dc,mc,nh,nq,eps1,niter)
	  endif
        do iq=1,nq
        aux2(if,iq)=mc(iq)
        enddo

1000    continue

c------ impose symmetry to the frequency axis
c       subr FX_symmetry (sfx,nx,nt,nfft,iflow,ifhigh)

        call FX_symmetry (aux2,nq,nt,nfft,kl,kh)

c------ transform the freq-q space to tau-q space

        call FX_go (nq,nt,nfft,m,aux2, 1)   


       
        return 
        end 

