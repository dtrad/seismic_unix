
      subroutine V_STACK_T0(d,nh,pos,nt,dt,m,nq,qmin,qmax,eps1,iter_end)
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
c        nh         - number of traces
c        pos		  - position of the traces in meters.
c        nt         - number of samples of each trace
c        dt         - time interval
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

 
        parameter    (nnx=128, nnt=2048)

        complex * 16 l(nnx,nnx),lh(nnx,nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx)
        real * 8     tmax,dt,dq,pos(nnx),qmin,qmax,eps1
        real * 8     dts ,f_low, f_high, fn, w ,pi

        pi=4.d0*datan(1.d0)
        
        Tmax=dt*(nt-1)
        tsmax=tmax**2
        dts=tsmax/dfloat(nt)        !sampling rate in t**2
        dq=(qmax-qmin)/dfloat(nq-1) !sampling interval of the variable q
        nts=nt

c------ compute the FFT number

        do 30 j=1,20
        nfft=2**j
        if(nfft.ge.nts) goto 40
30      continue
40      continue
        n2=nfft/2+1

c------ apply a t^2 transformation to have
c------ parabolic moveouts.
 
        call t_2 (nh,NT,Dt,nts,dts,d,1)

c------ Transform data to f-x 

        call FX_go (nh,nts,nfft,d,aux2,-1)

c------ Define the band in which the inversion is carried out

        fn=1./(2.d0*dts)
        f_low=0.
        f_high=fn
        kl=f_low*dfloat(nfft)*dts+1	!df=1/N/dt-->kl=f/df=f.N.dt
        kh=f_high*dfloat(nfft)*dts+1
        kl=kl+2
        kh=kh-2

c------ Make the inversion for each frequency

        do 1000 if=kl,kh
        print*,'freq index ',if
c------ dc contains the data at frequency if

        do ih=1,nh
        dc(ih)=aux2(if,ih)  
        enddo

        w=2.d0*pi*dfloat(if-1)/(dfloat(nfft)*dts) ! freq in rad/sec
												! (index to value)
c------ compute the linear operator L, that links the data with the velocity gather
c------ and LH, that links model to data.    m=LU,    d=LHm	

        call matrix_2(L,LH,qmin,qmax,pos,nq,nh,w)
					
c------ solve dc=L.mc using Cauchy-Gauss model
c------ subroutine CAUCHY_GAUSS_0(l,lh,d,model,nd,np,eps1,iter_end)
	  call GAUSS_GAUSS_0(l,lh,dc,mc,nh,nq,eps1)
  
c        call CAUCHY_GAUSS_T0(l,lh,dc,mc,nh,nq,eps1,iter_end)
	                           
c------ aux2(if,iq) contains the freq-q space 
 
        do iq=1,nq
        aux2(if,iq)=mc(iq)
        enddo

1000    continue

c------ impose symmetry to the frequency axis

        call FX_symmetry (aux2,nq,nts,nfft,kl,kh)

c------ transform the freq-q space to tau^2-q space

        call FX_go (nq,nts,nfft,m,aux2, 1)   

c------ Undo the t^2 transformation 
c------ m(nt,nq) is the Velocity gather in tau-q space 

        call t_2 (nq,nt,dt,nts,dts,m,-1)

       
        return 
        end 


