
        subroutine  P_STACK(method,pos,d,nh,nt,dt,
     #              m,nq,qmin,qmax,dq,eps1,eps2,iter_end)
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

        complex * 16 l(nnx,nnx),lh(nnx,nnx)
        complex * 16 mc(nnx),dc(nnx),aux2(nnt,nnx)
        real * 8     d(nnt,nnx),m(nnt,nnx),pos(nnx),dh(nnx)
        real * 8     tmax,dt,dq,near_offset,qmin,qmax,eps1
        real * 8     eps2,f_low, f_high, fn, w ,pi,Jtot(20),Jiter(20)
	  integer i, ji

        pi=4.d0*datan(1.d0)
        open(27,file='Jcost')
        Tmax=dt*(nt-1)
        dq=(qmax-qmin)/dfloat(nq-1) !sampling interval of the variable q

	  do i=1,nh	
		if (i.eq.1) then
			dh(i)=dabs(pos(i+1)-pos(i))
		elseif (i.eq.nh) then
			dh(i)=dabs(pos(i)-pos(i-1))
		else
			dh(i)=dabs((pos(i+1)-pos(i-1)))/2
		endif
		if (dh(i).lt.1) dh(i)=1.
	  enddo
	  

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

        do ji=1,iter_end
           Jiter(ji)=0   ! Total cost function per iteration
        enddo

        do 1000 if=kl,kh

c       print*,'Freq Index....',if
c------ dc contains the data at frequcny if

        do ih=1,nh
        dc(ih)=aux2(if,ih)  
        enddo

        w=2.d0*pi*dfloat(if-1)/(dfloat(nfft)*dt) ! freq in rad/sec

c------ compute the linear operator for non-equally spaced data
	
        call matrix_2(L,LH,qmin,qmax,pos,nq,nh,w,dh)
c****************************************************************
c------ solve dc=L.mc using Cauchy-Gauss model
c
c------sub GAUSS_GAUSS_0(l,lh,d,model,nd,np,eps1)
c------sub CAUCHY_GAUSS_0(l,lh,d,model,nd,np,eps1,iter_end)

        if(method.eq.2) then 
           call GAUSS_GAUSS_0(l,lh,dc,mc,nh,nq,eps1)
           call CAUCHY_GAUSS_0(l,lh,dc,mc,nh,nq,eps1,iter_end,Jtot)
           do ji=1,iter_end
           Jiter(ji)=Jiter(ji)+Jtot(ji)
           enddo
	  endif

        if(method.eq.1) then
        call GAUSS_GAUSS_T0(l,lh,dc,mc,nh,nq,eps1,dh)
        endif

        if(method.eq.3) then
           call graddesc(lh,l,dc,mc,nh,nq,eps1,iter_end,dh,dq,if,Jtot)
           do ji=1,iter_end
           Jiter(ji)=Jiter(ji)+Jtot(ji)
           enddo
        endif

        if(method.eq.4) then
           call conjgrad(lh,l,dc,mc,nh,nq,eps1,iter_end,dh,dq,if,Jtot)
           do ji=1,iter_end
           Jiter(ji)=Jiter(ji)+Jtot(ji)
           enddo
        endif

        if (method.eq.0) then
             call cAtimesx(mc,l,dc,nq,nh,nnx,nnx)
        endif
        
c****************************************************************
c------ aux2(if,iq) contains the freq-q space 
 
        do iq=1,nq
        aux2(if,iq)=mc(iq)
        enddo

1000    continue

c------ impose symmetry to the frequency axis

        call FX_symmetry (aux2,nq,nt,nfft,kl,kh)

c------ transform the freq-q space to tau-q space

        call FX_go (nq,nt,nfft,m,aux2, 1) 
        do ji=1,iter_end
           write(0,*) 'Iter=',ji,'Total J=',Jiter(ji)
           write(27,*) Jiter(ji) 
        enddo
        close(27)     
        return 
        end 
