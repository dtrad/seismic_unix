      subroutine hrrtf(pos,data,nh,nt,dt,method,eps1,iter_end,qmin
     #,model,q,dq,nq,freq,rtmethod)
c       subroutine: hrrt.f
c       HIGH RESOLUTION PARABOLIC TRANSFORMS
c       Mauricio D. Sacchi
c       E-mail: sacchi@geop.ubc.ca
c       *********************************************************
c       GENERAL DESCRIPTION:
c       Given the seimic data  computes a parabolic stack. The P.Stack 
c       is then maped back to the data space. Undesired components
c       can be masked in the velocity gather before the reconstruction 
c       of the data. Doing so, we can isolate and keep only
c       part of the information contain in the original data.
	  	  	          
	implicit none	
	integer nnt,nnx,nh,nh_r,nt,nq,i,ih,iq,nqt,freq, rtmethod
	integer method,iter_end,it
        parameter (nnt=2048,nnx=228)
	real * 8 d(nnt,nnx),m(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax 
        real * 8 dt,dh,t0,dq,dqt 
        real * 8 pos(nnx)
        real*8 dx_max, fmax, dx_av, qmaxt
        real*4 data(0:nnx*nnt-1),model(0:nnx*nnt-1),q(0:nnx-1)

        character*20 fileout_1,fileout_2,file_param

c------ read data
             
        fileout_1='vel_gather.out'   ! velocity gather is written here
        fileout_2='rec_data.out'     ! data computed from the velocity gather
        file_param='data1.par'

        eps2=eps1
        nh_r=nh
        t0=0
        
        open(30,file=fileout_1,status='unknown')
	open(40,file=fileout_2,status='unknown')  
	open(60,file=file_param,status='unknown')
        do ih=1,nh
           do it=1,nt    
              d(it,ih)=data((it-1)+(ih-1)*nt)
           enddo
	enddo
        
	call interval(pos,nh,dx_max,dx_av)
	write(0,*)  'dx_max=', dx_max, 'dx_av=', dx_av
        write(60,*) 'method', method
        write(60,*) 'eps1', eps1
        write(60,*)  'iter_end',iter_end 




c------ maximum and minimum q parameter to invert
        if (freq.eq.0) then
           fmax=1/(2*dt)
        else
           fmax=freq
        endif

   
	call radon_param(fmax,pos(1),pos(nh),dx_av,qmin,qmaxt,qmax,
     #  dq,nq,rtmethod)
        write(0,*) 'freq max', fmax
        write(0,*) 'qmax theor=',qmaxt,'qmax used',qmax
        write(0,*) 'dq=', dq

	do ih=1,nh
           write(60,*) ih,pos(ih)
        enddo
        if (rtmethod.eq.1) then

        call P_STACK(method,pos,d,nh,nt,dt,
     #              m,nq,qmin,qmax,dq,eps1,eps2,iter_end)

        elseif (rtmethod.eq.2) then

        call S_STACK_0(d,nh,nt,pos,dt,m,nq,qmin,qmax,eps1
     # ,iter_end,method,dq)   

        endif

	write(0,*) 'nq=',nq,'nt',nt

        do i=1,nq
           q(i-1)=sngl(qmin+(i-1)*dq)
           write(60,*) q(i-1)
        enddo
        call norma(m,nt,nq)

        do ih=1,nq
           do it=1,nt
              model((it-1)+(ih-1)*nt)=sngl(m(it,ih))
              write(30,*) model((it-1)+(ih-1)*nt)
           enddo
        enddo

        
        close(30)
        close(40)
        close(60)

        return
        end 
 










