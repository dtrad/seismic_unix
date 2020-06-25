      subroutine hrrt(pos,data,nh,nt,dt,method,eps1,iter_end,qmin)
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
	integer nnt,nnx,nh,nh_r,nt,nq,i,ih,iq
	integer method,iter_end,it
        parameter (nnt=2048,nnx=228)
	real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax 
        real * 8 dt,dh,t0,dq 
        real * 8 pos(nnx),pos_r(nnx)
        real*8 dx_max, fmax, dx_av, qmaxt
        real*4 qq, data(0:nnx*nnt-1)

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
	  print*, 'dx_max=', dx_max, 'dx_av=', dx_av
          print*,'method', method
          print*,'eps1', eps1
          print*, 'iter_end',iter_end
	  if (dh.lt.0) then
		nh_r=nh 
		do ih=1,nh
			pos_r(ih)=pos(ih)
		enddo
	  else 
		do ih=1,nh_r
			pos_r(ih)=pos(ih)
		enddo
	  endif


c------ maximum and minimum q parameter to invert
	  fmax=1/(2*dt)
	  call radon_param(fmax,pos(1),pos(nh),dx_av,qmin,qmaxt,qmax,dq,nq)
	  qmax=qmax
	  nq=nq
	  print*,'dq=',dq,'qmaxt=',qmaxt,'nq=', nq,'qmin',qmin

	  
        call P_STACK(method,pos,d,nh,nt,dt,
     #              m,nq,qmin,qmax,eps1,eps2,iter_end)

c        do i=1,nq
c           qq=1e10*(qmin+(i-1)*dq)
c          write(60,*) qq
c        enddo
   

        do iq=1,nq
           do it=1,nt
              write(30,*) m(it,iq)
           enddo
        enddo

        call  H_STACK(pos_r,d_rec,nh_r,nt,dt,m,nq,qmin,qmax)

        call norma(d_rec,nt,nh_r)

        do ih=1,nh_r
           do it=1,nt
              data((it-1)+(ih-1)*nt)=sngl(d_rec(it,ih))
              write(40,*) data((it-1)+(ih-1)*nt)
           enddo
        enddo

        
        close(30)
        close(40)
        close(60)

        return
        end 
 










