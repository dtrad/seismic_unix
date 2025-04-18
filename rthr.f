      subrouite
c
c       PROGRAM: hr_parab.for
c
c       HIGH RESOLUTION PARABOLIC TRANSFORMS
c
c       Jun/1994
c
c       Mauricio D. Sacchi, Geophysics and Astronomy, UBC
c
c       E-mail: sacchi@geop.ubc.ca
c
c       Tel (o): (604)-822-2267
c       Tel (h): (604)-224-0949
c       Fax:     (604)-822-6047
c
c
c       *********************************************************
c
c       GENERAL DESCRIPTION:
c
c       Given the seimic data  computes a parabolic stack. The P.Stack 
c       is then maped back to the data space. Undesired components
c       can be masked in the velocity gather before the reconstruction 
c       of the data. Doing so, we can isolate and keep only
c       part of the information contain in the original data.
c       
c
c       Input parameters:
c
c         The following parameters are read from unit 10 called
c         'hr_ss.in'
c
c         filein          - file that contains (unit 20):
c                              nh - number of traces
c                              nt - number of samples per trace
c                              dh - trace interval in meters
c                              dt - sampling interval in sec
c                     near_offset - position of the first trace in meters
c                              t0 - time where the window starts in sec
c                                   must start in t0=0sec.
c                        d(nt,nh) - CMP gather
c         fileout_1       - file where the velocity gather is written (unit 30).
c         fileout_2       - file where the reconstructed CMP is written (unit 40).
c         esp1            - prewhitening in % (1-5% is ok)         
c         qmin=nmo/h^2    - minimum velocity to invert
c         qmax=nmo/h^2    - maximum velocity to invert
c         nq              - number of traces of the velocity gather
c         iter_end        - maximum number of iterations (5 is 0.k.)
c         nw              - number of time-velocity windows that we desire
c                           to mask in the velocity gather (mask = filter out). 
c                           If nw=0 the velocity gather is not filtered.
c                           The following 4 parameters define each window:
c         tw1(nw) tw2(nw) - time range where the velocity is masked.
c         pw1(nw) pw2(nw) - range of slopes to mask.    

c        Output parameters:

c       m(nt,nq)       - Velocity gather. This is written in unit 30 (fileout_1)
c       d_rec(nt,nq)   - reconstructed CMP from the velocity gather.
c                        This is written in unit 40 (fileout_2)
c
c        Notes:
c      
c        The output is organized as follows:
c        m(it,iq), it=1,nt, iq=1,nq
c                    iq=1 is the first trace of the velocity
c                    gather corresponding to q parameter qmin. qmin=1./vmax^2.
c                    iq=nq is the last trace of the velocity
c                    gather corresponding to q parameter qmax. qmax=1./vmin^2.
c                    where vmin and vmax are the minimum and maximum stacking
c                    velocity seek by the procedure.
c                    it=1 corresponds to 0sec. and it=nt to the last sample.
c                 
c         The traces of the CMP (data) are considered equally spaced
c         the minumum  receiver-source distance is near_offset (meters)
c         The interval between traces of the CMP is dh. 
c  
c         Notice the each trace of the velocity gather correponds
c         to a parameter q=1/v^2. Therefore the panel gives the distribution
c         of energy in time-q, the varaiable q can then converted 
c         to velocity. In other words, you will see high 
c         velocities to the left (small q) and small velocities to
c         the right (big q).
c     
c         Everything is done by subr. 'V_STACK' the rest of the main is
c         for INPUT/OUTPUT and to compute the reconstructed CMP, this 
c         is done by subr. 'sum_velocity' 
c	  include "c:\daniel\radon\flibrary\cholesky.for"

c	  include "c:\daniel\radon\flibrary\cholesky2.for"
c	  include "c:\daniel\radon\flibrary\fx_go.for"
c	  include "c:\daniel\radon\flibrary\fx_symmetry.for"
c	  include "c:\daniel\radon\flibrary\fft.for"
c	  include "c:\daniel\radon\flibrary\nmoc.for"
c	  include "c:\daniel\radon\flibrary\norma.for"
c	  include "c:\daniel\radon\flibrary\matrix_2.for"
c	  include "c:\daniel\radon\flibrary\gauss_gauss_0.for"
c	  include "c:\daniel\radon\flibrary\cauchy_gauss_0.for"	 
c	  include "c:\daniel\radon\flibrary\p_stack.for"
c	  include "c:\daniel\radon\flibrary\h_stack.for"
c	  include "c:\daniel\radon\flibrary\sum_slope.for"
c	  include "c:\daniel\radon\flibrary\radon_param.for"
c	  include "c:\daniel\radon\flibrary\interval.for"
c	  include "c:\daniel\radon\flibrary\gauss_gauss_t0.for"
c	  include "c:\daniel\radon\flibrary\herm.for"
	  	  	          
	  implicit none	
	  integer nnt,nnx,nh,nh_r,nt,nv,nq,i,nw,ih
	  integer method,iter_end,index,iq,iq1,iq2,iw,it,it1,it2
        parameter (nnt=2048,nnx=228)
	  real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 d_nmoc(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax 
        real * 8 near_offset,dt,dh,t0,dq,pmin,pmax 
        real * 8 tw1(10),tw2(10),pw1(10),pw2(10)
        real * 8 tv(10),v(10),pos(nnx),pos_r(nnx)
        real*8 dx_max, fmax, dx_av, qmaxt
        real*4 singd

        character*20 filein,fileout_1,fileout_2,file_offset,file_param

        open(10,file='rtparab0.in',status='unknown')

c------ read data
 
        read(10,*) filein      ! data
        read(10,*) fileout_1   ! velocity gather is written here
        read(10,*) fileout_2   ! data computed from the velocity gather
        read(10,*) file_offset  ! offset data
        read(10,*) file_param

        read(10,*) eps1,eps2
        read(10,*) nh_r,nt
	  read(10,*) dh,dt
        read(10,*) near_offset,t0,qmin
	  read(10,*) iter_end
        print*, iter_end  
        read(10,*) nv
        if (nv.ne.0) read(10,*) (tv(i),v(i),i=1,nv)
        read(10,*) nw 
        if (nw.ne.0) read(10,*) (tw1(i),tw2(i),pw1(i),pw2(i),i=1,nw)

        open(20,file=filein,status='unknown')
        open(30,file=fileout_1,status='unknown')
	open(40,file=fileout_2,status='unknown')
	open(43,file='datanmo.out',status='unknown') 
	open(50,file=file_offset,status='unknown') 
	open(60,file=file_param,status='unknown')
         
	  do ih=1,1000
		read(50,*,end=1) pos(ih)
	  enddo
1	  continue
          print*,pos(ih-1)
	  nh=ih-1
	  print*,'nh=',nh

	  do ih=1,nh
		do it=1,nt    
                   read(20,*) d(it,ih)
		enddo
	  enddo


        
	  call interval(pos,nh,dx_max,dx_av)
	  print*, 'dx_max=', dx_max, 'dx_av=', dx_av

	  if (dh.lt.0) then
		nh_r=nh 
		do ih=1,nh
			pos_r(ih)=pos(ih)
		enddo
	  else 
		do ih=1,nh_r
			pos_r(ih)=(ih-1)*dh+near_offset
		enddo
	  endif
	  	
c------ apply a NMOC to correct the traces. 

        index=1
        call nmoc(d,nt,nh,pos,dt,t0,tv,v,nv,d_nmoc,index)
c        do  ih=1,nh
c		do it=1,nt
c			d_nmoc(it,ih)=d(it,ih)
c		enddo
c        enddo
	          
        write(60,*) nh,nt
        write(60,*) dh,dt
        write(60,*) near_offset,t0
        do  ih=1,nh
           do it=1,nt
             write(43,*) d_nmoc(it,ih)
          enddo
        enddo


c------ maximum and minimum q parameter to invert
	  fmax=1/(2*dt)
	  call radon_param(fmax,pos(1),pos(nh),dx_av,qmin,qmaxt,qmax,dq,nq)
	  qmax=qmax
	  nq=nq
	  print*,'dq=',dq,'qmaxt=', qmaxt, 'nq=', nq
	  pause	
	  
c-----  compute the slant stack

         write(*,2)  
2        format(3x,'method 1=ls 2=cg:  ',$)
         read(*,*) method

c      *************************************************
c-----subr.P_STACK(method,pos,d,nh,nt,dt,m,nq,qmin,qmax,eps1,eps2,iter_end)

        call P_STACK(method,pos,d_nmoc,nh,nt,dt,
     #              m,nq,qmin,qmax,eps1,eps2,iter_end)
c

c      *************************************************

        write(60,*) nq,nt
        write(60,*) dq,dt
        write(60,*) qmin,t0
        do iq=1,nq
           do it=1,nt
              write(30,*) m(it,iq)
           enddo
        enddo


 
c-----   THE UNDESIRE PARTS OF THE VELOCITY GATHER ARE
c-----   FILTERED HERE. THEN IT IS TRANSFORM BACK TO
c-----   THE OFFSET SPACE. 

        if (nw.ne.0) then
		do iw=1,nw
			it1=idint(tw1(iw)/dt+1)
			it2=idint(tw2(iw)/dt+1)
			iq1= idint((pw1(iw)-qmin)/dq)
			iq2= idint((pw2(iw)-qmin)/dq)
			do it=it1,it2 
				do iq=iq1,iq2
					m(it,iq)=0.d0
				enddo
			enddo
		enddo
	  endif

c-----  the velocity gather is used to recostruct the CMP gather
c-----  you can change the parameters of the recostruction. 
c-----  i.e., you can use a different near_offset and 
c-----  and a different nunber of traces. This is used 
c-----  to extrapolate near and far offset traces
c-----  sub sum_slope(d,nh,nt,pos,v,nq,qmin,qmax,dt,ipower).
c       call sum_slope(d_rec,nh_r,nt,pos_r,m,nq,qmin,qmax,dt,2)
c       sub  H_STACK(pos,d,nh,nt,dt,m,nq,qmin,qmax)
        call  H_STACK(pos_r,d_rec,nh_r,nt,dt,m,nq,qmin,qmax)

        index=-1
        call nmoc(d,nt,nh_r,pos_r,dt,t0,tv,v,nv,d_rec,index)                       
c------ write the reconstructed CMP gather. d_rec(nt,nh)


        call norma(d,nt,nh_r)

        write(60,*) nh_r,nt
        write(60,*) dh,dt
        write(60,*) near_offset,0.d0
        do ih=1,nh_r
           do it=1,nt
              write(40,*) d(it,ih)
           enddo
        enddo
	close(10)
        close(20)
        close(30)
        close(40)
        close(43)
        close(50)

        stop
        end 
 










