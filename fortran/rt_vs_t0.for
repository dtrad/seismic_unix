c       *********************************************************
c
c       PROGRAM: hr_vs.for
c
c       HIGH RESOLUTION VELOCITY GATHERS 
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
c       Given a CMP computes a velocity gather. The velocity gather 
c       is then maped back to the CMP space. Undesired components
c       can be masked in the velocity gather before the reconstruction 
c       of the CMP gather. Doing so, we can isolate and keep only
c       part of the information contain in the CMP gather.
c       
c
c       Input parameters:
c
c         The following parameters are read from unit 10 called
c         'hr_vs.in'
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
c         vmin            - minimum velocity to invert
c         vmax            - maximum velocity to invert
c         nq              - number of traces of the velocity gather
c         iter_end        - maximum number of iterations (5 is 0.k.)
c         nw              - number of time-velocity windows that we desire
c                           to mask in the velocity gather (mask = filter out). 
c                           If nw=0 the velocity gather is not filtered.
c                           The following 4 parameters define each window:
c         tw1(nw) tw2(nw) - time range where the velocity is masked.
c         vw1(nw) vw2(nw) - range of velocities to mask.    

c        Output parameters:

c          m(nt,nq)       - Velocity gather. This is written in unit 30 (fileout_1)
c          d_rec(nt,nq)   - reconstructed CMP from the velocity gather.
c                           This is written in unit 40 (fileout_2)
c
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

	  include "c:\daniel\radon\flibrary\cholesky.for"
	  include "c:\daniel\radon\flibrary\herm.for"
	  include "c:\daniel\radon\flibrary\fx_go.for"
	  include "c:\daniel\radon\flibrary\fx_symmetry.for"
	  include "c:\daniel\radon\flibrary\fft.for"
	  include "c:\daniel\radon\flibrary\t_2.for"
c	  include "c:\daniel\radon\flibrary\forward.for"
	  include "c:\daniel\radon\flibrary\matrix_2.for"
c	  include "c:\daniel\radon\flibrary\cauchy_gauss_T0.for"
	  include "c:\daniel\radon\flibrary\gauss_gauss_0.for"
	  include "c:\daniel\radon\flibrary\sum_velocity.for"  
	  include "c:\daniel\radon\flibrary\v_stack_T0.for" 

        parameter (nnt=2048,nnx=128)

        real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,qmin,qmax, pos(nnx)
        real * 8 near_offset,dt,dh,t0,dq,vmin,vmax 
        real * 8 tw1(10),tw2(10),vw1(10),vw2(10)
        real * 8 offset_r,hmax_r,dh_r
         
        character * 20 filein,fileout_1,fileout_2, fileout_3

        open(10,file='hr_vs_1.in',status='unknown')

c------ read data
 
        read(10,*) filein      ! data
        read(10,*) fileout_1   ! velocity gather is written here.
        read(10,*) fileout_2   ! data computed from the velocity gather.
        read(10,*) fileout_3   ! data computed from the velocity gather.
	  read(10,*) eps1
        read(10,*) vmin,vmax
        read(10,*) nq,nh,nh_r,nt
	  read(10,*) dh,dt
	  read(10,*) near_offset,t0
        read(10,*) iter_end 
        read(10,*) nw 
        if (nw.ne.0) read(10,*) (tw1(i),tw2(i),vw1(i),vw2(i),i=1,nw)
     

        open(20,file=filein,form='binary',status='unknown')
        open(30,file=fileout_1,form='binary',status='unknown')
        open(40,file=fileout_2,form='binary',status='unknown')
        open(50,file=fileout_3,status='unknown')
       
	  do ih=1,nh
		do it=1,nt
			read(20) d(it,ih)
		enddo
		read(50,*) pos(ih)
	  enddo


c------ maximum and minimum q parameter to invert

        qmin=1./vmax**2
        qmax=1./vmin**2
        dq=(qmax-qmin)/(nq-1) ! q interval

c-----  compute the velocty stack 

c      *************************************************

        call V_STACK_T0(d,nh,pos,nt,dt,m,nq,qmin,qmax,eps1,iter_end)

c      *************************************************

       
c       write the  velocity gather in unit 30. dq*1000**2 is the
c       sampling interval of q in sec^2.km^2 = mscec^2/m^2


        write(30) nq,nt
        write(30) dq*1000**2,dt
        write(30) qmin*1000**2,t0
        do 12 iq=1,nq
12      write(30) (m(it,iq),it=1,nt)

 
c-----   THE UNDESIRE PARTS OF THE VELOCITY GATHER ARE
c-----   FILTERED HERE. THEN IT IS TRANSFORM BACK TO
c-----   THE OFFSET SPACE. 

        if (nw.ne.0) then 
	        do iw=1,nw
			  it1=idint(tw1(iw)/dt+1)
			  it2=idint(tw2(iw)/dt+1)
			  iq1=idint((1./vw2(iw)**2-qmin)/dq)
			  iq2=idint((1./vw1(iw)**2-qmin)/dq)
			  do it=it1,it2 
				do iq=iq1,iq2
					m(it,iq)=0.d0
				enddo
			  enddo
			enddo
	  endif 

c-----  the velocity gather is used to recostruct the CMP gather
c-----  you can change the parametrs of the recostruction. 
c-----  i.e., you can use a different near_offset and 
c-----  and a different nunber of traces. This is used 
c-----  to extrapolate near and far offset traces.


           write(*,*)nh_r

c       ********************************

        call sum_velocity(d_rec,nh_r,nt,m,nq,
     #     qmin,qmax,dt,dh,near_offset)

c       ********************************

c------ write the reconstructed CMP gather. d_rec(nt,nh)

        write(40) nh_r,nt
        write(40) dh,dt
        write(40) near_offset,0.d0
        do 14 ih=1,nh_r
14      write(40) (d_rec(it,ih),it=1,nt)


31      format(1f15.5)
32	  format(1f15.5)	
	  close(20)
        close(30)
	  close(40)
        close(50)
	  stop
        end 

c-----------------------------------------------------------------------------
