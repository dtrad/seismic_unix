c       *********************************************************
c
c       PROGRAM: hr_ss_2.for
c
c	  SLANT STAKCS
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
c	  Modified for non- evenly spaced data
c	  by Daniel Trad - December 98		
c       E-mail: dtrad@geop.ubc.ca
c       *********************************************************
c
c       GENERAL DESCRIPTION:
c
c       Given the seimic data  computes a slant stack. The S.Stack 
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
c						d(nt,nh) - CMP gather
c         fileout_1       - file where the velocity gather is written (unit 30).
c         fileout_2       - file where the reconstructed CMP is written (unit 40).
c         filein          - file that contains offset axis:
c         esp1            - prewhitening in % (1-5% is ok)         
c         pmin            - minimum 1/velocity to invert
c         pmax            - maximum 1/velocity to invert
c         nq              - number of traces of the velocity gather
c         nh				- number of traces
c         nt				- number of samples per trace
c         dt				- sampling interval in sec
c         t0				- time where the window starts in sec
c                                   must start in t0=0sec.
c         nw              - number of time-velocity windows that we desire
c                           to mask in the velocity gather (mask = filter out). 
c                           If nw=0 the velocity gather is not filtered.
c                           The following 4 parameters define each window:
c         tw1(nw) tw2(nw) - time range where the velocity is masked.
c         pw1(nw) pw2(nw) - range of slopes to mask.    
c
c        Output parameters:
c
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
!	  include "c:\daniel\radon\flibrary\cholesky.for"
!	  include "c:\daniel\radon\flibrary\fx_go.for"
!	  include "c:\daniel\radon\flibrary\fx_symmetry.for"
!	  include "c:\daniel\radon\flibrary\fft.for"
!	  include "c:\daniel\radon\flibrary\nmoc.for"
!	  include "c:\daniel\radon\flibrary\norma.for"
!	  include "c:\daniel\radon\flibrary\matrix_1.for"
!	  include "c:\daniel\radon\flibrary\cauchy_gauss_0.for"
!	  include "c:\daniel\radon\flibrary\sum_slope.for"  
!	  include "c:\daniel\radon\flibrary\s_stack_0.for" 
!	  include "c:\daniel\radon\flibrary\ss_stack_forw.for" 	
	     
        parameter (nnt=2048,nnx=128)

        real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax, pos(nnx)
        real * 8 dt,t0,dq,pmin,pmax 
        real * 8 tw1(10),tw2(10),pw1(10),pw2(10),cn(nnx)
        integer niter
	   
        character * 20 filein,fileout_1,fileout_2,fileout_3, fileout_4

        open(10,file='hr_ss_2.in',status='unknown')

c------ read data
 
        read(10,*) filein      ! data
        read(10,*) fileout_1   ! velocity gather is written here.
	  read(10,*) fileout_4   ! velocity gather is written here.
        read(10,*) fileout_2   ! data computed from the velocity gather.
        read(10,*) fileout_3   ! offset file 
        read(10,*) eps1,eps2
        read(10,*) pmin,pmax
        read(10,*) nq,nh,nt
	  read(10,*) dt,t0
        read(10,*) niter
	  read(10,*) nw 
        if (nw.ne.0) read(10,*) (tw1(i),tw2(i),pw1(i),pw2(i),i=1,nw)


        open(20,file=filein,form='binary',status='unknown')
        open(30,file=fileout_1,form='binary',status='unknown')
        open(31,file=fileout_4,form='binary',status='unknown')
	  open(40,file=fileout_2,form='binary',status='unknown')
        open(50,file=fileout_3,status='unknown')      
	
        do ih=1,nh
	    do it=1,nt
			read(20) d(it,ih)
		enddo
		read(50,*) pos(ih)
		cn(ih)=1/eps2  ! assign a value to the inverse of covariance
	  enddo

c------ maximum and minimum q parameter to invert

        qmin=pmin
        qmax=pmax
        dq=(qmax-qmin)/(nq-1) ! q interval

c-----  compute the slant stack

c      *************************************************

        call S_STACK_0(d,nh,nt,pos,dt,
     #              m,nq,qmin,qmax,eps1,cn,niter)

c      *************************************************
       

        write(30) nq,nt
        write(30) dq,dt
        write(30) qmin,t0
        do 12 iq=1,nq
12      write(30) (m(it,iq),it=1,nt)

 
c-----   THE UNDESIRE PARTS OF THE VELOCITY GATHER ARE
c-----   FILTERED HERE. THEN IT IS TRANSFORM BACK TO
c-----   THE OFFSET SPACE. 

        if (nw.eq.0) goto  1

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

        write(31) nq,nt
        write(31) dq,dt
        write(31) qmin,t0
        
	  do iq=1,nq
		write(31) (m(it,iq),it=1,nt)
	  enddo	

c-----  the velocity gather is used to recostruct the CMP gather
c-----  you can change the parametrs of the recostruction. 
c-----  i.e., you can use a different near_offset and 
c-----  and a different nunber of traces. This is used 
c-----  to extrapolate near and far offset traces.


c       ********************************
c       subroutine sum_slope(d,nh,nt,pos,v,nq,qmin,qmax,dt,ipower)

c1       call sum_slope(d_rec,nh,nt,pos,m,nq,qmin,qmax,dt,1)
1        call  ss_stack_forw(pos,d_rec,nh,nt,dt,m,nq,qmin,qmax)
c       ********************************

c------ write the reconstructed CMP gather. d_rec(nt,nh)

        write(40) nh,nt
        write(40) 0.d0,dt
        write(40) 0.d0,0.d0
        do 14 ih=1,nh
14      write(40) (d_rec(it,ih),it=1,nt)


31      format(1f15.5)
32      format(1f15.5)
	  close(31)
	  close(30)
	  close(40)	
        stop
        end 
 

c-----------------------------------------------------------------------------



