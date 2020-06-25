!*********************************************************
!
!PROGRAM: hr_ss_2.for
!
!	  SLANT STAKCS
!
!Jun/1994
!
!Mauricio D. Sacchi, Geophysics and Astronomy, UBC
!
!E-mail: sacchi@geop.ubc.ca
!
!Tel (o): (604)-822-2267
!Tel (h): (604)-224-0949
!Fax:     (604)-822-6047
!
!	  Modified for non- evenly spaced data
!	  by Daniel Trad - December 98		
!E-mail: dtrad@geop.ubc.ca
!*********************************************************
!
!GENERAL DESCRIPTION:
!
!Given the slant stack computes the data after muting.
!Undesired components can be masked in the velocity gather before the reconstruction 
!of the data. Doing so, we can isolate and keep only
!part of the information contain in the original data.
!
!Input parameters:
!  The following parameters are read from unit 10 called
!  'hr_ss.in'
!
!  filein          - file that contains (unit 20):
!						d(nt,nh) - CMP gather
!  fileout_1       - file where the velocity gather is written (unit 30).
!  fileout_2       - file where the reconstructed CMP is written (unit 40).
!  filein          - file that contains offset axis:
!  esp1            - prewhitening in % (1-5% is ok)         
!  pmin            - minimum 1/velocity to invert
!  pmax            - maximum 1/velocity to invert
!  nq              - number of traces of the velocity gather
!  nh				- number of traces
!  nt				- number of samples per trace
!  dt				- sampling interval in sec
!  t0				- time where the window starts in sec
!                            must start in t0=0sec.
!  nw              - number of time-velocity windows that we desire
!                    to mask in the velocity gather (mask = filter out). 
!                    If nw=0 the velocity gather is not filtered.
!                    The following 4 parameters define each window:
!  tw1(nw) tw2(nw) - time range where the velocity is masked.
!  pw1(nw) pw2(nw) - range of slopes to mask.    
!
! Output parameters:
!
!   m(nt,nq)       - Velocity gather. This is written in unit 30 (fileout_1)
!   d_rec(nt,nq)   - reconstructed CMP from the velocity gather.
!                    This is written in unit 40 (fileout_2)
!
!            
! Notes:
!   
! The output is organized as follows:
! m(it,iq), it=1,nt, iq=1,nq
!             iq=1 is the first trace of the velocity
!             gather corresponding to q parameter qmin. qmin=1./vmax^2.
!             iq=nq is the last trace of the velocity
!             gather corresponding to q parameter qmax. qmax=1./vmin^2.
!             where vmin and vmax are the minimum and maximum stacking
!             velocity seek by the procedure.
!             it=1 corresponds to 0sec. and it=nt to the last sample.
!          

!	  include "c:\daniel\radon\flibrary\fx_go.for"
!	  include "c:\daniel\radon\flibrary\fx_symmetry.for"
!	  include "c:\daniel\radon\flibrary\fft.for"
!	  include "c:\daniel\radon\flibrary\matrix_1.for"
!	  include "c:\daniel\radon\flibrary\ss_stack_forw.for" 	
	     
        parameter (nnt=2048,nnx=228)

        real * 8 d(nnt,nnx),m(nnt,nnx),d_rec(nnt,nnx)
        real * 8 eps1,eps2,qmin,qmax, pos(nnx)
        real * 8 dt,t0,dq,pmin,pmax 
        real * 8 tw1(10),tw2(10),pw1(10),pw2(10),cn(nnx)
        integer niter
	   
        character * 20 filein,fileout_1,fileout_2,fileout_3,fileout_4

        open(10,file='hr_ss_2.in',status='unknown')

!--- read data
 
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
        
        open(30,file=fileout_1,action='read',form='binary',status='unknown')
        open(31,file=fileout_4,form='binary',status='unknown')
	    open(40,file=fileout_2,form='binary',status='unknown')
        open(50,file=fileout_3,status='unknown')      
	
        do ih=1,nh
			read(50,*) pos(ih)
		enddo

!--- maximum and minimum q parameter to invert

        qmin=pmin
        qmax=pmax
        dq=(qmax-qmin)/(nq-1) ! q interval
      

        read(30) nq,nt
        read(30) dq,dt
        read(30) qmin,t0
        do 12 iq=1,nq
12      read(30) (m(it,iq),it=1,nt)

 
!--   THE UNDESIRE PARTS OF THE VELOCITY GATHER ARE
!--   FILTERED HERE. THEN IT IS TRANSFORM BACK TO
!--   THE OFFSET SPACE. 

        if (nw.eq.0) goto  1

        do iw=1,nw
        it1=idint(tw1(iw)/dt+1)
        it2=idint(tw2(iw)/dt+1)
        iq1= idint((pw1(iw)-qmin)/dq)
        iq2= idint((pw2(iw)-qmin)/dq)
		iq2=min(iq2,nq)
		it2=min(it2,nt)
        do it=it1,it2 
        do iq=iq1,iq2
        m(it,iq)=m(it,iq)/2 !(abs(it2-it)*abs(iq2-iq)+1)
        enddo
        enddo
        enddo

        write(31) nq,nt
        write(31) dq,dt
        write(31) qmin,t0
        
	  do iq=1,nq
		write(31) (m(it,iq),it=1,nt)
	  enddo	

!--  the velocity gather is used to recostruct the CMP gather
!--  you can change the parametrs of the recostruction. 
!--  i.e., you can use a different near_offset and 
!--  and a different nunber of traces. This is used 
!--  to extrapolate near and far offset traces.


!********************************

1        call  ss_stack_forw(pos,d_rec,nh,nt,dt,m,nq,qmin,qmax)
!********************************

!--- write the reconstructed CMP gather. d_rec(nt,nh)

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
 

!--------------------------------------------------------------------------



