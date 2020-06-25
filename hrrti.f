      subroutine hrrti(pos,data,nh,nt,dt,model,qmin,qmax,nq,rtmethod)
c       subroutine: hrrti.f
c       HIGH RESOLUTION PARABOLIC TRANSFORMS
c       Mauricio D. Sacchi
c       Daniel Trad
c       E-mail: sacchi@geop.ubc.ca
c             : dtrad@geop.ubc.ca 
c       *********************************************************
c       GENERAL DESCRIPTION:
c       Given the seimic data  computes a parabolic stack. The P.Stack 
c       is then maped back to the data space. Undesired components
c       can be masked in the velocity gather before the reconstruction 
c       of the data. Doing so, we can isolate and keep only
c       part of the information contain in the original data.
	  	  	          
	implicit none	
	integer nnt,nnx,nh,nt,nq,i,ih,iq
	integer method,iter_end,it,rtmethod
        parameter (nnt=2048,nnx=228)
	real * 8 d(nnt,nnx),m(nnt,nnx)
        real * 8 pos(nnx),qmin,qmax,dt,dh,t0
        real*4 data(0:nnx*nnt-1),model(0:nnx*nnt-1)

        character*20 fileout_1,fileout_2,file_param

c------ read data
             
        fileout_2='rec_data.out'     ! data computed from the velocity gather
        file_param='data1.par'

        t0=0
c        print*,nh,nt,dt,qmin,qmax,nq

	open(40,file=fileout_2,status='unknown')  
	open(60,file=file_param,status='unknown')

        do ih=1,nq
           do it=1,nt    
              m(it,ih)=model((it-1)+(ih-1)*nt)
           enddo
	enddo
        
        if (rtmethod.eq.1) then
           call  H_STACK(pos,d,nh,nt,dt,m,nq,qmin,qmax)
        elseif (rtmethod.eq.2) then
           call  s_stacki(pos,d,nh,nt,dt,m,nq,qmin,qmax)
        endif   
        call norma(d,nt,nh)

        do ih=1,nh
           do it=1,nt
              data((it-1)+(ih-1)*nt)=sngl(d(it,ih))
              write(40,*) data((it-1)+(ih-1)*nt)
           enddo
        enddo
        
        close(40)
        close(60)
        return
        end 
 







