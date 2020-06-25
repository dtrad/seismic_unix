c*********************************************************************/
c                                                                    */  
c This code read cdps and put each cdp to cg_real.f to process, the  */
c output is primaries.bin, multiples.bin and velocities.bin.         */
c                                                                    */
c Input:                                                             */
c                                                                    */
c 1. Input SEGY binary file;                                         */
c                                                                    */
c 2. Output filename;                                                */
c                                                                    */
c 3. To time; (Double traveling time)                                */
c                                                                    */
c 4. Number of velocities, nv;                                       */
c                                                                    */ 
c 5. First velocity, v0;                                             */
c                                                                    */
c 6. CGLS or Huber, I_type;                                          */
c                                                                    */
c 7. Maximium CG iterations, max_iter_cg;                            */
c                                                                    */
c 8. Maximium times of updating Q;                                   */
c                                                                    */
c 9. \mu for sparseness;                                             */
c                                                                    */  
c 10.\sigma for weighting model;                                     */ 
c                                                                    */
c 11.a for Huber norm, a*\sigma is a tuning point;                   */
c                                                                    */
c 12.5 pair of number to divide the velocity panel.                  */
c                                                                    */
c                                                                    */
c Output:                                                            */
c                                                                    */
c 1. Priaries with headers in file 'primaries.bin';                  */
c                                                                    */
c 2. Multiples with headers in file 'multiples.bin';                 */
c                                                                    */
c 3. Velocities with headers in file 'velocities.bin'.               */
c                                                                    */
c                                                                    */
c How to use it,                                                     */
c                                                                    */
c  gcc -c readfuc.c -o segyread.o                                    */
c  gcc -c filelength.c -o filelength.o                               */
c  f77  segyread.o filelength.o master.f -o junk                     */
c                                                                    */
c  You can also use a script file called 'runcdp', then type         */
c  'junk', which is the excutable file.                              */
c                                                                    */
c                                                                    */
c  Author:  Qiang Li,  April, 2000                                   */
c           University of Alberta                                    */
c                                                                    */
c*********************************************************************/ 
      parameter (ntr=1000, ncdp=2000,n1x=500,n1t=3000,nnx=n1t*n1x)
      character*20 filename,
     #             filename1, 
     #             primary_out,
     #             multiple_out,
     #             velocity_out
      character*9  keyword,
     #             keyword1
      character    tail*2

      integer keys(ntr*ncdp), t_range(5000)

      real  D(nnx), H(n1x, 60), prim(nnx), N_pre(nnx), M(nnx), p(5,2),
     #       v(n1t,n1x), sep(n1t), Diag(nnx), offset(n1x)
  
      open (11, file='master.in' ) !put to 1 input file.    
        read(11,1) filename1
        read(11,1) primary_out
        read(11,1) multiple_out
        read(11,1) velocity_out
        read(11,*) to
        read(11,*) dv
        read(11,*) v0
        read(11,*) I_type
        read(11,*) max_iter_cg
        read(11,*) max_iter_update_q
        read(11,*) epsilon
        read(11,*) sigma
        read(11,*) a
        read(11,*) ((p(i,j),j=1,2),i=1,5)
c        write(*,*) p
 1      format(a20)
               
       keyword1='cdp'                       ! for cdp and trace counting
c--------------------------------------------
       tail='\0'
       i=1                                  !count the length of the filename
       do 10 while(filename1(i:i).ne.' ' )
          i=i+1
 10    continue

       j=1                                  !count the length of the keyword
       do 20 while(keyword1(j:j).ne.' ' )
          j=j+1
 20    continue
       filename=filename1(1:i-1)//tail       !append a '\0' to the end
       keyword=keyword1(1:j-1)//tail
c--------------------------------------------
       call flength(filename, ll)          !length(in bytes) of the file 
       call segyread(filename,'ns\0',output,1,1751)  ! get ns/trace
       l_trace=240+4*output                 ! calulate trace length
       ntrace=ll/l_trace                    ! calculate # of traces 
       ns=output
       
       do i=1,ntrace
       call segyread(filename,keyword,output,i,l_trace)  ! find keyword needed
       keys(i)=output
       enddo

c-------------------------------------------
       itrace=1
       num=2
       t_range(1)=0

       do 30 while((itrace+1).le.ntrace)
          if ((itrace+1).eq.ntrace) then
              t_range(num)=itrace+1 
          else
              if (keys(itrace).ne.keys(itrace+1)) then
                t_range(num)=itrace ! total  traces in each cdp
                num=num+1
              endif
          endif
       
          itrace=itrace+1
 30    continue
c                write(*,*) (t_range(i),i=1,13)
       n_cdp=num-1                           !total cdp number
         
      call segyread(filename,'dt\0',output,1,l_trace)
      dt=output/1000000.
c      write(*,*) n_cdp
c -------------------------------------------------------------------      
        do i_cdp=1,n_cdp                       ! cdp  loop 

c     get the offset vecotor for the cdp,and dh,nh,nv
       do i_offset=t_range(i_cdp)+1,t_range(i_cdp+1)
          call segyread(filename,'offset\0',output,i_offset,l_trace)
          offset(i_offset-t_range(i_cdp))=output
       enddo

       dh=offset(2)-offset(1)
       nh=t_range(2)-t_range(1)
       nv=nh
       nt=ns
            call read_bin_trace(filename1,ns,t_range(i_cdp)+1,
     #                             t_range(i_cdp+1),D,H)
                                                   ! 'D' is a cdp of data 
c       write(*,*)   ns,dt,dh,nh,nv,(offset(i),i=1,20)     
c ----------------------------------------------
c            main processing from cg.f
c ----------------------------------------------
c     I_type = 1 LS with CG
c            = 3 LS + Diag is updated to achieve HIGH RES. Huber norm

      if (I_type.eq.1) then       ! LS + CG
      do i=1,nv*ns
      M(i) = 0.
      Diag(i) = sqrt(epsilon)
      enddo
      call CG_LS(ns,dt,t0,nh,offset,nv,dv,v0,M,D,max_iter_cg,Diag)
c Write the velocity panel out for this cdp
      call write_bin_trace(velocity_out,ns,t_range(i_cdp)+1,
     #                             t_range(i_cdp+1),M,H)
      endif



      if (I_type.eq.3) then     ! Huber
       do i=1,nv*nt
       M(i) = 0.
       enddo

      do i=1,nv*ns
      M(i) = 0.
      Diag(i) = sqrt(epsilon)
      enddo
      call CG_LS(ns,dt,t0,nh,offset,nv,dv,v0,M,D,max_iter_cg,Diag)

      b=a*sigma
      do iter=1,max_iter_update_q
        do i=1,nv*ns
           if (abs(M(i)).le.b) then
             Diag(i)=sqrt(epsilon)/sigma
            else
             Diag(i) =sqrt(epsilon*a/(abs(M(i))*sigma))
           endif
        enddo
       call CG_LS(ns,dt,t0,nh,offset,nv,dv,v0,M,D,max_iter_cg,Diag)
      enddo
c Write the velocity panel out for this cdp
      call write_bin_trace(velocity_out,ns,t_range(i_cdp)+1,
     #                             t_range(i_cdp+1),M,H)
      endif

c Write the velocity panel out for this cdp

c        test='velo'   
c      call  write_bin(test,ns,nv,M) ! vel_gather2.bin  
c         test='velo' 
c      call write_bin_trace(test,ns,t_range(i_cdp)+1,
c     #                             t_range(i_cdp+1),M,H)                            
c ----------------------------------------------
c           main processing ends
c ---------------------------------------------- 
c put vel into a matrix form
      do iv=1,nv
      do it=1,ns
      k=(iv-1)*ns+it
      v(it,iv)=M(k)
      enddo
      enddo

c construct the seprator for the primary and multiples 
	do k=1,p(1,1)                    ! vertically from '0' to the 
	 sep(k)=p(1,2)                    ! first point 
	enddo

	do mm=2,5
	do k=p(mm-1,1)+1,p(mm,1)
       	 sep(k)=(p(mm,2)-p(mm-1,2))*(k-p(mm-1,1))/(p(mm,1)-p(mm-1,1))
     #          +p(mm-1,2)
	enddo
	enddo

	do k=p(5,1),ns
	 sep(k)=p(5,2)
	enddo
c        write(*,*) (sep(i),i=1,900)
c make the primary part zero
	do i=1,ns                       
	do j=1,nv
	if (j.gt.sep(i)) v(i,j)=0
	enddo
	enddo
c Put velocity into a vector form	
      do iv=1,nv
      do it=1,ns
      k=(iv-1)*ns+it
      M(k)=v(it,iv)
      enddo
      enddo

c      test='velo1' 
c       call write_bin(test,nt,nv,M)
c      call write_bin_trace(test,ns,t_range(i_cdp)+1,
c     #                             t_range(i_cdp+1),M,H) 

c Recover multiples from vel panel, and substracted from original data
      call  cvs(ns,dt,t0,nh,offset,nv,dv,v0,M,N_pre,-1)  ! recover multiples

	do i=1,nv*ns
	  prim(i)=D(i)-N_pre(i)
	enddo

c Write primary and multiples into two files.
            call write_bin_trace(primary_out,ns,t_range(i_cdp)+1,
     #                             t_range(i_cdp+1),prim,H)
         
            call write_bin_trace(multiple_out,ns,t_range(i_cdp)+1,
     #                             t_range(i_cdp+1),N_pre,H)

             
       enddo   ! cdp loop

       stop 
       end
c     -----------------------------------------------------
c                  Subs   starts   
c     -----------------------------------------------------
      subroutine write_bin_trace(filename,nt,nx1,nx2,x,H)
      parameter (n1x=500,n1t=3000,nnx=n1t*n1x) 
      real   x(nnx), H(n1x, 60)
      real   d(500,3000)
      character * 20 filename
      open(unit=9,file=filename,access='direct',recl=4*(nt+60))
      do ix=nx1,nx2
      do it=1,nt
      k=(ix-nx1)*nt+it
      d(ix,it)=real(x(k))
      enddo
      enddo
      do ix=nx1,nx2
      write(9,rec=ix,err=2000) (H(ix,ih),ih=1,60),(d(ix,it),it=1,nt)
      enddo
2000  continue
      close(9)
      return
      end 

c --------------------------------------------------
      SUBROUTINE  read_bin_trace(file,nt,nx1,nx2,s,H)
      parameter (nmax=500,ntt=3000,ns=nmax*ntt)
      real ss(nmax,ntt),H(nmax,60),s(ns)
      character*20 file

      open(unit=2,file=file,access='direct',recl=4*(nt+60))
      do ix=nx1,nx2
       read(2,rec=ix,err=2000)(H(ix,ih),ih=1,60),(ss(ix,it),it=1,nt)
      enddo 
2000  continue
      do i=nx1,nx2
         do j=1,nt
            s(j+(i-nx1)*nt)=ss(i,j)
         enddo
      enddo
      close(2) 
      return
      end  

c------------------------------------------------------------
       subroutine CG_LS(nt,dt,t0,nh,offset,nv,dv,v0,x,y,
     #                  max_iter_cg,Diag) 
c------------------------------------------------------------
       parameter (n1x=500,n1t=3000,nnx=n1t*n1x) 

       real   y0(nnx)
       real   x(nnx), y(nnx), r(nnx)
       real   g(nnx), s(nnx), ss(nnx) 
c             !,test_M1(nnx),test_y1(nnx)! for dot product
       real   alpha, beta, gamma, gammam, dot
       real   rms,e,den,e_old
       real   offset(n1x)
       real   Diag(nnx)

       nx = nv*nt 
       ny = nh*nt

       open(1,file='RMS')
       do i=1,nx
       y(i+ny) = 0. 
       enddo
c          call do_it(0,nx,x,ny+nx,y0)      ! y=Ax, get y   
       call cvs(nt,dt,t0,nh,offset,nv,dv,v0,x,y0,0)
       do i=1,nx
       y0(ny+i) = Diag(i)*x(i)
       enddo
       do i=1,ny+nx
       r(i)=y(i)-y0(i)
       enddo
c          call  do_it(1,nx,g,ny+nx,r)
       call cvs(nt,dt,t0,nh,offset,nv,dv,v0,g,r,1)
       do i=1,nx
       g(i) =  g(i) + Diag(i)*r(ny+i)
       enddo

       e_old=dot(ny+nx,r,r)

       do i=1,nx
       s(i)=g(i)                      
       enddo
 
       gammam=dot(nx,g,g)              

       do iter=1,max_iter_cg       

       call cvs(nt,dt,t0,nh,offset,nv,dv,v0,s,ss,0)

       do i=1,nx
       ss(ny+i) = Diag(i)*s(i)
       enddo
 
       den=dot(ny+nx,ss,ss)             
       alpha=gammam/den  
 
       do i=1,nx      
       x(i)=x(i)+alpha*s(i)  
       enddo
 
       do i=1,ny+nx
       r(i)=r(i)-alpha*ss(i)            
       enddo

       call cvs(nt,dt,t0,nh,offset,nv,dv,v0,g,r, 1)
       do i=1,nx
       g(i) =  g(i) + Diag(i)*r(ny+i)
       enddo
 
       gamma=dot(nx,g,g)
       beta=gamma/gammam   
       gammam=gamma
 
       do i=1,nx
       s(i)=g(i)+beta*s(i)   
       enddo
 
       e=dot(ny+nx,r,r)

c       ff= 200.0*(abs(e_old-e))/(e_old+e)
       
       r1=sqrt(e_old/(ny+nx))
       rms=sqrt(e/(ny+nx)) 

       if (abs(r1-rms).lt.1e-4) then   !if the error doesn't goo down..
          goto 3
       endif

       e_old = e

       write(*,2)  iter, rms
       write(1,2)  iter, rms
       enddo

 2     format(i5,3x,f15.4)
 3     continue

       return
       end

c     ----------------------------------------------------------
      subroutine cvs(nt,dt,t0,nh,offset,nv,dv,v0,M,D,iconj)
c     -----------------------------------------------------------
c     This subroutine uses both the time variable stack operator
c     and the componsating wavelet for the velocity inverse problem, 
c     
c               LwM=D, and the conjugate, M=w^T*L^T*D
c
c     for the first the code convolve the Model with the wavelet
c     first, then the stack operator is applied; for the second
c     the Data are processd first by inverse stacking operator,
c     then the result is correlated with wavelet, the equation are
c     adjount since convolution and correlation are conjugate 
c     operation in data processing.
c     The code passed the dot product test. which is done by
c     
c     A*X1=b1, A^T*b2=X2
c
c     then dot(X1,X2)=dot(b1,b2), where dot() means dot product.
c
c     Input and Output:
c
c     iconj=1 M=w^T*L^T*D
c     iconj=0 LwM=D
c     M       Model
c     D       Data
c
c     Qiang Li, Febuary, 2000
c     -------------------------------------------------------------

      parameter (n1x=500,n1t=3000,nnx=n1t*n1x) 

      real     M(nnx),D1(nnx),D(nnx)
      real     Dout(nnx),d_trace(nnx)
      real     offset(n1x),wavelet(n1x),f(4)
      real*8    tau, t, h,v

       f(1)=5.   !f1  (5)1
       f(2)=15.   !f2  (15)5
       f(3)=30.  !f3  (50)
       f(4)=40.  !f4  (70)

      call BP_wavelet(dt,f,wavelet,nw)

c      call ricker(28.,dt,wavelet,nw)
      if(iconj.eq.1) then
        do k=1,nv*nt
          M(k)=0.0
        enddo
c     correlate each trace of D with wavelet and truncation, w^T*D
        do ih=1,nh
          do it=1,nt
           k=(ih-1)*nt+it
            d_trace(it)=D(k)
              enddo
               call contran(1,0,nw,wavelet,nt,Dout,d_trace)
          
              do it=1,nt
             k=(ih-1)*nt+it  
            D1(k)=Dout(it)
          enddo
       enddo
      else
        do k=1,nh*nt
           D(k)=0.0
        enddo
      endif

c     Stacking operation
       do itau=1,nt
       tau=t0+(itau-1)*dt
        do ih=1,nh
         h=offset(ih)
          do iv=1,nv
           v=v0+(iv-1)*dv
           t=(sqrt(tau*tau+h*h/(v*v))-t0)/dt
           it=1+0.5+t
           j =(ih-1)*nt+it
           k =(iv-1)*nt+itau
c           aa=(D(j+1)-D(j))*(t-int(t))+D(j)
            if(it.le.nt) then
              if(iconj.eq.1) then
                M(k)=M(k)+D1(j)
              else
                D(j)=D(j)+M(k)
              endif
            endif
         enddo
        enddo
      enddo

c     get D by convolve  wLm=D
      if(iconj.eq.0) then
        do ih=1,nh
           do it=1,nt
            k=(ih-1)*nt+it
            d_trace(it)=D(k)
           enddo
              call contran(0,0,nw,wavelet,nt,d_trace,Dout)
           do it=1,nt
              k=(ih-1)*nt+it  
              D(k)=Dout(it)
           enddo
        enddo
      endif 
      return
      end
     
c ------------- conolution and correlation by Clearbout ----
      subroutine contran(conj,sum,nx,xx,nb,bb,yy)
c     
c     conj=0 convolution;
c     conj=1 correlation
c     this code passed the dot product, works with subroutine
c     'conjzero' below.
c
      parameter (nn=3000) 
      integer conj, sum,ib,nb,ix,nx,ny
      real xx(500)      !input signal
      real bb(nn)      !operaor
      real yy(nn)      !(nx+nb-1)
      real yyy(nn)     !(nx+nb-1)
      real bbb(nn)     !(nb)

      ns = nx/2+1
      ny=nx+nb-1
      call conjzero(conj,sum,nb,bb,ny,yy)

      if (conj .eq.0) then
c..convo...
         do ib=1,nb
            do ix=1,nx
               yy(ib+ix-1)=yy(ib+ix-1)+bb(ib)*xx(ix)
               enddo
               enddo
c shift
     
      do i=1,ny-ns
      yyy(i) = yy(i+ns)
      enddo
      do i=ny-ns+1,ny
      yyy(i) = 0.
      enddo
      do i=1,ny
      yy(i) = yyy(i)
      enddo

      else
c..corre...
       do ib=1,nb
           do ix=1,nx
              bb(ib)=bb(ib)+yy(ib+ix-1)*xx(ix)
              enddo
              enddo
c  shift    
      do i=ns+1,nb
      bbb(i) = bb(i-ns)
      enddo
      do i=1,ns
      bbb(i) =0.
      enddo
      do i=1,nb
      bb(i)=bbb(i)
      enddo

      endif
      return
      end

c -------------------------------------------------------
C  used by contran, this is a zeroing subroutine

      subroutine conjzero(conj,add,nx,x,ny,y)
      integer conj, add, ix,nx,iy,ny
      real x(nx),y(ny)
      if (add.eq.0) then
         if(conj.eq.0) then
            do iy=1,ny
               y(iy)=0.
            enddo   
          else
             do ix=1,nx
                x(ix)=0.
             enddo
          endif
       endif
       return
       end

c     -----------------------------------------------------
      function dot(n,x,y)
c     Compute the inner product  dot=(x,y)
      parameter (n1x=500,n1t=3000,nnx=n1t*n1x) 
      real  x(nnx), y(nnx)
      real  val ,dot
      val=0.0
      do i=1,n 
      val = val + x(i)*y(i)
      enddo
      dot=val
      return
      end


c     ------------------------------------------------
c
c     This subroutines are to construct the band-pass
c     wavelet 
c     ------------------------------------------------
      subroutine BP_wavelet(dt,f,w,nw)


      parameter  (nnt=500)
      real     w(nnt),f(4) 
      complex  aux(nnt)

      nf = 256
      if1 = f(1)*nf*dt+1.1 
      if2 = f(2)*nf*dt+1.1 
      if3 = f(3)*nf*dt+1.1 
      if4 = f(4)*nf*dt+1.1

      do i=1,if1-1
      w(i) = 0.
      enddo
      do i=if1,if2
      w(i) = float(i - if1)/float(if2-if1)
      enddo
      do i=if2+1,if3-1
      w(i) = 1. 
      enddo
      do i=if3,if4
      w(i) = 1.-float(i-if3)/float(if4-if3)
      enddo
      do i=if4+1,nf/2+1
      w(i) = 0.
      enddo
      do i = nf/2+2,nf
      w(i) = w(nf+2-i)
      enddo     
      do i = 1,nf
      aux(i) = cmplx(w(i),0.)
      enddo

      call fft(nf,aux,1)

      time = 2./(f(3)-f(2))
      L = ifix(time/dt) 


      do i=1,L+1
      w(i) = real(aux(L-i+2))
      enddo

      do i=L+2,2*L+1
      w(i) = real(aux(i-L))
      enddo

      nw = 2*L+1
      pi2 = 8.*atan(1.)
      do i = 1,nw
      w(i) = w(i)* (0.54-0.46*cos(pi2*(i-1)/float(nw-1)))
      enddo

      return
      end
 
c     -----------------------------------------------------
      subroutine fft(lx,x,isigni)
      parameter (nnt=500)
      complex    carg,cw,ctemp,x(nnt)
      real       signi,sc
      signi=float(isigni)
      j=1
      do 30 i=1,lx
      if(i.gt.j) go to 10
      ctemp=x(j)
      x(j)=x(i)
      x(i)=ctemp
10    m=lx/2
20    if(j.le.m) go to 30
      j=j-m
      m=m/2
      if(m.ge.1) go to 20
30    j=j+m
      l=1
40    istep=2*l
      do 50 m=1,l
      carg=(0.,1.)*(3.14159265*signi*(m-1))/l
      cw=cexp(carg)
      do 50 i=m,lx,istep
      ctemp=cw*x(i+l)
      x(i+l)=x(i)-ctemp
50    x(i)=x(i)+ctemp
      l=istep
      if(l.lt.lx) goto 40
      if(isigni.eq. 1 ) sc=1./float(lx)
      if(isigni.eq.-1 ) sc=1.
      do 70 i=1,lx
70    x(i)=x(i)*sc
60    return
      end

c     -----------------------------------------------------
      subroutine write_bin(filename,nt,nx,x)
      parameter (n1x=500,n1t=3000,nnx=n1t*n1x) 
      real   x(nnx)
      real   d(500,3000)
      character * 20 filename
      open(unit=9,file=filename,access='direct',recl=4*nt)
      do ix=1,nx
      do it=1,nt
      k=(ix-1)*nt+it
      d(ix,it)=real(x(k))
      enddo
      enddo
      do ix=1,nx
      write(9,rec=ix,err=2000) (d(ix,it),it=1,nt)
      enddo
2000  continue
      close(9)
      return
      end 

c --------------------------------------------------
      SUBROUTINE  read_bin(file,nt,nx,s,H)
      parameter (nmax=500,ntt=3000,ns=nmax*ntt)
      real ss(nmax,ntt),H(nmax,60),s(ns)
      character*20 file

      open(unit=2,file=file,access='direct',recl=4*(nt+60))
      do ix=1,nx
       read(2,rec=ix,err=2000)(H(ix,ih),ih=1,60),(ss(ix,it),it=1,nt)
      enddo 
2000  continue
      do i=1,nx
         do j=1,nt
            s(j+(i-1)*nt)=ss(i,j)
         enddo
      enddo
      close(2) 
      return
      end  

C------------------------------------------------------------     
c     This is a ricker wavelet with f central freq.
c     The length of it is wnw, since it is not only one cycle.
c     Sep. 10, 1999, Qiang Li
C------------------------------------------------------------

      subroutine ricker(f,dt,w,nl)
      real w(200) 
      wnw=7./f/dt
      nw=2*int(wnw/2.)+1
      nc=int(nw/2.)
      do i=1,nw
      alpha=(nc-i+1)*f*dt
      beta=alpha**2
      w(i)=(1.-beta*2)*exp(-beta)
      nl=nw
      enddo
     
      return
      end
