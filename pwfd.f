c---------------------------Rongfeng Zhang in July, in UBC------------
c
c                *** Physical Wavelet Frame Denoise ***
c
c                                             Copyright Rongfeng Zhang
c                                             Contact:zrfzrf@371.net
c---------------------------------------------------------------------
c     p0   :to store 2D wavelet
c     p00  :first get input data,then store all images in all scales
c     data :to store fft of input data
c     wt2  :fft of wavelet
c     buf  :work buffer, for example,buf=data*wt2
c
c     thresh:  factor that times the average of input in wavelet domain
c     headr,headi,heads: trace head of su data (240 bytes)
c     fname:name of the file that contains the data file names
c     tr   :one dimension bufer used to count number of traces in data
c     nt,nx:number of actual time samples and traces
c     n1,n2:they're power of 2,great than nt and nx
c
c     subroutines : rlft3(fft subroutine),fourn(fft subroutine)
c                   hpwavelet(double-side 2D hyperbolic wavelet)
c                   hpwavelet2(single-side 2D hyperbolic wavelet)
c     *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
c Note: if process double-side gather, please change all 
c           `call hpwavelet2'==> `call hpwavelet'
c---------------------------------------------------------------------
cIf modify,please write the date and your name just below:
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter(n1=2048,n2=128)
c      parameter(nt=512,nx=64,n1=512*2,n2=64*2)
        real*8 fac2
        complex cwt2(n1/2,n2),cdata(n1/2,n2),cbuf(n1/2,n2)
        complex wspeq(n2),dspeq(n2),bspeq(n2)
        real wt2(n1,n2),data(n1,n2),buf(n1,n2)
        equivalence (wt2,cwt2),(cdata,data),(cbuf,buf)
        real, allocatable :: p0(:,:),p00(:,:)

        real headr(60)
        integer headi(60),offset
        integer*2 heads(120)
        equivalence(headr,headi,heads)

        character*150 ifile,ofile
        logical ext,verbose
        real,pointer,dimension(:)::tr

        print*,'-------------------------------------------------------'
        print*,'           Denoise Program PWFD '
        print*,'                             writen by Rongfeng Zhang'
        print*,'                             July,2000'
        print*,' The file <fname> contains filenames needed processing'
        print*,'              fname should be located in current dir.'       
        print*,'-------------------------------------------------------'
        print*

        verbose=.false.
c        print*, ' Verbose?'
c        read*,verbose
c        mpart=2        ! the length of traces will be cut 1/mpart
 5      write(*,'(A$)') ' Input m (I will only process 1/m of the trace,
     + so m=2 means processing half of the trace) ->'
        read*,mpart
        if(mpart.lt.1.or.mpart.gt.8) goto 5 !reinput mpart
        write(*,'(A$)') ' Input threshold value(between 3 to 11,recomend 
     + 8) ->'
        read*,thresh
        if(thresh.gt.11 .or. thresh.lt.3)thresh=8

c////////////////////////////////////////////////////file
        open(22,file='fname')
        read(22,*) nf
        do ifname=1,nf
           read(22,'(a)')ifile
           write(ofile,'(A,A1)') trim(ifile),'z' !ofile=ifile+'z'
           print*,ifile
c////////////////////////////////////////////////////file
        open(11,status='scratch',form='binary')
        inquire(file=ifile,EXIST=ext)
        if(.not.ext) then
           print*,trim(ifile),' does not EXIST.'
           stop
        endif
        open(10,file=ifile,form='binary')
            read(10,end=98)heads
            nt=heads(58)
            dt=heads(59)
            allocate(tr(nt))
            read(10)tr
            nx=0
 10         nx=nx+1   !count length of input data file
            read(10,end=20)heads,tr
            goto 10
 20         rewind(10)
        nt=nt/mpart      !@get rid of the low part
            if(verbose)print*,'nt=',nt,'nx=',nx,' dt=',int(dt)
            allocate(p0(nt,nx),p00(nt,nx))
        offset=0
        do i=1,nx
           read(10)headr,tr
           p00(:,i)=tr(1:nt)
           dx=headi(10)-offset
           if(verbose)print*,'cdp=',headi(6),'   offset=',headi(10),
     +             '   step of offset=',dx
                heads(58)=nt   !@get rid of the low part
           write(11)headr
           offset=headi(10)
        enddo
        close(10)

        nt2=ceiling(nt/2.)
        nx2=ceiling(nx/2.)
        fac2=2./(n1*n2) !because fft sub rlft3 doesn't normalize
             fac2=fac2*fac2 !(nt2*nx2) is only for keep Amp.
             fac2=fac2/1000.!for stable
c test                  p00=0;        p00(nt2,nx2)=1
        data=0
        data(1:nt,1:nx)=p00

        sx=3400
        st=0.005
        dt=dt*0.000001  !0.004    
c        dt=0.002    !
        dx=100
c        dx=dx*0.3    !convert to meters  !100
        h=nt2*dt
        call rlft3(data,dspeq,n1,n2,1,1,1) ! data-->frq domain
c==------------------estimate threshold value-->cut
        call hpwavelet2(p0,sx,h,st,nt,nx,dt,dx)
        p0=p0/sqrt(sx*st)
        p0=cshift(p0,nt2,dim=1)
        p0=cshift(p0,nx2,dim=2)
        wt2=0
        wt2(1:nt2,1:nx2)=p0(1:nt2,1:nx2)
        wt2(n1-nt2+1:n1,1:nx2)=p0(nt2+1:nt,1:nx2)
c        wt2(1:nt2,n2-nx2+1:n2)=p0(1:nt2,nx2+1:nx)
c        wt2(n1-nt2+1:n1,n2-nx2+1:n2)=p0(nt2+1:nt,nx2+1:nx)
        wt2(1:nt2,n2-(nx-nx2-1):n2)=p0(1:nt2,nx2+1:nx)
        wt2(n1-nt2+1:n1,n2-(nx-nx2-1):n2)=p0(nt2+1:nt,nx2+1:nx)
c        print*,n2-n2+(nx-nx2-1),nx-nx2-1
        call rlft3(wt2,wspeq,n1,n2,1,1,1)
        cbuf=cdata*conjg(cwt2)
        bspeq=dspeq*conjg(wspeq)
        call rlft3(buf,bspeq,n1,n2,1,-1,-1)
             buf=buf*fac2
        cut=sum(abs(buf(1:nt,1:nx)))/(nt*nx)*thresh   !8
c==-----------------
        if(verbose)print*,'cut=',cut, '      Input cut->'
c             cut=0.385
        if(verbose)read*,cut
c        cut=0.04 !0.04        !depend on st

        open(10,file=ofile,form='binary')
c        write(10) p00
        p00=0
c-----------------------loop begin*****************************
        smax=0.015
        if(dt*nt.gt.4) smax=0.02
        do st=0.005,smax,0.005    
           st2=st*st*1e+5
        do sx=600,6000,800 !400
c        do sx=1400,6000,800 !400
           if(verbose)print*,'                       st=',st,'  sx=',sx
        call hpwavelet2(p0,sx,h,st,nt,nx,dt,dx)
        p0=p0/sqrt(sx*st)
c        write(10) p0
        p0=cshift(p0,nt2,dim=1)
        p0=cshift(p0,nx2,dim=2)
        wt2=0
        wt2(1:nt2,1:nx2)=p0(1:nt2,1:nx2)
        wt2(n1-nt2+1:n1,1:nx2)=p0(nt2+1:nt,1:nx2)
        wt2(1:nt2,n2-nx2+1:n2)=p0(1:nt2,nx2+1:nx)
        wt2(n1-nt2+1:n1,n2-nx2+1:n2)=p0(nt2+1:nt,nx2+1:nx)

        call rlft3(wt2,wspeq,n1,n2,1,1,1)
        cbuf=cdata*conjg(cwt2)
        bspeq=dspeq*conjg(wspeq)
        call rlft3(buf,bspeq,n1,n2,1,-1,-1)
             buf=buf*fac2
        where(abs(buf)<cut)buf=0
             buf(nt+1:n1,:)=0
             buf(1:nt,nx+1:nx)=0

c        do i=1,n1
c           do ii=1,n2
c              if(abs(buf(i,ii)).lt.cut) buf(i,ii)=0
c           enddo
c        enddo

        call rlft3(buf,bspeq,n1,n2,1,1,1)

        cbuf=cbuf*cwt2
        bspeq=bspeq*wspeq
        call rlft3(buf,bspeq,n1,n2,1,-1,-1)
             buf=buf*fac2
        p00=p00+(buf(1:nt,1:nx))/(st2) !sum each scale
c        print*,minval(buf(1:nt,1:nx)),maxval(buf(1:nt,1:nx))

c        write(10) (buf(1:nt,1:nx))
c        stop
        enddo !end of sx
        enddo !end of st
        p00=p00*(nt2*nx2)**2 !recover amplitude
c        write(10) real(wt2)
c        write(10) p00 !the file processed data

c su format output        
        rewind(11)
        do i=1,nx
           read(11) headr
           write(10)headr,p00(:,i)*1e+6
        enddo
        close(11)
        close(10)
        deallocate(p0,p00,tr)
        goto 99
 98   print*,trim(ifile), ' is an empty file. Process the next one.'
      close(11)
 99   enddo                     !end of loop of file
        close(22)
        stop

        end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of main program

c hiperbolic-Mexica Hat 2D wavelet
c p=hpwavelet(x,t,sx,h,st)
c   dx--space interval
c   dt--time interval
c   sx-- scale in x direction
c   h--determine the apix location in t direction of this wavelet    
c   st--scale in t directon
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine hpwavelet(p,sx,h,st,nt,nx,dt,dx)
      real p(nt,nx), sx,st,dx,dt
      double precision tmp,x,t
c      print*,sx,h,st,nt,nx,dt,dx
      p=0
      do ix=1,nx
         x=(ix-nx/2)*dx
         do i=1,nt
            t=(i-1)*dt
            tmp=((t-sqrt((x/sx)**2+h*h))/st)**2;
            p(i,ix)=((1-tmp)*exp(-tmp/2))
         enddo
      enddo
      end

c hiperbolic-Mexica Hat 2D wavelet with single side
c p=hpwavelet(x,t,sx,h,st)
c   dx--space interval
c   dt--time interval
c   sx-- scale in x direction
c   h--determine the apix location in t direction of this wavelet    
c   st--scale in t directon
      subroutine hpwavelet2(p,sx,h,st,nt,nx,dt,dx)
      real p(nt,nx), sx,st,dx,dt
      double precision tmp,x,t
c      print*,sx,h,st,nt,nx,dt,dx
      p=0
      do ix=nx/2,nx
         x=(ix-nx/2)*dx
         do i=1,nt
            t=(i-1)*dt
            tmp=((t-sqrt((x/sx)**2+h*h))/st)**2;
            p(i,ix)=((1-tmp)*exp(-tmp/2))
         enddo
      enddo
      end
     
ccccccccc debug
c        open(30,file='dtt.dat',form='binary')
c        write(30)abs(cwt2)
c        close(30)
c        stop
ccccccccc debug
      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX data(nn1/2,nn2,nn3),speq(nn2,nn3)
CU    USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX c1,c2,h1,h2,w
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=cmplx(sngl(wr),sngl(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5,29#(.

      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5,29#(.
