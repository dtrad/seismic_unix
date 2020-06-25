      program main
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
      parameter (n1 = 2048, n2 = 128)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 11:12:41 2/ 8/02
C...Switches:                     
c      parameter(nt=512,nx=64,n1=512*2,n2=64*2)
      double precision fac2
      complex  cwt2(1024,128), cdata(1024,128), cbuf(1024,128)
      complex  wspeq(128), dspeq(128), bspeq(128)
      real wt2(2048,128), data(2048,128), buf(2048,128)
      equivalence (wt2, cwt2), (cdata, data), (cbuf, buf)
 
      real headr(60)
      integer headi(60), offset
      integer*2 heads(120)
      equivalence (headr, headi, heads)
      character ifile*150, ofile*150
      logical ext, verbose
      double complex heap__(20)
      common /heap__f90/heap__
      integer tr1, tr2, tr3, p01, p02, p001, p002, ceiling_f90, heap__1,
     . j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17, j18, 
     .j19, j20, j21, j22, j25
      real r10, p0(9), p00(9), tr(9), r8(9), r9(9), r11(9), r12(9)
      equivalence (p0,heap__)
      integer p0_
      equivalence (p00,heap__)
      integer p00_
      equivalence (tr,heap__)
      integer tr_
      equivalence (r8,heap__)
      integer r8_
      equivalence (r9,heap__)
      integer r9_
      equivalence (r11,heap__)
      integer r11_
      equivalence (r12,heap__)
      integer r12_
      integer j3, j4
      data r12_/ 0/ 
      data r11_/ 0/ 
      data r9_/ 0/ 
      data r8_/ 0/ 
      heap__1 = %loc(heap__)
 
      data p00_/ 0/ 
      data p0_/ 0/ 
      print *, '-------------------------------------------------------'
      print *, '           Denoise Program PWFD '
      print *, '                             writen by Rongfeng Zhang'
      print *, '                             July,2000'
      print *, ' The file <fname> contains filenames needed processing'
      print *, '              fname should be located in current dir.'
      print *, '-------------------------------------------------------'
      print *
 
      verbose = .FALSE.
c        print*, ' Verbose?'
c        read*,verbose
c        mpart=2        ! the length of traces will be cut 1/mpart
    5 continue
      write (*, '(a$)') 
     1' Input m (I will only process 1/m of the trace, so m=2 means proc
     2essing half of the trace) ->'
      read * , mpart
      if (mpart.lt.1 .or. mpart.gt.8) go to 5    !reinput mpart
      write (*, '(a$)') 
     1   ' Input threshold value(between 3 to 11,recomend 8) ->'
      read * , thresh
      if (thresh.gt.11 .or. thresh.lt.3) thresh = 8
 
c////////////////////////////////////////////////////file
      open(22, file='fname')
      read (22, *) nf
      do ifname = 1, nf
         read (22, '(a)') ifile
         call trim_f90 (ifile, j21)
c                                                !ofile=ifile+'z'
         write (ofile, '(a,a1)') ifile(1:j21), 'z'
         print *, ifile
c////////////////////////////////////////////////////file
         open(11, status='scratch', form='binary')
         inquire(file=ifile, exist=ext)
         if (.not.ext) then
            call trim_f90 (ifile, j22)
            print *, ifile(1:j22), ' does not EXIST.'
            stop 
         endif
         open(10, file=ifile, form='binary')
         read (10, end=98) heads
         nt = heads(58)
         dt = heads(59)
         tr3 = 1
         tr1 = 1
         tr2 = nt
         tr_ = 0
         call allocate_f90 (4, tr2, 0, tr_, 0)
         tr_ = (tr_ - heap__1)/4
         read (10) (tr((j3-1)*tr3+1+tr_),j3=1,tr2 - tr1 + 1)
         nx = 0
   10    continue
         nx = nx + 1
c                                       !count length of input data file
         read(10,end=20)heads,(tr((j3-1)*tr3+1+tr_),j3=1,tr2-tr1+1)
         go to 10
   20    continue
         rewind (10)
         nt = nt/mpart                         !@get rid of the low part
         if (verbose) print *, 'nt=', nt, 'nx=', nx, ' dt=', int(dt)
         p01 = nt
         p02 = nx
         call allocate_f90 (4, p02*p01, 0, p0_, 0)
         p0_ = (p0_ - heap__1)/4
         p001 = nt
         p002 = nx
         call allocate_f90 (4, p002*p001, 0, p00_, 0)
         p00_ = (p00_ - heap__1)/4
         offset = 0
         do i = 1, nx
            read (10) headr, (tr((j3-1)*tr3+1+tr_),j3=1,tr2 - tr1 + 1)
            do j3 = 1, nt
               p00((i-1)*p001+j3+p00_) = tr(j3+tr_)
            end do
            dx = headi(10) - offset
            if (verbose) print *, 'cdp=', headi(6), '   offset=', headi(
     1         10), '   step of offset=', dx
            heads(58) = nt                     !@get rid of the low part
            write (11) headr
            offset = headi(10)
         end do
         close(10)
 
         nt2 = ceiling_f90(nt/2.)
         nx2 = ceiling_f90(nx/2.)
         fac2 = 2./(2048*128)   !because fft sub rlft3 doesn't normalize
         fac2 = fac2*fac2               !(nt2*nx2) is only for keep Amp.
         fac2 = fac2/1000.                       !for stable
c test                  p00=0;        p00(nt2,nx2)=1
         do j4 = 1, 128
            do j3 = 1, 2048
               data(j3,j4) = 0
            end do
         end do
         do j4 = 1, nx
            do j3 = 1, nt
               data(j3,j4) = p00((j4-1)*p001+j3+p00_)
            end do
         end do
 
         sx = 3400
         st = 0.005
         dt = dt*0.000001                        !0.004
c        dt=0.002    !
         dx = 100
c        dx=dx*0.3    !convert to meters  !100
         h = nt2*dt
c                                                ! data-->frq domain
         call rlft3 (data, dspeq, 2048, 128, 1, 1, 1)
c==------------------estimate threshold value-->cut
         call hpwavelet2 (p0(1+p0_), sx, h, st, nt, nx, dt, dx)
         do j4 = 1, p02
            do j3 = 1, p01
               p0((j4-1)*p01+j3+p0_) = p0((j4-1)*p01+j3+p0_)/sqrt(sx*st)
            end do
         end do
         j5 = p01
         j6 = p02
         call allocate_f90 (4, j6*j5, 0, r8_, 0)
         r8_ = (r8_ - heap__1)/4
         if (p01 .ne. 0) j8 = mod(nt2,p01)
         if (j8 .gt. 0) then
            do j4 = 1, p02
               do j3 = 1, p01 - j8
                  r8((j4-1)*j5+j3+r8_) = p0((j4-1)*p01+j3+j8+p0_)
               end do
               do j3 = p01 - j8 + 1, p01
                  r8((j4-1)*j5+j3+r8_) = p0((j4-2)*p01+j3+j8+p0_)
               end do
            end do
         else
            do j4 = 1, p02
               do j3 = 1 - j8, p01
                  r8((j4-1)*j5+j3+r8_) = p0((j4-1)*p01+j3+j8+p0_)
               end do
               do j3 = 1, -j8
                  r8((j4-1)*j5+j3+r8_) = p0(j4*p01+j3+j8+p0_)
               end do
            end do
         endif
         do j4 = 1, p02
            do j3 = 1, p01
               p0((j4-1)*p01+j3+p0_) = r8((j4-1)*j5+j3+r8_)
            end do
         end do
         r8_ = r8_*4 + heap__1
         call deallocate_f90 (4, j6*j5, 0, r8_, 0)
         j9 = p01
         j10 = p02
         call allocate_f90 (4, j10*j9, 0, r9_, 0)
         r9_ = (r9_ - heap__1)/4
         if (p02 .ne. 0) j12 = mod(nx2,p02)
         if (j12 .gt. 0) then
            do j3 = 1, p01
               do j4 = 1, p02 - j12
                  r9((j4-1)*j9+j3+r9_) = p0((j4+j12-1)*p01+j3+p0_)
               end do
               do j4 = p02 - j12 + 1, p02
                  r9((j4-1)*j9+j3+r9_) = p0((j4+j12-p02-1)*p01+j3+p0_)
               end do
            end do
         else
            do j3 = 1, p01
               do j4 = 1 - j12, p02
                  r9((j4-1)*j9+j3+r9_) = p0((j4+j12-1)*p01+j3+p0_)
               end do
               do j4 = 1, -j12
                  r9((j4-1)*j9+j3+r9_) = p0((j4+j12+p02-1)*p01+j3+p0_)
               end do
            end do
         endif
         do j4 = 1, p02
            do j3 = 1, p01
               p0((j4-1)*p01+j3+p0_) = r9((j4-1)*j9+j3+r9_)
            end do
         end do
         r9_ = r9_*4 + heap__1
         call deallocate_f90 (4, j10*j9, 0, r9_, 0)
         do j4 = 1, 128
            do j3 = 1, 2048
               wt2(j3,j4) = 0
            end do
         end do
         do j4 = 1, nx2
            do j3 = 1, nt2
               wt2(j3,j4) = p0((j4-1)*p01+j3+p0_)
            end do
            do j3 = 1, nt2
               wt2(j3+2048-nt2,j4) = p0((j4-1)*p01+j3+nt2+p0_)
            end do
         end do
         do j4 = 1, nx - nx2
            do j3 = 1, nt2
               wt2(j3,j4+128+nx2-nx) = p0((j4+nx2-1)*p01+j3+p0_)
            end do
            do j3 = 1, nt2
               wt2(j3+2048-nt2,j4+128+nx2-nx) = p0((j4+nx2-1)*p01+j3+nt2
     1            +p0_)
            end do
         end do
c        print*,n2-n2+(nx-nx2-1),nx-nx2-1
         call rlft3 (wt2, wspeq, 2048, 128, 1, 1, 1)
         do j4 = 1, 128
            do j3 = 1, 1024
               cbuf(j3,j4) = cdata(j3,j4)*conjg(cwt2(j3,j4))
            end do
         end do
         do j3 = 1, 128
            bspeq(j3) = dspeq(j3)*conjg(wspeq(j3))
         end do
         call rlft3 (buf, bspeq, 2048, 128, 1, -1, -1)
         do j4 = 1, 128
            do j3 = 1, 2048
               buf(j3,j4) = buf(j3,j4)*fac2
            end do
         end do
         r10 = 0
         do j4 = 1, nx
            do j3 = 1, nt
               r10 = r10 + abs(buf(j3,j4))
            end do
         end do
         cut = r10/(nt*nx)*thresh                !8
         if (verbose) then
c==-----------------
            print *, 'cut=', cut, '      Input cut->'
c             cut=0.385
            read * , cut
         endif
c        cut=0.04 !0.04        !depend on st
 
         open(10, file=ofile, form='binary')
c        write(10) p00
         do j4 = 1, p002
            do j3 = 1, p001
               p00((j4-1)*p001+j3+p00_) = 0
            end do
         end do
c-----------------------loop begin*****************************
         smax = 0.015
         if (dt*nt .gt. 4) smax = 0.02
         do st = 0.005, smax, 0.005
            st2 = st*st*1E+5
            do sx = 600, 6000, 800               !400
c        do sx=1400,6000,800 !400
               if (verbose) print *, '                       st=', st, 
     1            '  sx=', sx
               call hpwavelet2 (p0(1+p0_), sx, h, st, nt, nx, dt, dx)
               do j4 = 1, p02
                  do j3 = 1, p01
                     p0((j4-1)*p01+j3+p0_) = p0((j4-1)*p01+j3+p0_)/sqrt(
     1                  sx*st)
                  end do
               end do
c        write(10) p0
               j13 = p01
               j14 = p02
               call allocate_f90 (4, j14*j13, 0, r11_, 0)
               r11_ = (r11_ - heap__1)/4
               if (p01 .ne. 0) j16 = mod(nt2,p01)
               if (j16 .gt. 0) then
                  do j4 = 1, p02
                     do j3 = 1, p01 - j16
                        r11((j4-1)*j13+j3+r11_) = p0((j4-1)*p01+j3+j16+
     1                     p0_)
                     end do
                     do j3 = p01 - j16 + 1, p01
                        r11((j4-1)*j13+j3+r11_) = p0((j4-2)*p01+j3+j16+
     1                     p0_)
                     end do
                  end do
               else
                  do j4 = 1, p02
                     do j3 = 1 - j16, p01
                        r11((j4-1)*j13+j3+r11_) = p0((j4-1)*p01+j3+j16+
     1                     p0_)
                     end do
                     do j3 = 1, -j16
                        r11((j4-1)*j13+j3+r11_) = p0(j4*p01+j3+j16+p0_)
                     end do
                  end do
               endif
               do j4 = 1, p02
                  do j3 = 1, p01
                     p0((j4-1)*p01+j3+p0_) = r11((j4-1)*j13+j3+r11_)
                  end do
               end do
               r11_ = r11_*4 + heap__1
               call deallocate_f90 (4, j14*j13, 0, r11_, 0)
               j17 = p01
               j18 = p02
               call allocate_f90 (4, j18*j17, 0, r12_, 0)
               r12_ = (r12_ - heap__1)/4
               if (p02 .ne. 0) j20 = mod(nx2,p02)
               if (j20 .gt. 0) then
                  do j3 = 1, p01
                     do j4 = 1, p02 - j20
                        r12((j4-1)*j17+j3+r12_) = p0((j4+j20-1)*p01+j3+
     1                     p0_)
                     end do
                     do j4 = p02 - j20 + 1, p02
                        r12((j4-1)*j17+j3+r12_) = p0((j4+j20-p02-1)*p01+
     1                     j3+p0_)
                     end do
                  end do
               else
                  do j3 = 1, p01
                     do j4 = 1 - j20, p02
                        r12((j4-1)*j17+j3+r12_) = p0((j4+j20-1)*p01+j3+
     1                     p0_)
                     end do
                     do j4 = 1, -j20
                        r12((j4-1)*j17+j3+r12_) = p0((j4+j20+p02-1)*p01+
     1                     j3+p0_)
                     end do
                  end do
               endif
               do j4 = 1, p02
                  do j3 = 1, p01
                     p0((j4-1)*p01+j3+p0_) = r12((j4-1)*j17+j3+r12_)
                  end do
               end do
               r12_ = r12_*4 + heap__1
               call deallocate_f90 (4, j18*j17, 0, r12_, 0)
               do j4 = 1, 128
                  do j3 = 1, 2048
                     wt2(j3,j4) = 0
                  end do
               end do
               do j4 = 1, nx2
                  do j3 = 1, nt2
                     wt2(j3,j4) = p0((j4-1)*p01+j3+p0_)
                  end do
                  do j3 = 1, nt2
                     wt2(j3+2048-nt2,j4) = p0((j4-1)*p01+j3+nt2+p0_)
                  end do
               end do
               do j4 = 1, nx2
                  do j3 = 1, nt2
                     wt2(j3,j4+128-nx2) = p0((j4+nx2-1)*p01+j3+p0_)
                  end do
                  do j3 = 1, nt2
                     wt2(j3+2048-nt2,j4+128-nx2) = p0((j4+nx2-1)*p01+j3+
     1                  nt2+p0_)
                  end do
               end do
 
               call rlft3 (wt2, wspeq, 2048, 128, 1, 1, 1)
               do j4 = 1, 128
                  do j3 = 1, 1024
                     cbuf(j3,j4) = cdata(j3,j4)*conjg(cwt2(j3,j4))
                  end do
               end do
               do j3 = 1, 128
                  bspeq(j3) = dspeq(j3)*conjg(wspeq(j3))
               end do
               call rlft3 (buf, bspeq, 2048, 128, 1, -1, -1)
               do j4 = 1, 128
                  do j3 = 1, 2048
                     buf(j3,j4) = buf(j3,j4)*fac2
                  end do
               end do
               do j4 = 1, 128
                  do j3 = 1, 2048
                     if (abs(buf(j3,j4)) .lt. cut) buf(j3,j4) = 0
                  end do
               end do
               do j4 = 1, 128
                  do j3 = 1, 2048 - nt
                     buf(j3+nt,j4) = 0
                  end do
               end do
               do j4 = 1, 0
                  do j3 = 1, nt
                     buf(j3,j4+nx) = 0
                  end do
               end do
 
c        do i=1,n1
c           do ii=1,n2
c              if(abs(buf(i,ii)).lt.cut) buf(i,ii)=0
c           enddo
c        enddo
 
               call rlft3 (buf, bspeq, 2048, 128, 1, 1, 1)
 
               do j4 = 1, 128
                  do j3 = 1, 1024
                     cbuf(j3,j4) = cbuf(j3,j4)*cwt2(j3,j4)
                  end do
               end do
               do j3 = 1, 128
                  bspeq(j3) = bspeq(j3)*wspeq(j3)
               end do
               call rlft3 (buf, bspeq, 2048, 128, 1, -1, -1)
               do j4 = 1, 128
                  do j3 = 1, 2048
                     buf(j3,j4) = buf(j3,j4)*fac2
                  end do
               end do
               do j4 = 1, nx
                  do j3 = 1, nt
                     p00((j4-1)*p001+j3+p00_) = p00((j4-1)*p001+j3+p00_)
     1                   + buf(j3,j4)/st2        !sum each scale
                  end do
               end do
c        print*,minval(buf(1:nt,1:nx)),maxval(buf(1:nt,1:nx))
 
c        write(10) (buf(1:nt,1:nx))
c        stop
            end do                               !end of sx
         end do                                  !end of st
         do j4 = 1, p002
            do j3 = 1, p001
               p00((j4-1)*p001+j3+p00_) = p00((j4-1)*p001+j3+p00_)*(nt2*
     1            nx2)**2                        !recover amplitude
            end do
         end do
c        write(10) real(wt2)
c        write(10) p00 !the file processed data
 
c su format output
         rewind (11)
         do i = 1, nx
            read (11) headr
            write (10) headr, (p00((i-1)*p001+j3+p00_)*1E+6,j3=1,p001)
         end do
         close(11)
         close(10)
         p0_ = p0_*4 + heap__1
         call deallocate_f90 (4, p02*p01, 0, p0_, 0)
         p00_ = p00_*4 + heap__1
         call deallocate_f90 (4, p002*p001, 0, p00_, 0)
         tr_ = tr_*4 + heap__1
         call deallocate_f90 (4, tr2, 0, tr_, 0)
         go to 99
   98    continue
         call trim_f90 (ifile, j25)
         print *, ifile(1:j25), 
     1      ' is an empty file. Process the next one.'
         close(11)
   99    continue
      end do
      close(22)                                  !end of loop of file
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
      subroutine hpwavelet(p, sx, h, st, nt, nx, dt, dx)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 11:12:41 2/ 8/02
C...Switches:                     
      integer nt
      integer nx
      real p(nt,nx), sx, st, dx, dt
      double precision tmp, x, t
      integer j1, j2
      do j2 = 1, nx
         do j1 = 1, nt
            p(j1,j2) = 0
         end do
      end do
      do ix = 1, nx
         x = (ix - nx/2)*dx
         do i = 1, nt
            t = (i - 1)*dt
            tmp = ((t - sqrt((x/sx)**2 + h*h))/st)**2
            p(i,ix) = (1 - tmp)*exp((-tmp/2))
         end do
      end do
      end 
 
c hiperbolic-Mexica Hat 2D wavelet with single side
c p=hpwavelet(x,t,sx,h,st)
c   dx--space interval
c   dt--time interval
c   sx-- scale in x direction
c   h--determine the apix location in t direction of this wavelet
c   st--scale in t directon
      subroutine hpwavelet2(p, sx, h, st, nt, nx, dt, dx)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 11:12:41 2/ 8/02
C...Switches:                     
      integer nt
      integer nx
      real p(nt,nx), sx, st, dx, dt
      double precision tmp, x, t
      integer j1, j2
      do j2 = 1, nx
         do j1 = 1, nt
            p(j1,j2) = 0
         end do
      end do
      do ix = nx/2, nx
         x = (ix - nx/2)*dx
         do i = 1, nt
            t = (i - 1)*dt
            tmp = ((t - sqrt((x/sx)**2 + h*h))/st)**2
            p(i,ix) = (1 - tmp)*exp((-tmp/2))
         end do
      end do
      end 
 
ccccccccc debug
c        open(30,file='dtt.dat',form='binary')
c        write(30)abs(cwt2)
c        close(30)
c        stop
ccccccccc debug
      subroutine rlft3(data, speq, nn1, nn2, nn3, isign)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 11:12:41 2/ 8/02
C...Switches:                     
      integer isign, nn1, nn2, nn3
      complex  data(nn1/2,nn2,nn3), speq(nn2,nn3)
CU    USES fourn
      integer i1, i2, i3, j1, j2, j3, nn(3)
      double precision theta, wi, wpi, wpr, wr, wtemp
      complex  c1, c2, h1, h2, w
      c1 = cmplx(0.5,0.0)
      c2 = cmplx(0.0,(-0.5*isign))
      theta = 6.28318530717959D0/dble(isign*nn1)
      wpr = -2.0D0*sin(0.5D0*theta)**2
      wpi = sin(theta)
      nn(1) = nn1/2
      nn(2) = nn2
      nn(3) = nn3
      if (isign .eq. 1) then
         call fourn (data, nn, 3, isign)
         do i3 = 1, nn3
            do i2 = 1, nn2
               speq(i2,i3) = data(1,i2,i3)
            end do
         end do
      endif
      do i3 = 1, nn3
         j3 = 1
         if (i3 .ne. 1) j3 = nn3 - i3 + 2
         wr = 1.0D0
         wi = 0.0D0
         do i1 = 1, nn1/4 + 1
            j1 = nn1/2 - i1 + 2
            do i2 = 1, nn2
               j2 = 1
               if (i2 .ne. 1) j2 = nn2 - i2 + 2
               if (i1 .eq. 1) then
                  h1 = c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
                  h2 = c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
                  data(1,i2,i3) = h1 + h2
                  speq(j2,j3) = conjg(h1 - h2)
               else
                  h1 = c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
                  h2 = c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
                  data(i1,i2,i3) = h1 + w*h2
                  data(j1,j2,j3) = conjg(h1 - w*h2)
               endif
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
            w = cmplx(sngl(wr),sngl(wi))
         end do
      end do
      if (isign .eq. (-1)) call fourn (data, nn, 3, isign)
      return 
      end 
C  (C) Copr. 1986-92 Numerical Recipes Software 5,29#(.
 
      subroutine fourn(data, nn, ndim, isign)
C...Translated by Pacific-Sierra Research vf90 Personal 3.4N5 11:12:41 2/ 8/02
C...Switches:                     
      integer ndim
      integer isign, nn(ndim)
      real data(*)
      integer i1, i2, i2rev, i3, i3rev, ibit, idim, ifp1, ifp2, ip1, ip2
     1   , ip3, k1, k2, n, nprev, nrem, ntot
      real tempi, tempr
      double precision theta, wi, wpi, wpr, wr, wtemp
      ntot = 1
      do idim = 1, ndim
         ntot = ntot*nn(idim)
      end do
      nprev = 1
      do idim = 1, ndim
         n = nn(idim)
         nrem = ntot/(n*nprev)
         ip1 = 2*nprev
         ip2 = ip1*n
         ip3 = ip2*nrem
         i2rev = 1
         do i2 = 1, ip2, ip1
            if (i2 .lt. i2rev) then
               do i1 = i2, i2 + ip1 - 2, 2
                  do i3 = i1, ip3, ip2
                     i3rev = i2rev + i3 - i2
                     tempr = data(i3)
                     tempi = data(i3+1)
                     data(i3) = data(i3rev)
                     data(i3+1) = data(i3rev+1)
                     data(i3rev) = tempr
                     data(i3rev+1) = tempi
                  end do
               end do
            endif
            ibit = ip2/2
    1       continue
            if (ibit.ge.ip1 .and. i2rev.gt.ibit) then
               i2rev = i2rev - ibit
               ibit = ibit/2
               go to 1
            endif
            i2rev = i2rev + ibit
         end do
         ifp1 = ip1
    2    continue
         if (ifp1 .lt. ip2) then
            ifp2 = 2*ifp1
            theta = isign*6.28318530717959D0/(ifp2/ip1)
            wpr = -2.D0*sin(0.5D0*theta)**2
            wpi = sin(theta)
            wr = 1.D0
            wi = 0.D0
            do i3 = 1, ifp1, ip1
               do i1 = i3, i3 + ip1 - 2, 2
                  do i2 = i1, ip3, ifp2
                     k1 = i2
                     k2 = k1 + ifp1
                     tempr = sngl(wr)*data(k2) - sngl(wi)*data(k2+1)
                     tempi = sngl(wr)*data(k2+1) + sngl(wi)*data(k2)
                     data(k2) = data(k1) - tempr
                     data(k2+1) = data(k1+1) - tempi
                     data(k1) = data(k1) + tempr
                     data(k1+1) = data(k1+1) + tempi
                  end do
               end do
               wtemp = wr
               wr = wr*wpr - wi*wpi + wr
               wi = wi*wpr + wtemp*wpi + wi
            end do
            ifp1 = ifp2
            go to 2
         endif
         nprev = n*nprev
      end do
      return 
      end 
C  (C) Copr. 1986-92 Numerical Recipes Software 5,29#(.
