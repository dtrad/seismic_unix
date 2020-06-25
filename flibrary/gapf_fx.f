      subroutine gapf_fx(data,nt,nx,dt,dx,n1,n2,ip)
c     call to
c			SUBROUTINE BURG (N,IP,X,P,A,ISTAT) 


      parameter  (nnt=2048,nnx=228)
      real       data(nnt,nnx),mu
      complex    data_fx(nnt,nnx)
      complex    x(nnx),a(nnx)
      complex    x2(nnx),a2(nnx)


      nf=2048
      if(nf.lt.nt) pause 'Check NF'

 
      call sFX_go (nx,nt,nf,data,data_fx,-1)

c Do it for all freqs in the band

      iflow = 5
      ifhigh = nf/2+1-5

      do if=iflow,ifhigh
        do ix=1,nx
           x(ix)=data_fx(if,ix)
        enddo

        call BURG (N1-1,IP,X,P,A,ISTAT)
        do ix=1,nx-n2
           x2(ix)=data_fx(if,ix+n2)
        enddo
        call BURG (Nx-n2,IP,X2,P,A2,ISTAT)
        do i=1,ip
          a(i)=0.5*(a(i)+a2(i))
        enddo
        call gapf(x,nx,n1,n2,ip,a)   
        do ix=1,nx
          data_fx(if,ix)=x(ix)
        enddo
      enddo

      call  sFX_symmetry (data_fx,nx,nt,nf,iflow,ifhigh)
      call  sFX_go (nx,nt,nf,data,data_fx,1)

      return
      end
