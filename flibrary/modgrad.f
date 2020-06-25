      subroutine modgrad(u,eps2,np,sigma,norm,eps,Qp,nnp)
c     It computes the term Qp for the total gradient g=(L'L+Qp)u-L'd
c     This term corresponds to the probability model, so that
c     the distribution parameters sigma and norm are passed.
c     Qp is a vector of dimension np, but in fact is the diagonal 
c     of the Qp matrix of size np x np.
c     Input 
c          u: model
c          eps2: probability constant
c          np: number of model traces
c          sigma: standard deviation of the model
c          norm: implemented 1 Huber, 10 Cauchy
c          eps: minimum small number for avoiding division by zero
c     Ouput 
c          Qp: Vector such Matrix(diag(Qp)) defines the
c          model norm as Jmod= u'Qpu
c     Daniel Trad- 21 March 99. UBC- Canada
c     Based in Sacchi, 1996. phD thesis. UBC. Canada

      implicit none
      integer norm,i,np, nnp
      complex*16 u(nnp), Qp(nnp)
      real*8 eps2,sigma,power,power2(nnp),pmax,pmin,eps

      power=0.d0
      do i=1,np
           power2(i)=dreal(dconjg(u(i))*u(i))
           power=power+power2(i)
      enddo
      pmax=power2(1)
      do i=1,np
           if (power2(i).gt.pmax) pmax=power2(i)
      enddo
      pmin=power2(1)
      do i=1,np
           if (power2(i).lt.pmin) pmin=power2(i)
      enddo
c      if (pmin.gt.1.d-6) eps2=pmin
c      if (pmin.lt.1.d-6) eps2=1d-6

        power=power/np
        eps2=power
c        write(0,*) 'eps2',eps2
      do i=1,np
          if (norm.eq.10) then
             Qp(i)=dcmplx(1.d0/(eps2+power2(i)),0.d0)
          elseif (norm.eq.1) then   
             Qp(i)=dcmplx(1.d0/eps2/dsqrt(power2(i)),0.d0)
          endif   
      enddo
      return
      end

      
