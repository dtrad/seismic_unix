c-----------------------------------------------------------------------------


      subroutine FX_symmetry(sfx,nx,nt,nfft,iflow,ifhigh)


c       Imposes symmetry to the FX domain
c
c       Input parameters:
c
c         sfx(nfft/2+1,nx) - signal in the FX. Only positive freqs.
c
c
c       Output parameters:

c         sfx(nfft,nx)   - signal in the FX with conjugate symmetric
c                          freq. axis.



        parameter        (nnt=2048, nnx=228)
        
	  	
	  complex   * 16   sfx(nnt,nnx)

c-------makes the f-h or f-q symmetric

        do 200 ix=1,nx
        do 180 if=1,iflow
180     sfx(if,ix)=dcmplx(0.d0, 0.d0)
        do 190 if=ifhigh+1,nfft/2+1
190     sfx(if,ix)=dcmplx(0.d0, 0.d0)
200     continue

        do 210 ix=1,nx
        do 210 if=nfft/2+2,nfft
210     sfx(if,ix) = dconjg(sfx(nfft-if+2,ix))

        

        return
        end

c-----------------------------------------------------------------------------
