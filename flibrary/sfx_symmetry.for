
        subroutine sFX_symmetry (sfx,nx,nt,nfft,iflow,ifhigh)

        parameter       (nnt=2048,nnx=228)
        complex          sfx(nnt,nnx)
        do 200 ix=1,nx
        do 180 if=1,iflow
180     sfx(if,ix)=cmplx(0.0, 0.0)
        do 190 if=ifhigh+1,nfft/2+1
190     sfx(if,ix)=cmplx(0.0, 0.0)
200     continue
        do 210 ix=1,nx
        do 210 if=nfft/2+2,nfft
210     sfx(if,ix) = conjg(sfx(nfft-if+2,ix))
        return
        end
 

