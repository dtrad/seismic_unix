c-----------------------------------------------------------------------------

        subroutine FX_go (nx,nt,nfft,stx,sfx,index)


c       Transform the data to the FX domain
c       and viceversa depending on index
c
c       Input parameters:

c         stx(nt,nx)      -  if index.eq.-1  time-offset  data
c         stf(nfft,nx)    -  if index.eq. 1  freq-offset data
c         nfft            -  lenght of the transform

c       Output parameters:

c         stf(nfft,nx)    - if index.eq.-1  freq-offset
c         stx(nt,nx)      - if index.eq. 1  time-offset
c         nfft            - lenght of the transform
c
c       Notes:
c
c        The input/output changes according to index:
c
c          index = -1 TX ----> FX
c          index =  1 FX ----> TX



        parameter   (nnt=2048, nnx=228)
        complex * 16    sfx(nnt,nnx), aux(nnt)
        real    * 8     stx(nnt,nnx)

        if(index.eq.-1)  then

        do 140 ix=1,nx
        do 130 it=1,nt
130     aux(it)=dcmplx(stx(it,ix),0.d0)
        do 120 it=nt+1,nfft
120     aux(it)=dcmplx(0.d0, 0.d0)

        call fft(nfft,aux,-1)

        do 110 if=1,nfft
110     sfx(if,ix)=aux(if)
140     continue
         
        endif

        if(index.eq.1) then
        do 240 ix=1,nx
        do 220 if=1,nfft
220     aux(if)=sfx(if,ix)

        call fft(nfft,aux, 1)

        do 230 it=1,nfft
230     stx(it,ix)=dreal(aux(it))
240     continue
        
        endif

        return
        end


c-----------------------------------------------------------------------------
