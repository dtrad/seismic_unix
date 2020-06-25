
        subroutine sFX_go (nx,nt,nfft,stx,sfx,index)

c------ index=-1 TX ---> FX
c------ index= 1 FX----> TX

        parameter  (nnt=2048,nnx=228)
        complex         sfx(nnt,nnx), aux(nnt)
        real            stx(nnt,nnx)
        if(index.eq.-1)  then
        do 140 ix=1,nx
        do 130 it=1,nt
130     aux(it)=cmplx(stx(it,ix),0.0)
        do 120 it=nt+1,nfft
120     aux(it)=cmplx(0.0, 0.0)
        call fork(nfft,aux,-1.)
        do 110 if=1,nfft
110     sfx(if,ix)=aux(if)
140     continue
        endif
        if(index.eq.1) then
        do 240 ix=1,nx
        do 220 if=1,nfft
 
220     aux(if)=sfx(if,ix)
        call fork(nfft,aux, 1.)
        do 230 it=1,nfft
230     stx(it,ix)=real(aux(it))
240     continue
        endif
        return
        end
 
