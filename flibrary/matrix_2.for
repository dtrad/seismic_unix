c-----------------------------------------------------------------------------

      subroutine matrix_2(L,lh,pmin,pmax,pos,np,nh,w,dh)

c       Transformation matrix.
c       This matrix relates the cmp gather and the velocity
c       gather in the f-x space.

c       Input parameters:
c
c         np   - number of parameters= number of traces of the velocity gather

c         nh   - number of traces of the CMP
c         pmin - minimum parameter seek by the transform
c         pmax - maximum parameter seek by the transform
c            w - the normalized freq. at which the transform is evaluated
c   near_offet - near offset trace in meters.
c

c       Out parameter:
c
c       Notes:
c
c       The parameter p in the velocity gather is the slowness
c	  L=F.WU has size np x nh such that m=L.u
c	  LH=FH.WV has size nh x np such that u=LH.m
c
c       Note: I have eliminated dh from L because I use a formulation that allows cancelation 
c       of the dh terms. This allows us to use Cholesky instead of LU. 
c       If Levinson is used (gauss_gauss_To.f) Wu is multiplied by F to get L.
c		Daniel Trad- 22-02-99
 

        parameter (nnx=228)

        complex * 16 l(nnx,nnx), lh(nnx,nnx)
        real * 8     dp,dh(nnx),arg,w,pmin,pmax,p,h
        real * 8     pos(nnx)

        dp=(pmax-pmin)/dfloat(np-1)
        do 100 j=1,np
		p=pmin+(pmax-pmin)*dfloat(j-1)/dfloat(np-1)
		do 200 i=1,nh
			h=pos(i)
			arg=w*h*h*p
			l(j,i)=cdexp(dcmplx(0.d0,arg))
			lh(i,j)=dp*cdexp(dcmplx(0.d0,-arg))
200		continue
100     continue 	  	
        return
        end

c---------------------------------------------------
