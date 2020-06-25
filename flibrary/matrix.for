
        subroutine matrix(L,pmin,pmax,np,nh,dh,w,near_offset)

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
c         near_offet - near offset trace in meters.

c       Output units:
c
c         L(nh,np): matrix
c
c       Notes:
c
c       The parameter p in the velocity gather is the square solowness
c       p= 1/(velocity)^2


        parameter (nnx=128)

        complex * 16 l(nnx,nnx) 
        real * 8     dp,arg,dh,w,pmin,pmax,p,h,near_offset,sn

        dp=(pmax-pmin)/dfloat(np-1)
        sn=dsqrt(dp*dh)

        do 100 j=1,np
        p=pmin+(pmax-pmin)*dfloat(j-1)/dfloat(np-1)
        do 100 i=1,nh
        h=dfloat(i-1)*dh+near_offset 
        arg=w*h*h*p
100     l(i,j)=cdexp(dcmplx(0.d0,-arg))*sn

        return
        end
