I will try but if a problem is not simple for you perhaps it
does not have solution at all!

Let me see. I am confuse with the vector v_i, but since there is not any
other indication it seems to be a vector with constant coefficients that
somehow define the surface of integration.

Now, it seems that you need to compute the operator as usual
(hopefully is sparse) an 

for im
  for iv_i  ------> here I have a ??
        itheta=v_i^T m 

        L[0][im]=itheta
        L[1][im]=im
  end
end


And then the forward should be

for im=0:nsparse_elem
        b[L[0][im]]+=f[L[1][im]


The adjoint would be
   
for im=0:nsparse_elem
        f[L[1][im]]+=b[L[0][im]
