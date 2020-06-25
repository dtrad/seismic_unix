#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "cwp.h"
#include "clibrary.h"
void  circ_mult(int n,complex *RC,complex *f,complex *g,int nf,
		complex *fc, complex *gc)
{

  /*
     Perform the multiplication of a Toeplitz form with a vector.
     The procedure is done using a circulant matrix derived from
     the first row of the Toeplitz form. The Toeplitz form, in fact,
     is not needed. Everything you need is the firt row of the matrix.

     INPUT  r(n): first row of the autocorrelation matrix (Toeplitz)
            f(n): a filter
              nf: number of freq. samples to perform the circular conv.
            fc(n): a pointer to an already allocated vector
            gc(n): a pointer to an already allocated vector

     OUTPUT g(n): the product of the toeplitz form with f

     NOTE   that T(1) must be real in order to  have an Hermitian 
           form T

     Mauricio Sachi 

     adapted for SU and translated to C: Daniel Trad
  */
      
  complex czero,aux;
  int i;
  czero.r=czero.i=0.;

  for (i=0;i<n;i++) fc[i]=f[i];  
  for (i=n;i<nf;i++) fc[i]=czero;    // pad with zeros the autocorrelation
 
  pfacc(-1,nf,fc);
  
  for (i=0;i<nf;i++) gc[i]=RC[i]*fc[i];
  
  pfacc (1,nf,gc);
  
  for (i=0;i<n;i++) g[i]=gc[i];

  for (i=0;i<n;i++) g[i]/=nf;            

  return;
}










