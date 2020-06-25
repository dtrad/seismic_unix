
#include "f2c.h"

/* Subroutine */ int cgebrd_(integer *m, integer *n, complex *a, integer *lda,
	 real *d, real *e, complex *tauq, complex *taup, complex *work, 
	integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGEBRD reduces a general complex M-by-N matrix A to upper or lower   
    bidiagonal form B by a unitary transformation: Q**H * A * P = B.   

    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows in the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns in the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the M-by-N general matrix to be reduced.   
            On exit,   
            if m >= n, the diagonal and the first superdiagonal are   
              overwritten with the upper bidiagonal matrix B; the   
              elements below the diagonal, with the array TAUQ, represent 
  
              the unitary matrix Q as a product of elementary   
              reflectors, and the elements above the first superdiagonal, 
  
              with the array TAUP, represent the unitary matrix P as   
              a product of elementary reflectors;   
            if m < n, the diagonal and the first subdiagonal are   
              overwritten with the lower bidiagonal matrix B; the   
              elements below the first subdiagonal, with the array TAUQ, 
  
              represent the unitary matrix Q as a product of   
              elementary reflectors, and the elements above the diagonal, 
  
              with the array TAUP, represent the unitary matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    D       (output) REAL array, dimension (min(M,N))   
            The diagonal elements of the bidiagonal matrix B:   
            D(i) = A(i,i).   

    E       (output) REAL array, dimension (min(M,N)-1)   
            The off-diagonal elements of the bidiagonal matrix B:   
            if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;   
            if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.   

    TAUQ    (output) COMPLEX array dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix Q. See Further Details.   

    TAUP    (output) COMPLEX array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix P. See Further Details.   

    WORK    (workspace/output) COMPLEX array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= max(1,M,N).   
            For optimum performance LWORK >= (M+N)*NB, where NB   
            is the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

    If m >= n,   

       Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, and v and u are complex   
    vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in   
    A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in 
  
    A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n,   

       Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, and v and u are complex   
    vectors; v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in   
    A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in 
  
    A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The contents of A on exit are illustrated by the following examples: 
  

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )   
      (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )   
      (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )   
      (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )   
      (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )   
      (  v1  v2  v3  v4  v5 )   

    where d and e denote diagonal and off-diagonal elements of B, vi   
    denotes an element of the vector defining H(i), and ui an element of 
  
    the vector defining G(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {1.f,0.f};
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1;
    /* Local variables */
    static integer i, j;
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, integer *);
    static integer nbmin, iinfo, minmn;
    extern /* Subroutine */ int cgebd2_(integer *, integer *, complex *, 
	    integer *, real *, real *, complex *, complex *, complex *, 
	    integer *);
    static integer nb;
    extern /* Subroutine */ int clabrd_(integer *, integer *, integer *, 
	    complex *, integer *, real *, real *, complex *, complex *, 
	    complex *, integer *, complex *, integer *);
    static integer nx;
    static real ws;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwrkx, ldwrky;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAUQ(I) tauq[(I)-1]
#define TAUP(I) taup[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = max(1,*m);
	if (*lwork < max(i__1,*n)) {
	    *info = -10;
	}
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("CGEBRD", &i__1);
	return 0;
    }

/*     Quick return if possible */

    minmn = min(*m,*n);
    if (minmn == 0) {
	WORK(1).r = 1.f, WORK(1).i = 0.f;
	return 0;
    }

    ws = (real) max(*m,*n);
    ldwrkx = *m;
    ldwrky = *n;

/*     Set the block size NB and the crossover point NX.   

   Computing MAX */
    i__1 = 1, i__2 = ilaenv_(&c__1, "CGEBRD", " ", m, n, &c_n1, &c_n1, 6L, 1L)
	    ;
    nb = max(i__1,i__2);

    if (nb > 1 && nb < minmn) {

/*        Determine when to switch from blocked to unblocked code.   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "CGEBRD", " ", m, n, &c_n1, &c_n1, 
		6L, 1L);
	nx = max(i__1,i__2);
	if (nx < minmn) {
	    ws = (real) ((*m + *n) * nb);
	    if ((real) (*lwork) < ws) {

/*              Not enough work space for the optimal NB, cons
ider using   
                a smaller block size. */

		nbmin = ilaenv_(&c__2, "CGEBRD", " ", m, n, &c_n1, &c_n1, 6L, 
			1L);
		if (*lwork >= (*m + *n) * nbmin) {
		    nb = *lwork / (*m + *n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }

    i__1 = minmn - nx;
    i__2 = nb;
    for (i = 1; nb < 0 ? i >= minmn-nx : i <= minmn-nx; i += nb) {

/*        Reduce rows and columns i:i+ib-1 to bidiagonal form and retu
rn   
          the matrices X and Y which are needed to update the unreduce
d   
          part of the matrix */

	i__3 = *m - i + 1;
	i__4 = *n - i + 1;
	clabrd_(&i__3, &i__4, &nb, &A(i,i), lda, &D(i), &E(i), &
		TAUQ(i), &TAUP(i), &WORK(1), &ldwrkx, &WORK(ldwrkx * nb + 1), 
		&ldwrky);

/*        Update the trailing submatrix A(i+ib:m,i+ib:n), using   
          an update of the form  A := A - V*Y' - X*U' */

	i__3 = *m - i - nb + 1;
	i__4 = *n - i - nb + 1;
	q__1.r = -1.f, q__1.i = 0.f;
	cgemm_("No transpose", "Conjugate transpose", &i__3, &i__4, &nb, &
		q__1, &A(i+nb,i), lda, &WORK(ldwrkx * nb + nb + 
		1), &ldwrky, &c_b1, &A(i+nb,i+nb), lda);
	i__3 = *m - i - nb + 1;
	i__4 = *n - i - nb + 1;
	q__1.r = -1.f, q__1.i = 0.f;
	cgemm_("No transpose", "No transpose", &i__3, &i__4, &nb, &q__1, &
		WORK(nb + 1), &ldwrkx, &A(i,i+nb), lda, &c_b1, 
		&A(i+nb,i+nb), lda);

/*        Copy diagonal and off-diagonal elements of B back into A */

	if (*m >= *n) {
	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		i__4 = j + j * a_dim1;
		i__5 = j;
		A(j,j).r = D(j), A(j,j).i = 0.f;
		i__4 = j + (j + 1) * a_dim1;
		i__5 = j;
		A(j,j+1).r = E(j), A(j,j+1).i = 0.f;
/* L10: */
	    }
	} else {
	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		i__4 = j + j * a_dim1;
		i__5 = j;
		A(j,j).r = D(j), A(j,j).i = 0.f;
		i__4 = j + 1 + j * a_dim1;
		i__5 = j;
		A(j+1,j).r = E(j), A(j+1,j).i = 0.f;
/* L20: */
	    }
	}
/* L30: */
    }

/*     Use unblocked code to reduce the remainder of the matrix */

    i__2 = *m - i + 1;
    i__1 = *n - i + 1;
    cgebd2_(&i__2, &i__1, &A(i,i), lda, &D(i), &E(i), &TAUQ(i), &
	    TAUP(i), &WORK(1), &iinfo);
    WORK(1).r = ws, WORK(1).i = 0.f;
    return 0;

/*     End of CGEBRD */

} /* cgebrd_ */

#include "f2c.h"

/* Subroutine */ int cgebd2_(integer *m, integer *n, complex *a, integer *lda,
	 real *d, real *e, complex *tauq, complex *taup, complex *work, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CGEBD2 reduces a complex general m by n matrix A to upper or lower   
    real bidiagonal form B by a unitary transformation: Q' * A * P = B.   

    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows in the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns in the matrix A.  N >= 0.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the m by n general matrix to be reduced.   
            On exit,   
            if m >= n, the diagonal and the first superdiagonal are   
              overwritten with the upper bidiagonal matrix B; the   
              elements below the diagonal, with the array TAUQ, represent 
  
              the unitary matrix Q as a product of elementary   
              reflectors, and the elements above the first superdiagonal, 
  
              with the array TAUP, represent the unitary matrix P as   
              a product of elementary reflectors;   
            if m < n, the diagonal and the first subdiagonal are   
              overwritten with the lower bidiagonal matrix B; the   
              elements below the first subdiagonal, with the array TAUQ, 
  
              represent the unitary matrix Q as a product of   
              elementary reflectors, and the elements above the diagonal, 
  
              with the array TAUP, represent the unitary matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    D       (output) REAL array, dimension (min(M,N))   
            The diagonal elements of the bidiagonal matrix B:   
            D(i) = A(i,i).   

    E       (output) REAL array, dimension (min(M,N)-1)   
            The off-diagonal elements of the bidiagonal matrix B:   
            if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;   
            if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.   

    TAUQ    (output) COMPLEX array dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix Q. See Further Details.   

    TAUP    (output) COMPLEX array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix P. See Further Details.   

    WORK    (workspace) COMPLEX array, dimension (max(M,N))   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

    If m >= n,   

       Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, and v and u are complex   
    vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in   
    A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in 
  
    A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n,   

       Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, v and u are complex vectors; 
  
    v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i); 
  
    u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n); 
  
    tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The contents of A on exit are illustrated by the following examples: 
  

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )   
      (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )   
      (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )   
      (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )   
      (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )   
      (  v1  v2  v3  v4  v5 )   

    where d and e denote diagonal and off-diagonal elements of B, vi   
    denotes an element of the vector defining H(i), and ui an element of 
  
    the vector defining G(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static integer i;
    static complex alpha;
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, complex *
	    , integer *, complex *, complex *, integer *, complex *), 
	    clarfg_(integer *, complex *, complex *, integer *, complex *), 
	    clacgv_(integer *, complex *, integer *), xerbla_(char *, integer 
	    *);



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAUQ(I) tauq[(I)-1]
#define TAUP(I) taup[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("CGEBD2", &i__1);
	return 0;
    }

    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {

/*           Generate elementary reflector H(i) to annihilate A(i+
1:m,i) */

	    i__2 = i + i * a_dim1;
	    alpha.r = A(i,i).r, alpha.i = A(i,i).i;
	    i__2 = *m - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    clarfg_(&i__2, &alpha, &A(min(i+1,*m),i), &c__1, &
		    TAUQ(i));
	    i__2 = i;
	    D(i) = alpha.r;
	    i__2 = i + i * a_dim1;
	    A(i,i).r = 1.f, A(i,i).i = 0.f;

/*           Apply H(i)' to A(i:m,i+1:n) from the left */

	    i__2 = *m - i + 1;
	    i__3 = *n - i;
	    r_cnjg(&q__1, &TAUQ(i));
	    clarf_("Left", &i__2, &i__3, &A(i,i), &c__1, &q__1, &A(i,i+1), lda, &WORK(1));
	    i__2 = i + i * a_dim1;
	    i__3 = i;
	    A(i,i).r = D(i), A(i,i).i = 0.f;

	    if (i < *n) {

/*              Generate elementary reflector G(i) to annihila
te   
                A(i,i+2:n) */

		i__2 = *n - i;
		clacgv_(&i__2, &A(i,i+1), lda);
		i__2 = i + (i + 1) * a_dim1;
		alpha.r = A(i,i+1).r, alpha.i = A(i,i+1).i;
		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;
		clarfg_(&i__2, &alpha, &A(i,min(i+2,*n)), lda, &
			TAUP(i));
		i__2 = i;
		E(i) = alpha.r;
		i__2 = i + (i + 1) * a_dim1;
		A(i,i+1).r = 1.f, A(i,i+1).i = 0.f;

/*              Apply G(i) to A(i+1:m,i+1:n) from the right */

		i__2 = *m - i;
		i__3 = *n - i;
		clarf_("Right", &i__2, &i__3, &A(i,i+1), lda, &
			TAUP(i), &A(i+1,i+1), lda, &WORK(1));
		i__2 = *n - i;
		clacgv_(&i__2, &A(i,i+1), lda);
		i__2 = i + (i + 1) * a_dim1;
		i__3 = i;
		A(i,i+1).r = E(i), A(i,i+1).i = 0.f;
	    } else {
		i__2 = i;
		TAUP(i).r = 0.f, TAUP(i).i = 0.f;
	    }
/* L10: */
	}
    } else {

/*        Reduce to lower bidiagonal form */

	i__1 = *m;
	for (i = 1; i <= *m; ++i) {

/*           Generate elementary reflector G(i) to annihilate A(i,
i+1:n) */

	    i__2 = *n - i + 1;
	    clacgv_(&i__2, &A(i,i), lda);
	    i__2 = i + i * a_dim1;
	    alpha.r = A(i,i).r, alpha.i = A(i,i).i;
	    i__2 = *n - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    clarfg_(&i__2, &alpha, &A(i,min(i+1,*n)), lda, &TAUP(
		    i));
	    i__2 = i;
	    D(i) = alpha.r;
	    i__2 = i + i * a_dim1;
	    A(i,i).r = 1.f, A(i,i).i = 0.f;

/*           Apply G(i) to A(i+1:m,i:n) from the right */

	    i__2 = *m - i;
	    i__3 = *n - i + 1;
/* Computing MIN */
	    i__4 = i + 1;
	    clarf_("Right", &i__2, &i__3, &A(i,i), lda, &TAUP(i), &
		    A(min(i+1,*m),i), lda, &WORK(1));
	    i__2 = *n - i + 1;
	    clacgv_(&i__2, &A(i,i), lda);
	    i__2 = i + i * a_dim1;
	    i__3 = i;
	    A(i,i).r = D(i), A(i,i).i = 0.f;

	    if (i < *m) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(i+2:m,i) */

		i__2 = i + 1 + i * a_dim1;
		alpha.r = A(i+1,i).r, alpha.i = A(i+1,i).i;
		i__2 = *m - i;
/* Computing MIN */
		i__3 = i + 2;
		clarfg_(&i__2, &alpha, &A(min(i+2,*m),i), &c__1, &
			TAUQ(i));
		i__2 = i;
		E(i) = alpha.r;
		i__2 = i + 1 + i * a_dim1;
		A(i+1,i).r = 1.f, A(i+1,i).i = 0.f;

/*              Apply H(i)' to A(i+1:m,i+1:n) from the left */

		i__2 = *m - i;
		i__3 = *n - i;
		r_cnjg(&q__1, &TAUQ(i));
		clarf_("Left", &i__2, &i__3, &A(i+1,i), &c__1, &
			q__1, &A(i+1,i+1), lda, &WORK(1))
			;
		i__2 = i + 1 + i * a_dim1;
		i__3 = i;
		A(i+1,i).r = E(i), A(i+1,i).i = 0.f;
	    } else {
		i__2 = i;
		TAUQ(i).r = 0.f, TAUQ(i).i = 0.f;
	    }
/* L20: */
	}
    }
    return 0;

/*     End of CGEBD2 */

} /* cgebd2_ */

#include "f2c.h"

/* Subroutine */ int clabrd_(integer *m, integer *n, integer *nb, complex *a, 
	integer *lda, real *d, real *e, complex *tauq, complex *taup, complex 
	*x, integer *ldx, complex *y, integer *ldy)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLABRD reduces the first NB rows and columns of a complex general   
    m by n matrix A to upper or lower real bidiagonal form by a unitary   
    transformation Q' * A * P, and returns the matrices X and Y which   
    are needed to apply the transformation to the unreduced part of A.   

    If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower 
  
    bidiagonal form.   

    This is an auxiliary routine called by CGEBRD   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows in the matrix A.   

    N       (input) INTEGER   
            The number of columns in the matrix A.   

    NB      (input) INTEGER   
            The number of leading rows and columns of A to be reduced.   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            On entry, the m by n general matrix to be reduced.   
            On exit, the first NB rows and columns of the matrix are   
            overwritten; the rest of the array is unchanged.   
            If m >= n, elements on and below the diagonal in the first NB 
  
              columns, with the array TAUQ, represent the unitary   
              matrix Q as a product of elementary reflectors; and   
              elements above the diagonal in the first NB rows, with the 
  
              array TAUP, represent the unitary matrix P as a product   
              of elementary reflectors.   
            If m < n, elements below the diagonal in the first NB   
              columns, with the array TAUQ, represent the unitary   
              matrix Q as a product of elementary reflectors, and   
              elements on and above the diagonal in the first NB rows,   
              with the array TAUP, represent the unitary matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    D       (output) REAL array, dimension (NB)   
            The diagonal elements of the first NB rows and columns of   
            the reduced matrix.  D(i) = A(i,i).   

    E       (output) REAL array, dimension (NB)   
            The off-diagonal elements of the first NB rows and columns of 
  
            the reduced matrix.   

    TAUQ    (output) COMPLEX array dimension (NB)   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix Q. See Further Details.   

    TAUP    (output) COMPLEX array, dimension (NB)   
            The scalar factors of the elementary reflectors which   
            represent the unitary matrix P. See Further Details.   

    X       (output) COMPLEX array, dimension (LDX,NB)   
            The m-by-nb matrix X required to update the unreduced part   
            of A.   

    LDX     (input) INTEGER   
            The leading dimension of the array X. LDX >= max(1,M).   

    Y       (output) COMPLEX array, dimension (LDY,NB)   
            The n-by-nb matrix Y required to update the unreduced part   
            of A.   

    LDY     (output) INTEGER   
            The leading dimension of the array Y. LDY >= max(1,N).   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

       Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are complex scalars, and v and u are complex   
    vectors.   

    If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in   
    A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in   
    A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in   
    A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in   
    A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The elements of the vectors v and u together form the m-by-nb matrix 
  
    V and the nb-by-n matrix U' which are needed, with X and Y, to apply 
  
    the transformation to the unreduced part of the matrix, using a block 
  
    update of the form:  A := A - V*Y' - X*U'.   

    The contents of A on exit are illustrated by the following examples   
    with nb = 2:   

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )   
      (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )   
      (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )   
      (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )   
      (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )   
      (  v1  v2  a   a   a  )   

    where a denotes an element of the original matrix which is unchanged, 
  
    vi denotes an element of the vector defining H(i), and ui an element 
  
    of the vector defining G(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {0.f,0.f};
    static complex c_b2 = {1.f,0.f};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    complex q__1;
    /* Local variables */
    static integer i;
    static complex alpha;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), cgemv_(char *, integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, complex *, 
	    integer *), clarfg_(integer *, complex *, complex *, 
	    integer *, complex *), clacgv_(integer *, complex *, integer *);



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAUQ(I) tauq[(I)-1]
#define TAUP(I) taup[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]
#define Y(I,J) y[(I)-1 + ((J)-1)* ( *ldy)]

    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

	i__1 = *nb;
	for (i = 1; i <= *nb; ++i) {

/*           Update A(i:m,i) */

	    i__2 = i - 1;
	    clacgv_(&i__2, &Y(i,1), ldy);
	    i__2 = *m - i + 1;
	    i__3 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("No transpose", &i__2, &i__3, &q__1, &A(i,1), lda, &
		    Y(i,1), ldy, &c_b2, &A(i,i), &c__1)
		    ;
	    i__2 = i - 1;
	    clacgv_(&i__2, &Y(i,1), ldy);
	    i__2 = *m - i + 1;
	    i__3 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("No transpose", &i__2, &i__3, &q__1, &X(i,1), ldx, &
		    A(1,i), &c__1, &c_b2, &A(i,i), &
		    c__1);

/*           Generate reflection Q(i) to annihilate A(i+1:m,i) */

	    i__2 = i + i * a_dim1;
	    alpha.r = A(i,i).r, alpha.i = A(i,i).i;
	    i__2 = *m - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    clarfg_(&i__2, &alpha, &A(min(i+1,*m),i), &c__1, &
		    TAUQ(i));
	    i__2 = i;
	    D(i) = alpha.r;
	    if (i < *n) {
		i__2 = i + i * a_dim1;
		A(i,i).r = 1.f, A(i,i).i = 0.f;

/*              Compute Y(i+1:n,i) */

		i__2 = *m - i + 1;
		i__3 = *n - i;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(i,i+1), lda, &A(i,i), &c__1, &c_b1,
			 &Y(i+1,i), &c__1);
		i__2 = *m - i + 1;
		i__3 = i - 1;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(i,1), lda, &A(i,i), &c__1, &c_b1, &Y(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &Y(i+1,1)
			, ldy, &Y(1,i), &c__1, &c_b2, &Y(i+1,i), &c__1);
		i__2 = *m - i + 1;
		i__3 = i - 1;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &X(i,1), ldx, &A(i,i), &c__1, &c_b1, &Y(1,i), &c__1);
		i__2 = i - 1;
		i__3 = *n - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("Conjugate transpose", &i__2, &i__3, &q__1, &A(1,i+1), lda, &Y(1,i), &c__1, &c_b2, 
			&Y(i+1,i), &c__1);
		i__2 = *n - i;
		cscal_(&i__2, &TAUQ(i), &Y(i+1,i), &c__1);

/*              Update A(i,i+1:n) */

		i__2 = *n - i;
		clacgv_(&i__2, &A(i,i+1), lda);
		clacgv_(&i, &A(i,1), lda);
		i__2 = *n - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i, &q__1, &Y(i+1,1), 
			ldy, &A(i,1), lda, &c_b2, &A(i,i+1), lda);
		clacgv_(&i, &A(i,1), lda);
		i__2 = i - 1;
		clacgv_(&i__2, &X(i,1), ldx);
		i__2 = i - 1;
		i__3 = *n - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("Conjugate transpose", &i__2, &i__3, &q__1, &A(1,i+1), lda, &X(i,1), ldx, &c_b2, &A(i,i+1), lda);
		i__2 = i - 1;
		clacgv_(&i__2, &X(i,1), ldx);

/*              Generate reflection P(i) to annihilate A(i,i+2
:n) */

		i__2 = i + (i + 1) * a_dim1;
		alpha.r = A(i,i+1).r, alpha.i = A(i,i+1).i;
		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;
		clarfg_(&i__2, &alpha, &A(i,min(i+2,*n)), lda, &
			TAUP(i));
		i__2 = i;
		E(i) = alpha.r;
		i__2 = i + (i + 1) * a_dim1;
		A(i,i+1).r = 1.f, A(i,i+1).i = 0.f;

/*              Compute X(i+1:m,i) */

		i__2 = *m - i;
		i__3 = *n - i;
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &A(i+1,i+1), lda, &A(i,i+1), lda, &c_b1,
			 &X(i+1,i), &c__1);
		i__2 = *n - i;
		cgemv_("Conjugate transpose", &i__2, &i, &c_b2, &Y(i+1,1), ldy, &A(i,i+1), lda, &c_b1, &
			X(1,i), &c__1);
		i__2 = *m - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i, &q__1, &A(i+1,1), 
			lda, &X(1,i), &c__1, &c_b2, &X(i+1,i), &c__1);
		i__2 = i - 1;
		i__3 = *n - i;
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &A(1,i+1), lda, &A(i,i+1), lda, &
			c_b1, &X(1,i), &c__1);
		i__2 = *m - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &X(i+1,1)
			, ldx, &X(1,i), &c__1, &c_b2, &X(i+1,i), &c__1);
		i__2 = *m - i;
		cscal_(&i__2, &TAUP(i), &X(i+1,i), &c__1);
		i__2 = *n - i;
		clacgv_(&i__2, &A(i,i+1), lda);
	    }
/* L10: */
	}
    } else {

/*        Reduce to lower bidiagonal form */

	i__1 = *nb;
	for (i = 1; i <= *nb; ++i) {

/*           Update A(i,i:n) */

	    i__2 = *n - i + 1;
	    clacgv_(&i__2, &A(i,i), lda);
	    i__2 = i - 1;
	    clacgv_(&i__2, &A(i,1), lda);
	    i__2 = *n - i + 1;
	    i__3 = i - 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("No transpose", &i__2, &i__3, &q__1, &Y(i,1), ldy, &
		    A(i,1), lda, &c_b2, &A(i,i), lda);
	    i__2 = i - 1;
	    clacgv_(&i__2, &A(i,1), lda);
	    i__2 = i - 1;
	    clacgv_(&i__2, &X(i,1), ldx);
	    i__2 = i - 1;
	    i__3 = *n - i + 1;
	    q__1.r = -1.f, q__1.i = 0.f;
	    cgemv_("Conjugate transpose", &i__2, &i__3, &q__1, &A(1,i), lda, &X(i,1), ldx, &c_b2, &A(i,i),
		     lda);
	    i__2 = i - 1;
	    clacgv_(&i__2, &X(i,1), ldx);

/*           Generate reflection P(i) to annihilate A(i,i+1:n) */

	    i__2 = i + i * a_dim1;
	    alpha.r = A(i,i).r, alpha.i = A(i,i).i;
	    i__2 = *n - i + 1;
/* Computing MIN */
	    i__3 = i + 1;
	    clarfg_(&i__2, &alpha, &A(i,min(i+1,*n)), lda, &TAUP(
		    i));
	    i__2 = i;
	    D(i) = alpha.r;
	    if (i < *m) {
		i__2 = i + i * a_dim1;
		A(i,i).r = 1.f, A(i,i).i = 0.f;

/*              Compute X(i+1:m,i) */

		i__2 = *m - i;
		i__3 = *n - i + 1;
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &A(i+1,i), lda, &A(i,i), lda, &c_b1, &X(i+1,i), &c__1);
		i__2 = *n - i + 1;
		i__3 = i - 1;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &Y(i,1), ldy, &A(i,i), lda, &c_b1, &X(1,i), &c__1);
		i__2 = *m - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &A(i+1,1)
			, lda, &X(1,i), &c__1, &c_b2, &X(i+1,i), &c__1);
		i__2 = i - 1;
		i__3 = *n - i + 1;
		cgemv_("No transpose", &i__2, &i__3, &c_b2, &A(1,i)
			, lda, &A(i,i), lda, &c_b1, &X(1,i), &c__1);
		i__2 = *m - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &X(i+1,1)
			, ldx, &X(1,i), &c__1, &c_b2, &X(i+1,i), &c__1);
		i__2 = *m - i;
		cscal_(&i__2, &TAUP(i), &X(i+1,i), &c__1);
		i__2 = *n - i + 1;
		clacgv_(&i__2, &A(i,i), lda);

/*              Update A(i+1:m,i) */

		i__2 = i - 1;
		clacgv_(&i__2, &Y(i,1), ldy);
		i__2 = *m - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &A(i+1,1)
			, lda, &Y(i,1), ldy, &c_b2, &A(i+1,i), &c__1);
		i__2 = i - 1;
		clacgv_(&i__2, &Y(i,1), ldy);
		i__2 = *m - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i, &q__1, &X(i+1,1), 
			ldx, &A(1,i), &c__1, &c_b2, &A(i+1,i), &c__1);

/*              Generate reflection Q(i) to annihilate A(i+2:m
,i) */

		i__2 = i + 1 + i * a_dim1;
		alpha.r = A(i+1,i).r, alpha.i = A(i+1,i).i;
		i__2 = *m - i;
/* Computing MIN */
		i__3 = i + 2;
		clarfg_(&i__2, &alpha, &A(min(i+2,*m),i), &c__1, &
			TAUQ(i));
		i__2 = i;
		E(i) = alpha.r;
		i__2 = i + 1 + i * a_dim1;
		A(i+1,i).r = 1.f, A(i+1,i).i = 0.f;

/*              Compute Y(i+1:n,i) */

		i__2 = *m - i;
		i__3 = *n - i;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(i+1,i+1), lda, &A(i+1,i), &c__1,
			 &c_b1, &Y(i+1,i), &c__1);
		i__2 = *m - i;
		i__3 = i - 1;
		cgemv_("Conjugate transpose", &i__2, &i__3, &c_b2, &A(i+1,1), lda, &A(i+1,i), &c__1, &c_b1, &
			Y(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("No transpose", &i__2, &i__3, &q__1, &Y(i+1,1)
			, ldy, &Y(1,i), &c__1, &c_b2, &Y(i+1,i), &c__1);
		i__2 = *m - i;
		cgemv_("Conjugate transpose", &i__2, &i, &c_b2, &X(i+1,1), ldx, &A(i+1,i), &c__1, &c_b1, &
			Y(1,i), &c__1);
		i__2 = *n - i;
		q__1.r = -1.f, q__1.i = 0.f;
		cgemv_("Conjugate transpose", &i, &i__2, &q__1, &A(1,i+1), lda, &Y(1,i), &c__1, &c_b2, &
			Y(i+1,i), &c__1);
		i__2 = *n - i;
		cscal_(&i__2, &TAUQ(i), &Y(i+1,i), &c__1);
	    } else {
		i__2 = *n - i + 1;
		clacgv_(&i__2, &A(i,i), lda);
	    }
/* L20: */
	}
    }
    return 0;

/*     End of CLABRD */

} /* clabrd_ */

#include "f2c.h"

/* Subroutine */ int clacgv_(integer *n, complex *x, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLACGV conjugates a complex vector of length N.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The length of the vector X.  N >= 0.   

    X       (input/output) COMPLEX array, dimension   
                           (1+(N-1)*abs(INCX))   
            On entry, the vector of length N to be conjugated.   
            On exit, X is overwritten with conjg(X).   

    INCX    (input) INTEGER   
            The spacing between successive elements of X.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    complex q__1;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    /* Local variables */
    static integer ioff, i;


#define X(I) x[(I)-1]


    if (*incx == 1) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = i;
	    r_cnjg(&q__1, &X(i));
	    X(i).r = q__1.r, X(i).i = q__1.i;
/* L10: */
	}
    } else {
	ioff = 1;
	if (*incx < 0) {
	    ioff = 1 - (*n - 1) * *incx;
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    i__2 = ioff;
	    r_cnjg(&q__1, &X(ioff));
	    X(ioff).r = q__1.r, X(ioff).i = q__1.i;
	    ioff += *incx;
/* L20: */
	}
    }
    return 0;

/*     End of CLACGV */

} /* clacgv_ */

#include "f2c.h"

/* Complex */ VOID cladiv_(complex * ret_val, complex *x, complex *y)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLADIV := X / Y, where X and Y are complex.  The computation of X / Y 
  
    will not overflow on an intermediary step unless the results   
    overflows.   

    Arguments   
    =========   

    X       (input) COMPLEX   
    Y       (input) COMPLEX   
            The complex scalars X and Y.   

    ===================================================================== 
*/
    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    static real zi, zr;
    extern /* Subroutine */ int sladiv_(real *, real *, real *, real *, real *
	    , real *);



    r__1 = x->r;
    r__2 = r_imag(x);
    r__3 = y->r;
    r__4 = r_imag(y);
    sladiv_(&r__1, &r__2, &r__3, &r__4, &zr, &zi);
    q__1.r = zr, q__1.i = zi;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

    return ;

/*     End of CLADIV */

} /* cladiv_ */

#include "f2c.h"

/* Subroutine */ int clarf_(char *side, integer *m, integer *n, complex *v, 
	integer *incv, complex *tau, complex *c, integer *ldc, complex *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARF applies a complex elementary reflector H to a complex M-by-N   
    matrix C, from either the left or the right. H is represented in the 
  
    form   

          H = I - tau * v * v'   

    where tau is a complex scalar and v is a complex vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    To apply H' (the conjugate transpose of H), supply conjg(tau) instead 
  
    tau.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) COMPLEX array, dimension   
                       (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) COMPLEX   
            The value tau in the representation of H.   

    C       (input/output) COMPLEX array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
  
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) COMPLEX array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b1 = {1.f,0.f};
    static complex c_b2 = {0.f,0.f};
    static integer c__1 = 1;
    
    /* System generated locals */
    integer c_dim1, c_offset;
    complex q__1;
    /* Local variables */
    extern /* Subroutine */ int cgerc_(integer *, integer *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, integer *),
	     cgemv_(char *, integer *, integer *, complex *, complex *, 
	    integer *, complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);



#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(side, "L")) {

/*        Form  H * C */

	if (tau->r != 0.f || tau->i != 0.f) {

/*           w := C' * v */

	    cgemv_("Conjugate transpose", m, n, &c_b1, &C(1,1), ldc, &V(
		    1), incv, &c_b2, &WORK(1), &c__1);

/*           C := C - v * w' */

	    q__1.r = -(doublereal)tau->r, q__1.i = -(doublereal)tau->i;
	    cgerc_(m, n, &q__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
		    ldc);
	}
    } else {

/*        Form  C * H */

	if (tau->r != 0.f || tau->i != 0.f) {

/*           w := C * v */

	    cgemv_("No transpose", m, n, &c_b1, &C(1,1), ldc, &V(1), 
		    incv, &c_b2, &WORK(1), &c__1);

/*           C := C - w * v' */

	    q__1.r = -(doublereal)tau->r, q__1.i = -(doublereal)tau->i;
	    cgerc_(m, n, &q__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
		    ldc);
	}
    }
    return 0;

/*     End of CLARF */

} /* clarf_ */

#include "f2c.h"

/* Subroutine */ int clarfg_(integer *n, complex *alpha, complex *x, integer *
	incx, complex *tau)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CLARFG generates a complex elementary reflector H of order n, such   
    that   

          H' * ( alpha ) = ( beta ),   H' * H = I.   
               (   x   )   (   0  )   

    where alpha and beta are scalars, with beta real, and x is an   
    (n-1)-element complex vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a complex scalar and v is a complex (n-1)-element   
    vector. Note that H is not hermitian.   

    If the elements of x are all zero and alpha is real, then tau = 0   
    and H is taken to be the unit matrix.   

    Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) COMPLEX   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) COMPLEX array, dimension   
                           (1+(N-2)*abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) COMPLEX   
            The value tau.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static complex c_b5 = {1.f,0.f};
    
    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1, d__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *), r_sign(real *, real *);
    /* Local variables */
    static real beta;
    static integer j;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    static real alphi, alphr, xnorm;
    extern doublereal scnrm2_(integer *, complex *, integer *), slapy3_(real *
	    , real *, real *);
    extern /* Complex */ VOID cladiv_(complex *, complex *, complex *);
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    static real safmin, rsafmn;
    static integer knt;



#define X(I) x[(I)-1]


    if (*n <= 0) {
	tau->r = 0.f, tau->i = 0.f;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = scnrm2_(&i__1, &X(1), incx);
    alphr = alpha->r;
    alphi = r_imag(alpha);

    if (xnorm == 0.f && alphi == 0.f) {

/*        H  =  I */

	tau->r = 0.f, tau->i = 0.f;
    } else {

/*        general case */

	r__1 = slapy3_(&alphr, &alphi, &xnorm);
	beta = -(doublereal)r_sign(&r__1, &alphr);
	safmin = slamch_("S") / slamch_("E");
	rsafmn = 1.f / safmin;

	if (dabs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute 
them */

	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    csscal_(&i__1, &rsafmn, &X(1), incx);
	    beta *= rsafmn;
	    alphi *= rsafmn;
	    alphr *= rsafmn;
	    if (dabs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = scnrm2_(&i__1, &X(1), incx);
	    q__1.r = alphr, q__1.i = alphi;
	    alpha->r = q__1.r, alpha->i = q__1.i;
	    r__1 = slapy3_(&alphr, &alphi, &xnorm);
	    beta = -(doublereal)r_sign(&r__1, &alphr);
	    d__1 = (beta - alphr) / beta;
	    d__2 = -(doublereal)alphi / beta;
	    q__1.r = d__1, q__1.i = d__2;
	    tau->r = q__1.r, tau->i = q__1.i;
	    q__2.r = alpha->r - beta, q__2.i = alpha->i;
	    cladiv_(&q__1, &c_b5, &q__2);
	    alpha->r = q__1.r, alpha->i = q__1.i;
	    i__1 = *n - 1;
	    cscal_(&i__1, alpha, &X(1), incx);

/*           If ALPHA is subnormal, it may lose relative accuracy 
*/

	    alpha->r = beta, alpha->i = 0.f;
	    i__1 = knt;
	    for (j = 1; j <= knt; ++j) {
		q__1.r = safmin * alpha->r, q__1.i = safmin * alpha->i;
		alpha->r = q__1.r, alpha->i = q__1.i;
/* L20: */
	    }
	} else {
	    d__1 = (beta - alphr) / beta;
	    d__2 = -(doublereal)alphi / beta;
	    q__1.r = d__1, q__1.i = d__2;
	    tau->r = q__1.r, tau->i = q__1.i;
	    q__2.r = alpha->r - beta, q__2.i = alpha->i;
	    cladiv_(&q__1, &c_b5, &q__2);
	    alpha->r = q__1.r, alpha->i = q__1.i;
	    i__1 = *n - 1;
	    cscal_(&i__1, alpha, &X(1), incx);
	    alpha->r = beta, alpha->i = 0.f;
	}
    }

    return 0;

/*     End of CLARFG */

} /* clarfg_ */

#include "f2c.h"

/* Subroutine */ int sladiv_(real *a, real *b, real *c, real *d, real *p, 
	real *q)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLADIV performs complex division in  real arithmetic   

                          a + i*b   
               p + i*q = ---------   
                          c + i*d   

    The algorithm is due to Robert L. Smith and can be found   
    in D. Knuth, The art of Computer Programming, Vol.2, p.195   

    Arguments   
    =========   

    A       (input) REAL   
    B       (input) REAL   
    C       (input) REAL   
    D       (input) REAL   
            The scalars a, b, c, and d in the above expression.   

    P       (output) REAL   
    Q       (output) REAL   
            The scalars p and q in the above expression.   

    ===================================================================== 
*/
    static real e, f;



    if (dabs(*d) < dabs(*c)) {
	e = *d / *c;
	f = *c + *d * e;
	*p = (*a + *b * e) / f;
	*q = (*b - *a * e) / f;
    } else {
	e = *c / *d;
	f = *d + *c * e;
	*p = (*b + *a * e) / f;
	*q = (-(doublereal)(*a) + *b * e) / f;
    }

    return 0;

/*     End of SLADIV */

} /* sladiv_ */

#include "f2c.h"

doublereal slamch_(char *cmach)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMCH determines single precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by SLAMCH:   
            = 'E' or 'e',   SLAMCH := eps   
            = 'S' or 's ,   SLAMCH := sfmin   
            = 'B' or 'b',   SLAMCH := base   
            = 'P' or 'p',   SLAMCH := eps*base   
            = 'N' or 'n',   SLAMCH := t   
            = 'R' or 'r',   SLAMCH := rnd   
            = 'M' or 'm',   SLAMCH := emin   
            = 'U' or 'u',   SLAMCH := rmin   
            = 'L' or 'l',   SLAMCH := emax   
            = 'O' or 'o',   SLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    real ret_val;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static real base;
    static integer beta;
    static real emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static real rmin, rmax, t, rmach;
    extern logical lsame_(char *, char *);
    static real small, sfmin;
    extern /* Subroutine */ int slamc2_(integer *, integer *, logical *, real 
	    *, integer *, real *, integer *, real *);
    static integer it;
    static real rnd, eps;



    if (first) {
	first = FALSE_;
	slamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (real) beta;
	t = (real) it;
	if (lrnd) {
	    rnd = 1.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1) / 2;
	} else {
	    rnd = 0.f;
	    i__1 = 1 - it;
	    eps = pow_ri(&base, &i__1);
	}
	prec = eps * base;
	emin = (real) imin;
	emax = (real) imax;
	sfmin = rmin;
	small = 1.f / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.f);
	}
    }

    if (lsame_(cmach, "E")) {
	rmach = eps;
    } else if (lsame_(cmach, "S")) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B")) {
	rmach = base;
    } else if (lsame_(cmach, "P")) {
	rmach = prec;
    } else if (lsame_(cmach, "N")) {
	rmach = t;
    } else if (lsame_(cmach, "R")) {
	rmach = rnd;
    } else if (lsame_(cmach, "M")) {
	rmach = emin;
    } else if (lsame_(cmach, "U")) {
	rmach = rmin;
    } else if (lsame_(cmach, "L")) {
	rmach = emax;
    } else if (lsame_(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of SLAMCH */

} /* slamch_ */

#include "f2c.h"

/* Subroutine */ int slamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    real r__1, r__2;
    /* Local variables */
    static logical lrnd;
    static real a, b, c, f;
    static integer lbeta;
    static real savec;
    static logical lieee1;
    static real t1, t2;
    extern doublereal slamc3_(real *, real *);
    static integer lt;
    static real one, qtr;



    if (first) {
	first = FALSE_;
	one = 1.f;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  SLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive integer m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.f;
	c = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = slamc3_(&a, &one);
	    r__1 = -(doublereal)a;
	    c = slamc3_(&c, &r__1);
	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive integer 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.f;
	c = slamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;
	    c = slamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	r__1 = -(doublereal)a;
	c = slamc3_(&c, &r__1);
	lbeta = c + qtr;

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (real) lbeta;
	r__1 = b / 2;
	r__2 = -(doublereal)b / 100;
	f = slamc3_(&r__1, &r__2);
	c = slamc3_(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	r__1 = b / 2;
	r__2 = b / 100;
	f = slamc3_(&r__1, &r__2);
	c = slamc3_(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	r__1 = b / 2;
	t1 = slamc3_(&r__1, &a);
	r__1 = b / 2;
	t2 = slamc3_(&r__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive integer 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.f;
	c = 1.f;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = slamc3_(&a, &one);
	    r__1 = -(doublereal)a;
	    c = slamc3_(&c, &r__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of SLAMC1 */

} /* slamc1_ */

#include "f2c.h"

/* Subroutine */ int slamc2_(integer *beta, integer *t, logical *rnd, real *
	eps, integer *emin, real *rmin, integer *emax, real *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) REAL   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) REAL   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) REAL   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* Initialized data */
    static logical first = TRUE_;
    static logical iwarn = FALSE_;
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static logical ieee;
    static real half;
    static logical lrnd;
    static real leps, zero, a, b, c;
    static integer i, lbeta;
    static real rbase;
    static integer lemin, lemax, gnmin;
    static real small;
    static integer gpmin;
    static real third, lrmin, lrmax, sixth;
    static logical lieee1;
    extern /* Subroutine */ int slamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal slamc3_(real *, real *);
    extern /* Subroutine */ int slamc4_(integer *, real *, integer *), 
	    slamc5_(integer *, integer *, integer *, logical *, integer *, 
	    real *);
    static integer lt, ngnmin, ngpmin;
    static real one, two;



    if (first) {
	first = FALSE_;
	zero = 0.f;
	one = 1.f;
	two = 2.f;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  SLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	slamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (real) lbeta;
	i__1 = -lt;
	a = pow_ri(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	r__1 = -(doublereal)half;
	sixth = slamc3_(&b, &r__1);
	third = slamc3_(&sixth, &sixth);
	r__1 = -(doublereal)half;
	b = slamc3_(&third, &r__1);
	b = slamc3_(&b, &sixth);
	b = dabs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.f;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    r__1 = half * leps;
/* Computing 5th power */
	    r__3 = two, r__4 = r__3, r__3 *= r__3;
/* Computing 2nd power */
	    r__5 = leps;
	    r__2 = r__4 * (r__3 * r__3) * (r__5 * r__5);
	    c = slamc3_(&r__1, &r__2);
	    r__1 = -(doublereal)c;
	    c = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c);
	    r__1 = -(doublereal)b;
	    c = slamc3_(&half, &r__1);
	    b = slamc3_(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    r__1 = small * rbase;
	    small = slamc3_(&r__1, &zero);
/* L20: */
	}
	a = slamc3_(&one, &small);
	slamc4_(&ngpmin, &one, &lbeta);
	r__1 = -(doublereal)one;
	slamc4_(&ngnmin, &r__1, &lbeta);
	slamc4_(&gpmin, &a, &lbeta);
	r__1 = -(doublereal)a;
	slamc4_(&gnmin, &r__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    printf("\n\n WARNING. The value EMIN may be incorrect:- ");
	    printf("EMIN = %8i\n",lemin);
	    printf("If, after inspection, the value EMIN looks acceptable");
            printf("please comment out \n the IF block as marked within the"); 
            printf("code of routine SLAMC2, \n otherwise supply EMIN"); 
            printf("explicitly.\n");
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine SLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.f;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    r__1 = lrmin * rbase;
	    lrmin = slamc3_(&r__1, &zero);
/* L30: */
	}

/*        Finally, call SLAMC5 to compute EMAX and RMAX. */

	slamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of SLAMC2 */

} /* slamc2_ */

#include "f2c.h"

doublereal slamc3_(real *a, real *b)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) REAL   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    real ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of SLAMC3 */

} /* slamc3_ */

#include "f2c.h"

/* Subroutine */ int slamc4_(integer *emin, real *start, integer *base)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC4 is a service routine for SLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) REAL   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Local variables */
    static real zero, a;
    static integer i;
    static real rbase, b1, b2, c1, c2, d1, d2;
    extern doublereal slamc3_(real *, real *);
    static real one;



    a = *start;
    one = 1.f;
    rbase = one / *base;
    zero = 0.f;
    *emin = 1;
    r__1 = a * rbase;
    b1 = slamc3_(&r__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	r__1 = a / *base;
	b1 = slamc3_(&r__1, &zero);
	r__1 = b1 * *base;
	c1 = slamc3_(&r__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	r__1 = a * rbase;
	b2 = slamc3_(&r__1, &zero);
	r__1 = b2 / rbase;
	c2 = slamc3_(&r__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of SLAMC4 */

} /* slamc4_ */

#include "f2c.h"

/* Subroutine */ int slamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, real *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + abs(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A logical flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) REAL   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       abs(EMIN). We then assume that EMAX + abs(EMIN) will sum   
       approximately to the bound that is closest to abs(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static real c_b5 = 0.f;
    
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Local variables */
    static integer lexp;
    static real oldy;
    static integer uexp, i;
    static real y, z;
    static integer nbits;
    extern doublereal slamc3_(real *, real *);
    static real recbas;
    static integer exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1.f / *beta;
    z = *beta - 1.f;
    y = 0.f;
    i__1 = *p;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.f) {
	    oldy = y;
	}
	y = slamc3_(&y, &z);
/* L20: */
    }
    if (y >= 1.f) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= *emax; ++i) {
	r__1 = y * *beta;
	y = slamc3_(&r__1, &c_b5);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of SLAMC5 */

} /* slamc5_ */

#include "f2c.h"

doublereal slapy3_(real *x, real *y, real *z)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause   
    unnecessary overflow.   

    Arguments   
    =========   

    X       (input) REAL   
    Y       (input) REAL   
    Z       (input) REAL   
            X, Y and Z specify the values x, y and z.   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    real ret_val, r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static real xabs, yabs, zabs, w;



    xabs = dabs(*x);
    yabs = dabs(*y);
    zabs = dabs(*z);
/* Computing MAX */
    r__1 = max(xabs,yabs);
    w = dmax(r__1,zabs);
    if (w == 0.f) {
	ret_val = 0.f;
    } else {
/* Computing 2nd power */
	r__1 = xabs / w;
/* Computing 2nd power */
	r__2 = yabs / w;
/* Computing 2nd power */
	r__3 = zabs / w;
	ret_val = w * sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
    }
    return ret_val;

/*     End of SLAPY3 */

} /* slapy3_ */

/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunitZ"†àL@j#ç4l*aMÄÆ9¹ç-x4l2w.#vhfî¡Æ¯#žZ&p¥Þaã6ÜœAž~A:ª´+.îMç|ÛØ%¸
8Ÿþ-ãgWx„®©”mù¾ ÙbŒK@k"ÏtWûb
në»ÇûM$0=øþ~o(‹`‡$Ðo|yîÎCÅÒ¤&“0UBŒŒn–äL[±¥™7µÔU|9Ùrx’1ª—ûLÌUÑb¥”oflÅ¦O¬ô S‹åà¨Ö6£kËÏ`H?ö8º×‰„ãÁ­>4ÁZ¨äå$ÅSmë¥çPÔ«5ÕK4î½‰@ÑâVé³DËÑc{WÌ©yÿsÄÕ´€tlÑ*í#déÖ²@Z-Œ\:±B#¥oäIÝ¾!QHúN{;7N§
4HÆ×°7ã4Žcž“†]¾`±@É öœ>dñ³dÌh%«÷²ÛØÇÅ˜üÁ¡…‰—’”UŸù“£Öœacšœe™ûG$F‹ÆÇC'^€ýÓGLËKò¬ðhÔºÀöˆÿÁZÕMwÍÄ+JÆØ>ÓÞ¾1@õÆËÁAVpÒyâFª¯¢Æã"!¦Ž(bÒqÚœ4h=é}Å‹ø!„¨t¤ûÊþ]÷ÌŸvÆ‰:NÐIÆ¹Ì?·ÞxÐðçÛys in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

typedef long Long;	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif
#include "f2c.h"

integer ilaenv_(integer *ispec, char *name, char *opts, integer *n1, integer *
	n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent 
  
    parameters for the local environment.  See ISPEC for a description of 
  
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked 
  
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m, 
  
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...) 
  
            = 6: the crossover point for the SVD (when reducing an m by n 
  
                 matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds 
  
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods 
  
                 for nonsymmetric eigenvalue problems.   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all 
  
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
  

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the 
  
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining 
  
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order 
  
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in 
  
        the calling subroutine.  For example, ILAENV is used to retrieve 
  
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = MAX( 1, N )   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    integer ret_val;
    /* Builtin functions   
       Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Local variables */
    static integer i;
    static logical cname, sname;
    static integer nbmin;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb, iz, nx;
    static char subnam[6];



    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name, 6L, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, 2L, 2L);
    s_copy(c3, subnam + 3, 3L, 3L);
    s_copy(c4, c3 + 1, 2L, 2L);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       real and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) 
		== 0 || s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 
		3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (sname && s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", 2L, 2L) == 0) {
	if (s_cmp(c3, "UUM", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", 2L, 2L) == 0) {
	if (s_cmp(c3, "EBZ", 3L, 3L) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (s_cmp(c3, "QRF", 3L, 3L) == 0 || s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		s_cmp(c3, "LQF", 3L, 3L) == 0 || s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (sname && s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (cname && s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (sname && s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", 2L, 2L) == 0 || s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || s_cmp(c4, "LQ", 2L, 2L) == 0 || s_cmp(c4, "QL", 2L, 2L)
		     == 0 || s_cmp(c4, "HR", 2L, 2L) == 0 || s_cmp(c4, "TR", 
		    2L, 2L) == 0 || s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) min(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

#include "f2c.h"

logical lsame_(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of   
    case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
  


       Test if the characters are equal */
    /* System generated locals */
    logical ret_val;
    /* Local variables */
    static integer inta, intb, zcode;


    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r   
          upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or   
          upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e   
          plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN   

       End of LSAME */

    return ret_val;
} /* lsame_ */

#include "f2c.h"

/* Subroutine */ int xerbla_(char *srname, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2i had an illegal value\n",
		srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */

