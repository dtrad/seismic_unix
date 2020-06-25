
#include "f2c.h"

/* Subroutine */ int cbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, real *d, real *e, complex *vt, integer *ldvt, 
	complex *u, integer *ldu, complex *c, integer *ldc, real *rwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    CBDSQR computes the singular value decomposition (SVD) of a real   
    N-by-N (upper or lower) bidiagonal matrix B:  B = Q * S * P' (P'   
    denotes the transpose of P), where S is a diagonal matrix with   
    non-negative diagonal elements (the singular values of B), and Q   
    and P are orthogonal matrices.   

    The routine computes S, and optionally computes U * Q, P' * VT,   
    or Q' * C, for given complex input matrices U, VT, and C.   

    See "Computing  Small Singular Values of Bidiagonal Matrices With   
    Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,   
    LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,   
    no. 5, pp. 873-912, Sept 1990) and   
    "Accurate singular values and differential qd algorithms," by   
    B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics   
    Department, University of California at Berkeley, July 1992   
    for a detailed description of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  B is upper bidiagonal;   
            = 'L':  B is lower bidiagonal.   

    N       (input) INTEGER   
            The order of the matrix B.  N >= 0.   

    NCVT    (input) INTEGER   
            The number of columns of the matrix VT. NCVT >= 0.   

    NRU     (input) INTEGER   
            The number of rows of the matrix U. NRU >= 0.   

    NCC     (input) INTEGER   
            The number of columns of the matrix C. NCC >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the n diagonal elements of the bidiagonal matrix B. 
  
            On exit, if INFO=0, the singular values of B in decreasing   
            order.   

    E       (input/output) REAL array, dimension (N)   
            On entry, the elements of E contain the   
            offdiagonal elements of of the bidiagonal matrix whose SVD   
            is desired. On normal exit (INFO = 0), E is destroyed.   
            If the algorithm does not converge (INFO > 0), D and E   
            will contain the diagonal and superdiagonal elements of a   
            bidiagonal matrix orthogonally equivalent to the one given   
            as input. E(N) is used for workspace.   

    VT      (input/output) COMPLEX array, dimension (LDVT, NCVT)   
            On entry, an N-by-NCVT matrix VT.   
            On exit, VT is overwritten by P' * VT.   
            VT is not referenced if NCVT = 0.   

    LDVT    (input) INTEGER   
            The leading dimension of the array VT.   
            LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.   

    U       (input/output) COMPLEX array, dimension (LDU, N)   
            On entry, an NRU-by-N matrix U.   
            On exit, U is overwritten by U * Q.   
            U is not referenced if NRU = 0.   

    LDU     (input) INTEGER   
            The leading dimension of the array U.  LDU >= max(1,NRU).   

    C       (input/output) COMPLEX array, dimension (LDC, NCC)   
            On entry, an N-by-NCC matrix C.   
            On exit, C is overwritten by Q' * C.   
            C is not referenced if NCC = 0.   

    LDC     (input) INTEGER   
            The leading dimension of the array C.   
            LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.   

    RWORK   (workspace) REAL array, dimension   
              2*N  if only singular values wanted (NCVT = NRU = NCC = 0) 
  
              max( 1, 4*N-4 ) otherwise   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  If INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm did not converge; D and E contain the   
                  elements of a bidiagonal matrix which is orthogonally   
                  similar to the input matrix B;  if INFO = i, i   
                  elements of E have not converged to zero.   

    Internal Parameters   
    ===================   

    TOLMUL  REAL, default = max(10,min(100,EPS**(-1/8)))   
            TOLMUL controls the convergence criterion of the QR loop.   
            If it is positive, TOLMUL*EPS is the desired relative   
               precision in the computed singular values.   
            If it is negative, abs(TOLMUL*EPS*sigma_max) is the   
               desired absolute accuracy in the computed singular   
               values (corresponds to relative accuracy   
               abs(TOLMUL*EPS) in the largest singular value.   
            abs(TOLMUL) should be between 1 and 1/EPS, and preferably   
               between 10 (for fast convergence) and .1/EPS   
               (for there to be some accuracy in the results).   
            Default is to lose at either one eighth or 2 of the   
               available decimal digits in each computed singular value   
               (whichever is smaller).   

    MAXITR  INTEGER, default = 6   
            MAXITR controls the maximum number of passes of the   
            algorithm through its inner loop. The algorithms stops   
            (and so fails to converge) if the number of passes   
            through the inner loop exceeds MAXITR*N**2.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b15 = -.125;
    static integer c__1 = 1;
    static real c_b48 = 1.f;
    static real c_b71 = -1.f;
    
    /* System generated locals */
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, 
	    i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), r_sign(real *
	    , real *);
    /* Local variables */
    static real abse;
    static integer idir;
    static real abss;
    static integer oldm;
    static real cosl;
    static integer isub, iter;
    static real unfl, sinl, cosr, smin, smax, sinr;
    static integer irot;
    extern /* Subroutine */ int slas2_(real *, real *, real *, real *, real *)
	    ;
    static real f, g, h;
    static integer i, j, m;
    static real r;
    extern logical lsame_(char *, char *);
    static real oldcs;
    extern /* Subroutine */ int clasr_(char *, char *, char *, integer *, 
	    integer *, real *, real *, complex *, integer *);
    static integer oldll;
    static real shift, sigmn, oldsn;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    static integer maxit;
    static real sminl, sigmx;
    static integer iuplo;
    extern /* Subroutine */ int csrot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, real *), slasq1_(integer *, real *, 
	    real *, real *, integer *), slasv2_(real *, real *, real *, real *
	    , real *, real *, real *, real *, real *);
    static real cs;
    static integer ll;
    static real sn, mu;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *), xerbla_(char *, integer *);
    static real sminoa;
    extern /* Subroutine */ int slartg_(real *, real *, real *, real *, real *
	    );
    static real thresh;
    static logical rotate;
    static real sminlo;
    static integer nm1;
    static real tolmul;
    static integer nm12, nm13, lll;
    static real eps, sll, tol;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define RWORK(I) rwork[(I)-1]

#define VT(I,J) vt[(I)-1 + ((J)-1)* ( *ldvt)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    iuplo = 0;
    if (lsame_(uplo, "U")) {
	iuplo = 1;
    }
    if (lsame_(uplo, "L")) {
	iuplo = 2;
    }
    if (iuplo == 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ncvt < 0) {
	*info = -3;
    } else if (*nru < 0) {
	*info = -4;
    } else if (*ncc < 0) {
	*info = -5;
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1,*n)) {
	*info = -9;
    } else if (*ldu < max(1,*nru)) {
	*info = -11;
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1,*n)) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CBDSQR", &i__1);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    if (*n == 1) {
	goto L150;
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

    if (! rotate) {
	slasq1_(n, &D(1), &E(1), &RWORK(1), info);
	return 0;
    }

    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;

/*     Get machine constants */

    eps = slamch_("Epsilon");
    unfl = slamch_("Safe minimum");

/*     If matrix lower bidiagonal, rotate to be upper bidiagonal   
       by applying Givens rotations on the left */

    if (iuplo == 2) {
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    slartg_(&D(i), &E(i), &cs, &sn, &r);
	    D(i) = r;
	    E(i) = sn * D(i + 1);
	    D(i + 1) = cs * D(i + 1);
	    RWORK(i) = cs;
	    RWORK(nm1 + i) = sn;
/* L10: */
	}

/*        Update singular vectors if desired */

	if (*nru > 0) {
	    clasr_("R", "V", "F", nru, n, &RWORK(1), &RWORK(*n), &U(1,1),
		     ldu);
	}
	if (*ncc > 0) {
	    clasr_("L", "V", "F", n, ncc, &RWORK(1), &RWORK(*n), &C(1,1),
		     ldc);
	}
    }

/*     Compute singular values to relative accuracy TOL   
       (By setting TOL to be negative, algorithm will compute   
       singular values to absolute accuracy ABS(TOL)*norm(input matrix)) 
  

   Computing MAX   
   Computing MIN */
    d__1 = (doublereal) eps;
    r__3 = 100.f, r__4 = pow_dd(&d__1, &c_b15);
    r__1 = 10.f, r__2 = dmin(r__3,r__4);
    tolmul = dmax(r__1,r__2);
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

    smax = (r__1 = D(*n), dabs(r__1));
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	r__3 = smax, r__4 = (r__1 = D(i), dabs(r__1)), r__3 = max(r__3,r__4), 
		r__4 = (r__2 = E(i), dabs(r__2));
	smax = dmax(r__3,r__4);
/* L20: */
    }
    sminl = 0.f;
    if (tol >= 0.f) {

/*        Relative accuracy desired */

	sminoa = dabs(D(1));
	if (sminoa == 0.f) {
	    goto L40;
	}
	mu = sminoa;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    mu = (r__1 = D(i), dabs(r__1)) * (mu / (mu + (r__2 = E(i - 1), 
		    dabs(r__2))));
	    sminoa = dmin(sminoa,mu);
	    if (sminoa == 0.f) {
		goto L40;
	    }
/* L30: */
	}
L40:
	sminoa /= sqrt((real) (*n));
/* Computing MAX */
	r__1 = tol * sminoa, r__2 = *n * 6 * *n * unfl;
	thresh = dmax(r__1,r__2);
    } else {

/*        Absolute accuracy desired   

   Computing MAX */
	r__1 = dabs(tol) * smax, r__2 = *n * 6 * *n * unfl;
	thresh = dmax(r__1,r__2);
    }

/*     Prepare for main iteration loop for the singular values   
       (MAXIT is the maximum number of passes through the inner   
       loop permitted before nonconvergence signalled.) */

    maxit = *n * 6 * *n;
    iter = 0;
    oldll = -1;
    oldm = -1;

/*     M points to last element of unconverged part of matrix */

    m = *n;

/*     Begin main iteration loop */

L50:

/*     Check for convergence or exceeding iteration count */

    if (m <= 1) {
	goto L150;
    }
    if (iter > maxit) {
	goto L190;
    }

/*     Find diagonal block of matrix to work on */

    if (tol < 0.f && (r__1 = D(m), dabs(r__1)) <= thresh) {
	D(m) = 0.f;
    }
    smax = (r__1 = D(m), dabs(r__1));
    smin = smax;
    i__1 = m;
    for (lll = 1; lll <= m; ++lll) {
	ll = m - lll;
	if (ll == 0) {
	    goto L80;
	}
	abss = (r__1 = D(ll), dabs(r__1));
	abse = (r__1 = E(ll), dabs(r__1));
	if (tol < 0.f && abss <= thresh) {
	    D(ll) = 0.f;
	}
	if (abse <= thresh) {
	    goto L70;
	}
	smin = dmin(smin,abss);
/* Computing MAX */
	r__1 = max(smax,abss);
	smax = dmax(r__1,abse);
/* L60: */
    }
L70:
    E(ll) = 0.f;

/*     Matrix splits since E(LL) = 0 */

    if (ll == m - 1) {

/*        Convergence of bottom singular value, return to top of loop 
*/

	--m;
	goto L50;
    }
L80:
    ++ll;

/*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero */

    if (ll == m - 1) {

/*        2 by 2 block, handle separately */

	slasv2_(&D(m - 1), &E(m - 1), &D(m), &sigmn, &sigmx, &sinr, &cosr, &
		sinl, &cosl);
	D(m - 1) = sigmx;
	E(m - 1) = 0.f;
	D(m) = sigmn;

/*        Compute singular vectors, if desired */

	if (*ncvt > 0) {
	    csrot_(ncvt, &VT(m-1,1), ldvt, &VT(m,1), ldvt, &
		    cosr, &sinr);
	}
	if (*nru > 0) {
	    csrot_(nru, &U(1,m-1), &c__1, &U(1,m), &
		    c__1, &cosl, &sinl);
	}
	if (*ncc > 0) {
	    csrot_(ncc, &C(m-1,1), ldc, &C(m,1), ldc, &cosl, &
		    sinl);
	}
	m += -2;
	goto L50;
    }

/*     If working on new submatrix, choose shift direction   
       (from larger end diagonal element towards smaller) */

    if (ll > oldm || m < oldll) {
	if ((r__1 = D(ll), dabs(r__1)) >= (r__2 = D(m), dabs(r__2))) {

/*           Chase bulge from top (big end) to bottom (small end) 
*/

	    idir = 1;
	} else {

/*           Chase bulge from bottom (big end) to top (small end) 
*/

	    idir = 2;
	}
    }

/*     Apply convergence tests */

    if (idir == 1) {

/*        Run convergence test in forward direction   
          First apply standard test to bottom of matrix */

	if ((r__1 = E(m - 1), dabs(r__1)) <= dabs(tol) * (r__2 = D(m), dabs(
		r__2)) || tol < 0.f && (r__3 = E(m - 1), dabs(r__3)) <= 
		thresh) {
	    E(m - 1) = 0.f;
	    goto L50;
	}

	if (tol >= 0.f) {

/*           If relative accuracy desired,   
             apply convergence criterion forward */

	    mu = (r__1 = D(ll), dabs(r__1));
	    sminl = mu;
	    i__1 = m - 1;
	    for (lll = ll; lll <= m-1; ++lll) {
		if ((r__1 = E(lll), dabs(r__1)) <= tol * mu) {
		    E(lll) = 0.f;
		    goto L50;
		}
		sminlo = sminl;
		mu = (r__1 = D(lll + 1), dabs(r__1)) * (mu / (mu + (r__2 = E(
			lll), dabs(r__2))));
		sminl = dmin(sminl,mu);
/* L90: */
	    }
	}

    } else {

/*        Run convergence test in backward direction   
          First apply standard test to top of matrix */

	if ((r__1 = E(ll), dabs(r__1)) <= dabs(tol) * (r__2 = D(ll), dabs(
		r__2)) || tol < 0.f && (r__3 = E(ll), dabs(r__3)) <= thresh) {
	    E(ll) = 0.f;
	    goto L50;
	}

	if (tol >= 0.f) {

/*           If relative accuracy desired,   
             apply convergence criterion backward */

	    mu = (r__1 = D(m), dabs(r__1));
	    sminl = mu;
	    i__1 = ll;
	    for (lll = m - 1; lll >= ll; --lll) {
		if ((r__1 = E(lll), dabs(r__1)) <= tol * mu) {
		    E(lll) = 0.f;
		    goto L50;
		}
		sminlo = sminl;
		mu = (r__1 = D(lll), dabs(r__1)) * (mu / (mu + (r__2 = E(lll),
			 dabs(r__2))));
		sminl = dmin(sminl,mu);
/* L100: */
	    }
	}
    }
    oldll = ll;
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative   
       accuracy, and if so set the shift to zero.   

   Computing MAX */
    r__1 = eps, r__2 = tol * .01f;
    if (tol >= 0.f && *n * tol * (sminl / smax) <= dmax(r__1,r__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

	shift = 0.f;
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

	if (idir == 1) {
	    sll = (r__1 = D(ll), dabs(r__1));
	    slas2_(&D(m - 1), &E(m - 1), &D(m), &shift, &r);
	} else {
	    sll = (r__1 = D(m), dabs(r__1));
	    slas2_(&D(ll), &E(ll), &D(ll + 1), &shift, &r);
	}

/*        Test if shift negligible, and if so set to zero */

	if (sll > 0.f) {
/* Computing 2nd power */
	    r__1 = shift / sll;
	    if (r__1 * r__1 < eps) {
		shift = 0.f;
	    }
	}
    }

/*     Increment iteration count */

    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

    if (shift == 0.f) {
	if (idir == 1) {

/*           Chase bulge from top to bottom   
             Save cosines and sines for later singular vector upda
tes */

	    cs = 1.f;
	    oldcs = 1.f;
	    r__1 = D(ll) * cs;
	    slartg_(&r__1, &E(ll), &cs, &sn, &r);
	    r__1 = oldcs * r;
	    r__2 = D(ll + 1) * sn;
	    slartg_(&r__1, &r__2, &oldcs, &oldsn, &D(ll));
	    RWORK(1) = cs;
	    RWORK(nm1 + 1) = sn;
	    RWORK(nm12 + 1) = oldcs;
	    RWORK(nm13 + 1) = oldsn;
	    irot = 1;
	    i__1 = m - 1;
	    for (i = ll + 1; i <= m-1; ++i) {
		r__1 = D(i) * cs;
		slartg_(&r__1, &E(i), &cs, &sn, &r);
		E(i - 1) = oldsn * r;
		r__1 = oldcs * r;
		r__2 = D(i + 1) * sn;
		slartg_(&r__1, &r__2, &oldcs, &oldsn, &D(i));
		++irot;
		RWORK(irot) = cs;
		RWORK(irot + nm1) = sn;
		RWORK(irot + nm12) = oldcs;
		RWORK(irot + nm13) = oldsn;
/* L110: */
	    }
	    h = D(m) * cs;
	    D(m) = h * oldcs;
	    E(m - 1) = h * oldsn;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "F", &i__1, ncvt, &RWORK(1), &RWORK(*n), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		clasr_("R", "V", "F", nru, &i__1, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "F", &i__1, ncc, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &C(ll,1), ldc);
	    }

/*           Test convergence */

	    if ((r__1 = E(m - 1), dabs(r__1)) <= thresh) {
		E(m - 1) = 0.f;
	    }

	} else {

/*           Chase bulge from bottom to top   
             Save cosines and sines for later singular vector upda
tes */

	    cs = 1.f;
	    oldcs = 1.f;
	    r__1 = D(m) * cs;
	    slartg_(&r__1, &E(m - 1), &cs, &sn, &r);
	    r__1 = oldcs * r;
	    r__2 = D(m - 1) * sn;
	    slartg_(&r__1, &r__2, &oldcs, &oldsn, &D(m));
	    RWORK(m - ll) = cs;
	    RWORK(m - ll + nm1) = -(doublereal)sn;
	    RWORK(m - ll + nm12) = oldcs;
	    RWORK(m - ll + nm13) = -(doublereal)oldsn;
	    irot = m - ll;
	    i__1 = ll + 1;
	    for (i = m - 1; i >= ll+1; --i) {
		r__1 = D(i) * cs;
		slartg_(&r__1, &E(i - 1), &cs, &sn, &r);
		E(i) = oldsn * r;
		r__1 = oldcs * r;
		r__2 = D(i - 1) * sn;
		slartg_(&r__1, &r__2, &oldcs, &oldsn, &D(i));
		--irot;
		RWORK(irot) = cs;
		RWORK(irot + nm1) = -(doublereal)sn;
		RWORK(irot + nm12) = oldcs;
		RWORK(irot + nm13) = -(doublereal)oldsn;
/* L120: */
	    }
	    h = D(ll) * cs;
	    D(ll) = h * oldcs;
	    E(ll) = h * oldsn;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "B", &i__1, ncvt, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		clasr_("R", "V", "B", nru, &i__1, &RWORK(1), &RWORK(*n), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "B", &i__1, ncc, &RWORK(1), &RWORK(*n), &C(ll,1), ldc);
	    }

/*           Test convergence */

	    if ((r__1 = E(ll), dabs(r__1)) <= thresh) {
		E(ll) = 0.f;
	    }
	}
    } else {

/*        Use nonzero shift */

	if (idir == 1) {

/*           Chase bulge from top to bottom   
             Save cosines and sines for later singular vector upda
tes */

	    f = ((r__1 = D(ll), dabs(r__1)) - shift) * (r_sign(&c_b48, &D(ll))
		     + shift / D(ll));
	    g = E(ll);
	    slartg_(&f, &g, &cosr, &sinr, &r);
	    f = cosr * D(ll) + sinr * E(ll);
	    E(ll) = cosr * E(ll) - sinr * D(ll);
	    g = sinr * D(ll + 1);
	    D(ll + 1) = cosr * D(ll + 1);
	    slartg_(&f, &g, &cosl, &sinl, &r);
	    D(ll) = r;
	    f = cosl * E(ll) + sinl * D(ll + 1);
	    D(ll + 1) = cosl * D(ll + 1) - sinl * E(ll);
	    g = sinl * E(ll + 1);
	    E(ll + 1) = cosl * E(ll + 1);
	    RWORK(1) = cosr;
	    RWORK(nm1 + 1) = sinr;
	    RWORK(nm12 + 1) = cosl;
	    RWORK(nm13 + 1) = sinl;
	    irot = 1;
	    i__1 = m - 2;
	    for (i = ll + 1; i <= m-2; ++i) {
		slartg_(&f, &g, &cosr, &sinr, &r);
		E(i - 1) = r;
		f = cosr * D(i) + sinr * E(i);
		E(i) = cosr * E(i) - sinr * D(i);
		g = sinr * D(i + 1);
		D(i + 1) = cosr * D(i + 1);
		slartg_(&f, &g, &cosl, &sinl, &r);
		D(i) = r;
		f = cosl * E(i) + sinl * D(i + 1);
		D(i + 1) = cosl * D(i + 1) - sinl * E(i);
		g = sinl * E(i + 1);
		E(i + 1) = cosl * E(i + 1);
		++irot;
		RWORK(irot) = cosr;
		RWORK(irot + nm1) = sinr;
		RWORK(irot + nm12) = cosl;
		RWORK(irot + nm13) = sinl;
/* L130: */
	    }
	    slartg_(&f, &g, &cosr, &sinr, &r);
	    E(m - 2) = r;
	    f = cosr * D(m - 1) + sinr * E(m - 1);
	    E(m - 1) = cosr * E(m - 1) - sinr * D(m - 1);
	    g = sinr * D(m);
	    D(m) = cosr * D(m);
	    slartg_(&f, &g, &cosl, &sinl, &r);
	    D(m - 1) = r;
	    f = cosl * E(m - 1) + sinl * D(m);
	    D(m) = cosl * D(m) - sinl * E(m - 1);
	    ++irot;
	    RWORK(irot) = cosr;
	    RWORK(irot + nm1) = sinr;
	    RWORK(irot + nm12) = cosl;
	    RWORK(irot + nm13) = sinl;
	    E(m - 1) = f;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "F", &i__1, ncvt, &RWORK(1), &RWORK(*n), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		clasr_("R", "V", "F", nru, &i__1, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "F", &i__1, ncc, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &C(ll,1), ldc);
	    }

/*           Test convergence */

	    if ((r__1 = E(m - 1), dabs(r__1)) <= thresh) {
		E(m - 1) = 0.f;
	    }

	} else {

/*           Chase bulge from bottom to top   
             Save cosines and sines for later singular vector upda
tes */

	    f = ((r__1 = D(m), dabs(r__1)) - shift) * (r_sign(&c_b48, &D(m)) 
		    + shift / D(m));
	    g = E(m - 1);
	    slartg_(&f, &g, &cosr, &sinr, &r);
	    f = cosr * D(m) + sinr * E(m - 1);
	    E(m - 1) = cosr * E(m - 1) - sinr * D(m);
	    g = sinr * D(m - 1);
	    D(m - 1) = cosr * D(m - 1);
	    slartg_(&f, &g, &cosl, &sinl, &r);
	    D(m) = r;
	    f = cosl * E(m - 1) + sinl * D(m - 1);
	    D(m - 1) = cosl * D(m - 1) - sinl * E(m - 1);
	    g = sinl * E(m - 2);
	    E(m - 2) = cosl * E(m - 2);
	    RWORK(m - ll) = cosr;
	    RWORK(m - ll + nm1) = -(doublereal)sinr;
	    RWORK(m - ll + nm12) = cosl;
	    RWORK(m - ll + nm13) = -(doublereal)sinl;
	    irot = m - ll;
	    i__1 = ll + 2;
	    for (i = m - 1; i >= ll+2; --i) {
		slartg_(&f, &g, &cosr, &sinr, &r);
		E(i) = r;
		f = cosr * D(i) + sinr * E(i - 1);
		E(i - 1) = cosr * E(i - 1) - sinr * D(i);
		g = sinr * D(i - 1);
		D(i - 1) = cosr * D(i - 1);
		slartg_(&f, &g, &cosl, &sinl, &r);
		D(i) = r;
		f = cosl * E(i - 1) + sinl * D(i - 1);
		D(i - 1) = cosl * D(i - 1) - sinl * E(i - 1);
		g = sinl * E(i - 2);
		E(i - 2) = cosl * E(i - 2);
		--irot;
		RWORK(irot) = cosr;
		RWORK(irot + nm1) = -(doublereal)sinr;
		RWORK(irot + nm12) = cosl;
		RWORK(irot + nm13) = -(doublereal)sinl;
/* L140: */
	    }
	    slartg_(&f, &g, &cosr, &sinr, &r);
	    E(ll + 1) = r;
	    f = cosr * D(ll + 1) + sinr * E(ll);
	    E(ll) = cosr * E(ll) - sinr * D(ll + 1);
	    g = sinr * D(ll);
	    D(ll) = cosr * D(ll);
	    slartg_(&f, &g, &cosl, &sinl, &r);
	    D(ll + 1) = r;
	    f = cosl * E(ll) + sinl * D(ll);
	    D(ll) = cosl * D(ll) - sinl * E(ll);
	    --irot;
	    RWORK(irot) = cosr;
	    RWORK(irot + nm1) = -(doublereal)sinr;
	    RWORK(irot + nm12) = cosl;
	    RWORK(irot + nm13) = -(doublereal)sinl;
	    E(ll) = f;

/*           Test convergence */

	    if ((r__1 = E(ll), dabs(r__1)) <= thresh) {
		E(ll) = 0.f;
	    }

/*           Update singular vectors if desired */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "B", &i__1, ncvt, &RWORK(nm12 + 1), &RWORK(
			nm13 + 1), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;
		clasr_("R", "V", "B", nru, &i__1, &RWORK(1), &RWORK(*n), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;
		clasr_("L", "V", "B", &i__1, ncc, &RWORK(1), &RWORK(*n), &C(ll,1), ldc);
	    }
	}
    }

/*     QR iteration finished, go back and check convergence */

    goto L50;

/*     All singular values converged, so make them positive */

L150:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (D(i) < 0.f) {
	    D(i) = -(doublereal)D(i);

/*           Change sign of singular vectors, if desired */

	    if (*ncvt > 0) {
		csscal_(ncvt, &c_b71, &VT(i,1), ldvt);
	    }
	}
/* L160: */
    }

/*     Sort the singular values into decreasing order (insertion sort on 
  
       singular values, but only one transposition per singular vector) */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {

/*        Scan for smallest D(I) */

	isub = 1;
	smin = D(1);
	i__2 = *n + 1 - i;
	for (j = 2; j <= *n+1-i; ++j) {
	    if (D(j) <= smin) {
		isub = j;
		smin = D(j);
	    }
/* L170: */
	}
	if (isub != *n + 1 - i) {

/*           Swap singular values and vectors */

	    D(isub) = D(*n + 1 - i);
	    D(*n + 1 - i) = smin;
	    if (*ncvt > 0) {
		cswap_(ncvt, &VT(isub,1), ldvt, &VT(*n+1-i,1), ldvt);
	    }
	    if (*nru > 0) {
		cswap_(nru, &U(1,isub), &c__1, &U(1,*n+1-i), &c__1);
	    }
	    if (*ncc > 0) {
		cswap_(ncc, &C(isub,1), ldc, &C(*n+1-i,1), 
			ldc);
	    }
	}
/* L180: */
    }
    goto L210;

/*     Maximum number of iterations exceeded, failure to converge */

L190:
    *info = 0;
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	if (E(i) != 0.f) {
	    ++(*info);
	}
/* L200: */
    }
L210:
    return 0;

/*     End of CBDSQR */

} /* cbdsqr_ */

#include "f2c.h"

/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, real *c, real *s, complex *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    CLASR   performs the transformation   

       A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )   

       A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )   

    where A is an m by n complex matrix and P is an orthogonal matrix,   
    consisting of a sequence of plane rotations determined by the   
    parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l' 
  
    and z = n when SIDE = 'R' or 'r' ):   

    When  DIRECT = 'F' or 'f'  ( Forward sequence ) then   

       P = P( z - 1 )*...*P( 2 )*P( 1 ),   

    and when DIRECT = 'B' or 'b'  ( Backward sequence ) then   

       P = P( 1 )*P( 2 )*...*P( z - 1 ),   

    where  P( k ) is a plane rotation matrix for the following planes:   

       when  PIVOT = 'V' or 'v'  ( Variable pivot ),   
          the plane ( k, k + 1 )   

       when  PIVOT = 'T' or 't'  ( Top pivot ),   
          the plane ( 1, k + 1 )   

       when  PIVOT = 'B' or 'b'  ( Bottom pivot ),   
          the plane ( k, z )   

    c( k ) and s( k )  must contain the  cosine and sine that define the 
  
    matrix  P( k ).  The two by two plane rotation part of the matrix   
    P( k ), R( k ), is assumed to be of the form   

       R( k ) = (  c( k )  s( k ) ).   
                ( -s( k )  c( k ) )   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            Specifies whether the plane rotation matrix P is applied to   
            A on the left or the right.   
            = 'L':  Left, compute A := P*A   
            = 'R':  Right, compute A:= A*P'   

    DIRECT  (input) CHARACTER*1   
            Specifies whether P is a forward or backward sequence of   
            plane rotations.   
            = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )   
            = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )   

    PIVOT   (input) CHARACTER*1   
            Specifies the plane for which P(k) is a plane rotation   
            matrix.   
            = 'V':  Variable pivot, the plane (k,k+1)   
            = 'T':  Top pivot, the plane (1,k+1)   
            = 'B':  Bottom pivot, the plane (k,z)   

    M       (input) INTEGER   
            The number of rows of the matrix A.  If m <= 1, an immediate 
  
            return is effected.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  If n <= 1, an   
            immediate return is effected.   

    C, S    (input) REAL arrays, dimension   
                    (M-1) if SIDE = 'L'   
                    (N-1) if SIDE = 'R'   
            c(k) and s(k) contain the cosine and sine that define the   
            matrix P(k).  The two by two plane rotation part of the   
            matrix P(k), R(k), is assumed to be of the form   
            R( k ) = (  c( k )  s( k ) ).   
                     ( -s( k )  c( k ) )   

    A       (input/output) COMPLEX array, dimension (LDA,N)   
            The m by n matrix A.  On exit, A is overwritten by P*A if   
            SIDE = 'R' or by A*P' if SIDE = 'L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;
    /* Local variables */
    static integer info;
    static complex temp;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static real ctemp, stemp;
    extern /* Subroutine */ int xerbla_(char *, integer *);


#define C(I) c[(I)-1]
#define S(I) s[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! (lsame_(side, "L") || lsame_(side, "R"))) {
	info = 1;
    } else if (! (lsame_(pivot, "V") || lsame_(pivot, "T") || 
	    lsame_(pivot, "B"))) {
	info = 2;
    } else if (! (lsame_(direct, "F") || lsame_(direct, "B")))
	     {
	info = 3;
    } else if (*m < 0) {
	info = 4;
    } else if (*n < 0) {
	info = 5;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("CLASR ", &info);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }
    if (lsame_(side, "L")) {

/*        Form  P * A */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + 1 + i * a_dim1;
			    temp.r = A(j+1,i).r, temp.i = A(j+1,i).i;
			    i__3 = j + 1 + i * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__4 = j + i * a_dim1;
			    q__3.r = stemp * A(j,i).r, q__3.i = stemp * A(j,i).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(j+1,i).r = q__1.r, A(j+1,i).i = q__1.i;
			    i__3 = j + i * a_dim1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__4 = j + i * a_dim1;
			    q__3.r = ctemp * A(j,i).r, q__3.i = ctemp * A(j,i).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + 1 + i * a_dim1;
			    temp.r = A(j+1,i).r, temp.i = A(j+1,i).i;
			    i__2 = j + 1 + i * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__3 = j + i * a_dim1;
			    q__3.r = stemp * A(j,i).r, q__3.i = stemp * A(j,i).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(j+1,i).r = q__1.r, A(j+1,i).i = q__1.i;
			    i__2 = j + i * a_dim1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__3 = j + i * a_dim1;
			    q__3.r = ctemp * A(j,i).r, q__3.i = ctemp * A(j,i).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
/* L30: */
			}
		    }
/* L40: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m;
		for (j = 2; j <= *m; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__3 = j + i * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__4 = i * a_dim1 + 1;
			    q__3.r = stemp * A(1,i).r, q__3.i = stemp * A(1,i).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
			    i__3 = i * a_dim1 + 1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__4 = i * a_dim1 + 1;
			    q__3.r = ctemp * A(1,i).r, q__3.i = ctemp * A(1,i).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(1,i).r = q__1.r, A(1,i).i = q__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m; j >= 2; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__2 = j + i * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__3 = i * a_dim1 + 1;
			    q__3.r = stemp * A(1,i).r, q__3.i = stemp * A(1,i).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
			    i__2 = i * a_dim1 + 1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__3 = i * a_dim1 + 1;
			    q__3.r = ctemp * A(1,i).r, q__3.i = ctemp * A(1,i).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(1,i).r = q__1.r, A(1,i).i = q__1.i;
/* L70: */
			}
		    }
/* L80: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *m - 1;
		for (j = 1; j <= *m-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__3 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__3 = j + i * a_dim1;
			    i__4 = *m + i * a_dim1;
			    q__2.r = stemp * A(*m,i).r, q__2.i = stemp * A(*m,i).i;
			    q__3.r = ctemp * temp.r, q__3.i = ctemp * temp.i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
			    i__3 = *m + i * a_dim1;
			    i__4 = *m + i * a_dim1;
			    q__2.r = ctemp * A(*m,i).r, q__2.i = ctemp * A(*m,i).i;
			    q__3.r = stemp * temp.r, q__3.i = stemp * temp.i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(*m,i).r = q__1.r, A(*m,i).i = q__1.i;
/* L90: */
			}
		    }
/* L100: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *m - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *n;
			for (i = 1; i <= *n; ++i) {
			    i__2 = j + i * a_dim1;
			    temp.r = A(j,i).r, temp.i = A(j,i).i;
			    i__2 = j + i * a_dim1;
			    i__3 = *m + i * a_dim1;
			    q__2.r = stemp * A(*m,i).r, q__2.i = stemp * A(*m,i).i;
			    q__3.r = ctemp * temp.r, q__3.i = ctemp * temp.i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(j,i).r = q__1.r, A(j,i).i = q__1.i;
			    i__2 = *m + i * a_dim1;
			    i__3 = *m + i * a_dim1;
			    q__2.r = ctemp * A(*m,i).r, q__2.i = ctemp * A(*m,i).i;
			    q__3.r = stemp * temp.r, q__3.i = stemp * temp.i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(*m,i).r = q__1.r, A(*m,i).i = q__1.i;
/* L110: */
			}
		    }
/* L120: */
		}
	    }
	}
    } else if (lsame_(side, "R")) {

/*        Form A * P' */

	if (lsame_(pivot, "V")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + (j + 1) * a_dim1;
			    temp.r = A(i,j+1).r, temp.i = A(i,j+1).i;
			    i__3 = i + (j + 1) * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__4 = i + j * a_dim1;
			    q__3.r = stemp * A(i,j).r, q__3.i = stemp * A(i,j).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(i,j+1).r = q__1.r, A(i,j+1).i = q__1.i;
			    i__3 = i + j * a_dim1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__4 = i + j * a_dim1;
			    q__3.r = ctemp * A(i,j).r, q__3.i = ctemp * A(i,j).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L130: */
			}
		    }
/* L140: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + (j + 1) * a_dim1;
			    temp.r = A(i,j+1).r, temp.i = A(i,j+1).i;
			    i__2 = i + (j + 1) * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__3 = i + j * a_dim1;
			    q__3.r = stemp * A(i,j).r, q__3.i = stemp * A(i,j).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(i,j+1).r = q__1.r, A(i,j+1).i = q__1.i;
			    i__2 = i + j * a_dim1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__3 = i + j * a_dim1;
			    q__3.r = ctemp * A(i,j).r, q__3.i = ctemp * A(i,j).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
/* L150: */
			}
		    }
/* L160: */
		}
	    }
	} else if (lsame_(pivot, "T")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__3 = i + j * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__4 = i + a_dim1;
			    q__3.r = stemp * A(i,1).r, q__3.i = steöd\¥ˇTíÜé‡≥LâÏe_óh¶j40K°	…pˇ+5‹üÏ“¿–q2ºõ—mDvînfUZÅÓnAƒŒÿ/OÜÕ˚ù2Ùßˇ⁄ŸJk"@™£Ñ≤ë:Rz¬ÃY|ÇôÁƒıúÚ‡Ö§:Pà—/?†¡f`ªG«
AÏOS|{%¯.58µ}}@d;öÇ¥*ÆM®≈‡U%±Ûwˇ¿4ÓìyöÀ⁄ô™Ï˘,∆óÁ0©+‹Ö¢ˇZ±>#%f'ÍõÑÉ=Œ“∂‰X’·&•R›/f“∂Á‡o·K •f_÷x‚êz	G∑-8ÁÕV2ÆÌËFÏB•∂"‹¨Ä+;J⁄XE˛\Ÿ™ˆj‹¿Mg¶èÙ™ÊÉ/ØõCƒıs’•Z^ìYGÉx—h©Œ[æU?±˙≤˙è¬îÊp’’∫ 9kMOi][(Ö$cÙbç(#F°˝®¸sïx®&∆¢∆…]Ó(◊|µá∂£ˇ65pø:ráÛ˙å≈FHˆW\ı3?oÂÛo◊…˘©π∏ci@\æn°/˚¨qóæàJ‰æ—jœ0&ìS	Bñ4ÕÉÃªΩô{q¸èUêŒn/’	/Ê÷jµíΩ≠é˜ÁsíhÀåÏj	˛ EGsPcá+Û£el>P3Ø.∑Øï˛$5; --j) {
		    ctemp = C(j - 1);
		    stemp = S(j - 1);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__2 = i + j * a_dim1;
			    q__2.r = ctemp * temp.r, q__2.i = ctemp * temp.i;
			    i__3 = i + a_dim1;
			    q__3.r = stemp * A(i,1).r, q__3.i = stemp * A(i,1).i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
			    i__2 = i + a_dim1;
			    q__2.r = stemp * temp.r, q__2.i = stemp * temp.i;
			    i__3 = i + a_dim1;
			    q__3.r = ctemp * A(i,1).r, q__3.i = ctemp * A(i,1).i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(i,1).r = q__1.r, A(i,1).i = q__1.i;
/* L190: */
			}
		    }
/* L200: */
		}
	    }
	} else if (lsame_(pivot, "B")) {
	    if (lsame_(direct, "F")) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__2 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__3 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__3 = i + j * a_dim1;
			    i__4 = i + *n * a_dim1;
			    q__2.r = stemp * A(i,*n).r, q__2.i = stemp * A(i,*n).i;
			    q__3.r = ctemp * temp.r, q__3.i = ctemp * temp.i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
			    i__3 = i + *n * a_dim1;
			    i__4 = i + *n * a_dim1;
			    q__2.r = ctemp * A(i,*n).r, q__2.i = ctemp * A(i,*n).i;
			    q__3.r = stemp * temp.r, q__3.i = stemp * temp.i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(i,*n).r = q__1.r, A(i,*n).i = q__1.i;
/* L210: */
			}
		    }
/* L220: */
		}
	    } else if (lsame_(direct, "B")) {
		for (j = *n - 1; j >= 1; --j) {
		    ctemp = C(j);
		    stemp = S(j);
		    if (ctemp != 1.f || stemp != 0.f) {
			i__1 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__2 = i + j * a_dim1;
			    temp.r = A(i,j).r, temp.i = A(i,j).i;
			    i__2 = i + j * a_dim1;
			    i__3 = i + *n * a_dim1;
			    q__2.r = stemp * A(i,*n).r, q__2.i = stemp * A(i,*n).i;
			    q__3.r = ctemp * temp.r, q__3.i = ctemp * temp.i;
			    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + 
				    q__3.i;
			    A(i,j).r = q__1.r, A(i,j).i = q__1.i;
			    i__2 = i + *n * a_dim1;
			    i__3 = i + *n * a_dim1;
			    q__2.r = ctemp * A(i,*n).r, q__2.i = ctemp * A(i,*n).i;
			    q__3.r = stemp * temp.r, q__3.i = stemp * temp.i;
			    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - 
				    q__3.i;
			    A(i,*n).r = q__1.r, A(i,*n).i = q__1.i;
/* L230: */
			}
		    }
/* L240: */
		}
	    }
	}
    }

    return 0;

/*     End of CLASR */

} /* clasr_ */

#include "f2c.h"

/* Subroutine */ int csrot_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, real *c, real *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;
    /* Local variables */
    static integer i;
    static complex ctemp;
    static integer ix, iy;
/*     applies a plane rotation, where the cos and sin (c and s) are real 
  
       and the vectors cx and cy are complex.   
       jack dongarra, linpack, 3/11/78.   
    ===================================================================== 
  
    
   Parameter adjustments   
       Function Body */
#define CY(I) cy[(I)-1]
#define CX(I) cx[(I)-1]
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments not equal   
            to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = ix;
	q__2.r = *c * CX(ix).r, q__2.i = *c * CX(ix).i;
	i__3 = iy;
	q__3.r = *s * CY(iy).r, q__3.i = *s * CY(iy).i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
	i__2 = iy;
	i__3 = iy;
	q__2.r = *c * CY(iy).r, q__2.i = *c * CY(iy).i;
	i__4 = ix;
	q__3.r = *s * CX(ix).r, q__3.i = *s * CX(ix).i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	CY(iy).r = q__1.r, CY(iy).i = q__1.i;
	i__2 = ix;
	CX(ix).r = ctemp.r, CX(ix).i = ctemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	i__2 = i;
	q__2.r = *c * CX(i).r, q__2.i = *c * CX(i).i;
	i__3 = i;
	q__3.r = *s * CY(i).r, q__3.i = *s * CY(i).i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
	i__2 = i;
	i__3 = i;
	q__2.r = *c * CY(i).r, q__2.i = *c * CY(i).i;
	i__4 = i;
	q__3.r = *s * CX(i).r, q__3.i = *s * CX(i).i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	CY(i).r = q__1.r, CY(i).i = q__1.i;
	i__2 = i;
	CX(i).r = ctemp.r, CX(i).i = ctemp.i;
/* L30: */
    }
    return 0;
} /* csrot_ */


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

/* Subroutine */ int slartg_(real *f, real *g, real *cs, real *sn, real *r)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine SROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in SBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) REAL   
            The first component of vector to be rotated.   

    G       (input) REAL   
            The second component of vector to be rotated.   

    CS      (output) REAL   
            The cosine of the rotation.   

    SN      (output) REAL   
            The sine of the rotation.   

    R       (output) REAL   
            The nonzero component of the rotated vector.   

    ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal), pow_ri(real *, integer *), sqrt(doublereal);
    /* Local variables */
    static integer i;
    static real scale;
    static integer count;
    static real f1, g1, safmn2, safmx2;
    extern doublereal slamch_(char *);
    static real safmin, eps;



    if (first) {
	first = FALSE_;
	safmin = slamch_("S");
	eps = slamch_("E");
	r__1 = slamch_("B");
	i__1 = (integer) (log(safmin / eps) / log(slamch_("B")) / 2.f);
	safmn2 = pow_ri(&r__1, &i__1);
	safmx2 = 1.f / safmn2;
    }
    if (*g == 0.f) {
	*cs = 1.f;
	*sn = 0.f;
	*r = *f;
    } else if (*f == 0.f) {
	*cs = 0.f;
	*sn = 1.f;
	*r = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	r__1 = dabs(f1), r__2 = dabs(g1);
	scale = dmax(r__1,r__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    r__1 = dabs(f1), r__2 = dabs(g1);
	    scale = dmax(r__1,r__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    r__1 = dabs(f1), r__2 = dabs(g1);
	    scale = dmax(r__1,r__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    r__1 = f1;
/* Computing 2nd power */
	    r__2 = g1;
	    *r = sqrt(r__1 * r__1 + r__2 * r__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	}
	if (dabs(*f) > dabs(*g) && *cs < 0.f) {
	    *cs = -(doublereal)(*cs);
	    *sn = -(doublereal)(*sn);
	    *r = -(doublereal)(*r);
	}
    }
    return 0;

/*     End of SLARTG */

} /* slartg_ */

#include "f2c.h"

/* Subroutine */ int slas2_(real *f, real *g, real *h, real *ssmin, real *
	ssmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    SLAS2  computes the singular values of the 2-by-2 matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, SSMIN is the smaller singular value and SSMAX is the   
    larger singular value.   

    Arguments   
    =========   

    F       (input) REAL   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) REAL   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) REAL   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) REAL   
            The smaller singular value.   

    SSMAX   (output) REAL   
            The larger singular value.   

    Further Details   
    ===============   

    Barring over/underflow, all output quantities are correct to within   
    a few units in the last place (ulps), even in the absence of a guard 
  
    digit in addition/subtraction.   

    In IEEE arithmetic, the code works correctly if one matrix element is 
  
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows, or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

    ==================================================================== 
*/
    /* System generated locals */
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static real fhmn, fhmx, c, fa, ga, ha, as, at, au;



    fa = dabs(*f);
    ga = dabs(*g);
    ha = dabs(*h);
    fhmn = dmin(fa,ha);
    fhmx = dmax(fa,ha);
    if (fhmn == 0.f) {
	*ssmin = 0.f;
	if (fhmx == 0.f) {
	    *ssmax = ga;
	} else {
/* Computing 2nd power */
	    r__1 = dmin(fhmx,ga) / dmax(fhmx,ga);
	    *ssmax = dmax(fhmx,ga) * sqrt(r__1 * r__1 + 1.f);
	}
    } else {
	if (ga < fhmx) {
	    as = fhmn / fhmx + 1.f;
	    at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
	    r__1 = ga / fhmx;
	    au = r__1 * r__1;
	    c = 2.f / (sqrt(as * as + au) + sqrt(at * at + au));
	    *ssmin = fhmn * c;
	    *ssmax = fhmx / c;
	} else {
	    au = fhmx / ga;
	    if (au == 0.f) {

/*              Avoid possible harmful underflow if exponent r
ange   
                asymmetric (true SSMIN may not underflow even 
if   
                AU underflows) */

		*ssmin = fhmn * fhmx / ga;
		*ssmax = ga;
	    } else {
		as = fhmn / fhmx + 1.f;
		at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
		r__1 = as * au;
/* Computing 2nd power */
		r__2 = at * au;
		c = 1.f / (sqrt(r__1 * r__1 + 1.f) + sqrt(r__2 * r__2 + 1.f));
		*ssmin = fhmn * c * au;
		*ssmin += *ssmin;
		*ssmax = ga / (c + c);
	    }
	}
    }
    return 0;

/*     End of SLAS2 */

} /* slas2_ */

#include "f2c.h"

/* Subroutine */ int slascl_(char *type, integer *kl, integer *ku, real *
	cfrom, real *cto, integer *m, integer *n, real *a, integer *lda, 
	integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    SLASCL multiplies the M by N real matrix A by the real scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final 
  
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that 
  
    A may be full, upper triangular, lower triangular, upper Hessenberg, 
  
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper 
  
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) REAL   
    CTO     (input) REAL   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static logical done;
    static real ctoc;
    static integer i, j;
    extern logical lsame_(char *, char *);
    static integer itype, k1, k2, k3, k4;
    static real cfrom1;
    extern doublereal slamch_(char *);
    static real cfromc;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static real bignum, smlnum, mul, cto1;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;

    if (lsame_(type, "G")) {
	itype = 0;
    } else if (lsame_(type, "L")) {
	itype = 1;
    } else if (lsame_(type, "U")) {
	itype = 2;
    } else if (lsame_(type, "H")) {
	itype = 3;
    } else if (lsame_(type, "B")) {
	itype = 4;
    } else if (lsame_(type, "Q")) {
	itype = 5;
    } else if (lsame_(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.f) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
	*info = -7;
    } else if (itype <= 3 && *lda < max(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > max(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
		*info = -3;
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLASCL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = slamch_("S");
    bignum = 1.f / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (dabs(cfrom1) > dabs(ctoc) && ctoc != 0.f) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (dabs(cto1) > dabs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) *= mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = j; i <= *m; ++i) {
		A(i,j) *= mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = min(j,*m);
	    for (i = 1; i <= min(j,*m); ++i) {
		A(i,j) *= mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = min(i__3,*m);
	    for (i = 1; i <= min(j+1,*m); ++i) {
		A(i,j) *= mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = min(i__3,i__4);
	    for (i = 1; i <= min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i = max(k1-j,1); i <= k3; ++i) {
		A(i,j) *= mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = min(i__4,i__5);
	    for (i = max(k1-j,k2); i <= min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of SLASCL */

} /* slascl_ */

#include "f2c.h"

/* Subroutine */ int slasq1_(integer *n, real *d, real *e, real *work, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       SLASQ1 computes the singular values of a real N-by-N bidiagonal   
       matrix with diagonal D and off-diagonal E. The singular values are 
  
       computed to high relative accuracy, barring over/underflow or   
       denormalization. The algorithm is described in   

       "Accurate singular values and differential qd algorithms," by   
       K. V. Fernando and B. N. Parlett,   
       Numer. Math., Vol-67, No. 2, pp. 191-230,1994.   

       See also   
       "Implementation of differential qd algorithms," by   
       K. V. Fernando and B. N. Parlett, Technical Report,   
       Department of Mathematics, University of California at Berkeley,   
       1994 (Under preparation).   

       Arguments   
       =========   

    N       (input) INTEGER   
            The number of rows and columns in the matrix. N >= 0.   

    D       (input/output) REAL array, dimension (N)   
            On entry, D contains the diagonal elements of the   
            bidiagonal matrix whose SVD is desired. On normal exit,   
            D contains the singular values in decreasing order.   

    E       (input/output) REAL array, dimension (N)   
            On entry, elements E(1:N-1) contain the off-diagonal elements 
  
            of the bidiagonal matrix whose SVD is desired.   
            On exit, E is overwritten.   

    WORK    (workspace) REAL array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm did not converge;  i   
                  specifies how many superdiagonals did not converge.   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b8 = .125;
    static integer c__1 = 1;
    static integer c__0 = 0;
    
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    /* Local variables */
    static integer kend, ierr;
    extern /* Subroutine */ int slas2_(real *, real *, real *, real *, real *)
	    ;
    static integer i, j, m;
    static real sfmin, sigmn, sigmx;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static real small2;
    extern /* Subroutine */ int slasq2_(integer *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, integer *);
    static integer ke;
    static real dm, dx;
    static integer ny;
    extern doublereal slamch_(char *);
    extern /* Subroutine */ int xerbla_(char *, integer *), slascl_(
	    char *, integer *, integer *, real *, real *, integer *, integer *
	    , real *, integer *, integer *);
    static real thresh;
    extern /* Subroutine */ int slasrt_(char *, integer *, real *, integer *);
    static real tolmul;
    static logical restrt;
    static real scl, eps, tol, sig1, sig2, tol2;



#define WORK(I) work[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -2;
	i__1 = -(*info);
	xerbla_("SLASQ1", &i__1);
	return 0;
    } else if (*n == 0) {
	return 0;
    } else if (*n == 1) {
	D(1) = dabs(D(1));
	return 0;
    } else if (*n == 2) {
	slas2_(&D(1), &E(1), &D(2), &sigmn, &sigmx);
	D(1) = sigmx;
	D(2) = sigmn;
	return 0;
    }

/*     Estimate the largest singular value */

    sigmx = 0.f;
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	r__2 = sigmx, r__3 = (r__1 = E(i), dabs(r__1));
	sigmx = dmax(r__2,r__3);
/* L10: */
    }

/*     Early return if sigmx is zero (matrix is already diagonal) */

    if (sigmx == 0.f) {
	goto L70;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	D(i) = (r__1 = D(i), dabs(r__1));
/* Computing MAX */
	r__1 = sigmx, r__2 = D(i);
	sigmx = dmax(r__1,r__2);
/* L20: */
    }

/*     Get machine parameters */

    eps = slamch_("EPSILON");
    sfmin = slamch_("SAFE MINIMUM");

/*     Compute singular values to relative accuracy TOL   
       It is assumed that tol**2 does not underflow.   

   Computing MAX   
   Computing MIN */
    d__1 = (doublereal) eps;
    r__3 = 100.f, r__4 = pow_dd(&d__1, &c_b8);
    r__1 = 10.f, r__2 = dmin(r__3,r__4);
    tolmul = dmax(r__1,r__2);
    tol = tolmul * eps;
/* Computing 2nd power */
    r__1 = tol;
    tol2 = r__1 * r__1;

    thresh = sigmx * sqrt(sfmin) * tol;

/*     Scale matrix so the square of the largest element is   
       1 / ( 256 * SFMIN ) */

    scl = sqrt(1.f / (sfmin * 256.f));
/* Computing 2nd power */
    r__1 = tolmul;
    small2 = 1.f / (r__1 * r__1 * 256.f);
    scopy_(n, &D(1), &c__1, &WORK(1), &c__1);
    i__1 = *n - 1;
    scopy_(&i__1, &E(1), &c__1, &WORK(*n + 1), &c__1);
    slascl_("G", &c__0, &c__0, &sigmx, &scl, n, &c__1, &WORK(1), n, &ierr)
	    ;
    i__1 = *n - 1;
    i__2 = *n - 1;
    slascl_("G", &c__0, &c__0, &sigmx, &scl, &i__1, &c__1, &WORK(*n + 1), &
	    i__2, &ierr);

/*     Square D and E (the input for the qd algorithm) */

    i__1 = (*n << 1) - 1;
    for (j = 1; j <= (*n<<1)-1; ++j) {
/* Computing 2nd power */
	r__1 = WORK(j);
	WORK(j) = r__1 * r__1;
/* L30: */
    }

/*     Apply qd algorithm */

    m = 0;
    E(*n) = 0.f;
    dx = WORK(1);
    dm = dx;
    ke = 0;
    restrt = FALSE_;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if ((r__1 = E(i), dabs(r__1)) <= thresh || WORK(*n + i) <= tol2 * (dm 
		/ (real) (i - m))) {
	    ny = i - m;
	    if (ny == 1) {
		goto L50;
	    } else if (ny == 2) {
		slas2_(&D(m + 1), &E(m + 1), &D(m + 2), &sig1, &sig2);
		D(m + 1) = sig1;
		D(m + 2) = sig2;
	    } else {
		kend = ke + 1 - m;
		slasq2_(&ny, &D(m + 1), &E(m + 1), &WORK(m + 1), &WORK(m + *n 
			+ 1), &eps, &tol2, &small2, &dm, &kend, info);

/*                 Return, INFO = number of unconverged superd
iagonals */

		if (*info != 0) {
		    *info += i;
		    return 0;
		}

/*                 Undo scaling */

		i__2 = m + ny;
		for (j = m + 1; j <= m+ny; ++j) {
		    D(j) = sqrt(D(j));
/* L40: */
		}
		slascl_("G", &c__0, &c__0, &scl, &sigmx, &ny, &c__1, &D(m + 1)
			, &ny, &ierr);
	    }
L50:
	    m = i;
	    if (i != *n) {
		dx = WORK(i + 1);
		dm = dx;
		ke = i;
		restrt = TRUE_;
	    }
	}
	if (i != *n && ! restrt) {
	    dx = WORK(i + 1) * (dx / (dx + WORK(*n + i)));
	    if (dm > dx) {
		dm = dx;
		ke = i;
	    }
	}
	restrt = FALSE_;
/* L60: */
    }
    kend = ke + 1;

/*     Sort the singular values into decreasing order */

L70:
    slasrt_("D", n, &D(1), info);
    return 0;

/*     End of SLASQ1 */

} /* slasq1_ */

#include "f2c.h"

/* Subroutine */ int slasq2_(integer *m, real *q, real *e, real *qq, real *ee,
	 real *eps, real *tol2, real *small2, real *sup, integer *kend, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       SLASQ2 computes the singular values of a real N-by-N unreduced   
       bidiagonal matrix with squared diagonal elements in Q and   
       squared off-diagonal elements in E. The singular values are   
       computed to relative accuracy TOL, barring over/underflow or   
       denormalization.   

       Arguments   
       =========   

    M       (input) INTEGER   
            The number of rows and columns in the matrix. M >= 0.   

    Q       (output) REAL array, dimension (M)   
            On normal exit, contains the squared singular values.   

    E       (workspace) REAL array, dimension (M)   

    QQ      (input/output) REAL array, dimension (M)   
            On entry, QQ contains the squared diagonal elements of the   
            bidiagonal matrix whose SVD is desired.   
            On exit, QQ is overwritten.   

    EE      (input/output) REAL array, dimension (M)   
            On entry, EE(1:N-1) contains the squared off-diagonal   
            elements of the bidiagonal matrix whose SVD is desired.   
            On exit, EE is overwritten.   

    EPS     (input) REAL   
            Machine epsilon.   

    TOL2    (input) REAL   
            Desired relative accuracy of computed eigenvalues   
            as defined in SLASQ1.   

    SMALL2  (input) REAL   
            A threshold value as defined in SLASQ1.   

    SUP     (input/output) REAL   
            Upper bound for the smallest eigenvalue.   

    KEND    (input/output) INTEGER   
            Index where minimum d occurs.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm did not converge;  i   
                  specifies how many superdiagonals did not converge.   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    integer i_nint(real *);
    /* Local variables */
    static real xinf;
    static integer n;
    static real sigma, qemax;
    static integer iconv;
    extern /* Subroutine */ int slasq3_(integer *, real *, real *, real *, 
	    real *, real *, real *, integer *, integer *, integer *, integer *
	    , real *, real *, real *);
    static integer iphase;
    static real xx, yy;
    static integer off, isp, off1;


#define EE(I) ee[(I)-1]
#define QQ(I) qq[(I)-1]
#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    n = *m;

/*     Set the default maximum number of iterations */

    off = 0;
    off1 = off + 1;
    sigma = 0.f;
    xinf = 0.f;
    iconv = 0;
    iphase = 2;

/*     Try deflation at the bottom   

       1x1 deflation */

L10:
    if (n <= 2) {
	goto L20;
    }
/* Computing MAX */
    r__1 = QQ(n), r__1 = max(r__1,xinf);
    if (EE(n - 1) <= dmax(r__1,*small2) * *tol2) {
	Q(n) = QQ(n);
	--n;
	if (*kend > n) {
	    *kend = n;
	}
/* Computing MIN */
	r__1 = QQ(n), r__2 = QQ(n - 1);
	*sup = dmin(r__1,r__2);
	goto L10;
    }

/*     2x2 deflation   

   Computing MAX */
    r__1 = max(xinf,*small2), r__2 = QQ(n) / (QQ(n) + EE(n - 1) + QQ(n - 1)) *
	     QQ(n - 1);
    if (EE(n - 2) <= dmax(r__1,r__2) * *tol2) {
/* Computing MAX */
	r__1 = QQ(n), r__2 = QQ(n - 1), r__1 = max(r__1,r__2), r__2 = EE(n - 
		1);
	qemax = dmax(r__1,r__2);
	if (qemax != 0.f) {
	    if (qemax == QQ(n - 1)) {
/* Computing 2nd power */
		r__1 = (QQ(n) - QQ(n - 1) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(r__1 * 
			r__1 + EE(n - 1) * 4.f / qemax)) * .5f;
	    } else if (qemax == QQ(n)) {
/* Computing 2nd power */
		r__1 = (QQ(n - 1) - QQ(n) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(r__1 * 
			r__1 + EE(n - 1) * 4.f / qemax)) * .5f;
	    } else {
/* Computing 2nd power */
		r__1 = (QQ(n) - QQ(n - 1) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(r__1 * 
			r__1 + QQ(n - 1) * 4.f / qemax)) * .5f;
	    }
/* Computing MAX */
	    r__1 = QQ(n), r__2 = QQ(n - 1);
/* Computing MIN */
	    r__3 = QQ(n), r__4 = QQ(n - 1);
	    yy = dmax(r__1,r__2) / xx * dmin(r__3,r__4);
	} else {
	    xx = 0.f;
	    yy = 0.f;
	}
	Q(n - 1) = xx;
	Q(n) = yy;
	n += -2;
	if (*kend > n) {
	    *kend = n;
	}
	*sup = QQ(n);
	goto L10;
    }

L20:
    if (n == 0) {

/*         The lower branch is finished */

	if (off == 0) {

/*         No upper branch; return to SLASQ1 */

	    return 0;
	} else {

/*         Going back to upper branch */

	    xinf = 0.f;
	    if (EE(off) > 0.f) {
		isp = i_nint(&EE(off));
		iphase = 1;
	    } else {
		isp = -i_nint(&EE(off));
		iphase = 2;
	    }
	    sigma = E(off);
	    n = off - isp + 1;
	    off1 = isp;
	    off = off1 - 1;
	    if (n <= 2) {
		goto L20;
	    }
	    if (iphase == 1) {
/* Computing MIN */
		r__1 = Q(n + off), r__2 = Q(n - 1 + off), r__1 = min(r__1,
			r__2), r__2 = Q(n - 2 + off);
		*sup = dmin(r__1,r__2);
	    } else {
/* Computing MIN */
		r__1 = QQ(n + off), r__2 = QQ(n - 1 + off), r__1 = min(r__1,
			r__2), r__2 = QQ(n - 2 + off);
		*sup = dmin(r__1,r__2);
	    }
	    *kend = 0;
	    iconv = -3;
	}
    } else if (n == 1) {

/*     1x1 Solver */

	if (iphase == 1) {
	    Q(off1) += sigma;
	} else {
	    Q(off1) = QQ(off1) + sigma;
	}
	n = 0;
	goto L20;

/*     2x2 Solver */

    } else if (n == 2) {
	if (iphase == 2) {
/* Computing MAX */
	    r__1 = QQ(n + off), r__2 = QQ(n - 1 + off), r__1 = max(r__1,r__2),
		     r__2 = EE(n - 1 + off);
	    qemax = dmax(r__1,r__2);
	    if (qemax != 0.f) {
		if (qemax == QQ(n - 1 + off)) {
/* Computing 2nd power */
		    r__1 = (QQ(n + off) - QQ(n - 1 + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + EE(off + n - 1) * 4.f /
			     qemax)) * .5f;
		} else if (qemax == QQ(n + off)) {
/* Computing 2nd power */
		    r__1 = (QQ(n - 1 + off) - QQ(n + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + EE(n - 1 + off) * 4.f /
			     qemax)) * .5f;
		} else {
/* Computing 2nd power */
		    r__1 = (QQ(n + off) - QQ(n - 1 + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + QQ(n - 1 + off) * 4.f /
			     qemax)) * .5f;
		}
/* Computing MAX */
		r__1 = QQ(n + off), r__2 = QQ(n - 1 + off);
/* Computing MIN */
		r__3 = QQ(n + off), r__4 = QQ(n - 1 + off);
		yy = dmax(r__1,r__2) / xx * dmin(r__3,r__4);
	    } else {
		xx = 0.f;
		yy = 0.f;
	    }
	} else {
/* Computing MAX */
	    r__1 = Q(n + off), r__2 = Q(n - 1 + off), r__1 = max(r__1,r__2), 
		    r__2 = E(n - 1 + off);
	    qemax = dmax(r__1,r__2);
	    if (qemax != 0.f) {
		if (qemax == Q(n - 1 + off)) {
/* Computing 2nd power */
		    r__1 = (Q(n + off) - Q(n - 1 + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + E(n - 1 + off) * 4.f / 
			    qemax)) * .5f;
		} else if (qemax == Q(n + off)) {
/* Computing 2nd power */
		    r__1 = (Q(n - 1 + off) - Q(n + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + E(n - 1 + off) * 4.f / 
			    qemax)) * .5f;
		} else {
/* Computing 2nd power */
		    r__1 = (Q(n + off) - Q(n - 1 + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(r__1 * r__1 + Q(n - 1 + off) * 4.f / 
			    qemax)) * .5f;
		}
/* Computing MAX */
		r__1 = Q(n + off), r__2 = Q(n - 1 + off);
/* Computing MIN */
		r__3 = Q(n + off), r__4 = Q(n - 1 + off);
		yy = dmax(r__1,r__2) / xx * dmin(r__3,r__4);
	    } else {
		xx = 0.f;
		yy = 0.f;
	    }
	}
	Q(n - 1 + off) = sigma + xx;
	Q(n + off) = yy + sigma;
	n = 0;
	goto L20;
    }
    slasq3_(&n, &Q(off1), &E(off1), &QQ(off1), &EE(off1), sup, &sigma, kend, &
	    off, &iphase, &iconv, eps, tol2, small2);
    if (*sup < 0.f) {
	*info = n + off;
	return 0;
    }
    off1 = off + 1;
    goto L20;

/*     End of SLASQ2 */

} /* slasq2_ */

#include "f2c.h"

/* Subroutine */ int slasq3_(integer *n, real *q, real *e, real *qq, real *ee,
	 real *sup, real *sigma, integer *kend, integer *off, integer *iphase,
	 integer *iconv, real *eps, real *tol2, real *small2)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       SLASQ3 is the workhorse of the whole bidiagonal SVD algorithm.   
       This can be described as the differential qd with shifts.   

       Arguments   
       =========   

    N       (input/output) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 3.   
            On exit N is non-negative and less than the input value.   

    Q       (input/output) REAL array, dimension (N)   
            Q array in ping (see IPHASE below)   

    E       (input/output) REAL array, dimension (N)   
            E array in ping (see IPHASE below)   

    QQ      (input/output) REAL array, dimension (N)   
            Q array in pong (see IPHASE below)   

    EE      (input/output) REAL array, dimension (N)   
            E array in pong (see IPHASE below)   

    SUP     (input/output) REAL   
            Upper bound for the smallest eigenvalue   

    SIGMA   (input/output) REAL   
            Accumulated shift for the present submatrix   

    KEND    (input/output) INTEGER   
            Index where minimum D(i) occurs in recurrence for   
            splitting criterion   

    OFF     (input/output) INTEGER   
            Offset for arrays   

    IPHASE  (input/output) INTEGER   
            If IPHASE = 1 (ping) then data is in Q and E arrays   
            If IPHASE = 2 (pong) then data is in QQ and EE arrays   

    ICONV   (input) INTEGER   
            If ICONV = 0 a bottom part of a matrix (with a split)   
            If ICONV =-3 a top part of a matrix (with a split)   

    EPS     (input) REAL   
            Machine epsilon   

    TOL2    (input) REAL   
            Square of the relative tolerance TOL as defined in SLASQ1   

    SMALL2  (input) REAL   
            A threshold value as defined in SLASQ1   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static logical ldef;
    static integer icnt;
    static real tolx, toly, tolz;
    static integer k1end, k2end;
    static real d;
    static integer i;
    static real qemax;
    static integer maxit;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static integer n1, n2;
    static real t1;
    extern /* Subroutine */ int slasq4_(integer *, real *, real *, real *, 
	    real *);
    static integer ic, ke;
    static real dm;
    static integer ip, ks;
    static real xx, yy;
    static logical lsplit;
    static integer ifl;
    static real tau;
    static integer isp;



#define EE(I) ee[(I)-1]
#define QQ(I) qq[(I)-1]
#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    icnt = 0;
    tau = 0.f;
    dm = *sup;
    tolx = *sigma * *tol2;
    tolz = dmax(*small2,*sigma) * *tol2;

/*     Set maximum number of iterations */

    maxit = *n * 100;

/*     Flipping */

    ic = 2;
    if (*n > 3) {
	if (*iphase == 1) {
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {
		if (Q(i) > Q(i + 1)) {
		    ++ic;
		}
		if (E(i) > E(i + 1)) {
		    ++ic;
		}
/* L10: */
	    }
	    if (Q(*n - 1) > Q(*n)) {
		++ic;
	    }
	    if (ic < *n) {
		scopy_(n, &Q(1), &c__1, &QQ(1), &c_n1);
		i__1 = *n - 1;
		scopy_(&i__1, &E(1), &c__1, &EE(1), &c_n1);
		if (*kend != 0) {
		    *kend = *n - *kend + 1;
		}
		*iphase = 2;
	    }
	} else {
	    i__1 = *n - 2;
	    for (i = 1; i <= *n-2; ++i) {
		if (QQ(i) > QQ(i + 1)) {
		    ++ic;
		}
		if (EE(i) > EE(i + 1)) {
		    ++ic;
		}
/* L20: */
	    }
	    if (QQ(*n - 1) > QQ(*n)) {
		++ic;
	    }
	    if (ic < *n) {
		scopy_(n, &QQ(1), &c__1, &Q(1), &c_n1);
		i__1 = *n - 1;
		scopy_(&i__1, &EE(1), &c__1, &E(1), &c_n1);
		if (*kend != 0) {
		    *kend = *n - *kend + 1;
		}
		*iphase = 1;
	    }
	}
    }
    if (*iconv == -3) {
	if (*iphase == 1) {
	    goto L180;
	} else {
	    goto L80;
	}
    }
    if (*iphase == 2) {
	goto L130;
    }

/*     The ping section of the code */

L30:
    ifl = 0;

/*     Compute the shift */

    if (*kend == 0 || *sup == 0.f) {
	tau = 0.f;
    } else if (icnt > 0 && dm <= tolz) {
	tau = 0.f;
    } else {
/* Computing MAX */
	i__1 = 5, i__2 = *n / 32;
	ip = max(i__1,i__2);
	n2 = (ip << 1) + 1;
	if (n2 >= *n) {
	    n1 = 1;
	    n2 = *n;
	} else if (*kend + ip > *n) {
	    n1 = *n - (ip << 1);
	} else if (*kend - ip < 1) {
	    n1 = 1;
	} else {
	    n1 = *kend - ip;
	}
	slasq4_(&n2, &Q(n1), &E(n1), &tau, sup);
    }
L40:
    ++icnt;
    if (icnt > maxit) {
	*sup = -1.f;
	return 0;
    }
    if (tau == 0.f) {

/*     dqd algorithm */

	d = Q(1);
	dm = d;
	ke = 0;
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    QQ(i) = d + E(i);
	    d = d / QQ(i) * Q(i + 1);
	    if (dm > d) {
		dm = d;
		ke = i;
	    }
/* L50: */
	}
	++ke;

/*     Penultimate dqd step (in ping) */

	k2end = ke;
	QQ(*n - 2) = d + E(*n - 2);
	d = d / QQ(*n - 2) * Q(*n - 1);
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	}

/*     Final dqd step (in ping) */

	k1end = ke;
	QQ(*n - 1) = d + E(*n - 1);
	d = d / QQ(*n - 1) * Q(*n);
	if (dm > d) {
	    dm = d;
	    ke = *n;
	'”+£ºöÇÇù√¿∏¬î{XïÏ‹å_,†¯iW˛£yáÿ9fú˚õõYÉ&P]Íóƒ5ŸòyæÚ°ƒìWﬁx⁄aÄ°s›Ï
ú∏ÜÛ—[Î!ñÉ91∂*◊˙¸LI∑E»5ê˛≥,dÅ™˘s\ã…˙¢V!¿Y<±o_πÁ];.oñUûî˙¡¬åy´[À > Ó?{Âuz›€Ål€äê∆0IÛ≥ﬂBwõ#xN˘.x}0∫Ìó∆˚
{Ü|ˆ∂ì˚cÄ´∏ÎL®=’%]E– ∫ÅÃÕ±¶∆Ñ[´¬Û^ˇ$£—È5ª?%¸ØØÎ¯ª¡˘œNq-–!>NÁõkíÍye©1<ô	\#6I"‚æ˝*Àˇ√îø!¶5˝P@v∂Q*/[êÿ˛˝F#Í·ÂUÛ&π∫úª˙‡ {üâÜ ºIö#Ks‰q¬Ç
„7ˇÕÂNƒÚ∑3µæ◊ölwéc0‘8ﬁh1ËQ±ƒÏπv∆r1œ:¢;Áﬂÿ;'fƒö–ÿówqf¡ üM)Ô¥ãå•µ6Â·_©¨.P∂≥/(Üœç,˝ÁÂu(âÄXë†ò~ó’√˛ÀMy„(ÖmhCQ·{cd%»Pÿ„¿_ÒÉÓﬁÂÅ@[pï¥ëﬁá†´˛#	ƒ:”î¢^KÃÕˆÙ›∞LÔ3˛= *n - 1;
	    if (d < 0.f) {
		goto L120;
	    }
	}

/*     Final dqds step (in ping) */

	k1end = ke;
	QQ(*n - 1) = d + E(*n - 1);
	d = d / QQ(*n - 1) * Q(*n) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	QQ(*n) = d;
    }

/*        Convergence when QQ(N) is small (in ping) */

    if ((r__1 = QQ(*n), dabs(r__1)) <= *sigma * *tol2) {
	QQ(*n) = 0.f;
	dm = 0.f;
	ke = *n;
    }
    if (QQ(*n) < 0.f) {
	goto L120;
    }

/*     Non-negative qd array: Update the e's */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	EE(i) = E(i) / QQ(i) * Q(i + 1);
/* L70: */
    }

/*     Updating sigma and iphase in ping */

    *sigma += tau;
    *iphase = 2;
L80:
    tolx = *sigma * *tol2;
    toly = *sigma * *eps;
    tolz = dmax(*sigma,*small2) * *tol2;

/*     Checking for deflation and convergence (in ping) */

L90:
    if (*n <= 2) {
	return 0;
    }

/*        Deflation: bottom 1x1 (in ping) */

    ldef = FALSE_;
    if (EE(*n - 1) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.f) {
	if (EE(*n - 1) <= *eps * (*sigma + QQ(*n))) {
	    if (EE(*n - 1) * (QQ(*n) / (QQ(*n) + *sigma)) <= *tol2 * (QQ(*n) 
		    + *sigma)) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (EE(*n - 1) <= QQ(*n) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
	Q(*n) = QQ(*n) + *sigma;
	--(*n);
	++(*iconv);
	goto L90;
    }

/*        Deflation: bottom 2x2 (in ping) */

    ldef = FALSE_;
    if (EE(*n - 2) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.f) {
	t1 = *sigma + EE(*n - 1) * (*sigma / (*sigma + QQ(*n)));
	if (EE(*n - 2) * (t1 / (QQ(*n - 1) + t1)) <= toly) {
	    if (EE(*n - 2) * (QQ(*n - 1) / (QQ(*n - 1) + t1)) <= tolx) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (EE(*n - 2) <= QQ(*n) / (QQ(*n) + EE(*n - 1) + QQ(*n - 1)) * QQ(*n 
		- 1) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
/* Computing MAX */
	r__1 = QQ(*n), r__2 = QQ(*n - 1), r__1 = max(r__1,r__2), r__2 = EE(*n 
		- 1);
	qemax = dmax(r__1,r__2);
	if (qemax != 0.f) {
	    if (qemax == QQ(*n - 1)) {
/* Computing 2nd power */
		r__1 = (QQ(*n) - QQ(*n - 1) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + EE(*n - 1) * 4.f / qemax)) * .5f;
	    } else if (qemax == QQ(*n)) {
/* Computing 2nd power */
		r__1 = (QQ(*n - 1) - QQ(*n) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + EE(*n - 1) * 4.f / qemax)) * .5f;
	    } else {
/* Computing 2nd power */
		r__1 = (QQ(*n) - QQ(*n - 1) + EE(*n - 1)) / qemax;
		xx = (QQ(*n) + QQ(*n - 1) + EE(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + QQ(*n - 1) * 4.f / qemax)) * .5f;
	    }
/* Computing MAX */
	    r__1 = QQ(*n), r__2 = QQ(*n - 1);
/* Computing MIN */
	    r__3 = QQ(*n), r__4 = QQ(*n - 1);
	    yy = dmax(r__1,r__2) / xx * dmin(r__3,r__4);
	} else {
	    xx = 0.f;
	    yy = 0.f;
	}
	Q(*n - 1) = *sigma + xx;
	Q(*n) = yy + *sigma;
	*n += -2;
	*iconv += 2;
	goto L90;
    }

/*     Updating bounds before going to pong */

    if (*iconv == 0) {
	*kend = ke;
/* Computing MIN */
	r__1 = dm, r__2 = *sup - tau;
	*sup = dmin(r__1,r__2);
    } else if (*iconv > 0) {
/* Computing MIN */
	r__1 = QQ(*n), r__2 = QQ(*n - 1), r__1 = min(r__1,r__2), r__2 = QQ(*n 
		- 2), r__1 = min(r__1,r__2), r__1 = min(r__1,QQ(1)), r__1 = 
		min(r__1,QQ(2));
	*sup = dmin(r__1,QQ(3));
	if (*iconv == 1) {
	    *kend = k1end;
	} else if (*iconv == 2) {
	    *kend = k2end;
	} else {
	    *kend = *n;
	}
	icnt = 0;
	maxit = *n * 100;
    }

/*     Checking for splitting in ping */

    lsplit = FALSE_;
    for (ks = *n - 3; ks >= 3; --ks) {
	if (EE(ks) <= toly) {
/* Computing MIN */
	    r__1 = QQ(ks + 1), r__2 = QQ(ks);
/* Computing MIN */
	    r__3 = QQ(ks + 1), r__4 = QQ(ks);
	    if (EE(ks) * (dmin(r__1,r__2) / (dmin(r__3,r__4) + *sigma)) <= 
		    tolx) {
		lsplit = TRUE_;
		goto L110;
	    }
	}
/* L100: */
    }

    ks = 2;
    if (EE(2) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.f) {
	t1 = *sigma + EE(1) * (*sigma / (*sigma + QQ(1)));
	if (EE(2) * (t1 / (QQ(1) + t1)) <= toly) {
	    if (EE(2) * (QQ(1) / (QQ(1) + t1)) <= tolx) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (EE(2) <= QQ(1) / (QQ(1) + EE(1) + QQ(2)) * QQ(2) * *tol2) {
	    lsplit = TRUE_;
	}
    }
    if (lsplit) {
	goto L110;
    }

    ks = 1;
    if (EE(1) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.f) {
	if (EE(1) <= *eps * (*sigma + QQ(1))) {
	    if (EE(1) * (QQ(1) / (QQ(1) + *sigma)) <= *tol2 * (QQ(1) + *sigma)
		    ) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (EE(1) <= QQ(1) * *tol2) {
	    lsplit = TRUE_;
	}
    }

L110:
    if (lsplit) {
/* Computing MIN */
	r__1 = QQ(*n), r__2 = QQ(*n - 1), r__1 = min(r__1,r__2), r__2 = QQ(*n 
		- 2);
	*sup = dmin(r__1,r__2);
	isp = -(*off + 1);
	*off += ks;
	*n -= ks;
/* Computing MAX */
	i__1 = 1, i__2 = *kend - ks;
	*kend = max(i__1,i__2);
	E(ks) = *sigma;
	EE(ks) = (real) isp;
	*iconv = 0;
	return 0;
    }

/*     Coincidence */

    if (tau == 0.f && dm <= tolz && *kend != *n && *iconv == 0 && icnt > 0) {
	i__1 = *n - ke;
	scopy_(&i__1, &E(ke), &c__1, &QQ(ke), &c__1);
	QQ(*n) = 0.f;
	i__1 = *n - ke;
	scopy_(&i__1, &Q(ke + 1), &c__1, &EE(ke), &c__1);
	*sup = 0.f;
    }
    *iconv = 0;
    goto L130;

/*     A new shift when the previous failed (in ping) */

L120:
    ++ifl;
    *sup = tau;

/*     SUP is small or   
       Too many bad shifts (ping) */

    if (*sup <= tolz || ifl >= 2) {
	tau = 0.f;
	goto L40;

/*     The asymptotic shift (in ping) */

    } else {
/* Computing MAX */
	r__1 = tau + d;
	tau = dmax(r__1,0.f);
	if (tau <= tolz) {
	    tau = 0.f;
	}
	goto L40;
    }

/*     the pong section of the code */

L130:
    ifl = 0;

/*     Compute the shift (in pong) */

    if (*kend == 0 && *sup == 0.f) {
	tau = 0.f;
    } else if (icnt > 0 && dm <= tolz) {
	tau = 0.f;
    } else {
/* Computing MAX */
	i__1 = 5, i__2 = *n / 32;
	ip = max(i__1,i__2);
	n2 = (ip << 1) + 1;
	if (n2 >= *n) {
	    n1 = 1;
	    n2 = *n;
	} else if (*kend + ip > *n) {
	    n1 = *n - (ip << 1);
	} else if (*kend - ip < 1) {
	    n1 = 1;
	} else {
	    n1 = *kend - ip;
	}
	slasq4_(&n2, &QQ(n1), &EE(n1), &tau, sup);
    }
L140:
    ++icnt;
    if (icnt > maxit) {
	*sup = -(doublereal)(*sup);
	return 0;
    }
    if (tau == 0.f) {

/*     The dqd algorithm (in pong) */

	d = QQ(1);
	dm = d;
	ke = 0;
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    Q(i) = d + EE(i);
	    d = d / Q(i) * QQ(i + 1);
	    if (dm > d) {
		dm = d;
		ke = i;
	    }
/* L150: */
	}
	++ke;

/*     Penultimate dqd step (in pong) */

	k2end = ke;
	Q(*n - 2) = d + EE(*n - 2);
	d = d / Q(*n - 2) * QQ(*n - 1);
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	}

/*     Final dqd step (in pong) */

	k1end = ke;
	Q(*n - 1) = d + EE(*n - 1);
	d = d / Q(*n - 1) * QQ(*n);
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	Q(*n) = d;
    } else {

/*     The dqds algorithm (in pong) */

	d = QQ(1) - tau;
	dm = d;
	ke = 0;
	if (d < 0.f) {
	    goto L220;
	}
	i__1 = *n - 3;
	for (i = 1; i <= *n-3; ++i) {
	    Q(i) = d + EE(i);
	    d = d / Q(i) * QQ(i + 1) - tau;
	    if (dm > d) {
		dm = d;
		ke = i;
		if (d < 0.f) {
		    goto L220;
		}
	    }
/* L160: */
	}
	++ke;

/*     Penultimate dqds step (in pong) */

	k2end = ke;
	Q(*n - 2) = d + EE(*n - 2);
	d = d / Q(*n - 2) * QQ(*n - 1) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n - 1;
	    if (d < 0.f) {
		goto L220;
	    }
	}

/*     Final dqds step (in pong) */

	k1end = ke;
	Q(*n - 1) = d + EE(*n - 1);
	d = d / Q(*n - 1) * QQ(*n) - tau;
	if (dm > d) {
	    dm = d;
	    ke = *n;
	}
	Q(*n) = d;
    }

/*        Convergence when is small (in pong) */

    if ((r__1 = Q(*n), dabs(r__1)) <= *sigma * *tol2) {
	Q(*n) = 0.f;
	dm = 0.f;
	ke = *n;
    }
    if (Q(*n) < 0.f) {
	goto L220;
    }

/*     Non-negative qd array: Update the e's */

    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	E(i) = EE(i) / Q(i) * QQ(i + 1);
/* L170: */
    }

/*     Updating sigma and iphase in pong */

    *sigma += tau;
L180:
    *iphase = 1;
    tolx = *sigma * *tol2;
    toly = *sigma * *eps;

/*     Checking for deflation and convergence (in pong) */

L190:
    if (*n <= 2) {
	return 0;
    }

/*        Deflation: bottom 1x1 (in pong) */

    ldef = FALSE_;
    if (E(*n - 1) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.f) {
	if (E(*n - 1) <= *eps * (*sigma + Q(*n))) {
	    if (E(*n - 1) * (Q(*n) / (Q(*n) + *sigma)) <= *tol2 * (Q(*n) + *
		    sigma)) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (E(*n - 1) <= Q(*n) * *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
	Q(*n) += *sigma;
	--(*n);
	++(*iconv);
	goto L190;
    }

/*        Deflation: bottom 2x2 (in pong) */

    ldef = FALSE_;
    if (E(*n - 2) <= tolz) {
	ldef = TRUE_;
    } else if (*sigma > 0.f) {
	t1 = *sigma + E(*n - 1) * (*sigma / (*sigma + Q(*n)));
	if (E(*n - 2) * (t1 / (Q(*n - 1) + t1)) <= toly) {
	    if (E(*n - 2) * (Q(*n - 1) / (Q(*n - 1) + t1)) <= tolx) {
		ldef = TRUE_;
	    }
	}
    } else {
	if (E(*n - 2) <= Q(*n) / (Q(*n) + EE(*n - 1) + Q(*n - 1)) * Q(*n - 1) 
		* *tol2) {
	    ldef = TRUE_;
	}
    }
    if (ldef) {
/* Computing MAX */
	r__1 = Q(*n), r__2 = Q(*n - 1), r__1 = max(r__1,r__2), r__2 = E(*n - 
		1);
	qemax = dmax(r__1,r__2);
	if (qemax != 0.f) {
	    if (qemax == Q(*n - 1)) {
/* Computing 2nd power */
		r__1 = (Q(*n) - Q(*n - 1) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + E(*n - 1) * 4.f / qemax)) * .5f;
	    } else if (qemax == Q(*n)) {
/* Computing 2nd power */
		r__1 = (Q(*n - 1) - Q(*n) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + E(*n - 1) * 4.f / qemax)) * .5f;
	    } else {
/* Computing 2nd power */
		r__1 = (Q(*n) - Q(*n - 1) + E(*n - 1)) / qemax;
		xx = (Q(*n) + Q(*n - 1) + E(*n - 1) + qemax * sqrt(r__1 * 
			r__1 + Q(*n - 1) * 4.f / qemax)) * .5f;
	    }
/* Computing MAX */
	    r__1 = Q(*n), r__2 = Q(*n - 1);
/* Computing MIN */
	    r__3 = Q(*n), r__4 = Q(*n - 1);
	    yy = dmax(r__1,r__2) / xx * dmin(r__3,r__4);
	} else {
	    xx = 0.f;
	    yy = 0.f;
	}
	Q(*n - 1) = *sigma + xx;
	Q(*n) = yy + *sigma;
	*n += -2;
	*iconv += 2;
	goto L190;
    }

/*     Updating bounds before going to pong */

    if (*iconv == 0) {
	*kend = ke;
/* Computing MIN */
	r__1 = dm, r__2 = *sup - tau;
	*sup = dmin(r__1,r__2);
    } else if (*iconv > 0) {
/* Computing MIN */
	r__1 = Q(*n), r__2 = Q(*n - 1), r__1 = min(r__1,r__2), r__2 = Q(*n - 
		2), r__1 = min(r__1,r__2), r__1 = min(r__1,Q(1)), r__1 = min(
		r__1,Q(2));
	*sup = dmin(r__1,Q(3));
	if (*iconv == 1) {
	    *kend = k1end;
	} else if (*iconv == 2) {
	    *kend = k2end;
	} else {
	    *kend = *n;
	}
	icnt = 0;
	maxit = *n * 100;
    }

/*     Checking for splitting in pong */

    lsplit = FALSE_;
    for (ks = *n - 3; ks >= 3; --ks) {
	if (E(ks) <= toly) {
/* Computing MIN */
	    r__1 = Q(ks + 1), r__2 = Q(ks);
/* Computing MIN */
	    r__3 = Q(ks + 1), r__4 = Q(ks);
	    if (E(ks) * (dmin(r__1,r__2) / (dmin(r__3,r__4) + *sigma)) <= 
		    tolx) {
		lsplit = TRUE_;
		goto L210;
	    }
	}
/* L200: */
    }

    ks = 2;
    if (E(2) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.f) {
	t1 = *sigma + E(1) * (*sigma / (*sigma + Q(1)));
	if (E(2) * (t1 / (Q(1) + t1)) <= toly) {
	    if (E(2) * (Q(1) / (Q(1) + t1)) <= tolx) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (E(2) <= Q(1) / (Q(1) + E(1) + Q(2)) * Q(2) * *tol2) {
	    lsplit = TRUE_;
	}
    }
    if (lsplit) {
	goto L210;
    }

    ks = 1;
    if (E(1) <= tolz) {
	lsplit = TRUE_;
    } else if (*sigma > 0.f) {
	if (E(1) <= *eps * (*sigma + Q(1))) {
	    if (E(1) * (Q(1) / (Q(1) + *sigma)) <= *tol2 * (Q(1) + *sigma)) {
		lsplit = TRUE_;
	    }
	}
    } else {
	if (E(1) <= Q(1) * *tol2) {
	    lsplit = TRUE_;
	}
    }

L210:
    if (lsplit) {
/* Computing MIN */
	r__1 = Q(*n), r__2 = Q(*n - 1), r__1 = min(r__1,r__2), r__2 = Q(*n - 
		2);
	*sup = dmin(r__1,r__2);
	isp = *off + 1;
	*off += ks;
/* Computing MAX */
	i__1 = 1, i__2 = *kend - ks;
	*kend = max(i__1,i__2);
	*n -= ks;
	E(ks) = *sigma;
	EE(ks) = (real) isp;
	*iconv = 0;
	return 0;
    }

/*     Coincidence */

    if (tau == 0.f && dm <= tolz && *kend != *n && *iconv == 0 && icnt > 0) {
	i__1 = *n - ke;
	scopy_(&i__1, &EE(ke), &c__1, &Q(ke), &c__1);
	Q(*n) = 0.f;
	i__1 = *n - ke;
	scopy_(&i__1, &QQ(ke + 1), &c__1, &E(ke), &c__1);
	*sup = 0.f;
    }
    *iconv = 0;
    goto L30;

/*     Computation of a new shift when the previous failed (in pong) */

L220:
    ++ifl;
    *sup = tau;

/*     SUP is small or   
       Too many bad shifts (in pong) */

    if (*sup <= tolz || ifl >= 2) {
	tau = 0.f;
	goto L140;

/*     The asymptotic shift (in pong) */

    } else {
/* Computing MAX */
	r__1 = tau + d;
	tau = dmax(r__1,0.f);
	if (tau <= tolz) {
	    tau = 0.f;
	}
	goto L140;
    }

/*     End of SLASQ3 */

    return 0;
} /* slasq3_ */

#include "f2c.h"

/* Subroutine */ int slasq4_(integer *n, real *q, real *e, real *tau, real *
	sup)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       SLASQ4 estimates TAU, the smallest eigenvalue of a matrix. This   
       routine improves the input value of SUP which is an upper bound   
       for the smallest eigenvalue for this matrix .   

       Arguments   
       =========   

    N       (input) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 0.   

    Q       (input) REAL array, dimension (N)   
            Q array   

    E       (input) REAL array, dimension (N)   
            E array   

    TAU     (output) REAL   
            Estimate of the shift   

    SUP     (input/output) REAL   
            Upper bound for the smallest singular value   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static real c_b4 = .7f;
    
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double pow_ri(real *, integer *);
    /* Local variables */
    static real xinf, d;
    static integer i;
    static real dm;
    static integer ifl;



#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    ifl = 1;
/* Computing MIN */
    r__1 = min(*sup,Q(1)), r__1 = min(r__1,Q(2)), r__1 = min(r__1,Q(3)), r__2 
	    = Q(*n), r__1 = min(r__1,r__2), r__2 = Q(*n - 1), r__1 = min(r__1,
	    r__2), r__2 = Q(*n - 2);
    *sup = dmin(r__1,r__2);
    *tau = *sup * .9999f;
    xinf = 0.f;
L10:
    if (ifl == 5) {
	*tau = xinf;
	return 0;
    }
    d = Q(1) - *tau;
    dm = d;
    i__1 = *n - 2;
    for (i = 1; i <= *n-2; ++i) {
	d = d / (d + E(i)) * Q(i + 1) - *tau;
	if (dm > d) {
	    dm = d;
	}
	if (d < 0.f) {
	    *sup = *tau;
/* Computing MAX */
	    r__1 = *sup * pow_ri(&c_b4, &ifl), r__2 = d + *tau;
	    *tau = dmax(r__1,r__2);
	    ++ifl;
	    goto L10;
	}
/* L20: */
    }
    d = d / (d + E(*n - 1)) * Q(*n) - *tau;
    if (dm > d) {
	dm = d;
    }
    if (d < 0.f) {
	*sup = *tau;
/* Computing MAX */
	r__1 = xinf, r__2 = d + *tau;
	xinf = dmax(r__1,r__2);
	if (*sup * pow_ri(&c_b4, &ifl) <= xinf) {
	    *tau = xinf;
	} else {
	    *tau = *sup * pow_ri(&c_b4, &ifl);
	    ++ifl;
	    goto L10;
	}
    } else {
/* Computing MIN */
	r__1 = *sup, r__2 = dm + *tau;
	*sup = dmin(r__1,r__2);
    }
    return 0;

/*     End of SLASQ4 */

} /* slasq4_ */

#include "f2c.h"

/* Subroutine */ int slasrt_(char *id, integer *n, real *d, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    Sort the numbers in D in increasing order (if ID = 'I') or   
    in decreasing order (if ID = 'D' ).   

    Use Quick Sort, reverting to Insertion sort on arrays of   
    size <= 20. Dimension of STACK limits N to about 2**32.   

    Arguments   
    =========   

    ID      (input) CHARACTER*1   
            = 'I': sort D in increasing order;   
            = 'D': sort D in decreasing order.   

    N       (input) INTEGER   
            The length of the array D.   

    D       (input/output) REAL array, dimension (N)   
            On entry, the array to be sorted.   
            On exit, D has been sorted into increasing order   
            (D(1) <= ... <= D(N) ) or into decreasing order   
            (D(1) >= ... >= D(N) ), depending on ID.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input paramters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer endd, i, j;
    extern logical lsame_(char *, char *);
    static integer stack[64]	/* was [2][32] */;
    static real dmnmx, d1, d2, d3;
    static integer start;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer stkpnt, dir;
    static real tmp;


#define STACK(I) stack[(I)]
#define WAS(I) was[(I)]
#define D(I) d[(I)-1]


    *info = 0;
    dir = -1;
    if (lsame_(id, "D")) {
	dir = 0;
    } else if (lsame_(id, "I")) {
	dir = 1;
    }
    if (dir == -1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLASRT", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return 0;
    }

    stkpnt = 1;
    STACK(0) = 1;
    STACK(1) = *n;
L10:
    start = STACK((stkpnt << 1) - 2);
    endd = STACK((stkpnt << 1) - 1);
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {

/*        Do Insertion sort on D( START:ENDD ) */

	if (dir == 0) {

/*           Sort into decreasing order */

	    i__1 = endd;
	    for (i = start + 1; i <= endd; ++i) {
		i__2 = start + 1;
		for (j = i; j >= start+1; --j) {
		    if (D(j) > D(j - 1)) {
			dmnmx = D(j);
			D(j) = D(j - 1);
			D(j - 1) = dmnmx;
		    } else {
			goto L30;
		    }
/* L20: */
		}
L30:
		;
	    }

	} else {

/*           Sort into increasing order */

	    i__1 = endd;
	    for (i = start + 1; i <= endd; ++i) {
		i__2 = start + 1;
		for (j = i; j >= start+1; --j) {
		    if (D(j) < D(j - 1)) {
			dmnmx = D(j);
			D(j) = D(j - 1);
			D(j - 1) = dmnmx;
		    } else {
			goto L50;
		    }
/* L40: */
		}
L50:
		;
	    }

	}

    } else if (endd - start > 20) {

/*        Partition D( START:ENDD ) and stack parts, largest one first
   

          Choose partition entry as median of 3 */

	d1 = D(start);
	d2 = D(endd);
	i = (start + endd) / 2;
	d3 = D(i);
	if (d1 < d2) {
	    if (d3 < d1) {
		dmnmx = d1;
	    } else if (d3 < d2) {
		dmnmx = d3;
	    } else {
		dmnmx = d2;
	    }
	} else {
	    if (d3 < d2) {
		dmnmx = d2;
	    } else if (d3 < d1) {
		dmnmx = d3;
	    } else {
		dmnmx = d1;
	    }
	}

	if (dir == 0) {

/*           Sort into decreasing order */

	    i = start - 1;
	    j = endd + 1;
L60:
L70:
	    --j;
	    if (D(j) < dmnmx) {
		goto L70;
	    }
L80:
	    ++i;
	    if (D(i) > dmnmx) {
		goto L80;
	    }
	    if (i < j) {
		tmp = D(i);
		D(i) = D(j);
		D(j) = tmp;
		goto L60;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
	    } else {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
	    }
	} else {

/*           Sort into increasing order */

	    i = start - 1;
	    j = endd + 1;
L90:
L100:
	    --j;
	    if (D(j) > dmnmx) {
		goto L100;
	    }
L110:
	    ++i;
	    if (D(i) < dmnmx) {
		goto L110;
	    }
	    if (i < j) {
		tmp = D(i);
		D(i) = D(j);
		D(j) = tmp;
		goto L90;
	    }
	    if (j - start > endd - j - 1) {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
	    } else {
		++stkpnt;
		STACK((stkpnt << 1) - 2) = j + 1;
		STACK((stkpnt << 1) - 1) = endd;
		++stkpnt;
		STACK((stkpnt << 1) - 2) = start;
		STACK((stkpnt << 1) - 1) = j;
	    }
	}
    }
    if (stkpnt > 0) {
	goto L10;
    }
    return 0;

/*     End of SLASRT */

} /* slasrt_ */

#include "f2c.h"

/* Subroutine */ int slasv2_(real *f, real *g, real *h, real *ssmin, real *
	ssmax, real *snr, real *csr, real *snl, real *csl)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SLASV2 computes the singular value decomposition of a 2-by-2   
    triangular matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the 
  
    smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and 
  
    right singular vectors for abs(SSMAX), giving the decomposition   

       [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]   
       [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].   

    Arguments   
    =========   

    F       (input) REAL   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) REAL   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) REAL   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) REAL   
            abs(SSMIN) is the smaller singular value.   

    SSMAX   (output) REAL   
            abs(SSMAX) is the larger singular value.   

    SNL     (output) REAL   
    CSL     (output) REAL   
            The vector (CSL, SNL) is a unit left singular vector for the 
  
            singular value abs(SSMAX).   

    SNR     (output) REAL   
    CSR     (output) REAL   
            The vector (CSR, SNR) is a unit right singular vector for the 
  
            singular value abs(SSMAX).   

    Further Details   
    ===============   

    Any input parameter may be aliased with any output parameter.   

    Barring over/underflow and assuming a guard digit in subtraction, all 
  
    output quantities are correct to within a few units in the last   
    place (ulps).   

    In IEEE arithmetic, the code works correctly if one matrix element is 
  
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

   ===================================================================== 
*/
    /* Table of constant values */
    static real c_b3 = 2.f;
    static real c_b4 = 1.f;
    
    /* System generated locals */
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);
    /* Local variables */
    static integer pmax;
    static real temp;
    static logical swap;
    static real a, d, l, m, r, s, t, tsign, fa, ga, ha, ft, gt, ht, mm;
    static logical gasmal;
    extern doublereal slamch_(char *);
    static real tt, clt, crt, slt, srt;




    ft = *f;
    fa = dabs(ft);
    ht = *h;
    ha = dabs(*h);

/*     PMAX points to the maximum absolute element of matrix   
         PMAX = 1 if F largest in absolute values   
         PMAX = 2 if G largest in absolute values   
         PMAX = 3 if H largest in absolute values */

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

/*        Now FA .ge. HA */

    }
    gt = *g;
    ga = dabs(gt);
    if (ga == 0.f) {

/*        Diagonal matrix */

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.f;
	crt = 1.f;
	slt = 0.f;
	srt = 0.f;
    } else {
	gasmal = TRUE_;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < slamch_("EPS")) {

/*              Case of very large GA */

		gasmal = FALSE_;
		*ssmax = ga;
		if (ha > 1.f) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.f;
		slt = ht / gt;
		srt = 1.f;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

/*           Normal case */

	    d = fa - ha;
	    if (d == fa) {

/*              Copes with infinite F or H */

		l = 1.f;
	    } else {
		l = d / fa;
	    }

/*           Note that 0 .le. L .le. 1 */

	    m = gt / ft;

/*           Note that abs(M) .le. 1/macheps */

	    t = 2.f - l;

/*           Note that T .ge. 1 */

	    mm = m * m;
	    tt = t * t;
	    s = sqrt(tt + mm);

/*           Note that 1 .le. S .le. 1 + 1/macheps */

	    if (l == 0.f) {
		r = dabs(m);
	    } else {
		r = sqrt(l * l + mm);
	    }

/*           Note that 0 .le. R .le. 1 + 1/macheps */

	    a = (s + r) * .5f;

/*           Note that 1 .le. A .le. 1 + abs(M) */

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if (mm == 0.f) {

/*              Note that M is very tiny */

		if (l == 0.f) {
		    t = r_sign(&c_b3, &ft) * r_sign(&c_b4, &gt);
		} else {
		    t = gt / r_sign(&d, &ft) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r + l)) * (a + 1.f);
	    }
	    l = sqrt(t * t + 4.f);
	    crt = 2.f / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

/*     Correct signs of SSMAX and SSMIN */

    if (pmax == 1) {
	tsign = r_sign(&c_b4, csr) * r_sign(&c_b4, csl) * r_sign(&c_b4, f);
    }
    if (pmax == 2) {
	tsign = r_sign(&c_b4, snr) * r_sign(&c_b4, csl) * r_sign(&c_b4, g);
    }
    if (pmax == 3) {
	tsign = r_sign(&c_b4, snr) * r_sign(&c_b4, snl) * r_sign(&c_b4, h);
    }
    *ssmax = r_sign(ssmax, &tsign);
    r__1 = tsign * r_sign(&c_b4, f) * r_sign(&c_b4, h);
    *ssmin = r_sign(ssmin, &r__1);
    return 0;

/*     End of SLASV2 */

} /* slasv2_ */

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
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
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

