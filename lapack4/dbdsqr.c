#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dbdsqr_(char *uplo, int *n, int *ncvt, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qbdsqr(char *uplo, int *n, int *ncvt, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qbdsqr_(char *uplo, int *n, int *ncvt, int *
#endif

	nru, int *ncc, LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *vt, 
	int *ldvt, LONG DOUBLE *u, int *ldu, LONG DOUBLE *c, int *
	ldc, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DBDSQR computes the singular value decomposition (SVD) of a real   
    N-by-N (upper or lower) bidiagonal matrix B:  B = Q * S * P' (P'   
    denotes the transpose of P), where S is a diagonal matrix with   
    non-negative diagonal elements (the singular values of B), and Q   
    and P are orthogonal matrices.   

    The routine computes S, and optionally computes U * Q, P' * VT,   
    or Q' * C, for given real input matrices U, VT, and C.   

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

    N       (input) INT   
            The order of the matrix B.  N >= 0.   

    NCVT    (input) INT   
            The number of columns of the matrix VT. NCVT >= 0.   

    NRU     (input) INT   
            The number of rows of the matrix U. NRU >= 0.   

    NCC     (input) INT   
            The number of columns of the matrix C. NCC >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the bidiagonal matrix B. 
  
            On exit, if INFO=0, the singular values of B in decreasing   
            order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the elements of E contain the   
            offdiagonal elements of the bidiagonal matrix whose SVD   
            is desired. On normal exit (INFO = 0), E is destroyed.   
            If the algorithm does not converge (INFO > 0), D and E   
            will contain the diagonal and superdiagonal elements of a   
            bidiagonal matrix orthogonally equivalent to the one given   
            as input. E(N) is used for workspace.   

    VT      (input/output) LONG DOUBLE PRECISION array, dimension (LDVT, NCVT) 
  
            On entry, an N-by-NCVT matrix VT.   
            On exit, VT is overwritten by P' * VT.   
            VT is not referenced if NCVT = 0.   

    LDVT    (input) INT   
            The leading dimension of the array VT.   
            LDVT >= MAX(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.   

    U       (input/output) LONG DOUBLE PRECISION array, dimension (LDU, N)   
            On entry, an NRU-by-N matrix U.   
            On exit, U is overwritten by U * Q.   
            U is not referenced if NRU = 0.   

    LDU     (input) INT   
            The leading dimension of the array U.  LDU >= MAX(1,NRU).   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC, NCC)   
            On entry, an N-by-NCC matrix C.   
            On exit, C is overwritten by Q' * C.   
            C is not referenced if NCC = 0.   

    LDC     (input) INT   
            The leading dimension of the array C.   
            LDC >= MAX(1,N) if NCC > 0; LDC >=1 if NCC = 0.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
              2*N  if only singular values wanted (NCVT = NRU = NCC = 0) 
  
              MAX( 1, 4*N-4 ) otherwise   

    INFO    (output) INT   
            = 0:  successful exit   
            < 0:  If INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm did not converge; D and E contain the   
                  elements of a bidiagonal matrix which is orthogonally   
                  similar to the input matrix B;  if INFO = i, i   
                  elements of E have not converged to zero.   

    Internal Parameters   
    ===================   

    TOLMUL  LONG DOUBLE PRECISION, default = MAX(10,MIN(100,EPS**(-1/8)))   
            TOLMUL controls the convergence criterion of the QR loop.   
            If it is positive, TOLMUL*EPS is the desired relative   
               precision in the computed singular values.   
            If it is negative, ABS(TOLMUL*EPS*sigma_max) is the   
               desired absolute accuracy in the computed singular   
               values (corresponds to relative accuracy   
               ABS(TOLMUL*EPS) in the largest singular value.   
            ABS(TOLMUL) should be between 1 and 1/EPS, and preferably   
               between 10 (for fast convergence) and .1/EPS   
               (for there to be some accuracy in the results).   
            Default is to lose at either one eighth or 2 of the   
               available decimal digits in each computed singular value   
               (whichever is smaller).   

    MAXITR  INT, default = 6   
            MAXITR controls the maximum number of passes of the   
            algorithm through its inner loop. The algorithms stops   
            (and so fails to converge) if the number of passes   
            through the inner loop exceeds MAXITR*N**2.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b15 = -.125;
    static int c__1 = 1;
    static LONG DOUBLE c_b48 = 1.;
    static LONG DOUBLE c_b71 = -1.;
    
    /* System generated locals */
    int  i__1, 
	    i__2;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE abse;
    static int idir;
    static LONG DOUBLE abss;
    static int oldm;
    static LONG DOUBLE P_cosl;
    static int isub, iter;
    static LONG DOUBLE unfl, P_sinl, cosr, smin, smax, sinr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *);
    static int irot;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlas2_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlas2(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlas2_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE f, g, h;
    static int i, j, m;
    static LONG DOUBLE r;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);
    static LONG DOUBLE oldcs;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasr_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasr(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasr_(char *, char *, char *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int oldll;
    static LONG DOUBLE shift, sigmn, oldsn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dswap_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static int maxit;
    static LONG DOUBLE sminl, sigmx;
    static int iuplo;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasq1_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq1(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq1_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif


#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), dlasv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), qlasv2(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), qlasv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE cs;
    static int ll;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE sn, mu;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), xerbla_(char *, 
	    int *);
    static LONG DOUBLE sminoa, thresh;
    static long int rotate;
    static LONG DOUBLE sminlo;
    static int nm1;
    static LONG DOUBLE tolmul;
    static int nm12, nm13, lll;
    static LONG DOUBLE eps, sll, tol;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

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
    } else if ((*ncvt == 0 && *ldvt < 1) || (*ncvt > 0 && *ldvt < MAX(1,*n))) {
	*info = -9;
    } else if (*ldu < MAX(1,*nru)) {
	*info = -11;
    } else if ((*ncc == 0 && *ldc < 1) || (*ncc > 0 && *ldc < MAX(1,*n))) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DBDSQR", &i__1);
	return;
    }
    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	goto L150;
    }

/*     ROTATE is true if any singular vectors desired, false otherwise */

    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;

/*     If no singular vectors desired, use qd algorithm */

    if (! rotate) {

#ifdef PETSC_PREFIX_SUFFIX
	dlasq1_(n, &D(1), &E(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlasq1(n, &D(1), &E(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlasq1_(n, &D(1), &E(1), &WORK(1), info);
#endif

	return;
    }

    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;

/*     Get machine constants */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    unfl = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    unfl = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    unfl = qlamch_("Safe minimum");
#endif


/*     If matrix lower bidiagonal, rotate to be upper bidiagonal   
       by applying Givens rotations on the left */

    if (iuplo == 2) {
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&D(i), &E(i), &cs, &sn, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&D(i), &E(i), &cs, &sn, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&D(i), &E(i), &cs, &sn, &r);
#endif

	    D(i) = r;
	    E(i) = sn * D(i + 1);
	    D(i + 1) = cs * D(i + 1);
	    WORK(i) = cs;
	    WORK(nm1 + i) = sn;
/* L10: */
	}

/*        Update singular vectors if desired */

	if (*nru > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlasr_("R", "V", "F", nru, n, &WORK(1), &WORK(*n), &U(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasr("R", "V", "F", nru, n, &WORK(1), &WORK(*n), &U(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasr_("R", "V", "F", nru, n, &WORK(1), &WORK(*n), &U(1,1), 
#endif

		    ldu);
	}
	if (*ncc > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlasr_("L", "V", "F", n, ncc, &WORK(1), &WORK(*n), &C(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasr("L", "V", "F", n, ncc, &WORK(1), &WORK(*n), &C(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasr_("L", "V", "F", n, ncc, &WORK(1), &WORK(*n), &C(1,1), 
#endif

		    ldc);
	}
    }

/*     Compute singular values to relative accuracy TOL   
       (By setting TOL to be negative, algorithm will compute   
       singular values to absolute accuracy ABS(TOL)*norm(input matrix)) 
  

   Computing MAX   
   Computing MIN */
    d__3 = 100., d__4 = pow(eps, c_b15);
    d__1 = 10., d__2 = MIN(d__3,d__4);
    tolmul = MAX(d__1,d__2);
    tol = tolmul * eps;

/*     Compute approximate maximum, minimum singular values */

    smax = (d__1 = D(*n), ABS(d__1));
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	d__3 = smax, d__4 = (d__1 = D(i), ABS(d__1)), d__3 = MAX(d__3,d__4), 
		d__4 = (d__2 = E(i), ABS(d__2));
	smax = MAX(d__3,d__4);
/* L20: */
    }
    sminl = 0.;
    if (tol >= 0.) {

/*        Relative accuracy desired */

	sminoa = ABS(D(1));
	if (sminoa == 0.) {
	    goto L40;
	}
	mu = sminoa;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    mu = (d__1 = D(i), ABS(d__1)) * (mu / (mu + (d__2 = E(i - 1), ABS(
		    d__2))));
	    sminoa = MIN(sminoa,mu);
	    if (sminoa == 0.) {
		goto L40;
	    }
/* L30: */
	}
L40:
	sminoa /= sqrt((LONG DOUBLE) (*n));
/* Computing MAX */
	d__1 = tol * sminoa, d__2 = *n * 6 * *n * unfl;
	thresh = MAX(d__1,d__2);
    } else {

/*        Absolute accuracy desired   

   Computing MAX */
	d__1 = ABS(tol) * smax, d__2 = *n * 6 * *n * unfl;
	thresh = MAX(d__1,d__2);
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

    if (tol < 0. && (d__1 = D(m), ABS(d__1)) <= thresh) {
	D(m) = 0.;
    }
    smax = (d__1 = D(m), ABS(d__1));
    smin = smax;
    i__1 = m;
    for (lll = 1; lll <= m; ++lll) {
	ll = m - lll;
	if (ll == 0) {
	    goto L80;
	}
	abss = (d__1 = D(ll), ABS(d__1));
	abse = (d__1 = E(ll), ABS(d__1));
	if (tol < 0. && abss <= thresh) {
	    D(ll) = 0.;
	}
	if (abse <= thresh) {
	    goto L70;
	}
	smin = MIN(smin,abss);
/* Computing MAX */
	d__1 = MAX(smax,abss);
	smax = MAX(d__1,abse);
/* L60: */
    }
L70:
    E(ll) = 0.;

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


#ifdef PETSC_PREFIX_SUFFIX
	dlasv2_(&D(m - 1), &E(m - 1), &D(m), &sigmn, &sigmx, &sinr, &cosr, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlasv2(&D(m - 1), &E(m - 1), &D(m), &sigmn, &sigmx, &sinr, &cosr, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlasv2_(&D(m - 1), &E(m - 1), &D(m), &sigmn, &sigmx, &sinr, &cosr, &
#endif

		P_sinl, &P_cosl);
	D(m - 1) = sigmx;
	E(m - 1) = 0.;
	D(m) = sigmn;

/*        Compute singular vectors, if desired */

	if (*ncvt > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(ncvt, &VT(m-1,1), ldvt, &VT(m,1), ldvt, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(ncvt, &VT(m-1,1), ldvt, &VT(m,1), ldvt, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(ncvt, &VT(m-1,1), ldvt, &VT(m,1), ldvt, &
#endif

		    cosr, &sinr);
	}
	if (*nru > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(nru, &U(1,m-1), &c__1, &U(1,m), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(nru, &U(1,m-1), &c__1, &U(1,m), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(nru, &U(1,m-1), &c__1, &U(1,m), &
#endif

		    c__1, &P_cosl, &P_sinl);
	}
	if (*ncc > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(ncc, &C(m-1,1), ldc, &C(m,1), ldc, &P_cosl, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(ncc, &C(m-1,1), ldc, &C(m,1), ldc, &P_cosl, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(ncc, &C(m-1,1), ldc, &C(m,1), ldc, &P_cosl, &
#endif

		    P_sinl);
	}
	m += -2;
	goto L50;
    }

/*     If working on new submatrix, choose shift direction   
       (from larger end diagonal element towards smaller) */

    if (ll > oldm || m < oldll) {
	if ((d__1 = D(ll), ABS(d__1)) >= (d__2 = D(m), ABS(d__2))) {

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

	if ((d__1 = E(m - 1), ABS(d__1)) <= ABS(tol) * (d__2 = D(m), ABS(d__2)
		) || (tol < 0. && (d__3 = E(m - 1), ABS(d__3)) <= thresh)) {
	    E(m - 1) = 0.;
	    goto L50;
	}

	if (tol >= 0.) {

/*           If relative accuracy desired,   
             apply convergence criterion forward */

	    mu = (d__1 = D(ll), ABS(d__1));
	    sminl = mu;
	    i__1 = m - 1;
	    for (lll = ll; lll <= m-1; ++lll) {
		if ((d__1 = E(lll), ABS(d__1)) <= tol * mu) {
		    E(lll) = 0.;
		    goto L50;
		}
		sminlo = sminl;
		mu = (d__1 = D(lll + 1), ABS(d__1)) * (mu / (mu + (d__2 = E(
			lll), ABS(d__2))));
		sminl = MIN(sminl,mu);
/* L90: */
	    }
	}

    } else {

/*        Run convergence test in backward direction   
          First apply standard test to top of matrix */

	if ((d__1 = E(ll), ABS(d__1)) <= ABS(tol) * (d__2 = D(ll), ABS(d__2)) 
		|| (tol < 0. && (d__3 = E(ll), ABS(d__3)) <= thresh)) {
	    E(ll) = 0.;
	    goto L50;
	}

	if (tol >= 0.) {

/*           If relative accuracy desired,   
             apply convergence criterion backward */

	    mu = (d__1 = D(m), ABS(d__1));
	    sminl = mu;
	    i__1 = ll;
	    for (lll = m - 1; lll >= ll; --lll) {
		if ((d__1 = E(lll), ABS(d__1)) <= tol * mu) {
		    E(lll) = 0.;
		    goto L50;
		}
		sminlo = sminl;
		mu = (d__1 = D(lll), ABS(d__1)) * (mu / (mu + (d__2 = E(lll), 
			ABS(d__2))));
		sminl = MIN(sminl,mu);
/* L100: */
	    }
	}
    }
    oldll = ll;
    oldm = m;

/*     Compute shift.  First, test if shifting would ruin relative   
       accuracy, and if so set the shift to zero.   

   Computing MAX */
    d__1 = eps, d__2 = tol * .01;
    if (tol >= 0. && *n * tol * (sminl / smax) <= MAX(d__1,d__2)) {

/*        Use a zero shift to avoid loss of relative accuracy */

	shift = 0.;
    } else {

/*        Compute the shift from 2-by-2 block at end of matrix */

	if (idir == 1) {
	    sll = (d__1 = D(ll), ABS(d__1));

#ifdef PETSC_PREFIX_SUFFIX
	    dlas2_(&D(m - 1), &E(m - 1), &D(m), &shift, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlas2(&D(m - 1), &E(m - 1), &D(m), &shift, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlas2_(&D(m - 1), &E(m - 1), &D(m), &shift, &r);
#endif

	} else {
	    sll = (d__1 = D(m), ABS(d__1));

#ifdef PETSC_PREFIX_SUFFIX
	    dlas2_(&D(ll), &E(ll), &D(ll + 1), &shift, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlas2(&D(ll), &E(ll), &D(ll + 1), &shift, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlas2_(&D(ll), &E(ll), &D(ll + 1), &shift, &r);
#endif

	}

/*        Test if shift negligible, and if so set to zero */

	if (sll > 0.) {
/* Computing 2nd power */
	    d__1 = shift / sll;
	    if (d__1 * d__1 < eps) {
		shift = 0.;
	    }
	}
    }

/*     Increment iteration count */

    iter = iter + m - ll;

/*     If SHIFT = 0, do simplified QR iteration */

    if (shift == 0.) {
	if (idir == 1) {

/*           Chase bulge from top to bottom   
             Save cosines and sines for later singular vector upda
tes */

	    cs = 1.;
	    oldcs = 1.;
	    d__1 = D(ll) * cs;

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&d__1, &E(ll), &cs, &sn, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&d__1, &E(ll), &cs, &sn, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&d__1, &E(ll), &cs, &sn, &r);
#endif

	    d__1 = oldcs * r;
	    d__2 = D(ll + 1) * sn;

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(ll));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&d__1, &d__2, &oldcs, &oldsn, &D(ll));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(ll));
#endif

	    WORK(1) = cs;
	    WORK(nm1 + 1) = sn;
	    WORK(nm12 + 1) = oldcs;
	    WORK(nm13 + 1) = oldsn;
	    irot = 1;
	    i__1 = m - 1;
	    for (i = ll + 1; i <= m-1; ++i) {
		d__1 = D(i) * cs;

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&d__1, &E(i), &cs, &sn, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&d__1, &E(i), &cs, &sn, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&d__1, &E(i), &cs, &sn, &r);
#endif

		E(i - 1) = oldsn * r;
		d__1 = oldcs * r;
		d__2 = D(i + 1) * sn;

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif

		++irot;
		WORK(irot) = cs;
		WORK(irot + nm1) = sn;
		WORK(irot + nm12) = oldcs;
		WORK(irot + nm13) = oldsn;
/* L110: */
	    }
	    h = D(m) * cs;
	    D(m) = h * oldcs;
	    E(m - 1) = h * oldsn;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif

	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif

			+ 1), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif

			+ 1), &C(ll,1), ldc);
	    }

/*           Test convergence */

	    if ((d__1 = E(m - 1), ABS(d__1)) <= thresh) {
		E(m - 1) = 0.;
	    }

	} else {

/*           Chase bulge from bottom to top   
             Save cosines and sines for later singular vector upda
tes */

	    cs = 1.;
	    oldcs = 1.;
	    d__1 = D(m) * cs;

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&d__1, &E(m - 1), &cs, &sn, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&d__1, &E(m - 1), &cs, &sn, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&d__1, &E(m - 1), &cs, &sn, &r);
#endif

	    d__1 = oldcs * r;
	    d__2 = D(m - 1) * sn;

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(m));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&d__1, &d__2, &oldcs, &oldsn, &D(m));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(m));
#endif

	    WORK(m - ll) = cs;
	    WORK(m - ll + nm1) = -sn;
	    WORK(m - ll + nm12) = oldcs;
	    WORK(m - ll + nm13) = -oldsn;
	    irot = m - ll;
	    i__1 = ll + 1;
	    for (i = m - 1; i >= ll+1; --i) {
		d__1 = D(i) * cs;

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&d__1, &E(i - 1), &cs, &sn, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&d__1, &E(i - 1), &cs, &sn, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&d__1, &E(i - 1), &cs, &sn, &r);
#endif

		E(i) = oldsn * r;
		d__1 = oldcs * r;
		d__2 = D(i - 1) * sn;

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&d__1, &d__2, &oldcs, &oldsn, &D(i));
#endif

		--irot;
		WORK(irot) = cs;
		WORK(irot + nm1) = -sn;
		WORK(irot + nm12) = oldcs;
		WORK(irot + nm13) = -oldsn;
/* L120: */
	    }
	    h = D(ll) * cs;
	    D(ll) = h * oldcs;
	    E(ll) = h * oldsn;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif

			nm13 + 1), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif

	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif

	    }

/*           Test convergence */

	    if ((d__1 = E(ll), ABS(d__1)) <= thresh) {
		E(ll) = 0.;
	    }
	}
    } else {

/*        Use nonzero shift */

	if (idir == 1) {

/*           Chase bulge from top to bottom   
             Save cosines and sines for later singular vector upda
tes */

	    f = ((d__1 = D(ll), ABS(d__1)) - shift) * (SIGN(c_b48, D(ll)) 
		    + shift / D(ll));
	    g = E(ll);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

	    f = cosr * D(ll) + sinr * E(ll);
	    E(ll) = cosr * E(ll) - sinr * D(ll);
	    g = sinr * D(ll + 1);
	    D(ll + 1) = cosr * D(ll + 1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

	    D(ll) = r;
	    f = P_cosl * E(ll) + P_sinl * D(ll + 1);
	    D(ll + 1) = P_cosl * D(ll + 1) - P_sinl * E(ll);
	    g = P_sinl * E(ll + 1);
	    E(ll + 1) = P_cosl * E(ll + 1);
	    WORK(1) = cosr;
	    WORK(nm1 + 1) = sinr;
	    WORK(nm12 + 1) = P_cosl;
	    WORK(nm13 + 1) = P_sinl;
	    irot = 1;
	    i__1 = m - 2;
	    for (i = ll + 1; i <= m-2; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

		E(i - 1) = r;
		f = cosr * D(i) + sinr * E(i);
		E(i) = cosr * E(i) - sinr * D(i);
		g = sinr * D(i + 1);
		D(i + 1) = cosr * D(i + 1);

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

		D(i) = r;
		f = P_cosl * E(i) + P_sinl * D(i + 1);
		D(i + 1) = P_cosl * D(i + 1) - P_sinl * E(i);
		g = P_sinl * E(i + 1);
		E(i + 1) = P_cosl * E(i + 1);
		++irot;
		WORK(irot) = cosr;
		WORK(irot + nm1) = sinr;
		WORK(irot + nm12) = P_cosl;
		WORK(irot + nm13) = P_sinl;
/* L130: */
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

	    E(m - 2) = r;
	    f = cosr * D(m - 1) + sinr * E(m - 1);
	    E(m - 1) = cosr * E(m - 1) - sinr * D(m - 1);
	    g = sinr * D(m);
	    D(m) = cosr * D(m);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

	    D(m - 1) = r;
	    f = P_cosl * E(m - 1) + P_sinl * D(m);
	    D(m) = P_cosl * D(m) - P_sinl * E(m - 1);
	    ++irot;
	    WORK(irot) = cosr;
	    WORK(irot + nm1) = sinr;
	    WORK(irot + nm12) = P_cosl;
	    WORK(irot + nm13) = P_sinl;
	    E(m - 1) = f;

/*           Update singular vectors */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "F", &i__1, ncvt, &WORK(1), &WORK(*n), &VT(ll,1), ldvt);
#endif

	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "F", nru, &i__1, &WORK(nm12 + 1), &WORK(nm13 
#endif

			+ 1), &U(1,ll), ldu);
	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "F", &i__1, ncc, &WORK(nm12 + 1), &WORK(nm13 
#endif

			+ 1), &C(ll,1), ldc);
	    }

/*           Test convergence */

	    if ((d__1 = E(m - 1), ABS(d__1)) <= thresh) {
		E(m - 1) = 0.;
	    }

	} else {

/*           Chase bulge from bottom to top   
             Save cosines and sines for later singular vector upda
tes */

	    f = ((d__1 = D(m), ABS(d__1)) - shift) * (SIGN(c_b48, D(m)) + 
		    shift / D(m));
	    g = E(m - 1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

	    f = cosr * D(m) + sinr * E(m - 1);
	    E(m - 1) = cosr * E(m - 1) - sinr * D(m);
	    g = sinr * D(m - 1);
	    D(m - 1) = cosr * D(m - 1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

	    D(m) = r;
	    f = P_cosl * E(m - 1) + P_sinl * D(m - 1);
	    D(m - 1) = P_cosl * D(m - 1) - P_sinl * E(m - 1);
	    g = P_sinl * E(m - 2);
	    E(m - 2) = P_cosl * E(m - 2);
	    WORK(m - ll) = cosr;
	    WORK(m - ll + nm1) = -sinr;
	    WORK(m - ll + nm12) = P_cosl;
	    WORK(m - ll + nm13) = -P_sinl;
	    irot = m - ll;
	    i__1 = ll + 2;
	    for (i = m - 1; i >= ll+2; --i) {

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

		E(i) = r;
		f = cosr * D(i) + sinr * E(i - 1);
		E(i - 1) = cosr * E(i - 1) - sinr * D(i);
		g = sinr * D(i - 1);
		D(i - 1) = cosr * D(i - 1);

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

		D(i) = r;
		f = P_cosl * E(i - 1) + P_sinl * D(i - 1);
		D(i - 1) = P_cosl * D(i - 1) - P_sinl * E(i - 1);
		g = P_sinl * E(i - 2);
		E(i - 2) = P_cosl * E(i - 2);
		--irot;
		WORK(irot) = cosr;
		WORK(irot + nm1) = -sinr;
		WORK(irot + nm12) = P_cosl;
		WORK(irot + nm13) = -P_sinl;
/* L140: */
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &cosr, &sinr, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &cosr, &sinr, &r);
#endif

	    E(ll + 1) = r;
	    f = cosr * D(ll + 1) + sinr * E(ll);
	    E(ll) = cosr * E(ll) - sinr * D(ll + 1);
	    g = sinr * D(ll);
	    D(ll) = cosr * D(ll);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&f, &g, &P_cosl, &P_sinl, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&f, &g, &P_cosl, &P_sinl, &r);
#endif

	    D(ll + 1) = r;
	    f = P_cosl * E(ll) + P_sinl * D(ll);
	    D(ll) = P_cosl * D(ll) - P_sinl * E(ll);
	    --irot;
	    WORK(irot) = cosr;
	    WORK(irot + nm1) = -sinr;
	    WORK(irot + nm12) = P_cosl;
	    WORK(irot + nm13) = -P_sinl;
	    E(ll) = f;

/*           Test convergence */

	    if ((d__1 = E(ll), ABS(d__1)) <= thresh) {
		E(ll) = 0.;
	    }

/*           Update singular vectors if desired */

	    if (*ncvt > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "B", &i__1, ncvt, &WORK(nm12 + 1), &WORK(
#endif

			nm13 + 1), &VT(ll,1), ldvt);
	    }
	    if (*nru > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "B", nru, &i__1, &WORK(1), &WORK(*n), &U(1,ll), ldu);
#endif

	    }
	    if (*ncc > 0) {
		i__1 = m - ll + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("L", "V", "B", &i__1, ncc, &WORK(1), &WORK(*n), &C(ll,1), ldc);
#endif

	    }
	}
    }

/*     QR iteration finished, go back and check convergence */

    goto L50;

/*     All singular values converged, so make them positive */

L150:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (D(i) < 0.) {
	    D(i) = -D(i);

/*           Change sign of singular vectors, if desired */

	    if (*ncvt > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(ncvt, &c_b71, &VT(i,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(ncvt, &c_b71, &VT(i,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(ncvt, &c_b71, &VT(i,1), ldvt);
#endif

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

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(ncvt, &VT(isub,1), ldvt, &VT(*n+1-i,1), ldvt);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(ncvt, &VT(isub,1), ldvt, &VT(*n+1-i,1), ldvt);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(ncvt, &VT(isub,1), ldvt, &VT(*n+1-i,1), ldvt);
#endif

	    }
	    if (*nru > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(nru, &U(1,isub), &c__1, &U(1,*n+1-i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(nru, &U(1,isub), &c__1, &U(1,*n+1-i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(nru, &U(1,isub), &c__1, &U(1,*n+1-i), &c__1);
#endif

	    }
	    if (*ncc > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(ncc, &C(isub,1), ldc, &C(*n+1-i,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(ncc, &C(isub,1), ldc, &C(*n+1-i,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(ncc, &C(isub,1), ldc, &C(*n+1-i,1), 
#endif

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
	if (E(i) != 0.) {
	    ++(*info);
	}
/* L200: */
    }
L210:
    return;

/*     End of DBDSQR */

} /* dbdsqr_ */

