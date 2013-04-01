#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlasq1_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlasq1(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlasq1_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif

	LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ1 computes the singular values of a real N-by-N bidiagonal   
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

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, D contains the diagonal elements of the   
            bidiagonal matrix whose SVD is desired. On normal exit,   
            D contains the singular values in decreasing order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, elements E(1:N-1) contain the off-diagonal elements 
  
            of the bidiagonal matrix whose SVD is desired.   
            On exit, E is overwritten.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm did not converge;  i   
                  specifies how many superdiagonals did not converge.   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b8 = .125;
    static int c__1 = 1;
    static int c__0 = 0;
    
    /* System generated locals */
    int i__1, i__2;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static int kend, ierr;

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
    static int i, j, m;
    static LONG DOUBLE sfmin, sigmn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static LONG DOUBLE sigmx;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasq2_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq2(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq2_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static LONG DOUBLE small2;
    static int ke;
    static LONG DOUBLE dm;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE dx;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlascl_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, int *);
    static int ny;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dlasrt_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlasrt(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlasrt_(
#endif

	    char *, int *, LONG DOUBLE *, int *);
    static LONG DOUBLE thresh, tolmul;
    static long int restrt;
    static LONG DOUBLE scl, eps, tol, sig1, sig2, tol2;



#define WORK(I) work[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -2;
	i__1 = -(*info);
	xerbla_("DLASQ1", &i__1);
	return;
    } else if (*n == 0) {
	return;
    } else if (*n == 1) {
	D(1) = ABS(D(1));
	return;
    } else if (*n == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlas2_(&D(1), &E(1), &D(2), &sigmn, &sigmx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlas2(&D(1), &E(1), &D(2), &sigmn, &sigmx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlas2_(&D(1), &E(1), &D(2), &sigmn, &sigmx);
#endif

	D(1) = sigmx;
	D(2) = sigmn;
	return;
    }

/*     Estimate the largest singular value */

    sigmx = 0.;
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	d__2 = sigmx, d__3 = (d__1 = E(i), ABS(d__1));
	sigmx = MAX(d__2,d__3);
/* L10: */
    }

/*     Early return if sigmx is zero (matrix is already diagonal) */

    if (sigmx == 0.) {
	goto L70;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	D(i) = (d__1 = D(i), ABS(d__1));
/* Computing MAX */
	d__1 = sigmx, d__2 = D(i);
	sigmx = MAX(d__1,d__2);
/* L20: */
    }

/*     Get machine parameters */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("EPSILON");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("EPSILON");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("EPSILON");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    sfmin = dlamch_("SAFE MINIMUM");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sfmin = qlamch("SAFE MINIMUM");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sfmin = qlamch_("SAFE MINIMUM");
#endif


/*     Compute singular values to relative accuracy TOL   
       It is assumed that tol**2 does not underflow.   

   Computing MAX   
   Computing MIN */
    d__3 = 100., d__4 = pow(eps, c_b8);
    d__1 = 10., d__2 = MIN(d__3,d__4);
    tolmul = MAX(d__1,d__2);
    tol = tolmul * eps;
/* Computing 2nd power */
    d__1 = tol;
    tol2 = d__1 * d__1;

    thresh = sigmx * sqrt(sfmin) * tol;

/*     Scale matrix so the square of the largest element is   
       1 / ( 256 * SFMIN ) */

    scl = sqrt(1. / (sfmin * 256.));
/* Computing 2nd power */
    d__1 = tolmul;
    small2 = 1. / (d__1 * d__1 * 256.);

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(n, &D(1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(n, &D(1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(n, &D(1), &c__1, &WORK(1), &c__1);
#endif

    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&i__1, &E(1), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&i__1, &E(1), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&i__1, &E(1), &c__1, &WORK(*n + 1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dlascl_("G", &c__0, &c__0, &sigmx, &scl, n, &c__1, &WORK(1), n, &ierr)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlascl("G", &c__0, &c__0, &sigmx, &scl, n, &c__1, &WORK(1), n, &ierr)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlascl_("G", &c__0, &c__0, &sigmx, &scl, n, &c__1, &WORK(1), n, &ierr)
#endif

	    ;
    i__1 = *n - 1;
    i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
    dlascl_("G", &c__0, &c__0, &sigmx, &scl, &i__1, &c__1, &WORK(*n + 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlascl("G", &c__0, &c__0, &sigmx, &scl, &i__1, &c__1, &WORK(*n + 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlascl_("G", &c__0, &c__0, &sigmx, &scl, &i__1, &c__1, &WORK(*n + 1), &
#endif

	    i__2, &ierr);

/*     Square D and E (the input for the qd algorithm) */

    i__1 = (*n << 1) - 1;
    for (j = 1; j <= (*n<<1)-1; ++j) {
/* Computing 2nd power */
	d__1 = WORK(j);
	WORK(j) = d__1 * d__1;
/* L30: */
    }

/*     Apply qd algorithm */

    m = 0;
    E(*n) = 0.;
    dx = WORK(1);
    dm = dx;
    ke = 0;
    restrt = 0;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if ((d__1 = E(i), ABS(d__1)) <= thresh || WORK(*n + i) <= tol2 * (dm /
		 (LONG DOUBLE) (i - m))) {
	    ny = i - m;
	    if (ny == 1) {
		goto L50;
	    } else if (ny == 2) {

#ifdef PETSC_PREFIX_SUFFIX
		dlas2_(&D(m + 1), &E(m + 1), &D(m + 2), &sig1, &sig2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlas2(&D(m + 1), &E(m + 1), &D(m + 2), &sig1, &sig2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlas2_(&D(m + 1), &E(m + 1), &D(m + 2), &sig1, &sig2);
#endif

		D(m + 1) = sig1;
		D(m + 2) = sig2;
	    } else {
		kend = ke + 1 - m;

#ifdef PETSC_PREFIX_SUFFIX
		dlasq2_(&ny, &D(m + 1), &E(m + 1), &WORK(m + 1), &WORK(m + *n 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasq2(&ny, &D(m + 1), &E(m + 1), &WORK(m + 1), &WORK(m + *n 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasq2_(&ny, &D(m + 1), &E(m + 1), &WORK(m + 1), &WORK(m + *n 
#endif

			+ 1), &eps, &tol2, &small2, &dm, &kend, info);

/*                 Return, INFO = number of unconverged superd
iagonals */

		if (*info != 0) {
		    *info += i;
		    return;
		}

/*                 Undo scaling */

		i__2 = m + ny;
		for (j = m + 1; j <= m+ny; ++j) {
		    D(j) = sqrt(D(j));
/* L40: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		dlascl_("G", &c__0, &c__0, &scl, &sigmx, &ny, &c__1, &D(m + 1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlascl("G", &c__0, &c__0, &scl, &sigmx, &ny, &c__1, &D(m + 1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlascl_("G", &c__0, &c__0, &scl, &sigmx, &ny, &c__1, &D(m + 1)
#endif

			, &ny, &ierr);
	    }
L50:
	    m = i;
	    if (i != *n) {
		dx = WORK(i + 1);
		dm = dx;
		ke = i;
		restrt = 1;
	    }
	}
	if (i != *n && ! restrt) {
	    dx = WORK(i + 1) * (dx / (dx + WORK(*n + i)));
	    if (dm > dx) {
		dm = dx;
		ke = i;
	    }
	}
	restrt = 0;
/* L60: */
    }
    kend = ke + 1;

/*     Sort the singular values into decreasing order */

L70:

#ifdef PETSC_PREFIX_SUFFIX
    dlasrt_("D", n, &D(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlasrt("D", n, &D(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlasrt_("D", n, &D(1), info);
#endif

    return;

/*     End of DLASQ1 */

} /* dlasq1_ */

