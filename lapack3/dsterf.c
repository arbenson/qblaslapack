#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsterf_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsterf(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsterf_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif

	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   
    using the Pal-Walker-Kahan variant of the QL or QR algorithm.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix. 
  
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm failed to find all of the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__0 = 0;
    static int c__1 = 1;
    static LONG DOUBLE c_b32 = 1.;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE oldc;
    static int lend, jtot;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlae2_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlae2(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlae2_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE c;
    static int i, l, m;
    static LONG DOUBLE p, P_gamma, r, s, alpha, sigma, anorm;
    static int l1, lendm1, lendp1;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE bb;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int iscale;

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
    static LONG DOUBLE oldgam, safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE safmax;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasrt_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasrt(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasrt_(char *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static int lendsv;
    static LONG DOUBLE ssfmin;
    static int nmaxit;
    static LONG DOUBLE ssfmax;
    static int lm1, mm1, nm1;
    static LONG DOUBLE rt1, rt2, eps, rte;
    static int lsv;
    static LONG DOUBLE tst, eps2;



#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;

/*     Quick return if possible */

    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DSTERF", &i__1);
	return;
    }
    if (*n <= 1) {
	return;
    }

/*     Determine the unit roundoff for this environment. */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("E");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("E");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("E");
#endif

/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;

#ifdef PETSC_PREFIX_SUFFIX
    safmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    safmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    safmin = qlamch_("S");
#endif

    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues of the tridiagonal matrix. */

    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L170;
    }
    if (l1 > 1) {
	E(l1 - 1) = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= nm1; ++m) {
	    tst = (d__1 = E(m), ABS(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    d__1 = D(m);d__2 = D(m + 1);
            if (tst <= sqrt(( ABS(d__1))) * sqrt((
		     ABS(d__2))) * eps) {
		E(m) = 0.;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlanst_("I", &i__1, &D(l), &E(l));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlanst("I", &i__1, &D(l), &E(l));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlanst_("I", &i__1, &D(l), &E(l));
#endif

    iscale = 0;
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
#endif

		info);
	i__1 = lend - l;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
#endif

		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
#endif

		info);
	i__1 = lend - l;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
#endif

		info);
    }

    i__1 = lend - 1;
    for (i = l; i <= lend-1; ++i) {
/* Computing 2nd power */
	d__1 = E(i);
	E(i) = d__1 * d__1;
/* L40: */
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), ABS(d__1)) < (d__2 = D(l), ABS(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend >= l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L50:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
		tst = (d__1 = E(m), ABS(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m + 1), ABS(d__1))) {
		    goto L70;
		}
/* L60: */
	    }
	}

	m = lend;

L70:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L90;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l + 1) {
	    rte = sqrt(E(l));

#ifdef PETSC_PREFIX_SUFFIX
	    dlae2_(&D(l), &rte, &D(l + 1), &rt1, &rt2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlae2(&D(l), &rte, &D(l + 1), &rt1, &rt2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlae2_(&D(l), &rte, &D(l + 1), &rt1, &rt2);
#endif

	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L50;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l));
	sigma = (D(l + 1) - p) / (rte * 2.);

#ifdef PETSC_PREFIX_SUFFIX
	r = dlapy2_(&sigma, &c_b32);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	r = qlapy2(&sigma, &c_b32);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	r = qlapy2_(&sigma, &c_b32);
#endif

	sigma = p - rte / (sigma + SIGN(r, sigma));

	c = 1.;
	s = 0.;
	P_gamma = D(m) - sigma;
	p = P_gamma * P_gamma;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m - 1) {
		E(i + 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = P_gamma;
	    alpha = D(i);
	    P_gamma = c * (alpha - sigma) - s * oldgam;
	    D(i + 1) = oldgam + (alpha - P_gamma);
	    if (c != 0.) {
		p = P_gamma * P_gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L80: */
	}

	E(l) = s * p;
	D(l) = sigma + P_gamma;
	goto L50;

/*        Eigenvalue found. */

L90:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L50;
	}
	goto L150;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L100:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
		tst = (d__1 = E(m - 1), ABS(d__1));
		if (tst <= eps2 * (d__1 = D(m) * D(m - 1), ABS(d__1))) {
		    goto L120;
		}
/* L110: */
	    }
	}

	m = lend;

L120:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L140;
	}

/*        If remaining matrix is 2 by 2, use DLAE2 to compute its   
          eigenvalues. */

	if (m == l - 1) {
	    rte = sqrt(E(l - 1));

#ifdef PETSC_PREFIX_SUFFIX
	    dlae2_(&D(l), &rte, &D(l - 1), &rt1, &rt2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlae2(&D(l), &rte, &D(l - 1), &rt1, &rt2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlae2_(&D(l), &rte, &D(l - 1), &rt1, &rt2);
#endif

	    D(l) = rt1;
	    D(l - 1) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L100;
	    }
	    goto L150;
	}

	if (jtot == nmaxit) {
	    goto L150;
	}
	++jtot;

/*        Form shift. */

	rte = sqrt(E(l - 1));
	sigma = (D(l - 1) - p) / (rte * 2.);

#ifdef PETSC_PREFIX_SUFFIX
	r = dlapy2_(&sigma, &c_b32);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	r = qlapy2(&sigma, &c_b32);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	r = qlapy2_(&sigma, &c_b32);
#endif

	sigma = p - rte / (sigma + SIGN(r, sigma));

	c = 1.;
	s = 0.;
	P_gamma = D(m) - sigma;
	p = P_gamma * P_gamma;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    bb = E(i);
	    r = p + bb;
	    if (i != m) {
		E(i - 1) = s * r;
	    }
	    oldc = c;
	    c = p / r;
	    s = bb / r;
	    oldgam = P_gamma;
	    alpha = D(i + 1);
	    P_gamma = c * (alpha - sigma) - s * oldgam;
	    D(i) = oldgam + (alpha - P_gamma);
	    if (c != 0.) {
		p = P_gamma * P_gamma / c;
	    } else {
		p = oldc * bb;
	    }
/* L130: */
	}

	E(lm1) = s * p;
	D(l) = sigma + P_gamma;
	goto L100;

/*        Eigenvalue found. */

L140:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L100;
	}
	goto L150;

    }

/*     Undo scaling if necessary */

L150:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif

		info);
    }
    if (iscale == 2) {
	i__1 = lendsv - lsv + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
#endif

		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot == nmaxit) {
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    if (E(i) != 0.) {
		++(*info);
	    }
/* L160: */
	}
	return;
    }
    goto L10;

/*     Sort eigenvalues in increasing order. */

L170:

#ifdef PETSC_PREFIX_SUFFIX
    dlasrt_("I", n, &D(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlasrt("I", n, &D(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlasrt_("I", n, &D(1), info);
#endif


    return;

/*     End of DSTERF */

} /* dsterf_ */

