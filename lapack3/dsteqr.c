#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsteqr_(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsteqr(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsteqr_(char *compz, int *n, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int 
	*info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the implicit QL or QR method.   
    The eigenvectors of a full or band symmetric matrix can also be found 
  
    if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to 
  
    tridiagonal form.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvalues and eigenvectors of the original 
  
                    symmetric matrix.  On entry, Z must contain the   
                    orthogonal matrix used to reduce the original matrix 
  
                    to tridiagonal form.   
            = 'I':  Compute eigenvalues and eigenvectors of the   
                    tridiagonal matrix.  Z is initialized to the identity 
  
                    matrix.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            On entry, if  COMPZ = 'V', then Z contains the orthogonal   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original symmetric matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            eigenvectors are desired, then  LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (MAX(1,2*N-2)) 
  
            If COMPZ = 'N', then WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm has failed to find all the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero; on exit, D   
                  and E contain the elements of a symmetric tridiagonal   
                  matrix which is orthogonally similar to the original   
                  matrix.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b9 = 0.;
    static LONG DOUBLE c_b10 = 1.;
    static int c__0 = 0;
    static int c__1 = 1;
    static int c__2 = 2;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
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
    static LONG DOUBLE b, c, f, g;
    static int i, j, k, l, m;
    static LONG DOUBLE p, r, s;
    extern long int lsame_(char *, char *);

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
    static LONG DOUBLE anorm;

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
    static int l1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaev2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaev2(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaev2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *);
    static int lendm1, lendp1;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static int ii;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int mm, iscale;

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

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dlaset_(char *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qlaset(char *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qlaset_(char *, int *, int 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE safmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE safmax;
    extern /* Subroutine */ void xerbla_(char *, int *);

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
    static int nmaxit, icompz;
    static LONG DOUBLE ssfmax;
    static int lm1, mm1, nm1;
    static LONG DOUBLE rt1, rt2, eps;
    static int lsv;
    static LONG DOUBLE tst, eps2;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < MAX(1,*n))) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEQR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    Z(1,1) = 1.;
	}
	return;
    }

/*     Determine the unit roundoff and over/underflow thresholds. */


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

/*     Compute the eigenvalues and eigenvectors of the tridiagonal   
       matrix. */

    if (icompz == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b9, &c_b10, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b9, &c_b10, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b9, &c_b10, &Z(1,1), ldz);
#endif

    }

    nmaxit = *n * 30;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
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
            d__1 = D(m); d__2 = D(m + 1);
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
    if (anorm == 0.) {
	goto L10;
    }
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

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), ABS(d__1)) < (d__2 = D(l), ABS(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m), ABS(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), ABS(d__1)) * (d__2 = D(m + 1),
			 ABS(d__2)) + safmin) {
		    goto L60;
		}
/* L50: */
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L80;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l + 1) {
	    if (icompz > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		dlaev2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaev2(&D(l), &E(l), &D(l + 1), &rt1, &rt2, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaev2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2, &c, &s);
#endif

		WORK(l) = c;
		WORK(*n - 1 + l) = s;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "B", n, &c__2, &WORK(l), &WORK(*n - 1 + l), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "B", n, &c__2, &WORK(l), &WORK(*n - 1 + l), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "B", n, &c__2, &WORK(l), &WORK(*n - 1 + l), &
#endif

			Z(1,l), ldz);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dlae2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlae2(&D(l), &E(l), &D(l + 1), &rt1, &rt2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlae2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2);
#endif

	    }
	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l + 1) - p) / (E(l) * 2.);

#ifdef PETSC_PREFIX_SUFFIX
	r = dlapy2_(&g, &c_b10);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	r = qlapy2(&g, &c_b10);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	r = qlapy2_(&g, &c_b10);
#endif

	g = D(m) - p + E(l) / (g + SIGN(r, g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    f = s * E(i);
	    b = c * E(i);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&g, &f, &c, &s, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&g, &f, &c, &s, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&g, &f, &c, &s, &r);
#endif

	    if (i != m - 1) {
		E(i + 1) = r;
	    }
	    g = D(i + 1) - p;
	    r = (D(i) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i + 1) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = -s;
	    }

/* L70: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = m - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlasr_("R", "V", "B", n, &mm, &WORK(l), &WORK(*n - 1 + l), &Z(1,l), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasr("R", "V", "B", n, &mm, &WORK(l), &WORK(*n - 1 + l), &Z(1,l), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasr_("R", "V", "B", n, &mm, &WORK(l), &WORK(*n - 1 + l), &Z(1,l), ldz);
#endif

	}

	D(l) -= p;
	E(l) = g;
	goto L40;

/*        Eigenvalue found. */

L80:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m - 1), ABS(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), ABS(d__1)) * (d__2 = D(m - 1),
			 ABS(d__2)) + safmin) {
		    goto L110;
		}
/* L100: */
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L130;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l - 1) {
	    if (icompz > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		dlaev2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaev2(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaev2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2, &c, &s);
#endif

		WORK(m) = c;
		WORK(*n - 1 + m) = s;

#ifdef PETSC_PREFIX_SUFFIX
		dlasr_("R", "V", "F", n, &c__2, &WORK(m), &WORK(*n - 1 + m), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlasr("R", "V", "F", n, &c__2, &WORK(m), &WORK(*n - 1 + m), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlasr_("R", "V", "F", n, &c__2, &WORK(m), &WORK(*n - 1 + m), &
#endif

			Z(1,l-1), ldz);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dlae2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlae2(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlae2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2);
#endif

	    }
	    D(l - 1) = rt1;
	    D(l) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l - 1) - p) / (E(l - 1) * 2.);

#ifdef PETSC_PREFIX_SUFFIX
	r = dlapy2_(&g, &c_b10);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	r = qlapy2(&g, &c_b10);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	r = qlapy2_(&g, &c_b10);
#endif

	g = D(m) - p + E(l - 1) / (g + SIGN(r, g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    f = s * E(i);
	    b = c * E(i);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&g, &f, &c, &s, &r);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&g, &f, &c, &s, &r);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&g, &f, &c, &s, &r);
#endif

	    if (i != m) {
		E(i - 1) = r;
	    }
	    g = D(i) - p;
	    r = (D(i + 1) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = s;
	    }

/* L120: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = l - m + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlasr_("R", "V", "F", n, &mm, &WORK(m), &WORK(*n - 1 + m), &Z(1,m), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasr("R", "V", "F", n, &mm, &WORK(m), &WORK(*n - 1 + m), &Z(1,m), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasr_("R", "V", "F", n, &mm, &WORK(m), &WORK(*n - 1 + m), &Z(1,m), ldz);
#endif

	}

	D(l) -= p;
	E(lm1) = g;
	goto L90;

/*        Eigenvalue found. */

L130:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

/*     Undo scaling if necessary */

L140:
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
	i__1 = lendsv - lsv;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif

		info);
    } else if (iscale == 2) {
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
	i__1 = lendsv - lsv;

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &E(lsv), n, 
#endif

		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	if (E(i) != 0.) {
	    ++(*info);
	}
/* L150: */
    }
    goto L190;

/*     Order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0) {

/*        Use Quick Sort */


#ifdef PETSC_PREFIX_SUFFIX
	dlasrt_("I", n, &D(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlasrt("I", n, &D(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlasrt_("I", n, &D(1), info);
#endif


    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= *n; ++ii) {
	    i = ii - 1;
	    k = i;
	    p = D(i);
	    i__2 = *n;
	    for (j = ii; j <= *n; ++j) {
		if (D(j) < p) {
		    k = j;
		    p = D(j);
		}
/* L170: */
	    }
	    if (k != i) {
		D(k) = D(i);
		D(i) = p;

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(n, &Z(1,i), &c__1, &Z(1,k), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(n, &Z(1,i), &c__1, &Z(1,k), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(n, &Z(1,i), &c__1, &Z(1,k), &
#endif

			c__1);
	    }
/* L180: */
	}
    }

L190:
    return;

/*     End of DSTEQR */

} /* dsteqr_ */

