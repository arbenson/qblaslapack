#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>
#define P_nint(x)      ((int)(x+0.5))



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlasq2_(int *m, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlasq2(int *m, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlasq2_(int *m, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif

	LONG DOUBLE *qq, LONG DOUBLE *ee, LONG DOUBLE *eps, LONG DOUBLE *tol2, 
	LONG DOUBLE *small2, LONG DOUBLE *sup, int *kend, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ2 computes the singular values of a real N-by-N unreduced   
       bidiagonal matrix with squared diagonal elements in Q and   
       squared off-diagonal elements in E. The singular values are   
       computed to relative accuracy TOL, barring over/underflow or   
       denormalization.   

       Arguments   
       =========   

    M       (input) INTEGER   
            The number of rows and columns in the matrix. M >= 0.   

    Q       (output) LONG DOUBLE PRECISION array, dimension (M)   
            On normal exit, contains the squared singular values.   

    E       (workspace) LONG DOUBLE PRECISION array, dimension (M)   

    QQ      (input/output) LONG DOUBLE PRECISION array, dimension (M)   
            On entry, QQ contains the squared diagonal elements of the   
            bidiagonal matrix whose SVD is desired.   
            On exit, QQ is overwritten.   

    EE      (input/output) LONG DOUBLE PRECISION array, dimension (M)   
            On entry, EE(1:N-1) contains the squared off-diagonal   
            elements of the bidiagonal matrix whose SVD is desired.   
            On exit, EE is overwritten.   

    EPS     (input) LONG DOUBLE PRECISION   
            Machine epsilon.   

    TOL2    (input) LONG DOUBLE PRECISION   
            Desired relative accuracy of computed eigenvalues   
            as defined in DLASQ1.   

    SMALL2  (input) LONG DOUBLE PRECISION   
            A threshold value as defined in DLASQ1.   

    SUP     (input/output) LONG DOUBLE PRECISION   
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
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE xinf;
    static int n;
    static LONG DOUBLE sigma, qemax;
    static int iconv;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasq3_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq3(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasq3_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *
	    , int *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *);
    static int iphase;
    static LONG DOUBLE xx, yy;
    static int off, isp, off1;


#define EE(I) ee[(I)-1]
#define QQ(I) qq[(I)-1]
#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    n = *m;

/*     Set the default maximum number of iterations */

    off = 0;
    off1 = off + 1;
    sigma = 0.;
    xinf = 0.;
    iconv = 0;
    iphase = 2;

/*     Try deflation at the bottom   

       1x1 deflation */

L10:
    if (n <= 2) {
	goto L20;
    }
/* Computing MAX */
    d__1 = QQ(n), d__1 = MAX(d__1,xinf);
    if (EE(n - 1) <= MAX(d__1,*small2) * *tol2) {
	Q(n) = QQ(n);
	--n;
	if (*kend > n) {
	    *kend = n;
	}
/* Computing MIN */
	d__1 = QQ(n), d__2 = QQ(n - 1);
	*sup = MIN(d__1,d__2);
	goto L10;
    }

/*     2x2 deflation   

   Computing MAX */
    d__1 = MAX(xinf,*small2), d__2 = QQ(n) / (QQ(n) + EE(n - 1) + QQ(n - 1)) *
	     QQ(n - 1);
    if (EE(n - 2) <= MAX(d__1,d__2) * *tol2) {
/* Computing MAX */
	d__1 = QQ(n), d__2 = QQ(n - 1), d__1 = MAX(d__1,d__2), d__2 = EE(n - 
		1);
	qemax = MAX(d__1,d__2);
	if (qemax != 0.) {
	    if (qemax == QQ(n - 1)) {
/* Computing 2nd power */
		d__1 = (QQ(n) - QQ(n - 1) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(d__1 * 
			d__1 + EE(n - 1) * 4. / qemax)) * .5;
	    } else if (qemax == QQ(n)) {
/* Computing 2nd power */
		d__1 = (QQ(n - 1) - QQ(n) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(d__1 * 
			d__1 + EE(n - 1) * 4. / qemax)) * .5;
	    } else {
/* Computing 2nd power */
		d__1 = (QQ(n) - QQ(n - 1) + EE(n - 1)) / qemax;
		xx = (QQ(n) + QQ(n - 1) + EE(n - 1) + qemax * sqrt(d__1 * 
			d__1 + QQ(n - 1) * 4. / qemax)) * .5;
	    }
/* Computing MAX */
	    d__1 = QQ(n), d__2 = QQ(n - 1);
/* Computing MIN */
	    d__3 = QQ(n), d__4 = QQ(n - 1);
	    yy = MAX(d__1,d__2) / xx * MIN(d__3,d__4);
	} else {
	    xx = 0.;
	    yy = 0.;
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

/*         No upper branch; return to DLASQ1 */

	    return;
	} else {

/*         Going back to upper branch */

	    xinf = 0.;
	    if (EE(off) > 0.) {
		isp = P_nint(EE(off));
		iphase = 1;
	    } else {
		isp = -P_nint(EE(off));
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
		d__1 = Q(n + off), d__2 = Q(n - 1 + off), d__1 = MIN(d__1,
			d__2), d__2 = Q(n - 2 + off);
		*sup = MIN(d__1,d__2);
	    } else {
/* Computing MIN */
		d__1 = QQ(n + off), d__2 = QQ(n - 1 + off), d__1 = MIN(d__1,
			d__2), d__2 = QQ(n - 2 + off);
		*sup = MIN(d__1,d__2);
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
	    d__1 = QQ(n + off), d__2 = QQ(n - 1 + off), d__1 = MAX(d__1,d__2),
		     d__2 = EE(n - 1 + off);
	    qemax = MAX(d__1,d__2);
	    if (qemax != 0.) {
		if (qemax == QQ(n - 1 + off)) {
/* Computing 2nd power */
		    d__1 = (QQ(n + off) - QQ(n - 1 + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + EE(off + n - 1) * 4. / 
			    qemax)) * .5;
		} else if (qemax == QQ(n + off)) {
/* Computing 2nd power */
		    d__1 = (QQ(n - 1 + off) - QQ(n + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + EE(n - 1 + off) * 4. / 
			    qemax)) * .5;
		} else {
/* Computing 2nd power */
		    d__1 = (QQ(n + off) - QQ(n - 1 + off) + EE(n - 1 + off)) /
			     qemax;
		    xx = (QQ(n + off) + QQ(n - 1 + off) + EE(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + QQ(n - 1 + off) * 4. / 
			    qemax)) * .5;
		}
/* Computing MAX */
		d__1 = QQ(n + off), d__2 = QQ(n - 1 + off);
/* Computing MIN */
		d__3 = QQ(n + off), d__4 = QQ(n - 1 + off);
		yy = MAX(d__1,d__2) / xx * MIN(d__3,d__4);
	    } else {
		xx = 0.;
		yy = 0.;
	    }
	} else {
/* Computing MAX */
	    d__1 = Q(n + off), d__2 = Q(n - 1 + off), d__1 = MAX(d__1,d__2), 
		    d__2 = E(n - 1 + off);
	    qemax = MAX(d__1,d__2);
	    if (qemax != 0.) {
		if (qemax == Q(n - 1 + off)) {
/* Computing 2nd power */
		    d__1 = (Q(n + off) - Q(n - 1 + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + E(n - 1 + off) * 4. / 
			    qemax)) * .5;
		} else if (qemax == Q(n + off)) {
/* Computing 2nd power */
		    d__1 = (Q(n - 1 + off) - Q(n + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + E(n - 1 + off) * 4. / 
			    qemax)) * .5;
		} else {
/* Computing 2nd power */
		    d__1 = (Q(n + off) - Q(n - 1 + off) + E(n - 1 + off)) / 
			    qemax;
		    xx = (Q(n + off) + Q(n - 1 + off) + E(n - 1 + off) + 
			    qemax * sqrt(d__1 * d__1 + Q(n - 1 + off) * 4. / 
			    qemax)) * .5;
		}
/* Computing MAX */
		d__1 = Q(n + off), d__2 = Q(n - 1 + off);
/* Computing MIN */
		d__3 = Q(n + off), d__4 = Q(n - 1 + off);
		yy = MAX(d__1,d__2) / xx * MIN(d__3,d__4);
	    } else {
		xx = 0.;
		yy = 0.;
	    }
	}
	Q(n - 1 + off) = sigma + xx;
	Q(n + off) = yy + sigma;
	n = 0;
	goto L20;
    }

#ifdef PETSC_PREFIX_SUFFIX
    dlasq3_(&n, &Q(off1), &E(off1), &QQ(off1), &EE(off1), sup, &sigma, kend, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlasq3(&n, &Q(off1), &E(off1), &QQ(off1), &EE(off1), sup, &sigma, kend, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlasq3_(&n, &Q(off1), &E(off1), &QQ(off1), &EE(off1), sup, &sigma, kend, &
#endif

	    off, &iphase, &iconv, eps, tol2, small2);
    if (*sup < 0.) {
	*info = n + off;
	return;
    }
    off1 = off + 1;
    goto L20;

/*     End of DLASQ2 */

} /* dlasq2_ */

