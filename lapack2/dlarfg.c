#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlarfg_(int *n, LONG DOUBLE *alpha, LONG DOUBLE *x, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlarfg(int *n, LONG DOUBLE *alpha, LONG DOUBLE *x, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlarfg_(int *n, LONG DOUBLE *alpha, LONG DOUBLE *x, 
#endif

	int *incx, LONG DOUBLE *tau)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARFG generates a real elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element real   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a real scalar and v is a real (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be 
  
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) LONG DOUBLE PRECISION   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) LONG DOUBLE PRECISION array, dimension   
                           (1+(N-2)*ABS(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) LONG DOUBLE PRECISION   
            The value tau.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE beta;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

    static int j;

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
    static LONG DOUBLE xnorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif

    static LONG DOUBLE safmin, rsafmn;
    static int knt;


#define X(I) x[(I)-1]


    if (*n <= 1) {
	*tau = 0.;
	return;
    }

    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
    xnorm = dnrm2_(&i__1, &X(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    xnorm = qnrm2(&i__1, &X(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    xnorm = qnrm2_(&i__1, &X(1), incx);
#endif


    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */


#ifdef PETSC_PREFIX_SUFFIX
	d__1 = dlapy2_(alpha, &xnorm);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	d__1 = qlapy2(alpha, &xnorm);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	d__1 = qlapy2_(alpha, &xnorm);
#endif

	beta = -SIGN(d__1, *alpha);

#ifdef PETSC_PREFIX_SUFFIX
	safmin = dlamch_("S") / dlamch_("E");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	safmin = qlamch("S") / dlamch_("E");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	safmin = qlamch_("S") / dlamch_("E");
#endif

	if (ABS(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute 
them */

	    rsafmn = 1. / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__1, &rsafmn, &X(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__1, &rsafmn, &X(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__1, &rsafmn, &X(1), incx);
#endif

	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (ABS(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    xnorm = dnrm2_(&i__1, &X(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    xnorm = qnrm2(&i__1, &X(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    xnorm = qnrm2_(&i__1, &X(1), incx);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    d__1 = dlapy2_(alpha, &xnorm);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    d__1 = qlapy2(alpha, &xnorm);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    d__1 = qlapy2_(alpha, &xnorm);
#endif

	    beta = -SIGN(d__1, *alpha);
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__1, &d__1, &X(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__1, &d__1, &X(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__1, &d__1, &X(1), incx);
#endif


/*           If ALPHA is subnormal, it may lose relative accuracy 
*/

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= knt; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__1, &d__1, &X(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__1, &d__1, &X(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__1, &d__1, &X(1), incx);
#endif

	    *alpha = beta;
	}
    }

    return;

/*     End of DLARFG */

} /* dlarfg_ */

