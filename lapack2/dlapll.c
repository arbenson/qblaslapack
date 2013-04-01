#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlapll_(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlapll(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlapll_(int *n, LONG DOUBLE *x, int *incx, 
#endif

	LONG DOUBLE *y, int *incy, LONG DOUBLE *ssmin)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    Given two column vectors X and Y, let   

                         A = ( X Y ).   

    The subroutine first computes the QR factorization of A = Q*R,   
    and then computes the SVD of the 2-by-2 upper triangular matrix R.   
    The smaller singular value of R is returned in SSMIN, which is used   
    as the measurement of the linear dependency of the vectors X and Y.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The length of the vectors X and Y.   

    X       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCX)   
            On entry, X contains the N-vector X.   
            On exit, X is overwritten.   

    INCX    (input) INTEGER   
            The increment between successive elements of X. INCX > 0.   

    Y       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCY)   
            On entry, Y contains the N-vector Y.   
            On exit, Y is overwritten.   

    INCY    (input) INTEGER   
            The increment between successive elements of Y. INCY > 0.   

    SSMIN   (output) LONG DOUBLE PRECISION   
            The smallest singular value of the N-by-2 matrix A = ( X Y ). 
  

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);

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
    static LONG DOUBLE c;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static LONG DOUBLE ssmax, a11, a12, a22;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *);
    static LONG DOUBLE tau;


#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    if (*n <= 1) {
	*ssmin = 0.;
	return;
    }

/*     Compute the QR factorization of the N-by-2 matrix ( X Y ) */


#ifdef PETSC_PREFIX_SUFFIX
    dlarfg_(n, &X(1), &X(*incx + 1), incx, &tau);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlarfg(n, &X(1), &X(*incx + 1), incx, &tau);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlarfg_(n, &X(1), &X(*incx + 1), incx, &tau);
#endif

    a11 = X(1);
    X(1) = 1.;


#ifdef PETSC_PREFIX_SUFFIX
    c = -tau * ddot_(n, &X(1), incx, &Y(1), incy);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    c = -tau * qdot(n, &X(1), incx, &Y(1), incy);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    c = -tau * qdot_(n, &X(1), incx, &Y(1), incy);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(n, &c, &X(1), incx, &Y(1), incy);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(n, &c, &X(1), incx, &Y(1), incy);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(n, &c, &X(1), incx, &Y(1), incy);
#endif


    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
    dlarfg_(&i__1, &Y(*incy + 1), &Y((*incy << 1) + 1), incy, &tau);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlarfg(&i__1, &Y(*incy + 1), &Y((*incy << 1) + 1), incy, &tau);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlarfg_(&i__1, &Y(*incy + 1), &Y((*incy << 1) + 1), incy, &tau);
#endif


    a12 = Y(1);
    a22 = Y(*incy + 1);

/*     Compute the SVD of 2-by-2 Upper triangular matrix. */


#ifdef PETSC_PREFIX_SUFFIX
    dlas2_(&a11, &a12, &a22, ssmin, &ssmax);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlas2(&a11, &a12, &a22, ssmin, &ssmax);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlas2_(&a11, &a12, &a22, ssmin, &ssmax);
#endif


    return;

/*     End of DLAPLL */

} /* dlapll_ */

