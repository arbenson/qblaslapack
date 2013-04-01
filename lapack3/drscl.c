#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void drscl_(int *n, LONG DOUBLE *sa, LONG DOUBLE *sx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qrscl(int *n, LONG DOUBLE *sa, LONG DOUBLE *sx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qrscl_(int *n, LONG DOUBLE *sa, LONG DOUBLE *sx, 
#endif

	int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DRSCL multiplies an n-element real vector x by the real scalar 1/a.   
    This is done without overflow or underflow as long as   
    the final result x/a does not overflow or underflow.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of components of the vector x.   

    SA      (input) LONG DOUBLE PRECISION   
            The scalar a which is used to divide each component of x.   
            SA must be >= 0, or the subroutine will divide by zero.   

    SX      (input/output) LONG DOUBLE PRECISION array, dimension   
                           (1+(N-1)*ABS(INCX))   
            The n-element vector x.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector SX.   
            > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n 
  

   ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    static LONG DOUBLE cden;
    static long int done;
    static LONG DOUBLE cnum, cden1, cnum1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlabad(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE bignum, smlnum, mul;


#define SX(I) sx[(I)-1]


    if (*n <= 0) {
	return;
    }

/*     Get machine parameters */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S");
#endif

    bignum = 1. / smlnum;

#ifdef PETSC_PREFIX_SUFFIX
    dlabad_(&smlnum, &bignum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlabad(&smlnum, &bignum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlabad_(&smlnum, &bignum);
#endif


/*     Initialize the denominator to SA and the numerator to 1. */

    cden = *sa;
    cnum = 1.;

L10:
    cden1 = cden * smlnum;
    cnum1 = cnum / bignum;
    if (ABS(cden1) > ABS(cnum) && cnum != 0.) {

/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. 
*/

	mul = smlnum;
	done = 0;
	cden = cden1;
    } else if (ABS(cnum1) > ABS(cden)) {

/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. 
*/

	mul = bignum;
	done = 0;
	cnum = cnum1;
    } else {

/*        Multiply X by CNUM / CDEN and return. */

	mul = cnum / cden;
	done = 1;
    }

/*     Scale the vector X by MUL */


#ifdef PETSC_PREFIX_SUFFIX
    dscal_(n, &mul, &SX(1), incx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qscal(n, &mul, &SX(1), incx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qscal_(n, &mul, &SX(1), incx);
#endif


    if (! done) {
	goto L10;
    }

    return;

/*     End of DRSCL */

} /* drscl_ */

