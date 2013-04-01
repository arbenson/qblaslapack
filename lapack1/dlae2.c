#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlae2_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlae2(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlae2_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif

	LONG DOUBLE *rt1, LONG DOUBLE *rt2)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix   
       [  A   B  ]   
       [  B   C  ].   
    On return, RT1 is the eigenvalue of larger absolute value, and RT2   
    is the eigenvalue of smaller absolute value.   

    Arguments   
    =========   

    A       (input) LONG DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    B       (input) LONG DOUBLE PRECISION   
            The (1,2) and (2,1) elements of the 2-by-2 matrix.   

    C       (input) LONG DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    RT1     (output) LONG DOUBLE PRECISION   
            The eigenvalue of larger absolute value.   

    RT2     (output) LONG DOUBLE PRECISION   
            The eigenvalue of smaller absolute value.   

    Further Details   
    ===============   

    RT1 is accurate to a few ulps barring over/underflow.   

    RT2 may be inaccurate if there is massive cancellation in the   
    determinant A*C-B*B; higher precision or correctly rounded or   
    correctly truncated arithmetic would be needed to compute RT2   
    accurately in all cases.   

    Overflow is possible only if RT1 is within a factor of 5 of overflow. 
  
    Underflow is harmless if the input data is 0 or exceeds   
       underflow_threshold / macheps.   

   ===================================================================== 
  


       Compute the eigenvalues */
    /* System generated locals */
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE acmn, acmx, ab, df, tb, sm, rt, adf;


    sm = *a + *c;
    df = *a - *c;
    adf = ABS(df);
    tb = *b + *b;
    ab = ABS(tb);
    if (ABS(*a) > ABS(*c)) {
	acmx = *a;
	acmn = *c;
    } else {
	acmx = *c;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
/* Computing 2nd power */
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;

/*        Order of execution important.   
          To get fully accurate smaller eigenvalue,   
          next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5;
	*rt2 = rt * -.5;
    }
    return;

/*     End of DLAE2 */

} /* dlae2_ */

