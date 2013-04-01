#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlasq4_(int *n, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlasq4(int *n, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlasq4_(int *n, LONG DOUBLE *q, LONG DOUBLE *e, 
#endif

	LONG DOUBLE *tau, LONG DOUBLE *sup)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


       Purpose   
       =======   

       DLASQ4 estimates TAU, the smallest eigenvalue of a matrix. This   
       routine improves the input value of SUP which is an upper bound   
       for the smallest eigenvalue for this matrix .   

       Arguments   
       =========   

    N       (input) INTEGER   
            On entry, N specifies the number of rows and columns   
            in the matrix. N must be at least 0.   

    Q       (input) LONG DOUBLE PRECISION array, dimension (N)   
            Q array   

    E       (input) LONG DOUBLE PRECISION array, dimension (N)   
            E array   

    TAU     (output) LONG DOUBLE PRECISION   
            Estimate of the shift   

    SUP     (input/output) LONG DOUBLE PRECISION   
            Upper bound for the smallest singular value   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b4 = .7;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE xinf, d;
    static int i;
    static LONG DOUBLE dm;
    static int ifl;



#define E(I) e[(I)-1]
#define Q(I) q[(I)-1]


    ifl = 1;
/* Computing MIN */
    d__1 = MIN(*sup,Q(1)), d__1 = MIN(d__1,Q(2)), d__1 = MIN(d__1,Q(3)), d__2 
	    = Q(*n), d__1 = MIN(d__1,d__2), d__2 = Q(*n - 1), d__1 = MIN(d__1,
	    d__2), d__2 = Q(*n - 2);
    *sup = MIN(d__1,d__2);
    *tau = *sup * .9999;
    xinf = 0.;
L10:
    if (ifl == 5) {
	*tau = xinf;
	return;
    }
    d = Q(1) - *tau;
    dm = d;
    i__1 = *n - 2;
    for (i = 1; i <= *n-2; ++i) {
	d = d / (d + E(i)) * Q(i + 1) - *tau;
	if (dm > d) {
	    dm = d;
	}
	if (d < 0.) {
	    *sup = *tau;
/* Computing MAX */
	    d__1 = *sup * pow(c_b4, (LONG DOUBLE)ifl), d__2 = d + *tau;
	    *tau = MAX(d__1,d__2);
	    ++ifl;
	    goto L10;
	}
/* L20: */
    }
    d = d / (d + E(*n - 1)) * Q(*n) - *tau;
    if (dm > d) {
	dm = d;
    }
    if (d < 0.) {
	*sup = *tau;
/* Computing MAX */
	d__1 = xinf, d__2 = d + *tau;
	xinf = MAX(d__1,d__2);
	if (*sup * pow(c_b4, (LONG DOUBLE)ifl) <= xinf) {
	    *tau = xinf;
	} else {
	    *tau = *sup * pow(c_b4, (LONG DOUBLE)ifl);
	    ++ifl;
	    goto L10;
	}
    } else {
/* Computing MIN */
	d__1 = *sup, d__2 = dm + *tau;
	*sup = MIN(d__1,d__2);
    }
    return;

/*     End of DLASQ4 */

} /* dlasq4_ */

