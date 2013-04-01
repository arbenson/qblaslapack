#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlas2_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlas2(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlas2_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif

	LONG DOUBLE *ssmin, LONG DOUBLE *ssmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAS2  computes the singular values of the 2-by-2 matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, SSMIN is the smaller singular value and SSMAX is the   
    larger singular value.   

    Arguments   
    =========   

    F       (input) LONG DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) LONG DOUBLE PRECISION   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) LONG DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) LONG DOUBLE PRECISION   
            The smaller singular value.   

    SSMAX   (output) LONG DOUBLE PRECISION   
            The larger singular value.   

    Further Details   
    ===============   

    Barring over/underflow, all output quantities are correct to within   
    a few units in the last place (ulps), even in the absence of a guard 
  
    digit in addition/subtraction.   

    In IEEE arithmetic, the code works correctly if one matrix element is 
  
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows, or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

    ==================================================================== 
*/
    /* System generated locals */
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */

    /* Local variables */
    static LONG DOUBLE fhmn, fhmx, c, fa, ga, ha, as, at, au;



    fa = ABS(*f);
    ga = ABS(*g);
    ha = ABS(*h);
    fhmn = MIN(fa,ha);
    fhmx = MAX(fa,ha);
    if (fhmn == 0.) {
	*ssmin = 0.;
	if (fhmx == 0.) {
	    *ssmax = ga;
	} else {
/* Computing 2nd power */
	    d__1 = MIN(fhmx,ga) / MAX(fhmx,ga);
	    *ssmax = MAX(fhmx,ga) * sqrt(d__1 * d__1 + 1.);
	}
    } else {
	if (ga < fhmx) {
	    as = fhmn / fhmx + 1.;
	    at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
	    d__1 = ga / fhmx;
	    au = d__1 * d__1;
	    c = 2. / (sqrt(as * as + au) + sqrt(at * at + au));
	    *ssmin = fhmn * c;
	    *ssmax = fhmx / c;
	} else {
	    au = fhmx / ga;
	    if (au == 0.) {

/*              Avoid possible harmful underflow if exponent r
ange   
                asymmetric (true SSMIN may not underflow even 
if   
                AU underflows) */

		*ssmin = fhmn * fhmx / ga;
		*ssmax = ga;
	    } else {
		as = fhmn / fhmx + 1.;
		at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
		d__1 = as * au;
/* Computing 2nd power */
		d__2 = at * au;
		c = 1. / (sqrt(d__1 * d__1 + 1.) + sqrt(d__2 * d__2 + 1.));
		*ssmin = fhmn * c * au;
		*ssmin += *ssmin;
		*ssmax = ga / (c + c);
	    }
	}
    }
    return;

/*     End of DLAS2 */

} /* dlas2_ */

