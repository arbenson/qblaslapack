#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dladiv_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qladiv(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qladiv_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif

	LONG DOUBLE *d, LONG DOUBLE *p, LONG DOUBLE *q)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLADIV performs complex division in  real arithmetic   

                          a + i*b   
               p + i*q = ---------   
                          c + i*d   

    The algorithm is due to Robert L. Smith and can be found   
    in D. Knuth, The art of Computer Programming, Vol.2, p.195   

    Arguments   
    =========   

    A       (input) LONG DOUBLE PRECISION   
    B       (input) LONG DOUBLE PRECISION   
    C       (input) LONG DOUBLE PRECISION   
    D       (input) LONG DOUBLE PRECISION   
            The scalars a, b, c, and d in the above expression.   

    P       (output) LONG DOUBLE PRECISION   
    Q       (output) LONG DOUBLE PRECISION   
            The scalars p and q in the above expression.   

    ===================================================================== 
*/
    static LONG DOUBLE e, f;



    if (ABS(*d) < ABS(*c)) {
	e = *d / *c;
	f = *c + *d * e;
	*p = (*a + *b * e) / f;
	*q = (*b - *a * e) / f;
    } else {
	e = *c / *d;
	f = *d + *c * e;
	*p = (*b + *a * e) / f;
	*q = (-(*a) + *b * e) / f;
    }

    return;

/*     End of DLADIV */

} /* dladiv_ */

