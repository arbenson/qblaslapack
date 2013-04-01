#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlartg_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *cs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlartg(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *cs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlartg_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *cs, 
#endif

	LONG DOUBLE *sn, LONG DOUBLE *r)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine DROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in DBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) LONG DOUBLE PRECISION   
            The first component of vector to be rotated.   

    G       (input) LONG DOUBLE PRECISION   
            The second component of vector to be rotated.   

    CS      (output) LONG DOUBLE PRECISION   
            The cosine of the rotation.   

    SN      (output) LONG DOUBLE PRECISION   
            The sine of the rotation.   

    R       (output) LONG DOUBLE PRECISION   
            The nonzero component of the rotated vector.   

    ===================================================================== 
*/
    /* Initialized data */
    static long int first = 1;
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static int i;
    static LONG DOUBLE scale;
    static int count;
    static LONG DOUBLE f1, g1, safmn2, safmx2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE safmin, eps;



    if (first) {
	first = 0;

#ifdef PETSC_PREFIX_SUFFIX
	safmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	safmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	safmin = qlamch_("S");
#endif


#ifdef PETSC_PREFIX_SUFFIX
	eps = dlamch_("E");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	eps = qlamch("E");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	eps = qlamch_("E");
#endif


#ifdef PETSC_PREFIX_SUFFIX
	d__1 = dlamch_("B");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	d__1 = qlamch("B");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	d__1 = qlamch_("B");
#endif


#ifdef PETSC_PREFIX_SUFFIX
	i__1 = (int) (log(safmin / eps) / log(dlamch_("B")) / 2.);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	i__1 = (int) (log(safmin / eps) / log(qlamch("B")) / 2.);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	i__1 = (int) (log(safmin / eps) / log(qlamch_("B")) / 2.);
#endif

	safmn2 = pow(d__1, (LONG DOUBLE) i__1);
	safmx2 = 1. / safmn2;
    }
    if (*g == 0.) {
	*cs = 1.;
	*sn = 0.;
	*r = *f;
    } else if (*f == 0.) {
	*cs = 0.;
	*sn = 1.;
	*r = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	d__1 = ABS(f1), d__2 = ABS(g1);
	scale = MAX(d__1,d__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    d__1 = ABS(f1), d__2 = ABS(g1);
	    scale = MAX(d__1,d__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    d__1 = ABS(f1), d__2 = ABS(g1);
	    scale = MAX(d__1,d__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	    i__1 = count;
	    for (i = 1; i <= count; ++i) {
		*r *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r;
	    *sn = g1 / *r;
	}
	if (ABS(*f) > ABS(*g) && *cs < 0.) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r = -(*r);
	}
    }
    return;

/*     End of DLARTG */

} /* dlartg_ */

