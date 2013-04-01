#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>

/* Table of constant values */



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void drotg_(LONG DOUBLE *da, LONG DOUBLE *db, LONG DOUBLE *c, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qrotg(LONG DOUBLE *da, LONG DOUBLE *db, LONG DOUBLE *c, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qrotg_(LONG DOUBLE *da, LONG DOUBLE *db, LONG DOUBLE *c, 
#endif

	LONG DOUBLE *s)
{


    /* System generated locals */
    LONG DOUBLE d__1, d__2;

    /* Builtin functions */

    /* Local variables */
    static LONG DOUBLE r, scale, z, roe;


/*     construct givens plane rotation.   
       jack dongarra, linpack, 3/11/78. */


    roe = *db;
    if (ABS(*da) > ABS(*db)) {
	roe = *da;
    }
    scale = ABS(*da) + ABS(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c = 1.;
    *s = 0.;
    r = 0.;
    z = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r = (( roe < 0 ) ? -1.0 : 1.0) * r;
    *c = *da / r;
    *s = *db / r;
    z = 1.;
    if (ABS(*da) > ABS(*db)) {
	z = *s;
    }
    if (ABS(*db) >= ABS(*da) && *c != 0.) {
	z = 1. / *c;
    }
L20:
    *da = r;
    *db = z;
    return;
} /* drotg_ */

