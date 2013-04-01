#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

typedef struct { LONG DOUBLE r, i; } doublecomplex;

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dcabs1_(doublecomplex *z)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qcabs1(doublecomplex *z)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qcabs1_(doublecomplex *z)
#endif

{
/* >>Start of File<<   

       System generated locals */
    LONG DOUBLE ret_val;
    static doublecomplex equiv_0[1];

    /* Local variables */
#define t ((LONG DOUBLE *)equiv_0)
#define zz (equiv_0)

    zz->r = z->r, zz->i = z->i;
    ret_val = ABS(t[0]) + ABS(t[1]);
    return ret_val;
} /* dcabs1_ */

#undef zz
#undef t


