#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>
#define P_nint(x)      ((int)(x+0.5))


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlacon_(int *n, LONG DOUBLE *v, LONG DOUBLE *x, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlacon(int *n, LONG DOUBLE *v, LONG DOUBLE *x, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlacon_(int *n, LONG DOUBLE *v, LONG DOUBLE *x, 
#endif

	int *isgn, LONG DOUBLE *est, int *kase)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLACON estimates the 1-norm of a square, real matrix A.   
    Reverse communication is used for evaluating matrix-vector products. 
  

    Arguments   
    =========   

    N      (input) INTEGER   
           The order of the matrix.  N >= 1.   

    V      (workspace) LONG DOUBLE PRECISION array, dimension (N)   
           On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
           (W is not returned).   

    X      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
           On an intermediate return, X should be overwritten by   
                 A * X,   if KASE=1,   
                 A' * X,  if KASE=2,   
           and DLACON must be re-called with all the other parameters   
           unchanged.   

    ISGN   (workspace) INTEGER array, dimension (N)   

    EST    (output) LONG DOUBLE PRECISION   
           An estimate (a lower bound) for norm(A).   

    KASE   (input/output) INTEGER   
           On the initial call to DLACON, KASE should be 0.   
           On an intermediate return, KASE will be 1 or 2, indicating   
           whether X should be overwritten by A * X  or A' * X.   
           On the final return from DLACON, KASE will again be 0.   

    Further Details   
    ======= =======   

    Contributed by Nick Higham, University of Manchester.   
    Originally named SONEST, dated March 16, 1988.   

    Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of 
  
    a real or complex matrix, with applications to condition estimation", 
  
    ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b11 = 1.;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int iter;
    static LONG DOUBLE temp;
    static int jump, i, j;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dasum_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qasum(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qasum_(int *, LONG DOUBLE *, int *);
#endif

    static int jlast;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE altsgn, estold;



#define ISGN(I) isgn[(I)-1]
#define X(I) x[(I)-1]
#define V(I) v[(I)-1]


    if (*kase == 0) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    X(i) = 1. / (LONG DOUBLE) (*n);
/* L10: */
	}
	*kase = 1;
	jump = 1;
	return;
    }

    switch (jump) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

/*     ................ ENTRY   (JUMP = 1)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

L20:
    if (*n == 1) {
	V(1) = X(1);
	*est = ABS(V(1));
/*        ... QUIT */
	goto L150;
    }

#ifdef PETSC_PREFIX_SUFFIX
    *est = dasum_(n, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    *est = qasum(n, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    *est = qasum_(n, &X(1), &c__1);
#endif


    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	X(i) = SIGN(c_b11, X(i));
	ISGN(i) = P_nint(X(i));
/* L30: */
    }
    *kase = 2;
    jump = 2;
    return;

/*     ................ ENTRY   (JUMP = 2)   
       FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L40:

#ifdef PETSC_PREFIX_SUFFIX
    j = idamax_(n, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    j = iqamax(n, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    j = iqamax_(n, &X(1), &c__1);
#endif

    iter = 2;

/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

L50:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	X(i) = 0.;
/* L60: */
    }
    X(j) = 1.;
    *kase = 1;
    jump = 3;
    return;

/*     ................ ENTRY   (JUMP = 3)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L70:

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(n, &X(1), &c__1, &V(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(n, &X(1), &c__1, &V(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(n, &X(1), &c__1, &V(1), &c__1);
#endif

    estold = *est;

#ifdef PETSC_PREFIX_SUFFIX
    *est = dasum_(n, &V(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    *est = qasum(n, &V(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    *est = qasum_(n, &V(1), &c__1);
#endif

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	d__1 = SIGN(c_b11, X(i));
	if (P_nint(d__1) != ISGN(i)) {
	    goto L90;
	}
/* L80: */
    }
/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
/*     TEST FOR CYCLING. */
    if (*est <= estold) {
	goto L120;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	X(i) = SIGN(c_b11, X(i));
	ISGN(i) = P_nint(X(i));
/* L100: */
    }
    *kase = 2;
    jump = 4;
    return;

/*     ................ ENTRY   (JUMP = 4)   
       X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

L110:
    jlast = j;

#ifdef PETSC_PREFIX_SUFFIX
    j = idamax_(n, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    j = iqamax(n, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    j = iqamax_(n, &X(1), &c__1);
#endif

    if (X(jlast) != (d__1 = X(j), ABS(d__1)) && iter < 5) {
	++iter;
	goto L50;
    }

/*     ITERATION COMPLETE.  FINAL STAGE. */

L120:
    altsgn = 1.;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	X(i) = altsgn * ((LONG DOUBLE) (i - 1) / (LONG DOUBLE) (*n - 1) + 1.);
	altsgn = -altsgn;
/* L130: */
    }
    *kase = 1;
    jump = 5;
    return;

/*     ................ ENTRY   (JUMP = 5)   
       X HAS BEEN OVERWRITTEN BY A*X. */

L140:

#ifdef PETSC_PREFIX_SUFFIX
    temp = dasum_(n, &X(1), &c__1) / (LONG DOUBLE) (*n * 3) * 2.;
#endif
#ifdef Q_C_PREFIX_SUFFIX
    temp = qasum(n, &X(1), &c__1) / (LONG DOUBLE) (*n * 3) * 2.;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    temp = qasum_(n, &X(1), &c__1) / (LONG DOUBLE) (*n * 3) * 2.;
#endif

    if (temp > *est) {

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &X(1), &c__1, &V(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &X(1), &c__1, &V(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &X(1), &c__1, &V(1), &c__1);
#endif

	*est = temp;
    }

L150:
    *kase = 0;
    return;

/*     End of DLACON */

} /* dlacon_ */

