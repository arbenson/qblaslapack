#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlagts_(int *job, int *n, LONG DOUBLE *a, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlagts(int *job, int *n, LONG DOUBLE *a, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlagts_(int *job, int *n, LONG DOUBLE *a, 
#endif

	LONG DOUBLE *b, LONG DOUBLE *c, LONG DOUBLE *d, int *in, LONG DOUBLE *
	y, LONG DOUBLE *tol, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAGTS may be used to solve one of the systems of equations   

       (T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,   

    where T is an n by n tridiagonal matrix, for x, following the   
    factorization of (T - lambda*I) as   

       (T - lambda*I) = P*L*U ,   

    by routine DLAGTF. The choice of equation to be solved is   
    controlled by the argument JOB, and in each case there is an option   
    to perturb zero or very small diagonal elements of U, this option   
    being intended for use in applications such as inverse iteration.   

    Arguments   
    =========   

    JOB     (input) INTEGER   
            Specifies the job to be performed by DLAGTS as follows:   
            =  1: The equations  (T - lambda*I)x = y  are to be solved,   
                  but diagonal elements of U are not to be perturbed.   
            = -1: The equations  (T - lambda*I)x = y  are to be solved   
                  and, if overflow would otherwise occur, the diagonal   
                  elements of U are to be perturbed. See argument TOL   
                  below.   
            =  2: The equations  (T - lambda*I)'x = y  are to be solved, 
  
                  but diagonal elements of U are not to be perturbed.   
            = -2: The equations  (T - lambda*I)'x = y  are to be solved   
                  and, if overflow would otherwise occur, the diagonal   
                  elements of U are to be perturbed. See argument TOL   
                  below.   

    N       (input) INTEGER   
            The order of the matrix T.   

    A       (input) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, A must contain the diagonal elements of U as   
            returned from DLAGTF.   

    B       (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, B must contain the first super-diagonal elements of 
  
            U as returned from DLAGTF.   

    C       (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, C must contain the sub-diagonal elements of L as   
            returned from DLAGTF.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N-2)   
            On entry, D must contain the second super-diagonal elements   
            of U as returned from DLAGTF.   

    IN      (input) INTEGER array, dimension (N)   
            On entry, IN must contain details of the matrix P as returned 
  
            from DLAGTF.   

    Y       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the right hand side vector y.   
            On exit, Y is overwritten by the solution vector x.   

    TOL     (input/output) LONG DOUBLE PRECISION   
            On entry, with  JOB .lt. 0, TOL should be the minimum   
            perturbation to be made to very small diagonal elements of U. 
  
            TOL should normally be chosen as about eps*norm(U), where eps 
  
            is the relative machine precision, but if TOL is supplied as 
  
            non-positive, then it is reset to eps*MAX( ABS( u(i,j) ) ).   
            If  JOB .gt. 0  then TOL is not referenced.   

            On exit, TOL is changed as described above, only if TOL is   
            non-positive on entry. Otherwise TOL is unchanged.   

    INFO    (output) INTEGER   
            = 0   : successful exit   
            .lt. 0: if INFO = -i, the i-th argument had an illegal value 
  
            .gt. 0: overflow would occur when computing the INFO(th)   
                    element of the solution vector x. This can only occur 
  
                    when JOB is supplied as positive and either means   
                    that a diagonal element of U is very small, or that   
                    the elements of the right-hand side vector y are very 
  
                    large.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp, pert;
    static int k;
    static LONG DOUBLE absak, sfmin, ak;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum, eps;


#define Y(I) y[(I)-1]
#define IN(I) in[(I)-1]
#define D(I) d[(I)-1]
#define C(I) c[(I)-1]
#define B(I) b[(I)-1]
#define A(I) a[(I)-1]


    *info = 0;
    if (ABS(*job) > 2 || *job == 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAGTS", &i__1);
	return;
    }

    if (*n == 0) {
	return;
    }


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    sfmin = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sfmin = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sfmin = qlamch_("Safe minimum");
#endif

    bignum = 1. / sfmin;

    if (*job < 0) {
	if (*tol <= 0.) {
	    *tol = ABS(A(1));
	    if (*n > 1) {
/* Computing MAX */
		d__1 = *tol, d__2 = ABS(A(2)), d__1 = MAX(d__1,d__2), d__2 = 
			ABS(B(1));
		*tol = MAX(d__1,d__2);
	    }
	    i__1 = *n;
	    for (k = 3; k <= *n; ++k) {
/* Computing MAX */
		d__4 = *tol, d__5 = (d__1 = A(k), ABS(d__1)), d__4 = MAX(d__4,
			d__5), d__5 = (d__2 = B(k - 1), ABS(d__2)), d__4 = 
			MAX(d__4,d__5), d__5 = (d__3 = D(k - 2), ABS(d__3));
		*tol = MAX(d__4,d__5);
/* L10: */
	    }
	    *tol *= eps;
	    if (*tol == 0.) {
		*tol = eps;
	    }
	}
    }

    if (ABS(*job) == 1) {
	i__1 = *n;
	for (k = 2; k <= *n; ++k) {
	    if (IN(k - 1) == 0) {
		Y(k) -= C(k - 1) * Y(k - 1);
	    } else {
		temp = Y(k - 1);
		Y(k - 1) = Y(k);
		Y(k) = temp - C(k - 1) * Y(k);
	    }
/* L20: */
	}
	if (*job == 1) {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = Y(k) - B(k) * Y(k + 1) - D(k) * Y(k + 2);
		} else if (k == *n - 1) {
		    temp = Y(k) - B(k) * Y(k + 1);
		} else {
		    temp = Y(k);
		}
		ak = A(k);
		absak = ABS(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (absak == 0. || ABS(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (ABS(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		Y(k) = temp / ak;
/* L30: */
	    }
	} else {
	    for (k = *n; k >= 1; --k) {
		if (k <= *n - 2) {
		    temp = Y(k) - B(k) * Y(k + 1) - D(k) * Y(k + 2);
		} else if (k == *n - 1) {
		    temp = Y(k) - B(k) * Y(k + 1);
		} else {
		    temp = Y(k);
		}
		ak = A(k);
		pert = SIGN(*tol, ak);
L40:
		absak = ABS(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (absak == 0. || ABS(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L40;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (ABS(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L40;
		    }
		}
		Y(k) = temp / ak;
/* L50: */
	    }
	}
    } else {

/*        Come to here if  JOB = 2 or -2 */

	if (*job == 2) {
	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {
		if (k >= 3) {
		    temp = Y(k) - B(k - 1) * Y(k - 1) - D(k - 2) * Y(k - 2);
		} else if (k == 2) {
		    temp = Y(k) - B(k - 1) * Y(k - 1);
		} else {
		    temp = Y(k);
		}
		ak = A(k);
		absak = ABS(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (absak == 0. || ABS(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (ABS(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		Y(k) = temp / ak;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {
		if (k >= 3) {
		    temp = Y(k) - B(k - 1) * Y(k - 1) - D(k - 2) * Y(k - 2);
		} else if (k == 2) {
		    temp = Y(k) - B(k - 1) * Y(k - 1);
		} else {
		    temp = Y(k);
		}
		ak = A(k);
		pert = SIGN(*tol, ak);
L70:
		absak = ABS(ak);
		if (absak < 1.) {
		    if (absak < sfmin) {
			if (absak == 0. || ABS(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L70;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (ABS(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L70;
		    }
		}
		Y(k) = temp / ak;
/* L80: */
	    }
	}

	for (k = *n; k >= 2; --k) {
	    if (IN(k - 1) == 0) {
		Y(k - 1) -= C(k - 1) * Y(k);
	    } else {
		temp = Y(k - 1);
		Y(k - 1) = Y(k);
		Y(k) = temp - C(k - 1) * Y(k);
	    }
/* L90: */
	}
    }

/*     End of DLAGTS */

    return;
} /* dlagts_ */

