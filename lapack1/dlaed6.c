#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed6_(int *kniter, long int *orgati, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed6(int *kniter, long int *orgati, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed6_(int *kniter, long int *orgati, LONG DOUBLE *
#endif

	rho, LONG DOUBLE *d, LONG DOUBLE *z, LONG DOUBLE *finit, LONG DOUBLE *tau,
	 int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED6 computes the positive or negative root (closest to the origin) 
  
    of   
                     z(1)        z(2)        z(3)   
    f(x) =   rho + --------- + ---------- + ---------   
                    d(1)-x      d(2)-x      d(3)-x   

    It is assumed that   

          if ORGATI = .true. the root is between d(2) and d(3);   
          otherwise it is between d(1) and d(2)   

    This routine will be called by DLAED4 when necessary. In most cases, 
  
    the root sought is the smallest in magnitude, though it might not be 
  
    in some extremely rare situations.   

    Arguments   
    =========   

    KNITER       (input) INTEGER   
                 Refer to DLAED4 for its significance.   

    ORGATI       (input) LOGICAL   
                 If ORGATI is true, the needed root is between d(2) and   
                 d(3); otherwise it is between d(1) and d(2).  See   
                 DLAED4 for further details.   

    RHO          (input) LONG DOUBLE PRECISION   
                 Refer to the equation f(x) above.   

    D            (input) LONG DOUBLE PRECISION array, dimension (3)   
                 D satisfies d(1) < d(2) < d(3).   

    Z            (input) LONG DOUBLE PRECISION array, dimension (3)   
                 Each of the elements in z must be positive.   

    FINIT        (input) LONG DOUBLE PRECISION   
                 The value of f at 0. It is more accurate than the one   
                 evaluated inside this routine (if someone wants to do   
                 so).   

    TAU          (output) LONG DOUBLE PRECISION   
                 The root of the equation f(x).   

    INFO         (output) INTEGER   
                 = 0: successful exit   
                 > 0: if INFO = 1, failure to converge   

    ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Initialized data */
    static long int first = 1;
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
   /* Local variables */
    static LONG DOUBLE base;
    static int iter;
    static LONG DOUBLE temp, temp1, temp2, temp3, a, b, c, f;
    static int i;
    static long int scale;
    static int niter;
    static LONG DOUBLE small1, small2, fc, df, sminv1, sminv2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE dscale[3], sclfac, zscale[3], erretm, sclinv, pretau, 
	    ddf, eta, eps;


#define DSCALE(I) dscale[(I)]
#define ZSCALE(I) zscale[(I)]
#define Z(I) z[(I)-1]
#define D(I) d[(I)-1]



    *info = 0;

    niter = 1;
    *tau = 0.;
    if (*kniter == 2) {
	if (*orgati) {
	    temp = (D(3) - D(2)) / 2.;
	    c = *rho + Z(1) / (D(1) - D(2) - temp);
	    a = c * (D(2) + D(3)) + Z(2) + Z(3);
	    b = c * D(2) * D(3) + Z(2) * D(3) + Z(3) * D(2);
	} else {
	    temp = (D(1) - D(2)) / 2.;
	    c = *rho + Z(3) / (D(3) - D(2) - temp);
	    a = c * (D(1) + D(2)) + Z(1) + Z(2);
	    b = c * D(1) * D(2) + Z(1) * D(2) + Z(2) * D(1);
	}
	if (c == 0.) {
	    *tau = b / a;
	} else if (a <= 0.) {
	    d__1 = a * a - b * 4. * c;
            *tau = (a - sqrt(( ABS(d__1)))) / (c * 
		    2.);
	} else {
	    d__1 = a * a - b * 4. * c;
            *tau = b * 2. / (a + sqrt(( ABS(d__1))))
		    ;
	}
	temp = *rho + Z(1) / (D(1) - *tau) + Z(2) / (D(2) - *tau) + Z(3) / (D(
		3) - *tau);
	if (ABS(*finit) <= ABS(temp)) {
	    *tau = 0.;
	}
    }

/*     On first call to routine, get machine parameters for   
       possible scaling to avoid overflow */

    if (first) {

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
	base = dlamch_("Base");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	base = qlamch("Base");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	base = qlamch_("Base");
#endif


#ifdef PETSC_PREFIX_SUFFIX
	i__1 = (int) (log(dlamch_("SafMin")) / log(base) / 3.);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	i__1 = (int) (log(qlamch("SafMin")) / log(base) / 3.);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	i__1 = (int) (log(qlamch_("SafMin")) / log(base) / 3.);
#endif

	small1 = pow(base, (LONG DOUBLE)i__1);
	sminv1 = 1. / small1;
	small2 = small1 * small1;
	sminv2 = sminv1 * sminv1;
	first = 0;
    }

/*     Determine if scaling of inputs necessary to avoid overflow   
       when computing 1/TEMP**3 */

    if (*orgati) {
/* Computing MIN */
	d__3 = (d__1 = D(2) - *tau, ABS(d__1)), d__4 = (d__2 = D(3) - *tau, 
		ABS(d__2));
	temp = MIN(d__3,d__4);
    } else {
/* Computing MIN */
	d__3 = (d__1 = D(1) - *tau, ABS(d__1)), d__4 = (d__2 = D(2) - *tau, 
		ABS(d__2));
	temp = MIN(d__3,d__4);
    }
    scale = 0;
    if (temp <= small1) {
	scale = 1;
	if (temp <= small2) {

/*        Scale up by power of radix nearest 1/SAFMIN**(2/3) */

	    sclfac = sminv2;
	    sclinv = small2;
	} else {

/*        Scale up by power of radix nearest 1/SAFMIN**(1/3) */

	    sclfac = sminv1;
	    sclinv = small1;
	}

/*        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1
) */

	for (i = 1; i <= 3; ++i) {
	    DSCALE(i - 1) = D(i) * sclfac;
	    ZSCALE(i - 1) = Z(i) * sclfac;
/* L10: */
	}
	*tau *= sclfac;
    } else {

/*        Copy D and Z to DSCALE and ZSCALE */

	for (i = 1; i <= 3; ++i) {
	    DSCALE(i - 1) = D(i);
	    ZSCALE(i - 1) = Z(i);
/* L20: */
	}
    }

    fc = 0.;
    df = 0.;
    ddf = 0.;
    for (i = 1; i <= 3; ++i) {
	temp = 1. / (DSCALE(i - 1) - *tau);
	temp1 = ZSCALE(i - 1) * temp;
	temp2 = temp1 * temp;
	temp3 = temp2 * temp;
	fc += temp1 / DSCALE(i - 1);
	df += temp2;
	ddf += temp3;
/* L30: */
    }
    f = *finit + *tau * fc;

    if (ABS(f) <= 0.) {
	goto L60;
    }

/*        Iteration begins   

       It is not hard to see that   

             1) Iterations will go up monotonically   
                if FINIT < 0;   

             2) Iterations will go down monotonically   
                if FINIT > 0. */

    iter = niter + 1;

    for (niter = iter; niter <= 20; ++niter) {

	if (*orgati) {
	    temp1 = DSCALE(1) - *tau;
	    temp2 = DSCALE(2) - *tau;
	} else {
	    temp1 = DSCALE(0) - *tau;
	    temp2 = DSCALE(1) - *tau;
	}
	a = (temp1 + temp2) * f - temp1 * temp2 * df;
	b = temp1 * temp2 * f;
	c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
	if (c == 0.) {
	    eta = b / a;
	} else if (a <= 0.) {
	    d__1 = a * a - b * 4. * c;
            eta = (a - sqrt(( ABS(d__1)))) / (c * 
		    2.);
	} else {
	    d__1 = a * a - b * 4. * c;
            eta = b * 2. / (a + sqrt(( ABS(d__1))));
	}
	if (f * eta >= 0.) {
	    eta = -f / df;
	}

	temp = eta + *tau;
	if (*orgati) {
	    if (eta > 0. && temp >= DSCALE(2)) {
		eta = (DSCALE(2) - *tau) / 2.;
	    }
	    if (eta < 0. && temp <= DSCALE(1)) {
		eta = (DSCALE(1) - *tau) / 2.;
	    }
	} else {
	    if (eta > 0. && temp >= DSCALE(1)) {
		eta = (DSCALE(1) - *tau) / 2.;
	    }
	    if (eta < 0. && temp <= DSCALE(0)) {
		eta = (DSCALE(0) - *tau) / 2.;
	    }
	}
	pretau = *tau;
	*tau += eta;

	fc = 0.;
	erretm = ABS(*rho);
	df = 0.;
	ddf = 0.;
	for (i = 1; i <= 3; ++i) {
	    temp = 1. / (DSCALE(i - 1) - *tau);
	    temp1 = ZSCALE(i - 1) * temp;
	    temp2 = temp1 * temp;
	    temp3 = temp2 * temp;
	    fc += temp1 / (DSCALE(i - 1) - pretau);
	    erretm += ABS(temp1);
	    df += temp2;
	    ddf += temp3;
/* L40: */
	}
	f += eta * fc;
	erretm = erretm * 8. + ABS(*tau) * df;
	if (ABS(f) <= eps * erretm) {
	    goto L60;
	}
/* L50: */
    }
    *info = 1;
L60:

/*     Undo scaling */

    if (scale) {
	*tau *= sclinv;
    }
    return;

/*     End of DLAED6 */

} /* dlaed6_ */

