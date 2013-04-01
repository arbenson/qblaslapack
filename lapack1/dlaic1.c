#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaic1_(int *job, int *j, LONG DOUBLE *x, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaic1(int *job, int *j, LONG DOUBLE *x, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaic1_(int *job, int *j, LONG DOUBLE *x, 
#endif

	LONG DOUBLE *sest, LONG DOUBLE *w, LONG DOUBLE *P_gamma, LONG DOUBLE *
	sestpr, LONG DOUBLE *s, LONG DOUBLE *c)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAIC1 applies one step of incremental condition estimation in   
    its simplest version:   

    Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j 
  
    lower triangular matrix L, such that   
             twonorm(L*x) = sest   
    Then DLAIC1 computes sestpr, s, c such that   
    the vector   
                    [ s*x ]   
             xhat = [  c  ]   
    is an approximate singular vector of   
                    [ L     0  ]   
             Lhat = [ w' P_gamma ]   
    in the sense that   
             twonorm(Lhat*xhat) = sestpr.   

    Depending on JOB, an estimate for the largest or smallest singular   
    value is computed.   

    Note that [s c]' and sestpr**2 is an eigenpair of the system   

        diag(sest*sest, 0) + [alpha  P_gamma] * [ alpha ]   
                                              [ P_gamma ]   

    where  alpha =  x'*w.   

    Arguments   
    =========   

    JOB     (input) INTEGER   
            = 1: an estimate for the largest singular value is computed. 
  
            = 2: an estimate for the smallest singular value is computed. 
  

    J       (input) INTEGER   
            Length of X and W   

    X       (input) LONG DOUBLE PRECISION array, dimension (J)   
            The j-vector x.   

    SEST    (input) LONG DOUBLE PRECISION   
            Estimated singular value of j by j matrix L   

    W       (input) LONG DOUBLE PRECISION array, dimension (J)   
            The j-vector w.   

    P_GAMMA   (input) LONG DOUBLE PRECISION   
            The diagonal element P_gamma.   

    SEDTPR  (output) LONG DOUBLE PRECISION   
            Estimated singular value of (j+1) by (j+1) matrix Lhat.   

    S       (output) LONG DOUBLE PRECISION   
            Sine needed in forming xhat.   

    C       (output) LONG DOUBLE PRECISION   
            Cosine needed in forming xhat.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b5 = 1.;
    
    /* System generated locals */
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE sine, test, zeta1, zeta2, b, t, alpha, norma, s1, s2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE absgam, absalp, cosine, absest, eps, tmp;



#define W(I) w[(I)-1]
#define X(I) x[(I)-1]



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
    alpha = ddot_(j, &X(1), &c__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    alpha = qdot(j, &X(1), &c__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    alpha = qdot_(j, &X(1), &c__1, &W(1), &c__1);
#endif


    absalp = ABS(alpha);
    absgam = ABS(*P_gamma);
    absest = ABS(*sest);

    if (*job == 1) {

/*        Estimating largest singular value   

          special cases */

	if (*sest == 0.) {
	    s1 = MAX(absgam,absalp);
	    if (s1 == 0.) {
		*s = 0.;
		*c = 1.;
		*sestpr = 0.;
	    } else {
		*s = alpha / s1;
		*c = *P_gamma / s1;
		tmp = sqrt(*s * *s + *c * *c);
		*s /= tmp;
		*c /= tmp;
		*sestpr = s1 * tmp;
	    }
	    return;
	} else if (absgam <= eps * absest) {
	    *s = 1.;
	    *c = 0.;
	    tmp = MAX(absest,absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	    return;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = 1.;
		*c = 0.;
		*sestpr = s2;
	    } else {
		*s = 0.;
		*c = 1.;
		*sestpr = s1;
	    }
	    return;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		*s = sqrt(tmp * tmp + 1.);
		*sestpr = s2 * *s;
		*c = *P_gamma / s2 / *s;
		*s = SIGN(c_b5, alpha) / *s;
	    } else {
		tmp = s2 / s1;
		*c = sqrt(tmp * tmp + 1.);
		*sestpr = s1 * *c;
		*s = alpha / s1 / *c;
		*c = SIGN(c_b5, *P_gamma) / *c;
	    }
	    return;
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *P_gamma / absest;

	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
	    *c = zeta1 * zeta1;
	    if (b > 0.) {
		t = *c / (b + sqrt(b * b + *c));
	    } else {
		t = sqrt(b * b + *c) - b;
	    }

	    sine = -zeta1 / t;
	    cosine = -zeta2 / (t + 1.);
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    *sestpr = sqrt(t + 1.) * absest;
	    return;
	}

    } else if (*job == 2) {

/*        Estimating smallest singular value   

          special cases */

	if (*sest == 0.) {
	    *sestpr = 0.;
	    if (MAX(absgam,absalp) == 0.) {
		sine = 1.;
		cosine = 0.;
	    } else {
		sine = -(*P_gamma);
		cosine = alpha;
	    }
/* Computing MAX */
	    d__1 = ABS(sine), d__2 = ABS(cosine);
	    s1 = MAX(d__1,d__2);
	    *s = sine / s1;
	    *c = cosine / s1;
	    tmp = sqrt(*s * *s + *c * *c);
	    *s /= tmp;
	    *c /= tmp;
	    return;
	} else if (absgam <= eps * absest) {
	    *s = 0.;
	    *c = 1.;
	    *sestpr = absgam;
	    return;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = 0.;
		*c = 1.;
		*sestpr = s1;
	    } else {
		*s = 1.;
		*c = 0.;
		*sestpr = s2;
	    }
	    return;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		*c = sqrt(tmp * tmp + 1.);
		*sestpr = absest * (tmp / *c);
		*s = -(*P_gamma / s2) / *c;
		*c = SIGN(c_b5, alpha) / *c;
	    } else {
		tmp = s2 / s1;
		*s = sqrt(tmp * tmp + 1.);
		*sestpr = absest / *s;
		*c = alpha / s1 / *s;
		*s = -SIGN(c_b5, *P_gamma) / *s;
	    }
	    return;
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *P_gamma / absest;

/* Computing MAX */
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, ABS(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, ABS(d__2)) + zeta2 * zeta2;
	    norma = MAX(d__3,d__4);

/*           See if root is closer to zero or to ONE */

	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
		*c = zeta2 * zeta2;
		d__1 = b * b - *c;
                t = *c / (b + sqrt(( ABS(d__1))));
		sine = zeta1 / (1. - t);
		cosine = -zeta2 / t;
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
	    } else {

/*              root is closer to ONE, shift by that amount */

		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
		*c = zeta1 * zeta1;
		if (b >= 0.) {
		    t = -(*c) / (b + sqrt(b * b + *c));
		} else {
		    t = b - sqrt(b * b + *c);
		}
		sine = -zeta1 / t;
		cosine = -zeta2 / (t + 1.);
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    return;

	}
    }
    return;

/*     End of DLAIC1 */

} /* dlaic1_ */

