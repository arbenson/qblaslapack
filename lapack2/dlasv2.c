#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlasv2_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlasv2(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlasv2_(LONG DOUBLE *f, LONG DOUBLE *g, LONG DOUBLE *h, 
#endif

	LONG DOUBLE *ssmin, LONG DOUBLE *ssmax, LONG DOUBLE *snr, LONG DOUBLE *
	csr, LONG DOUBLE *snl, LONG DOUBLE *csl)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASV2 computes the singular value decomposition of a 2-by-2   
    triangular matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, ABS(SSMAX) is the larger singular value, ABS(SSMIN) is the 
  
    smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and 
  
    right singular vectors for ABS(SSMAX), giving the decomposition   

       [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]   
       [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].   

    Arguments   
    =========   

    F       (input) LONG DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) LONG DOUBLE PRECISION   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) LONG DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) LONG DOUBLE PRECISION   
            ABS(SSMIN) is the smaller singular value.   

    SSMAX   (output) LONG DOUBLE PRECISION   
            ABS(SSMAX) is the larger singular value.   

    SNL     (output) LONG DOUBLE PRECISION   
    CSL     (output) LONG DOUBLE PRECISION   
            The vector (CSL, SNL) is a unit left singular vector for the 
  
            singular value ABS(SSMAX).   

    SNR     (output) LONG DOUBLE PRECISION   
    CSR     (output) LONG DOUBLE PRECISION   
            The vector (CSR, SNR) is a unit right singular vector for the 
  
            singular value ABS(SSMAX).   

    Further Details   
    ===============   

    Any input parameter may be aliased with any output parameter.   

    Barring over/underflow and assuming a guard digit in subtraction, all 
  
    output quantities are correct to within a few units in the last   
    place (ulps).   

    In IEEE arithmetic, the code works correctly if one matrix element is 
  
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

   ===================================================================== 
*/
    /* Table of constant values */
    static LONG DOUBLE c_b3 = 2.;
    static LONG DOUBLE c_b4 = 1.;
    
    /* System generated locals */
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int pmax;
    static LONG DOUBLE temp;
    static long int swap;
    static LONG DOUBLE a, d, l, m, r, s, t, tsign, fa, ga, ha;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE ft, gt, ht, mm;
    static long int gasmal;
    static LONG DOUBLE tt, clt, crt, slt, srt;




    ft = *f;
    fa = ABS(ft);
    ht = *h;
    ha = ABS(*h);

/*     PMAX points to the maximum absolute element of matrix   
         PMAX = 1 if F largest in absolute values   
         PMAX = 2 if G largest in absolute values   
         PMAX = 3 if H largest in absolute values */

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

/*        Now FA .ge. HA */

    }
    gt = *g;
    ga = ABS(gt);
    if (ga == 0.) {

/*        Diagonal matrix */

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = 1;
	if (ga > fa) {
	    pmax = 2;

#ifdef PETSC_PREFIX_SUFFIX
	    if (fa / ga < dlamch_("EPS")) {
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    if (fa / ga < qlamch("EPS")) {
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    if (fa / ga < qlamch_("EPS")) {
#endif


/*              Case of very large GA */

		gasmal = 0;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

/*           Normal case */

	    d = fa - ha;
	    if (d == fa) {

/*              Copes with infinite F or H */

		l = 1.;
	    } else {
		l = d / fa;
	    }

/*           Note that 0 .le. L .le. 1 */

	    m = gt / ft;

/*           Note that ABS(M) .le. 1/macheps */

	    t = 2. - l;

/*           Note that T .ge. 1 */

	    mm = m * m;
	    tt = t * t;
	    s = sqrt(tt + mm);

/*           Note that 1 .le. S .le. 1 + 1/macheps */

	    if (l == 0.) {
		r = ABS(m);
	    } else {
		r = sqrt(l * l + mm);
	    }

/*           Note that 0 .le. R .le. 1 + 1/macheps */

	    a = (s + r) * .5;

/*           Note that 1 .le. A .le. 1 + ABS(M) */

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if (mm == 0.) {

/*              Note that M is very tiny */

		if (l == 0.) {
		    t = SIGN(c_b3, ft) * SIGN(c_b4, gt);
		} else {
		    t = gt / SIGN(d, ft) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r + l)) * (a + 1.);
	    }
	    l = sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

/*     Correct signs of SSMAX and SSMIN */

    if (pmax == 1) {
	tsign = SIGN(c_b4, *csr) * SIGN(c_b4, *csl) * SIGN(c_b4, *f);
    }
    if (pmax == 2) {
	tsign = SIGN(c_b4, *snr) * SIGN(c_b4, *csl) * SIGN(c_b4, *g);
    }
    if (pmax == 3) {
	tsign = SIGN(c_b4, *snr) * SIGN(c_b4, *snl) * SIGN(c_b4, *h);
    }
    *ssmax = SIGN(*ssmax, tsign);
    d__1 = tsign * SIGN(c_b4, *f) * SIGN(c_b4, *h);
    *ssmin = SIGN(*ssmin, d__1);
    return;

/*     End of DLASV2 */

} /* dlasv2_ */

