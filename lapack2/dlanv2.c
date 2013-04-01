#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlanv2_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlanv2(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlanv2_(LONG DOUBLE *a, LONG DOUBLE *b, LONG DOUBLE *c, 
#endif

	LONG DOUBLE *d, LONG DOUBLE *rt1r, LONG DOUBLE *rt1i, LONG DOUBLE *rt2r, 
	LONG DOUBLE *rt2i, LONG DOUBLE *cs, LONG DOUBLE *sn)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric 
  
    matrix in standard form:   

         [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]   
         [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]   

    where either   
    1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or   
    2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex   
    conjugate eigenvalues.   

    Arguments   
    =========   

    A       (input/output) LONG DOUBLE PRECISION   
    B       (input/output) LONG DOUBLE PRECISION   
    C       (input/output) LONG DOUBLE PRECISION   
    D       (input/output) LONG DOUBLE PRECISION   
            On entry, the elements of the input matrix.   
            On exit, they are overwritten by the elements of the   
            standardised Schur form.   

    RT1R    (output) LONG DOUBLE PRECISION   
    RT1I    (output) LONG DOUBLE PRECISION   
    RT2R    (output) LONG DOUBLE PRECISION   
    RT2I    (output) LONG DOUBLE PRECISION   
            The real and imaginary parts of the eigenvalues. If the   
            eigenvalues are both real, ABS(RT1R) >= ABS(RT2R); if the   
            eigenvalues are a complex conjugate pair, RT1I > 0.   

    CS      (output) LONG DOUBLE PRECISION   
    SN      (output) LONG DOUBLE PRECISION   
            Parameters of the rotation matrix.   

    ===================================================================== 
  


       Initialize CS and SN */
    /* Table of constant values */
    static LONG DOUBLE c_b3 = 1.;
    
    /* System generated locals */
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp, p, sigma;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static LONG DOUBLE aa, bb, cc, dd, cs1, sn1, sab, sac, tau;



    *cs = 1.;
    *sn = 0.;

    if (*c == 0.) {
	goto L10;

    } else if (*b == 0.) {

/*        Swap rows and columns */

	*cs = 0.;
	*sn = 1.;
	temp = *d;
	*d = *a;
	*a = temp;
	*b = -(*c);
	*c = 0.;
	goto L10;
    } else if (*a - *d == 0. && SIGN(c_b3, *b) != SIGN(c_b3, *c)) {
	goto L10;
    } else {

/*        Make diagonal elements equal */

	temp = *a - *d;
	p = temp * .5;
	sigma = *b + *c;

#ifdef PETSC_PREFIX_SUFFIX
	tau = dlapy2_(&sigma, &temp);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	tau = qlapy2(&sigma, &temp);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	tau = qlapy2_(&sigma, &temp);
#endif

	cs1 = sqrt((ABS(sigma) / tau + 1.) * .5);
	sn1 = -(p / (tau * cs1)) * SIGN(c_b3, sigma);

/*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]   
                  [ CC  DD ]   [ C  D ] [ SN1  CS1 ] */

	aa = *a * cs1 + *b * sn1;
	bb = -(*a) * sn1 + *b * cs1;
	cc = *c * cs1 + *d * sn1;
	dd = -(*c) * sn1 + *d * cs1;

/*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]   
                  [ C  D ]   [-SN1  CS1 ] [ CC  DD ] */

	*a = aa * cs1 + cc * sn1;
	*b = bb * cs1 + dd * sn1;
	*c = -aa * sn1 + cc * cs1;
	*d = -bb * sn1 + dd * cs1;

/*        Accumulate transformation */

	temp = *cs * cs1 - *sn * sn1;
	*sn = *cs * sn1 + *sn * cs1;
	*cs = temp;

	temp = (*a + *d) * .5;
	*a = temp;
	*d = temp;

	if (*c != 0.) {
	    if (*b != 0.) {
		if (SIGN(c_b3, *b) == SIGN(c_b3, *c)) {

/*                 Real eigenvalues: reduce to upper trian
gular form */

		    sab = sqrt((ABS(*b)));
		    sac = sqrt((ABS(*c)));
		    d__1 = sab * sac;
		    p = SIGN(d__1, *c);
		    d__1 = *b + *c;
                    tau = 1. / sqrt(( ABS(d__1)));
		    *a = temp + p;
		    *d = temp - p;
		    *b -= *c;
		    *c = 0.;
		    cs1 = sab * tau;
		    sn1 = sac * tau;
		    temp = *cs * cs1 - *sn * sn1;
		    *sn = *cs * sn1 + *sn * cs1;
		    *cs = temp;
		}
	    } else {
		*b = -(*c);
		*c = 0.;
		temp = *cs;
		*cs = -(*sn);
		*sn = temp;
	    }
	}
    }

L10:

/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

    *rt1r = *a;
    *rt2r = *d;
    if (*c == 0.) {
	*rt1i = 0.;
	*rt2i = 0.;
    } else {
	*rt1i = sqrt((ABS(*b))) * sqrt((ABS(*c)));
	*rt2i = -(*rt1i);
    }
    return;

/*     End of DLANV2 */

} /* dlanv2_ */

