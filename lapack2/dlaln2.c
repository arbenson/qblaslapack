#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaln2_(long int *ltrans, int *na, int *nw, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaln2(long int *ltrans, int *na, int *nw, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaln2_(long int *ltrans, int *na, int *nw, 
#endif

	LONG DOUBLE *smin, LONG DOUBLE *ca, LONG DOUBLE *a, int *lda, 
	LONG DOUBLE *d1, LONG DOUBLE *d2, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *wr, LONG DOUBLE *wi, LONG DOUBLE *x, int *ldx, 
	LONG DOUBLE *scale, LONG DOUBLE *xnorm, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLALN2 solves a system of the form  (ca A - w D ) X = s B   
    or (ca A' - w D) X = s B   with possible scaling ("s") and   
    perturbation of A.  (A' means A-transpose.)   

    A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA   
    real diagonal matrix, w is a real or complex value, and X and B are   
    NA x 1 matrices -- real if w is real, complex if w is complex.  NA   
    may be 1 or 2.   

    If w is complex, X and B are represented as NA x 2 matrices,   
    the first column of each being the real part and the second   
    being the imaginary part.   

    "s" is a scaling factor (.LE. 1), computed by DLALN2, which is   
    so chosen that X can be computed without overflow.  X is further   
    scaled if necessary to assure that norm(ca A - w D)*norm(X) is less   
    than overflow.   

    If both singular values of (ca A - w D) are less than SMIN,   
    SMIN*identity will be used instead of (ca A - w D).  If only one   
    singular value is less than SMIN, one element of (ca A - w D) will be 
  
    perturbed enough to make the smallest singular value roughly SMIN.   
    If both singular values are at least SMIN, (ca A - w D) will not be   
    perturbed.  In any case, the perturbation will be at most some small 
  
    multiple of MAX( SMIN, ulp*norm(ca A - w D) ).  The singular values   
    are computed by infinity-norm approximations, and thus will only be   
    correct to a factor of 2 or so.   

    Note: all input quantities are assumed to be smaller than overflow   
    by a reasonable factor.  (See BIGNUM.)   

    Arguments   
    ==========   

    LTRANS  (input) LOGICAL   
            =.TRUE.:  A-transpose will be used.   
            =.FALSE.: A will be used (not transposed.)   

    NA      (input) INTEGER   
            The size of the matrix A.  It may (only) be 1 or 2.   

    NW      (input) INTEGER   
            1 if "w" is real, 2 if "w" is complex.  It may only be 1   
            or 2.   

    SMIN    (input) LONG DOUBLE PRECISION   
            The desired lower bound on the singular values of A.  This   
            should be a safe distance away from underflow or overflow,   
            say, between (underflow/machine precision) and  (machine   
            precision * overflow ).  (See BIGNUM and ULP.)   

    CA      (input) LONG DOUBLE PRECISION   
            The coefficient c, which A is multiplied by.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,NA)   
            The NA x NA matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of A.  It must be at least NA.   

    D1      (input) LONG DOUBLE PRECISION   
            The 1,1 element in the diagonal matrix D.   

    D2      (input) LONG DOUBLE PRECISION   
            The 2,2 element in the diagonal matrix D.  Not used if NW=1. 
  

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NW)   
            The NA x NW matrix B (right-hand side).  If NW=2 ("w" is   
            complex), column 1 contains the real part of B and column 2   
            contains the imaginary part.   

    LDB     (input) INTEGER   
            The leading dimension of B.  It must be at least NA.   

    WR      (input) LONG DOUBLE PRECISION   
            The real part of the scalar "w".   

    WI      (input) LONG DOUBLE PRECISION   
            The imaginary part of the scalar "w".  Not used if NW=1.   

    X       (output) LONG DOUBLE PRECISION array, dimension (LDX,NW)   
            The NA x NW matrix X (unknowns), as computed by DLALN2.   
            If NW=2 ("w" is complex), on exit, column 1 will contain   
            the real part of X and column 2 will contain the imaginary   
            part.   

    LDX     (input) INTEGER   
            The leading dimension of X.  It must be at least NA.   

    SCALE   (output) LONG DOUBLE PRECISION   
            The scale factor that B must be multiplied by to insure   
            that overflow does not occur when computing X.  Thus,   
            (ca A - w D) X  will be SCALE*B, not B (ignoring   
            perturbations of A.)  It will be at most 1.   

    XNORM   (output) LONG DOUBLE PRECISION   
            The infinity-norm of X, when X is regarded as an NA x NW   
            real matrix.   

    INFO    (output) INTEGER   
            An error flag.  It will be set to zero if no error occurs,   
            a negative number if an argument is in error, or a positive   
            number if  ca A - w D  had to be perturbed.   
            The possible values are:   
            = 0: No error occurred, and (ca A - w D) did not have to be   
                   perturbed.   
            = 1: (ca A - w D) had to be perturbed to make its smallest   
                 (or only) singular value greater than SMIN.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

   ===================================================================== 
  

    
   Parameter adjustments   
       Function Body */
    /* Initialized data */
    static long int zswap[4] = { 0,0,1,1 };
    static long int rswap[4] = { 0,1,0,1 };
    static int ipivot[16]	/* was [4][4] */ = { 1,2,3,4,2,1,4,3,3,4,1,2,
	    4,3,2,1 };
    /* System generated locals */
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5, d__6;
    static LONG DOUBLE equiv_0[4], equiv_1[4];
    /* Local variables */
    static LONG DOUBLE bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s;
    static int j;
    static LONG DOUBLE u22abs;
    static int icmax;
    static LONG DOUBLE bnorm, cnorm, smini;
#define ci (equiv_0)
#define cr (equiv_1)

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE bignum, bi1, bi2, br1, br2, smlnum, xi1, xi2, xr1, xr2, 
	    ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
    static LONG DOUBLE csr, ur11, ur12, ur22;
#define crv (equiv_1)


#define IPIVOT(I) ipivot[(I)]
#define WAS(I) was[(I)]
#define EQUIV_0(I) equiv_0[(I)]
#define EQUIV_1(I) equiv_1[(I)]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]


/*     Compute BIGNUM */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("Safe minimum") * 2.;
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("Safe minimum") * 2.;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("Safe minimum") * 2.;
#endif

    bignum = 1. / smlnum;
    smini = MAX(*smin,smlnum);

/*     Don't check for input errors */

    *info = 0;

/*     Standard Initializations */

    *scale = 1.;

    if (*na == 1) {

/*        1 x 1  (i.e., scalar) system   C X = B */

	if (*nw == 1) {

/*           Real 1x1 system.   

             C = ca A - w D */

	    csr = *ca * A(1,1) - *wr * *d1;
	    cnorm = ABS(csr);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = B(1,1), ABS(d__1));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    X(1,1) = B(1,1) * *scale / csr;
	    *xnorm = (d__1 = X(1,1), ABS(d__1));
	} else {

/*           Complex 1x1 system (w is complex)   

             C = ca A - w D */

	    csr = *ca * A(1,1) - *wr * *d1;
	    csi = -(*wi) * *d1;
	    cnorm = ABS(csr) + ABS(csi);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		csi = 0.;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = B(1,1), ABS(d__1)) + (d__2 = B(1,2), ABS(d__2));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    d__1 = *scale * B(1,1);
	    d__2 = *scale * B(1,2);

#ifdef PETSC_PREFIX_SUFFIX
	    dladiv_(&d__1, &d__2, &csr, &csi, &X(1,1), &X(1,2));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qladiv(&d__1, &d__2, &csr, &csi, &X(1,1), &X(1,2));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qladiv_(&d__1, &d__2, &csr, &csi, &X(1,1), &X(1,2));
#endif

	    *xnorm = (d__1 = X(1,1), ABS(d__1)) + (d__2 = X(1,2), ABS(d__2));
	}

    } else {

/*        2x2 System   

          Compute the real part of  C = ca A - w D  (or  ca A' - w D )
 */

	cr[0] = *ca * A(1,1) - *wr * *d1;
	cr[3] = *ca * A(2,2) - *wr * *d2;
	if (*ltrans) {
	    cr[2] = *ca * A(2,1);
	    cr[1] = *ca * A(1,2);
	} else {
	    cr[1] = *ca * A(2,1);
	    cr[2] = *ca * A(1,2);
	}

	if (*nw == 1) {

/*           Real 2x2 system  (w is real)   

             Find the largest element in C */

	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], ABS(d__1)) > cmax) {
		    cmax = (d__1 = crv[j - 1], ABS(d__1));
		    icmax = j;
		}
/* L10: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__3 = (d__1 = B(1,1), ABS(d__1)), d__4 = (d__2 = B(2,1), ABS(d__2));
		bnorm = MAX(d__3,d__4);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		X(1,1) = temp * B(1,1);
		X(2,1) = temp * B(2,1);
		*xnorm = temp * bnorm;
		*info = 1;
		return;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    cr21 = crv[IPIVOT((icmax << 2) - 3) - 1];
	    ur12 = crv[IPIVOT((icmax << 2) - 2) - 1];
	    cr22 = crv[IPIVOT((icmax << 2) - 1) - 1];
	    ur11r = 1. / ur11;
	    lr21 = ur11r * cr21;
	    ur22 = cr22 - ur12 * lr21;

/*           If smaller pivot < SMINI, use SMINI */

	    if (ABS(ur22) < smini) {
		ur22 = smini;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br1 = B(2,1);
		br2 = B(1,1);
	    } else {
		br1 = B(1,1);
		br2 = B(2,1);
	    }
	    br2 -= lr21 * br1;
/* Computing MAX */
	    d__2 = (d__1 = br1 * (ur22 * ur11r), ABS(d__1)), d__3 = ABS(br2);
	    bbnd = MAX(d__2,d__3);
	    if (bbnd > 1. && ABS(ur22) < 1.) {
		if (bbnd >= bignum * ABS(ur22)) {
		    *scale = 1. / bbnd;
		}
	    }

	    xr2 = br2 * *scale / ur22;
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
	    if (zswap[icmax - 1]) {
		X(1,1) = xr2;
		X(2,1) = xr1;
	    } else {
		X(1,1) = xr1;
		X(2,1) = xr2;
	    }
/* Computing MAX */
	    d__1 = ABS(xr1), d__2 = ABS(xr2);
	    *xnorm = MAX(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    X(1,1) = temp * X(1,1);
		    X(2,1) = temp * X(2,1);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	} else {

/*           Complex 2x2 system  (w is complex)   

             Find the largest element in C */

	    ci[0] = -(*wi) * *d1;
	    ci[1] = 0.;
	    ci[2] = 0.;
	    ci[3] = -(*wi) * *d2;
	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], ABS(d__1)) + (d__2 = civ[j - 1], ABS(
			d__2)) > cmax) {
		    cmax = (d__1 = crv[j - 1], ABS(d__1)) + (d__2 = civ[j - 1]
			    , ABS(d__2));
		    icmax = j;
		}
/* L20: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__5 = (d__1 = B(1,1), ABS(d__1)) + (d__2 = B(1,2), ABS(d__2)), d__6 = (d__3 = B(2,1), 
			ABS(d__3)) + (d__4 = B(2,2), ABS(d__4));
		bnorm = MAX(d__5,d__6);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		X(1,1) = temp * B(1,1);
		X(2,1) = temp * B(2,1);
		X(1,2) = temp * B(1,2);
		X(2,2) = temp * B(2,2);
		*xnorm = temp * bnorm;
		*info = 1;
		return;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    ui11 = civ[icmax - 1];
	    cr21 = crv[IPIVOT((icmax << 2) - 3) - 1];
	    ci21 = civ[IPIVOT((icmax << 2) - 3) - 1];
	    ur12 = crv[IPIVOT((icmax << 2) - 2) - 1];
	    ui12 = civ[IPIVOT((icmax << 2) - 2) - 1];
	    cr22 = crv[IPIVOT((icmax << 2) - 1) - 1];
	    ci22 = civ[IPIVOT((icmax << 2) - 1) - 1];
	    if (icmax == 1 || icmax == 4) {

/*              Code when off-diagonals of pivoted C are real 
*/

		if (ABS(ur11) > ABS(ui11)) {
		    temp = ui11 / ur11;
/* Computing 2nd power */
		    d__1 = temp;
		    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
		    ui11r = -temp * ur11r;
		} else {
		    temp = ur11 / ui11;
/* Computing 2nd power */
		    d__1 = temp;
		    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
		    ur11r = -temp * ui11r;
		}
		lr21 = cr21 * ur11r;
		li21 = cr21 * ui11r;
		ur12s = ur12 * ur11r;
		ui12s = ur12 * ui11r;
		ur22 = cr22 - ur12 * lr21;
		ui22 = ci22 - ur12 * li21;
	    } else {

/*              Code when diagonals of pivoted C are real */

		ur11r = 1. / ur11;
		ui11r = 0.;
		lr21 = cr21 * ur11r;
		li21 = ci21 * ur11r;
		ur12s = ur12 * ur11r;
		ui12s = ui12 * ur11r;
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
		ui22 = -ur12 * li21 - ui12 * lr21;
	    }
	    u22abs = ABS(ur22) + ABS(ui22);

/*           If smaller pivot < SMINI, use SMINI */

	    if (u22abs < smini) {
		ur22 = smini;
		ui22 = 0.;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br2 = B(1,1);
		br1 = B(2,1);
		bi2 = B(1,2);
		bi1 = B(2,2);
	    } else {
		br1 = B(1,1);
		br2 = B(2,1);
		bi1 = B(1,2);
		bi2 = B(2,2);
	    }
	    br2 = br2 - lr21 * br1 + li21 * bi1;
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
/* Computing MAX */
	    d__1 = (ABS(br1) + ABS(bi1)) * (u22abs * (ABS(ur11r) + ABS(ui11r))
		    ), d__2 = ABS(br2) + ABS(bi2);
	    bbnd = MAX(d__1,d__2);
	    if (bbnd > 1. && u22abs < 1.) {
		if (bbnd >= bignum * u22abs) {
		    *scale = 1. / bbnd;
		    br1 = *scale * br1;
		    bi1 = *scale * bi1;
		    br2 = *scale * br2;
		    bi2 = *scale * bi2;
		}
	    }


#ifdef PETSC_PREFIX_SUFFIX
	    dladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qladiv(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
#endif

	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
	    if (zswap[icmax - 1]) {
		X(1,1) = xr2;
		X(2,1) = xr1;
		X(1,2) = xi2;
		X(2,2) = xi1;
	    } else {
		X(1,1) = xr1;
		X(2,1) = xr2;
		X(1,2) = xi1;
		X(2,2) = xi2;
	    }
/* Computing MAX */
	    d__1 = ABS(xr1) + ABS(xi1), d__2 = ABS(xr2) + ABS(xi2);
	    *xnorm = MAX(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    X(1,1) = temp * X(1,1);
		    X(2,1) = temp * X(2,1);
		    X(1,2) = temp * X(1,2);
		    X(2,2) = temp * X(2,2);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	}
    }

    return;

/*     End of DLALN2 */

} /* dlaln2_ */

#undef crv
#undef civ
#undef cr
#undef ci


