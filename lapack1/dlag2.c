#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlag2_(LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlag2(LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlag2_(LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
#endif

	int *ldb, LONG DOUBLE *safmin, LONG DOUBLE *scale1, LONG DOUBLE *
	scale2, LONG DOUBLE *wr1, LONG DOUBLE *wr2, LONG DOUBLE *wi)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue   
    problem  A - w B, with scaling as necessary to avoid over-/underflow. 
  

    The scaling factor "s" results in a modified eigenvalue equation   

        s A - w B   

    where  s  is a non-negative scaling factor chosen so that  w,  w B,   
    and  s A  do not overflow and, if possible, do not underflow, either. 
  

    Arguments   
    =========   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA, 2)   
            On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm 
  
            is less than 1/SAFMIN.  Entries less than   
            sqrt(SAFMIN)*norm(A) are subject to being treated as zero.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= 2.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB, 2)   
            On entry, the 2 x 2 upper triangular matrix B.  It is   
            assumed that the one-norm of B is less than 1/SAFMIN.  The   
            diagonals should be at least sqrt(SAFMIN) times the largest   
            element of B (in absolute value); if a diagonal is smaller   
            than that, then  +/- sqrt(SAFMIN) will be used instead of   
            that diagonal.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= 2.   

    SAFMIN  (input) LONG DOUBLE PRECISION   
            The smallest positive number s.t. 1/SAFMIN does not   
            overflow.  (This should always be DLAMCH('S') -- it is an   
            argument in order to avoid having to call DLAMCH frequently.) 
  

    SCALE1  (output) LONG DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the first eigenvalue.  If   
            the eigenvalues are complex, then the eigenvalues are   
            ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the   
            exponent range of the machine), SCALE1=SCALE2, and SCALE1   
            will always be positive.  If the eigenvalues are real, then   
            the first (real) eigenvalue is  WR1 / SCALE1 , but this may   
            overflow or underflow, and in fact, SCALE1 may be zero or   
            less than the underflow threshhold if the exact eigenvalue   
            is sufficiently large.   

    SCALE2  (output) LONG DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the second eigenvalue.  If 
  
            the eigenvalues are complex, then SCALE2=SCALE1.  If the   
            eigenvalues are real, then the second (real) eigenvalue is   
            WR2 / SCALE2 , but this may overflow or underflow, and in   
            fact, SCALE2 may be zero or less than the underflow   
            threshhold if the exact eigenvalue is sufficiently large.   

    WR1     (output) LONG DOUBLE PRECISION   
            If the eigenvalue is real, then WR1 is SCALE1 times the   
            eigenvalue closest to the (2,2) element of A B**(-1).  If the 
  
            eigenvalue is complex, then WR1=WR2 is SCALE1 times the real 
  
            part of the eigenvalues.   

    WR2     (output) LONG DOUBLE PRECISION   
            If the eigenvalue is real, then WR2 is SCALE2 times the   
            other eigenvalue.  If the eigenvalue is complex, then   
            WR1=WR2 is SCALE1 times the real part of the eigenvalues.   

    WI      (output) LONG DOUBLE PRECISION   
            If the eigenvalue is real, then WI is zero.  If the   
            eigenvalue is complex, then WI is SCALE1 times the imaginary 
  
            part of the eigenvalues.  WI will always be non-negative.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5, d__6;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE diff, bmin, wbig, wabs, wdet, r, binv11, binv22, discr, 
	    anorm, bnorm, bsize, shift, c1, c2, c3, c4, c5, rtmin, rtmax, 
	    wsize, s1, s2, a11, a12, a21, a22, b11, b12, b22, ascale, bscale, 
	    pp, qq, ss, wscale, safmax, wsmall, as11, as12, as22, sum, abi22;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    rtmin = sqrt(*safmin);
    rtmax = 1. / rtmin;
    safmax = 1. / *safmin;

/*     Scale A   

   Computing MAX */
    d__5 = (d__1 = A(1,1), ABS(d__1)) + (d__2 = A(2,1), ABS(
	    d__2)), d__6 = (d__3 = A(1,2), ABS(d__3)) + (d__4 = 
	    A(2,2), ABS(d__4)), d__5 = MAX(d__5,d__6);
    anorm = MAX(d__5,*safmin);
    ascale = 1. / anorm;
    a11 = ascale * A(1,1);
    a21 = ascale * A(2,1);
    a12 = ascale * A(1,2);
    a22 = ascale * A(2,2);

/*     Perturb B if necessary to insure non-singularity */

    b11 = B(1,1);
    b12 = B(1,2);
    b22 = B(2,2);
/* Computing MAX */
    d__1 = ABS(b11), d__2 = ABS(b12), d__1 = MAX(d__1,d__2), d__2 = ABS(b22), 
	    d__1 = MAX(d__1,d__2);
    bmin = rtmin * MAX(d__1,rtmin);
    if (ABS(b11) < bmin) {
	b11 = SIGN(bmin, b11);
    }
    if (ABS(b22) < bmin) {
	b22 = SIGN(bmin, b22);
    }

/*     Scale B   

   Computing MAX */
    d__1 = ABS(b11), d__2 = ABS(b12) + ABS(b22), d__1 = MAX(d__1,d__2);
    bnorm = MAX(d__1,*safmin);
/* Computing MAX */
    d__1 = ABS(b11), d__2 = ABS(b22);
    bsize = MAX(d__1,d__2);
    bscale = 1. / bsize;
    b11 *= bscale;
    b12 *= bscale;
    b22 *= bscale;

/*     Compute larger eigenvalue by method described by C. van Loan   

       ( AS is A shifted by -SHIFT*B ) */

    binv11 = 1. / b11;
    binv22 = 1. / b22;
    s1 = a11 * binv11;
    s2 = a22 * binv22;
    if (ABS(s1) <= ABS(s2)) {
	as12 = a12 - s1 * b12;
	as22 = a22 - s1 * b22;
	ss = a21 * (binv11 * binv22);
	abi22 = as22 * binv22 - ss * b12;
	pp = abi22 * .5;
	shift = s1;
    } else {
	as12 = a12 - s2 * b12;
	as11 = a11 - s2 * b11;
	ss = a21 * (binv11 * binv22);
	abi22 = -ss * b12;
	pp = (as11 * binv11 + abi22) * .5;
	shift = s2;
    }
    qq = ss * as12;
    if ((d__1 = pp * rtmin, ABS(d__1)) >= 1.) {
/* Computing 2nd power */
	d__1 = rtmin * pp;
	discr = d__1 * d__1 + qq * *safmin;
	r = sqrt((ABS(discr))) * rtmax;
    } else {
/* Computing 2nd power */
	d__1 = pp;
	if (d__1 * d__1 + ABS(qq) <= *safmin) {
/* Computing 2nd power */
	    d__1 = rtmax * pp;
	    discr = d__1 * d__1 + qq * safmax;
	    r = sqrt((ABS(discr))) * rtmin;
	} else {
/* Computing 2nd power */
	    d__1 = pp;
	    discr = d__1 * d__1 + qq;
	    r = sqrt((ABS(discr)));
	}
    }

/*     Note: the test of R in the following IF is to cover the case when 
  
             DISCR is small and negative and is flushed to zero during   
             the calculation of R.  On machines which have a consistent   
             flush-to-zero threshhold and handle numbers above that   
             threshhold correctly, it would not be necessary. */

    if (discr >= 0. || r == 0.) {
	sum = pp + SIGN(r, pp);
	diff = pp - SIGN(r, pp);
	wbig = shift + sum;

/*        Compute smaller eigenvalue */

	wsmall = shift + diff;
/* Computing MAX */
	d__1 = ABS(wsmall);
	if (ABS(wbig) * .5 > MAX(d__1,*safmin)) {
	    wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
	    wsmall = wdet / wbig;
	}

/*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
   
          for WR1. */

	if (pp > abi22) {
	    *wr1 = MIN(wbig,wsmall);
	    *wr2 = MAX(wbig,wsmall);
	} else {
	    *wr1 = MAX(wbig,wsmall);
	    *wr2 = MIN(wbig,wsmall);
	}
	*wi = 0.;
    } else {

/*        Complex eigenvalues */

	*wr1 = shift + pp;
	*wr2 = *wr1;
	*wi = r;
    }

/*     Further scaling to avoid underflow and overflow in computing   
       SCALE1 and overflow in computing w*B.   

       This scale factor (WSCALE) is bounded from above using C1 and C2, 
  
       and from below using C3 and C4.   
          C1 implements the condition  s A  must never overflow.   
          C2 implements the condition  w B  must never overflow.   
          C3, with C2,   
             implement the condition that s A - w B must never overflow. 
  
          C4 implements the condition  s    should not underflow.   
          C5 implements the condition  MAX(s,|w|) should be at least 2. */

    c1 = bsize * (*safmin * MAX(1.,ascale));
    c2 = *safmin * MAX(1.,bnorm);
    c3 = bsize * *safmin;
    if (ascale <= 1. && bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale / *safmin * bsize;
	c4 = MIN(d__1,d__2);
    } else {
	c4 = 1.;
    }
    if (ascale <= 1. || bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale * bsize;
	c5 = MIN(d__1,d__2);
    } else {
	c5 = 1.;
    }

/*     Scale first eigenvalue */

    wabs = ABS(*wr1) + ABS(*wi);
/* Computing MAX   
   Computing MIN */
    d__3 = c4, d__4 = MAX(wabs,c5) * .5;
    d__1 = MAX(*safmin,c1), d__2 = (wabs * c2 + c3) * 1.0000100000000001, 
	    d__1 = MAX(d__1,d__2), d__2 = MIN(d__3,d__4);
    wsize = MAX(d__1,d__2);
    if (wsize != 1.) {
	wscale = 1. / wsize;
	if (wsize > 1.) {
	    *scale1 = MAX(ascale,bsize) * wscale * MIN(ascale,bsize);
	} else {
	    *scale1 = MIN(ascale,bsize) * wscale * MAX(ascale,bsize);
	}
	*wr1 *= wscale;
	if (*wi != 0.) {
	    *wi *= wscale;
	    *wr2 = *wr1;
	    *scale2 = *scale1;
	}
    } else {
	*scale1 = ascale * bsize;
	*scale2 = *scale1;
    }

/*     Scale second eigenvalue (if real) */

    if (*wi == 0.) {
/* Computing MAX   
   Computing MIN   
   Computing MAX */
	d__5 = ABS(*wr2);
	d__3 = c4, d__4 = MAX(d__5,c5) * .5;
	d__1 = MAX(*safmin,c1), d__2 = (ABS(*wr2) * c2 + c3) * 
		1.0000100000000001, d__1 = MAX(d__1,d__2), d__2 = MIN(d__3,
		d__4);
	wsize = MAX(d__1,d__2);
	if (wsize != 1.) {
	    wscale = 1. / wsize;
	    if (wsize > 1.) {
		*scale2 = MAX(ascale,bsize) * wscale * MIN(ascale,bsize);
	    } else {
		*scale2 = MIN(ascale,bsize) * wscale * MAX(ascale,bsize);
	    }
	    *wr2 *= wscale;
	} else {
	    *scale2 = ascale * bsize;
	}
    }

/*     End of DLAG2 */

    return;
} /* dlag2_ */

