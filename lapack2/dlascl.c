#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlascl_(char *type, int *kl, int *ku, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlascl(char *type, int *kl, int *ku, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlascl_(char *type, int *kl, int *ku, LONG DOUBLE 
#endif

	*cfrom, LONG DOUBLE *cto, int *m, int *n, LONG DOUBLE *a, 
	int *lda, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLASCL multiplies the M by N real matrix A by the real scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final 
  
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that 
  
    A may be full, upper triangular, lower triangular, upper Hessenberg, 
  
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper 
  
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) LONG DOUBLE PRECISION   
    CTO     (input) LONG DOUBLE PRECISION   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static long int done;
    static LONG DOUBLE ctoc;
    static int i, j;
    extern long int lsame_(char *, char *);
    static int itype, k1, k2, k3, k4;
    static LONG DOUBLE cfrom1;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE cfromc;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum, smlnum, mul, cto1;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;

    if (lsame_(type, "G")) {
	itype = 0;
    } else if (lsame_(type, "L")) {
	itype = 1;
    } else if (lsame_(type, "U")) {
	itype = 2;
    } else if (lsame_(type, "H")) {
	itype = 3;
    } else if (lsame_(type, "B")) {
	itype = 4;
    } else if (lsame_(type, "Q")) {
	itype = 5;
    } else if (lsame_(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || (itype == 4 && *n != *m) || (itype == 5 && *n != *m)) {
	*info = -7;
    } else if (itype <= 3 && *lda < MAX(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > MAX(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > MAX(i__1,0) || ((itype == 4 || itype == 5) && 
		    *kl != *ku)) {
		*info = -3;
	    } else if ((itype == 4 && *lda < *kl + 1) || (itype == 5 && *lda < *
		    ku + 1) || (itype == 6 && *lda < (*kl << 1) + *ku + 1)) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLASCL", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return;
    }

/*     Get machine parameters */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S");
#endif

    bignum = 1. / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (ABS(cfrom1) > ABS(ctoc) && ctoc != 0.) {
	mul = smlnum;
	done = 0;
	cfromc = cfrom1;
    } else if (ABS(cto1) > ABS(cfromc)) {
	mul = bignum;
	done = 0;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = 1;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) *= mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = j; i <= *m; ++i) {
		A(i,j) *= mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = MIN(j,*m);
	    for (i = 1; i <= MIN(j,*m); ++i) {
		A(i,j) *= mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = MIN(i__3,*m);
	    for (i = 1; i <= MIN(j+1,*m); ++i) {
		A(i,j) *= mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = MIN(i__3,i__4);
	    for (i = 1; i <= MIN(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i = MAX(k1-j,1); i <= k3; ++i) {
		A(i,j) *= mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = MIN(i__4,i__5);
	    for (i = MAX(k1-j,k2); i <= MIN(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return;

/*     End of DLASCL */

} /* dlascl_ */

