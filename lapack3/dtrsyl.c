#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrsyl_(char *trana, char *tranb, int *isgn, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrsyl(char *trana, char *tranb, int *isgn, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrsyl_(char *trana, char *tranb, int *isgn, int 
#endif

	*m, int *n, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *
	ldb, LONG DOUBLE *c, int *ldc, LONG DOUBLE *scale, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTRSYL solves the real Sylvester matrix equation:   

       op(A)*X + X*op(B) = scale*C or   
       op(A)*X - X*op(B) = scale*C,   

    where op(A) = A or A**T, and  A and B are both upper quasi-   
    triangular. A is M-by-M and B is N-by-N; the right hand side C and   
    the solution X are M-by-N; and scale is an output scale factor, set   
    <= 1 to avoid overflow in X.   

    A and B must be in Schur canonical form (as returned by DHSEQR), that 
  
    is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;   
    each 2-by-2 diagonal block has its diagonal elements equal and its   
    off-diagonal elements of opposite sign.   

    Arguments   
    =========   

    TRANA   (input) CHARACTER*1   
            Specifies the option op(A):   
            = 'N': op(A) = A    (No transpose)   
            = 'T': op(A) = A**T (Transpose)   
            = 'C': op(A) = A**H (Conjugate transpose = Transpose)   

    TRANB   (input) CHARACTER*1   
            Specifies the option op(B):   
            = 'N': op(B) = B    (No transpose)   
            = 'T': op(B) = B**T (Transpose)   
            = 'C': op(B) = B**H (Conjugate transpose = Transpose)   

    ISGN    (input) INTEGER   
            Specifies the sign in the equation:   
            = +1: solve op(A)*X + X*op(B) = scale*C   
            = -1: solve op(A)*X - X*op(B) = scale*C   

    M       (input) INTEGER   
            The order of the matrix A, and the number of rows in the   
            matrices X and C. M >= 0.   

    N       (input) INTEGER   
            The order of the matrix B, and the number of columns in the   
            matrices X and C. N >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,M)   
            The upper quasi-triangular matrix A, in Schur canonical form. 
  

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            The upper quasi-triangular matrix B, in Schur canonical form. 
  

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,N).   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N right hand side matrix C.   
            On exit, C is overwritten by the solution matrix X.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= MAX(1,M)   

    SCALE   (output) LONG DOUBLE PRECISION   
            The scale factor, scale, set <= 1 to avoid overflow in X.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1: A and B have common or very close eigenvalues; perturbed 
  
                 values were used to solve the equation (but the matrices 
  
                 A and B are unchanged).   

    ===================================================================== 
  


       Decode and Test input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static long int c_false = 0;
    static int c__2 = 2;
    static LONG DOUBLE c_b26 = 1.;
    static LONG DOUBLE c_b30 = 0.;
    static long int c_true = 1;
    
    /* System generated locals */
    int i__1, i__2, 
	    i__3, i__4;
    LONG DOUBLE d__1, d__2;
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
    static int ierr;
    static LONG DOUBLE smin, suml, sumr;
    static int j, k, l;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE x[4]	/* was [2][2] */;
    extern long int lsame_(char *, char *);
    static int knext, lnext, k1, k2, l1, l2;
    static LONG DOUBLE xnorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaln2_(long int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaln2(long int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaln2_(long int *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *
	    , LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *),

#ifdef PETSC_PREFIX_SUFFIX
	     dlasy2_(long int *, long int *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     qlasy2(long int *, long int *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     qlasy2_(long int *, long int *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static LONG DOUBLE a11, db;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *), dlange_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlange_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlange_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static LONG DOUBLE scaloc;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum;
    static long int notrna, notrnb;
    static LONG DOUBLE smlnum, da11, vec[4]	/* was [2][2] */, dum[1], eps,
	     sgn;



#define X(I) x[(I)]
#define WAS(I) was[(I)]
#define DUM(I) dum[(I)]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    notrna = lsame_(trana, "N");
    notrnb = lsame_(tranb, "N");

    *info = 0;
    if (! notrna && ! lsame_(trana, "T") && ! lsame_(trana, "C")) {
	*info = -1;
    } else if (! notrnb && ! lsame_(tranb, "T") && ! lsame_(tranb, 
	    "C")) {
	*info = -2;
    } else if (*isgn != 1 && *isgn != -1) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < MAX(1,*m)) {
	*info = -7;
    } else if (*ldb < MAX(1,*n)) {
	*info = -9;
    } else if (*ldc < MAX(1,*m)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTRSYL", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return;
    }

/*     Set constants to control overflow */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("P");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("P");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("P");
#endif


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

#ifdef PETSC_PREFIX_SUFFIX
    dlabad_(&smlnum, &bignum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlabad(&smlnum, &bignum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlabad_(&smlnum, &bignum);
#endif

    smlnum = smlnum * (LONG DOUBLE) (*m * *n) / eps;
    bignum = 1. / smlnum;

/* Computing MAX */

#ifdef PETSC_PREFIX_SUFFIX
    d__1 = smlnum, d__2 = eps * dlange_("M", m, m, &A(1,1), lda, dum)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    d__1 = smlnum, d__2 = eps * qlange("M", m, m, &A(1,1), lda, dum)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    d__1 = smlnum, d__2 = eps * qlange_("M", m, m, &A(1,1), lda, dum)
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    , d__1 = MAX(d__1,d__2), d__2 = eps * dlange_("M", n, n, &B(1,1), ldb, dum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    , d__1 = MAX(d__1,d__2), d__2 = eps * qlange("M", n, n, &B(1,1), ldb, dum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    , d__1 = MAX(d__1,d__2), d__2 = eps * qlange_("M", n, n, &B(1,1), ldb, dum);
#endif

    smin = MAX(d__1,d__2);

    *scale = 1.;
    sgn = (LONG DOUBLE) (*isgn);

    if (notrna && notrnb) {

/*        Solve    A*X + ISGN*X*B = scale*C.   

          The (K,L)th block of X is determined starting from   
          bottom-left corner column by column by   

           A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)   

          Where   
                    M                         L-1   
          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].   
                  I=K+1                       J=1   

          Start column loop (index = L)   
          L1 (L2) : column index of the first (first) row of X(K,L). 
*/

	lnext = 1;
	i__1 = *n;
	for (l = 1; l <= *n; ++l) {
	    if (l < lnext) {
		goto L60;
	    }
	    if (l == *n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l+1,l) != 0.) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L).
 */

	    knext = *m;
	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L50;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k,k-1) != 0.) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = ABS(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = ABS(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L10: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif

			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L20: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__2 = *m - k1;
/* Computing MIN */
		    i__3 = k1 + 1;
/* Computing MIN */
		    i__4 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l2), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif

			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L30: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l2), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__2 = *m - k2;
/* Computing MIN */
		    i__3 = k2 + 1;
/* Computing MIN */
		    i__4 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l2), &c__1);
		    i__2 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[3] = C(k2,l2) - (suml + sgn * sumr);


#ifdef PETSC_PREFIX_SUFFIX
		    dlasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlasy2(&c_false, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlasy2_(&c_false, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec,
#endif

			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L40: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L50:
		;
	    }

L60:
	    ;
	}

    } else if (! notrna && notrnb) {

/*        Solve    A' *X + ISGN*X*B = scale*C.   

          The (K,L)th block of X is determined starting from   
          upper-left corner column by column by   

            A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)   

          Where   
                     K-1                        L-1   
            R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]   
                     I=1                        J=1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = 1;
	i__1 = *n;
	for (l = 1; l <= *n; ++l) {
	    if (l < lnext) {
		goto L120;
	    }
	    if (l == *n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l+1,l) != 0.) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = 1;
	    i__2 = *m;
	    for (k = 1; k <= *m; ++k) {
		if (k < knext) {
		    goto L110;
		}
		if (k == *m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k+1,k) != 0.) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = ABS(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = ABS(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L70: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif

			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L80: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif

			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L90: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k1,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k2,1), ldc, &B(1,l1), &c__1);
#endif

		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__3 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__3, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__3, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__3, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif

		    i__3 = l1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__3, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__3, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__3, &C(k2,1), ldc, &B(1,l2), &c__1);
#endif

		    vec[3] = C(k2,l2) - (suml + sgn * sumr);


#ifdef PETSC_PREFIX_SUFFIX
		    dlasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlasy2(&c_true, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlasy2_(&c_true, &c_false, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif

			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__3 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L100: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L110:
		;
	    }
L120:
	    ;
	}

    } else if (! notrna && ! notrnb) {

/*        Solve    A'*X + ISGN*X*B' = scale*C.   

          The (K,L)th block of X is determined starting from   
          top-right corner column by column by   

             A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)   

          Where   
                       K-1                          N   
              R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'
].   
                       I=1                        J=L+1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = *n;
	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L180;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l,l-1) != 0.) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = 1;
	    i__1 = *m;
	    for (k = 1; k <= *m; ++k) {
		if (k < knext) {
		    goto L170;
		}
		if (k == *m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k+1,k) != 0.) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l1;
/* Computing MIN */
		    i__3 = l1 + 1;
/* Computing MIN */
		    i__4 = l1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l1+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l1+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l1+1,*n)), ldc, &
#endif

			    B(l1,MIN(l1+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = ABS(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = ABS(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L130: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1, 
#endif

			    &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L140: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif

			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L150: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k1), &c__1, &C(1,l2), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k2), &c__1, &C(1,l1), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__2 = k1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__2, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__2, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__2, &A(1,k2), &c__1, &C(1,l2), &c__1);
#endif

		    i__2 = *n - l2;
/* Computing MIN */
		    i__3 = l2 + 1;
/* Computing MIN */
		    i__4 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__2, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);


#ifdef PETSC_PREFIX_SUFFIX
		    dlasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlasy2(&c_true, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlasy2_(&c_true, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif

			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__2 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L160: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L170:
		;
	    }
L180:
	    ;
	}

    } else if (notrna && ! notrnb) {

/*        Solve    A*X + ISGN*X*B' = scale*C.   

          The (K,L)th block of X is determined starting from   
          bottom-right corner column by column by   

              A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)   

          Where   
                        M                          N   
              R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)']
.   
                      I=K+1                      J=L+1   

          Start column loop (index = L)   
          L1 (L2): column index of the first (last) row of X(K,L) */

	lnext = *n;
	for (l = *n; l >= 1; --l) {
	    if (l > lnext) {
		goto L240;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B(l,l-1) != 0.) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }

/*           Start row loop (index = K)   
             K1 (K2): row index of the first (last) row of X(K,L) 
*/

	    knext = *m;
	    for (k = *m; k >= 1; --k) {
		if (k > knext) {
		    goto L230;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A(k,k-1) != 0.) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l1), &c__1);
		    i__1 = *n - l1;
/* Computing MIN */
		    i__2 = l1 + 1;
/* Computing MIN */
		    i__3 = l1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l1+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l1+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l1+1,*n)), ldc, &
#endif

			    B(l1,MIN(l1+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);
		    scaloc = 1.;

		    a11 = A(k1,k1) + sgn * B(l1,l1);
		    da11 = ABS(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = ABS(vec[0]);
		    if (da11 < 1. && db > 1.) {
			if (db > bignum * da11) {
			    scaloc = 1. / db;
			}
		    }
		    X(0) = vec[0] * scaloc / a11;

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L190: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);

		} else if (l1 == l2 && k1 != k2) {

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    d__1 = -sgn * B(l1,l1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &A(k1,k1), lda, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif

			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L200: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k2,l1) = X(1);

		} else if (l1 != l2 && k1 == k2) {

		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = sgn * (C(k1,l1) - (suml + sgn * sumr))
			    ;

		    i__1 = *m - k1;
/* Computing MIN */
		    i__2 = k1 + 1;
/* Computing MIN */
		    i__3 = k1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k1+1,*m)), lda, &
#endif

			    C(MIN(k1+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[1] = sgn * (C(k1,l2) - (suml + sgn * sumr))
			    ;

		    d__1 = -sgn * A(k1,k1);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__1, &smin, &c_b26, &B(l1,l1), ldb, &c_b26, &c_b26, vec, &c__2, &d__1,
#endif

			     &c_b30, x, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L210: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(1);

		} else if (l1 != l2 && k1 != k2) {

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[0] = C(k1,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k1,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k1,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[2] = C(k1,l2) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l1), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l1,MIN(l2+1,*n)), ldb);
		    vec[1] = C(k2,l1) - (suml + sgn * sumr);

		    i__1 = *m - k2;
/* Computing MIN */
		    i__2 = k2 + 1;
/* Computing MIN */
		    i__3 = k2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    suml = ddot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    suml = qdot(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    suml = qdot_(&i__1, &A(k2,MIN(k2+1,*m)), lda, &
#endif

			    C(MIN(k2+1,*m),l2), &c__1);
		    i__1 = *n - l2;
/* Computing MIN */
		    i__2 = l2 + 1;
/* Computing MIN */
		    i__3 = l2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    sumr = ddot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    sumr = qdot(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    sumr = qdot_(&i__1, &C(k2,MIN(l2+1,*n)), ldc, &
#endif

			    B(l2,MIN(l2+1,*n)), ldb);
		    vec[3] = C(k2,l2) - (suml + sgn * sumr);


#ifdef PETSC_PREFIX_SUFFIX
		    dlasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlasy2(&c_false, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlasy2_(&c_false, &c_true, isgn, &c__2, &c__2, &A(k1,k1), lda, &B(l1,l1), ldb, vec, &
#endif

			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }

		    if (scaloc != 1.) {
			i__1 = *n;
			for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(m, &scaloc, &C(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(m, &scaloc, &C(1,j), &c__1);
#endif

/* L220: */
			}
			*scale *= scaloc;
		    }
		    C(k1,l1) = X(0);
		    C(k1,l2) = X(2);
		    C(k2,l1) = X(1);
		    C(k2,l2) = X(3);
		}

L230:
		;
	    }
L240:
	    ;
	}

    }

    return;

/*     End of DTRSYL */

} /* dtrsyl_ */

