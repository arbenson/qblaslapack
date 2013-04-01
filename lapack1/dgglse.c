#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgglse_(int *m, int *n, int *p, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgglse(int *m, int *n, int *p, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgglse_(int *m, int *n, int *p, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *c, 
	LONG DOUBLE *d, LONG DOUBLE *x, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGLSE solves the linear equality-constrained least squares (LSE)   
    problem:   

            minimize || c - A*x ||_2   subject to   B*x = d   

    where A is an M-by-N matrix, B is a P-by-N matrix, c is a given   
    M-vector, and d is a given P-vector. It is assumed that   
    P <= N <= M+P, and   

             rank(B) = P and  rank( ( A ) ) = N.   
                                  ( ( B ) )   

    These conditions ensure that the LSE problem has a unique solution,   
    which is obtained using a GRQ factorization of the matrices B and A. 
  

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B. N >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B. 0 <= P <= N <= M+P.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, B is destroyed.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,P).   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (M)   
            On entry, C contains the right hand side vector for the   
            least squares part of the LSE problem.   
            On exit, the residual sum of squares for the solution   
            is given by the sum of squares of elements N-P+1 to M of   
            vector C.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (P)   
            On entry, D contains the right hand side vector for the   
            constrained equation.   
            On exit, D is destroyed.   

    X       (output) LONG DOUBLE PRECISION array, dimension (N)   
            On exit, X is the solution of the LSE problem.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,M+N+P).   
            For optimum performance LWORK >= P+MIN(M,N)+MAX(M,N)*NB,   
            where NB is an upper bound for the optimal blocksizes for   
            DGEQRF, SGERQF, DORMQR and SORMRQ.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b11 = -1.;
    static LONG DOUBLE c_b13 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int lopt;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dcopy_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy_(int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), daxpy_(int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy(int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy_(int 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *)

#ifdef PETSC_PREFIX_SUFFIX
	    , dtrmv_(char *, char *, char *, int *, LONG DOUBLE *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    , qtrmv(char *, char *, char *, int *, LONG DOUBLE *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    , qtrmv_(char *, char *, char *, int *, LONG DOUBLE *, int 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *), dtrsv_(char *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *), qtrsv(char *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *), qtrsv_(char *
#endif

	    , char *, char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *
	    , int *);
    static int mn, nr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dggrqf_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qggrqf(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qggrqf_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), xerbla_(char *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *), dormqr_(char *, char *, int *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *), qormqr(char *, char *, int *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *), qormqr_(char *, char *, int *, int *,
#endif

	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dormrq_(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormrq(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormrq_(char *, char *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, int *);



#define C(I) c[(I)-1]
#define D(I) d[(I)-1]
#define X(I) x[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    mn = MIN(*m,*n);
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*p < 0 || *p > *n || *p < *n - *m) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*ldb < MAX(1,*p)) {
	*info = -7;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m + *n + *p;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGLSE", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Compute the GRQ factorization of matrices B and A:   

              B*Q' = (  0  T12 ) P   Z'*A*Q' = ( R11 R12 ) N-P   
                       N-P  P                  (  0  R22 ) M+P-N   
                                                 N-P  P   

       where T12 and R11 are upper triangular, and Q and Z are   
       orthogonal. */

    i__1 = *lwork - *p - mn;

#ifdef PETSC_PREFIX_SUFFIX
    dggrqf_(p, m, n, &B(1,1), ldb, &WORK(1), &A(1,1), lda, &WORK(*p 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qggrqf(p, m, n, &B(1,1), ldb, &WORK(1), &A(1,1), lda, &WORK(*p 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qggrqf_(p, m, n, &B(1,1), ldb, &WORK(1), &A(1,1), lda, &WORK(*p 
#endif

	    + 1), &WORK(*p + mn + 1), &i__1, info);
    lopt = (int) WORK(*p + mn + 1);

/*     Update c = Z'*c = ( c1 ) N-P   
                         ( c2 ) M+P-N */

    i__1 = MAX(1,*m);
    i__2 = *lwork - *p - mn;

#ifdef PETSC_PREFIX_SUFFIX
    dormqr_("Left", "Transpose", m, &c__1, &mn, &A(1,1), lda, &WORK(*p + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormqr("Left", "Transpose", m, &c__1, &mn, &A(1,1), lda, &WORK(*p + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormqr_("Left", "Transpose", m, &c__1, &mn, &A(1,1), lda, &WORK(*p + 
#endif

	    1), &C(1), &i__1, &WORK(*p + mn + 1), &i__2, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(*p + mn + 1);
    lopt = MAX(i__1,i__2);

/*     Solve T12*x2 = d for x2 */


#ifdef PETSC_PREFIX_SUFFIX
    dtrsv_("Upper", "No transpose", "Non unit", p, &B(1,*n-*p+1), ldb, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrsv("Upper", "No transpose", "Non unit", p, &B(1,*n-*p+1), ldb, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrsv_("Upper", "No transpose", "Non unit", p, &B(1,*n-*p+1), ldb, &D(1), &c__1);
#endif


/*     Update c1 */

    i__1 = *n - *p;

#ifdef PETSC_PREFIX_SUFFIX
    dgemv_("No transpose", &i__1, p, &c_b11, &A(1,*n-*p+1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgemv("No transpose", &i__1, p, &c_b11, &A(1,*n-*p+1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgemv_("No transpose", &i__1, p, &c_b11, &A(1,*n-*p+1), 
#endif

	    lda, &D(1), &c__1, &c_b13, &C(1), &c__1);

/*     Sovle R11*x1 = c1 for x1 */

    i__1 = *n - *p;

#ifdef PETSC_PREFIX_SUFFIX
    dtrsv_("Upper", "No transpose", "Non unit", &i__1, &A(1,1), lda, &C(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrsv("Upper", "No transpose", "Non unit", &i__1, &A(1,1), lda, &C(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrsv_("Upper", "No transpose", "Non unit", &i__1, &A(1,1), lda, &C(
#endif

	    1), &c__1);

/*     Put the solutions in X */

    i__1 = *n - *p;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&i__1, &C(1), &c__1, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&i__1, &C(1), &c__1, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&i__1, &C(1), &c__1, &X(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(p, &D(1), &c__1, &X(*n - *p + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(p, &D(1), &c__1, &X(*n - *p + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(p, &D(1), &c__1, &X(*n - *p + 1), &c__1);
#endif


/*     Compute the residual vector: */

    if (*m < *n) {
	nr = *m + *p - *n;
	i__1 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("No transpose", &nr, &i__1, &c_b11, &A(*n-*p+1,*m+1), lda, &D(nr + 1), &c__1, &c_b13, &C(*n - *p + 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("No transpose", &nr, &i__1, &c_b11, &A(*n-*p+1,*m+1), lda, &D(nr + 1), &c__1, &c_b13, &C(*n - *p + 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("No transpose", &nr, &i__1, &c_b11, &A(*n-*p+1,*m+1), lda, &D(nr + 1), &c__1, &c_b13, &C(*n - *p + 1), &
#endif

		c__1);
    } else {
	nr = *p;
    }

#ifdef PETSC_PREFIX_SUFFIX
    dtrmv_("Upper", "No transpose", "Non unit", &nr, &A(*n-*p+1,*n-*p+1), lda, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrmv("Upper", "No transpose", "Non unit", &nr, &A(*n-*p+1,*n-*p+1), lda, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrmv_("Upper", "No transpose", "Non unit", &nr, &A(*n-*p+1,*n-*p+1), lda, &D(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(&nr, &c_b11, &D(1), &c__1, &C(*n - *p + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(&nr, &c_b11, &D(1), &c__1, &C(*n - *p + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(&nr, &c_b11, &D(1), &c__1, &C(*n - *p + 1), &c__1);
#endif


/*     Backward transformation x = Q'*x */

    i__1 = *lwork - *p - mn;

#ifdef PETSC_PREFIX_SUFFIX
    dormrq_("Left", "Transpose", n, &c__1, p, &B(1,1), ldb, &WORK(1), &X(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormrq("Left", "Transpose", n, &c__1, p, &B(1,1), ldb, &WORK(1), &X(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormrq_("Left", "Transpose", n, &c__1, p, &B(1,1), ldb, &WORK(1), &X(
#endif

	    1), n, &WORK(*p + mn + 1), &i__1, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(*p + mn + 1);
    WORK(1) = (LONG DOUBLE) (*p + mn + MAX(i__1,i__2));

    return;

/*     End of DGGLSE */

} /* dgglse_ */

