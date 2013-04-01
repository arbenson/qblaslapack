#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggglm_(int *n, int *m, int *p, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggglm(int *n, int *m, int *p, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggglm_(int *n, int *m, int *p, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *d, 
	LONG DOUBLE *x, LONG DOUBLE *y, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGGLM solves a general Gauss-Markov linear model (GLM) problem:   

            minimize || y ||_2   subject to   d = A*x + B*y   
                x   

    where A is an N-by-M matrix, B is an N-by-P matrix, and d is a   
    given N-vector. It is assumed that M <= N <= M+P, and   

               rank(A) = M    and    rank( A B ) = N.   

    Under these assumptions, the constrained equation is always   
    consistent, and there is a unique solution x and a minimal 2-norm   
    solution y, which is obtained using a generalized QR factorization   
    of A and B.   

    In particular, if matrix B is square nonsingular, then the problem   
    GLM is equivalent to the following weighted linear least squares   
    problem   

                 minimize || inv(B)*(d-A*x) ||_2   
                     x   

    where inv(B) denotes the inverse of B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of rows of the matrices A and B.  N >= 0.   

    M       (input) INTEGER   
            The number of columns of the matrix A.  0 <= M <= N.   

    P       (input) INTEGER   
            The number of columns of the matrix B.  P >= N-M.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,M)   
            On entry, the N-by-M matrix A.   
            On exit, A is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,P)   
            On entry, the N-by-P matrix B.   
            On exit, B is destroyed.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,N).   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, D is the left hand side of the GLM equation.   
            On exit, D is destroyed.   

    X       (output) LONG DOUBLE PRECISION array, dimension (M)   
    Y       (output) LONG DOUBLE PRECISION array, dimension (P)   
            On exit, X and Y are the solutions of the GLM problem.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,N+M+P).   
            For optimum performance, LWORK >= M+MIN(N,P)+MAX(N,P)*NB,   
            where NB is an upper bound for the optimal blocksizes for   
            DGEQRF, SGERQF, DORMQR and SORMRQ.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================   


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b14 = -1.;
    static LONG DOUBLE c_b16 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int lopt, i;

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
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), dtrsv_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrsv(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrsv_(char *, 
#endif

	    char *, char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int np;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dggqrf_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qggqrf(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qggqrf_(int *, int *, int *, 
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



#define D(I) d[(I)-1]
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    np = MIN(*n,*p);
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0 || *m > *n) {
	*info = -2;
    } else if (*p < 0 || *p < *n - *m) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n + *m + *p;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGGLM", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Compute the GQR factorization of matrices A and B:   

              Q'*A = ( R11 ) M,    Q'*B*Z' = ( T11   T12 ) M   
                     (  0  ) N-M             (  0    T22 ) N-M   
                        M                     M+P-N  N-M   

       where R11 and T22 are upper triangular, and Q and Z are   
       orthogonal. */

    i__1 = *lwork - *m - np;

#ifdef PETSC_PREFIX_SUFFIX
    dggqrf_(n, m, p, &A(1,1), lda, &WORK(1), &B(1,1), ldb, &WORK(*m 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qggqrf(n, m, p, &A(1,1), lda, &WORK(1), &B(1,1), ldb, &WORK(*m 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qggqrf_(n, m, p, &A(1,1), lda, &WORK(1), &B(1,1), ldb, &WORK(*m 
#endif

	    + 1), &WORK(*m + np + 1), &i__1, info);
    lopt = (int) WORK(*m + np + 1);

/*     Update left-hand-side vector d = Q'*d = ( d1 ) M   
                                               ( d2 ) N-M */

    i__1 = MAX(1,*n);
    i__2 = *lwork - *m - np;

#ifdef PETSC_PREFIX_SUFFIX
    dormqr_("Left", "Transpose", n, &c__1, m, &A(1,1), lda, &WORK(1), &D(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormqr("Left", "Transpose", n, &c__1, m, &A(1,1), lda, &WORK(1), &D(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormqr_("Left", "Transpose", n, &c__1, m, &A(1,1), lda, &WORK(1), &D(
#endif

	    1), &i__1, &WORK(*m + np + 1), &i__2, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(*m + np + 1);
    lopt = MAX(i__1,i__2);

/*     Solve T22*y2 = d2 for y2 */

    i__1 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
    dtrsv_("Upper", "No transpose", "Non unit", &i__1, &B(*m+1,*m+*p-*n+1), ldb, &D(*m + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrsv("Upper", "No transpose", "Non unit", &i__1, &B(*m+1,*m+*p-*n+1), ldb, &D(*m + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrsv_("Upper", "No transpose", "Non unit", &i__1, &B(*m+1,*m+*p-*n+1), ldb, &D(*m + 1), &c__1);
#endif

    i__1 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(&i__1, &D(*m + 1), &c__1, &Y(*m + *p - *n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(&i__1, &D(*m + 1), &c__1, &Y(*m + *p - *n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(&i__1, &D(*m + 1), &c__1, &Y(*m + *p - *n + 1), &c__1);
#endif


/*     Set y1 = 0 */

    i__1 = *m + *p - *n;
    for (i = 1; i <= *m+*p-*n; ++i) {
	Y(i) = 0.;
/* L10: */
    }

/*     Update d1 = d1 - T12*y2 */

    i__1 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
    dgemv_("No transpose", m, &i__1, &c_b14, &B(1,*m+*p-*n+1), ldb, &Y(*m + *p - *n + 1), &c__1, &c_b16, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgemv("No transpose", m, &i__1, &c_b14, &B(1,*m+*p-*n+1), ldb, &Y(*m + *p - *n + 1), &c__1, &c_b16, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgemv_("No transpose", m, &i__1, &c_b14, &B(1,*m+*p-*n+1), ldb, &Y(*m + *p - *n + 1), &c__1, &c_b16, &D(1), &c__1);
#endif


/*     Solve triangular system: R11*x = d1 */


#ifdef PETSC_PREFIX_SUFFIX
    dtrsv_("Upper", "No Transpose", "Non unit", m, &A(1,1), lda, &D(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrsv("Upper", "No Transpose", "Non unit", m, &A(1,1), lda, &D(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrsv_("Upper", "No Transpose", "Non unit", m, &A(1,1), lda, &D(1), &
#endif

	    c__1);

/*     Copy D to X */


#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(m, &D(1), &c__1, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(m, &D(1), &c__1, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(m, &D(1), &c__1, &X(1), &c__1);
#endif


/*     Backward transformation y = Z'*y   

   Computing MAX */
    i__1 = 1, i__2 = *n - *p + 1;
    i__3 = MAX(1,*p);
    i__4 = *lwork - *m - np;

#ifdef PETSC_PREFIX_SUFFIX
    dormrq_("Left", "Transpose", p, &c__1, &np, &B(MAX(1,*n-*p+1),1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormrq("Left", "Transpose", p, &c__1, &np, &B(MAX(1,*n-*p+1),1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormrq_("Left", "Transpose", p, &c__1, &np, &B(MAX(1,*n-*p+1),1), 
#endif

	    ldb, &WORK(*m + 1), &Y(1), &i__3, &WORK(*m + np + 1), &i__4, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(*m + np + 1);
    WORK(1) = (LONG DOUBLE) MAX(i__1,i__2);

    return;

/*     End of DGGGLM */

} /* dggglm_ */

