#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgesvx_(char *fact, char *trans, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgesvx(char *fact, char *trans, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgesvx_(char *fact, char *trans, int *n, int *
#endif

	nrhs, LONG DOUBLE *a, int *lda, LONG DOUBLE *af, int *ldaf, 
	int *ipiv, char *equed, LONG DOUBLE *r, LONG DOUBLE *c, LONG DOUBLE *
	b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *rcond, 
	LONG DOUBLE *ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *iwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGESVX uses the LU factorization to compute the solution to a real   
    system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'E', real scaling factors are computed to equilibrate   
       the system:   
          TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B 
  
          TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B 
  
          TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B 
  
       Whether or not the system will be equilibrated depends on the   
       scaling of the matrix A, but if equilibration is used, A is   
       overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') 
  
       or diag(C)*B (if TRANS = 'T' or 'C').   

    2. If FACT = 'N' or 'E', the LU decomposition is used to factor the   
       matrix A (after equilibration if FACT = 'E') as   
          A = P * L * U,   
       where P is a permutation matrix, L is a unit lower triangular   
       matrix, and U is upper triangular.   

    3. The factored form of A is used to estimate the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 4-6 are skipped.   

    4. The system of equations is solved for X using the factored form   
       of A.   

    5. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    6. If equilibration was used, the matrix X is premultiplied by   
       diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so   
       that it solves the original system before equilibration.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of the matrix A is 
  
            supplied on entry, and if not, whether the matrix A should be 
  
            equilibrated before it is factored.   
            = 'F':  On entry, AF and IPIV contain the factored form of A. 
  
                    If EQUED is not 'N', the matrix A has been   
                    equilibrated with scaling factors given by R and C.   
                    A, AF, and IPIV are not modified.   
            = 'N':  The matrix A will be copied to AF and factored.   
            = 'E':  The matrix A will be equilibrated if necessary, then 
  
                    copied to AF and factored.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Transpose)   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is   
            not 'N', then A must have been equilibrated by the scaling   
            factors in R and/or C.  A is not modified if FACT = 'F' or   
            'N', or if FACT = 'E' and EQUED = 'N' on exit.   

            On exit, if EQUED .ne. 'N', A is scaled as follows:   
            EQUED = 'R':  A := diag(R) * A   
            EQUED = 'C':  A := A * diag(C)   
            EQUED = 'B':  A := diag(R) * A * diag(C).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    AF      (input or output) LONG DOUBLE PRECISION array, dimension (LDAF,N) 
  
            If FACT = 'F', then AF is an input argument and on entry   
            contains the factors L and U from the factorization   
            A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then   
            AF is the factored form of the equilibrated matrix A.   

            If FACT = 'N', then AF is an output argument and on exit   
            returns the factors L and U from the factorization A = P*L*U 
  
            of the original matrix A.   

            If FACT = 'E', then AF is an output argument and on exit   
            returns the factors L and U from the factorization A = P*L*U 
  
            of the equilibrated matrix A (see the description of A for   
            the form of the equilibrated matrix).   

    LDAF    (input) INTEGER   
            The leading dimension of the array AF.  LDAF >= MAX(1,N).   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains the pivot indices from the factorization A = P*L*U   
            as computed by DGETRF; row i of the matrix was interchanged   
            with row IPIV(i).   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = P*L*U   
            of the original matrix A.   

            If FACT = 'E', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = P*L*U   
            of the equilibrated matrix A.   

    EQUED   (input or output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration (always true if FACT = 'N').   
            = 'R':  Row equilibration, i.e., A has been premultiplied by 
  
                    diag(R).   
            = 'C':  Column equilibration, i.e., A has been postmultiplied 
  
                    by diag(C).   
            = 'B':  Both row and column equilibration, i.e., A has been   
                    replaced by diag(R) * A * diag(C).   
            EQUED is an input argument if FACT = 'F'; otherwise, it is an 
  
            output argument.   

    R       (input or output) LONG DOUBLE PRECISION array, dimension (N)   
            The row scale factors for A.  If EQUED = 'R' or 'B', A is   
            multiplied on the left by diag(R); if EQUED = 'N' or 'C', R   
            is not accessed.  R is an input argument if FACT = 'F';   
            otherwise, R is an output argument.  If FACT = 'F' and   
            EQUED = 'R' or 'B', each element of R must be positive.   

    C       (input or output) LONG DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.  If EQUED = 'C' or 'B', A is 
  
            multiplied on the right by diag(C); if EQUED = 'N' or 'R', C 
  
            is not accessed.  C is an input argument if FACT = 'F';   
            otherwise, C is an output argument.  If FACT = 'F' and   
            EQUED = 'C' or 'B', each element of C must be positive.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit,   
            if EQUED = 'N', B is not modified;   
            if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by   
            diag(R)*B;   
            if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is   
            overwritten by diag(C)*B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X to the original 
  
            system of equations.  Note that A and B are modified on exit 
  
            if EQUED .ne. 'N', and the solution to the equilibrated   
            system is inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or   
            'B', or inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' 
  
            or 'B'.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= MAX(1,N).   

    RCOND   (output) LONG DOUBLE PRECISION   
            The estimate of the reciprocal condition number of the matrix 
  
            A after equilibration (if done).  If RCOND is less than the   
            machine precision (in particular, if RCOND = 0), the matrix   
            is singular to working precision.  This condition is   
            indicated by a return code of INFO > 0, and the solution and 
  
            error bounds are not computed.   

    FERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (4*N)   
            On exit, WORK(1) contains the reciprocal pivot growth   
            factor norm(A)/norm(U). The "max absolute element" norm is   
            used. If WORK(1) is much less than 1, then the stability   
            of the LU factorization of the (equilibrated) matrix A   
            could be poor. This also means that the solution X, condition 
  
            estimator RCOND, and forward error bound FERR could be   
            unreliable. If factorization fails with 0<INFO<=N, then   
            WORK(1) contains the reciprocal pivot growth factor for the   
            leading INFO columns of A.   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= N:  U(i,i) is exactly zero.  The factorization has   
                         been completed, but the factor U is exactly   
                         singular, so the solution and error bounds   
                         could not be computed.   
                  = N+1: RCOND is less than machine precision.  The   
                         factorization has been completed, but the   
                         matrix is singular to working precision, and   
                         the solution and error bounds have not been   
                         computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static LONG DOUBLE amax;
    static char norm[1];
    static int i, j;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE rcmin, rcmax, anorm;
    static long int equil;

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

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaqge_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqge(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqge_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, char *), dgecon_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, char *), qgecon(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, char *), qgecon_(char *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *, int *);
    static LONG DOUBLE colcnd;
    static long int nofact;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgeequ_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeequ(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeequ_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), dgerfs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), qgerfs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, int *), qgerfs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dgetrf_(int *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgetrf(int *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgetrf_(int *, int *, LONG DOUBLE *, int *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlacpy_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlacpy(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlacpy_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *), xerbla_(char *, 
	    int *);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlantr_(char *, char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlantr(char *, char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlantr_(char *, char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *);
    static int infequ;
    static long int colequ;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgetrs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgetrs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgetrs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    int *);
    static LONG DOUBLE rowcnd;
    static long int notran;
    static LONG DOUBLE smlnum;
    static long int rowequ;
    static LONG DOUBLE rpvgrw;


#define IPIV(I) ipiv[(I)-1]
#define R(I) r[(I)-1]
#define C(I) c[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define AF(I,J) af[(I)-1 + ((J)-1)* ( *ldaf)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    notran = lsame_(trans, "N");
    if (nofact || equil) {
	*(unsigned char *)equed = 'N';
	rowequ = 0;
	colequ = 0;
    } else {
	rowequ = lsame_(equed, "R") || lsame_(equed, "B");
	colequ = lsame_(equed, "C") || lsame_(equed, "B");

#ifdef PETSC_PREFIX_SUFFIX
	smlnum = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
	smlnum = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	smlnum = qlamch_("Safe minimum");
#endif

	bignum = 1. / smlnum;
    }

/*     Test the input parameters. */

    if (! nofact && ! equil && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < MAX(1,*n)) {
	*info = -6;
    } else if (*ldaf < MAX(1,*n)) {
	*info = -8;
    } else if (lsame_(fact, "F") && ! (rowequ || colequ || lsame_(
	    equed, "N"))) {
	*info = -10;
    } else {
	if (rowequ) {
	    rcmin = bignum;
	    rcmax = 0.;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		d__1 = rcmin, d__2 = R(j);
		rcmin = MIN(d__1,d__2);
/* Computing MAX */
		d__1 = rcmax, d__2 = R(j);
		rcmax = MAX(d__1,d__2);
/* L10: */
	    }
	    if (rcmin <= 0.) {
		*info = -11;
	    } else if (*n > 0) {
		rowcnd = MAX(rcmin,smlnum) / MIN(rcmax,bignum);
	    } else {
		rowcnd = 1.;
	    }
	}
	if (colequ && *info == 0) {
	    rcmin = bignum;
	    rcmax = 0.;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		d__1 = rcmin, d__2 = C(j);
		rcmin = MIN(d__1,d__2);
/* Computing MAX */
		d__1 = rcmax, d__2 = C(j);
		rcmax = MAX(d__1,d__2);
/* L20: */
	    }
	    if (rcmin <= 0.) {
		*info = -12;
	    } else if (*n > 0) {
		colcnd = MAX(rcmin,smlnum) / MIN(rcmax,bignum);
	    } else {
		colcnd = 1.;
	    }
	}
	if (*info == 0) {
	    if (*ldb < MAX(1,*n)) {
		*info = -14;
	    } else if (*ldx < MAX(1,*n)) {
		*info = -16;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGESVX", &i__1);
	return;
    }

    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A.
 */


#ifdef PETSC_PREFIX_SUFFIX
	dgeequ_(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgeequ(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgeequ_(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif

		amax, &infequ);
	if (infequ == 0) {

/*           Equilibrate the matrix. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlaqge_(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaqge(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaqge_(n, n, &A(1,1), lda, &R(1), &C(1), &rowcnd, &colcnd, &
#endif

		    amax, equed);
	    rowequ = lsame_(equed, "R") || lsame_(equed, "B");
	    colequ = lsame_(equed, "C") || lsame_(equed, "B");
	}
    }

/*     Scale the right hand side. */

    if (notran) {
	if (rowequ) {
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    B(i,j) = R(i) * B(i,j);
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (colequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		B(i,j) = C(i) * B(i,j);
/* L50: */
	    }
/* L60: */
	}
    }

    if (nofact || equil) {

/*        Compute the LU factorization of A. */


#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("Full", n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("Full", n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("Full", n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dgetrf_(n, n, &AF(1,1), ldaf, &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgetrf(n, n, &AF(1,1), ldaf, &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgetrf_(n, n, &AF(1,1), ldaf, &IPIV(1), info);
#endif


/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {

/*              Compute the reciprocal pivot growth factor of 
the   
                leading rank-deficient INFO columns of A. */


#ifdef PETSC_PREFIX_SUFFIX
		rpvgrw = dlantr_("M", "U", "N", info, info, &AF(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		rpvgrw = qlantr("M", "U", "N", info, info, &AF(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		rpvgrw = qlantr_("M", "U", "N", info, info, &AF(1,1), 
#endif

			ldaf, &WORK(1));
		if (rpvgrw == 0.) {
		    rpvgrw = 1.;
		} else {

#ifdef PETSC_PREFIX_SUFFIX
		    rpvgrw = dlange_("M", n, info, &A(1,1), lda, &WORK(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    rpvgrw = qlange("M", n, info, &A(1,1), lda, &WORK(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    rpvgrw = qlange_("M", n, info, &A(1,1), lda, &WORK(1)
#endif

			    ) / rpvgrw;
		}
		WORK(1) = rpvgrw;
		*rcond = 0.;
	    }
	    return;
	}
    }

/*     Compute the norm of the matrix A and the   
       reciprocal pivot growth factor RPVGRW. */

    if (notran) {
	*(unsigned char *)norm = '1';
    } else {
	*(unsigned char *)norm = 'I';
    }

#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlange_(norm, n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlange(norm, n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlange_(norm, n, n, &A(1,1), lda, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
    rpvgrw = dlantr_("M", "U", "N", n, n, &AF(1,1), ldaf, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    rpvgrw = qlantr("M", "U", "N", n, n, &AF(1,1), ldaf, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    rpvgrw = qlantr_("M", "U", "N", n, n, &AF(1,1), ldaf, &WORK(1));
#endif

    if (rpvgrw == 0.) {
	rpvgrw = 1.;
    } else {

#ifdef PETSC_PREFIX_SUFFIX
	rpvgrw = dlange_("M", n, n, &A(1,1), lda, &WORK(1)) / rpvgrw;
#endif
#ifdef Q_C_PREFIX_SUFFIX
	rpvgrw = qlange("M", n, n, &A(1,1), lda, &WORK(1)) / rpvgrw;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	rpvgrw = qlange_("M", n, n, &A(1,1), lda, &WORK(1)) / rpvgrw;
#endif

    }

/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dgecon_(norm, n, &AF(1,1), ldaf, &anorm, rcond, &WORK(1), &IWORK(1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgecon(norm, n, &AF(1,1), ldaf, &anorm, rcond, &WORK(1), &IWORK(1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgecon_(norm, n, &AF(1,1), ldaf, &anorm, rcond, &WORK(1), &IWORK(1),
#endif

	     info);

/*     Return if the matrix is singular to working precision. */


#ifdef PETSC_PREFIX_SUFFIX
    if (*rcond < dlamch_("Epsilon")) {
#endif
#ifdef Q_C_PREFIX_SUFFIX
    if (*rcond < qlamch("Epsilon")) {
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    if (*rcond < qlamch_("Epsilon")) {
#endif

	WORK(1) = rpvgrw;
	*info = *n + 1;
	return;
    }

/*     Compute the solution matrix X. */


#ifdef PETSC_PREFIX_SUFFIX
    dlacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlacpy("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dgetrs_(trans, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgetrs(trans, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgetrs_(trans, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx,
#endif

	     info);

/*     Use iterative refinement to improve the computed solution and   
       compute error bounds and backward error estimates for it. */


#ifdef PETSC_PREFIX_SUFFIX
    dgerfs_(trans, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgerfs(trans, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgerfs_(trans, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1),
#endif

	     &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(
	    1), &IWORK(1), info);

/*     Transform the solution matrix X to a solution of the original   
       system. */

    if (notran) {
	if (colequ) {
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    X(i,j) = C(i) * X(i,j);
/* L70: */
		}
/* L80: */
	    }
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		FERR(j) /= colcnd;
/* L90: */
	    }
	}
    } else if (rowequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		X(i,j) = R(i) * X(i,j);
/* L100: */
	    }
/* L110: */
	}
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) /= rowcnd;
/* L120: */
	}
    }

    WORK(1) = rpvgrw;
    return;

/*     End of DGESVX */

} /* dgesvx_ */

