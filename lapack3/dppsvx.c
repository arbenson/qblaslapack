#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dppsvx_(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qppsvx(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qppsvx_(char *fact, char *uplo, int *n, int *
#endif

	nrhs, LONG DOUBLE *ap, LONG DOUBLE *afp, char *equed, LONG DOUBLE *s, 
	LONG DOUBLE *b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *
	rcond, LONG DOUBLE *ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *
	iwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPPSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to   
    compute the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric positive definite matrix stored in   
    packed format and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'E', real scaling factors are computed to equilibrate   
       the system:   
          diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B   
       Whether or not the system will be equilibrated depends on the   
       scaling of the matrix A, but if equilibration is used, A is   
       overwritten by diag(S)*A*diag(S) and B by diag(S)*B.   

    2. If FACT = 'N' or 'E', the Cholesky decomposition is used to   
       factor the matrix A (after equilibration if FACT = 'E') as   
          A = U**T* U,  if UPLO = 'U', or   
          A = L * L**T,  if UPLO = 'L',   
       where U is an upper triangular matrix and L is a lower triangular 
  
       matrix.   

    3. The factored form of A is used to estimate the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 4-6 are skipped.   

    4. The system of equations is solved for X using the factored form   
       of A.   

    5. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    6. If equilibration was used, the matrix X is premultiplied by   
       diag(S) so that it solves the original system before   
       equilibration.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of the matrix A is 
  
            supplied on entry, and if not, whether the matrix A should be 
  
            equilibrated before it is factored.   
            = 'F':  On entry, AFP contains the factored form of A.   
                    If EQUED = 'Y', the matrix A has been equilibrated   
                    with scaling factors given by S.  AP and AFP will not 
  
                    be modified.   
            = 'N':  The matrix A will be copied to AFP and factored.   
            = 'E':  The matrix A will be equilibrated if necessary, then 
  
                    copied to AFP and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array, except if FACT = 'F' 
  
            and EQUED = 'Y', then A must contain the equilibrated matrix 
  
            diag(S)*A*diag(S).  The j-th column of A is stored in the   
            array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.  A is not modified if   
            FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. 
  

            On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by   
            diag(S)*A*diag(S).   

    AFP     (input or output) LONG DOUBLE PRECISION array, dimension   
                              (N*(N+1)/2)   
            If FACT = 'F', then AFP is an input argument and on entry   
            contains the triangular factor U or L from the Cholesky   
            factorization A = U'*U or A = L*L', in the same storage   
            format as A.  If EQUED .ne. 'N', then AFP is the factored   
            form of the equilibrated matrix A.   

            If FACT = 'N', then AFP is an output argument and on exit   
            returns the triangular factor U or L from the Cholesky   
            factorization A = U'*U or A = L*L' of the original matrix A. 
  

            If FACT = 'E', then AFP is an output argument and on exit   
            returns the triangular factor U or L from the Cholesky   
            factorization A = U'*U or A = L*L' of the equilibrated   
            matrix A (see the description of AP for the form of the   
            equilibrated matrix).   

    EQUED   (input or output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration (always true if FACT = 'N').   
            = 'Y':  Equilibration was done, i.e., A has been replaced by 
  
                    diag(S) * A * diag(S).   
            EQUED is an input argument if FACT = 'F'; otherwise, it is an 
  
            output argument.   

    S       (input or output) LONG DOUBLE PRECISION array, dimension (N)   
            The scale factors for A; not accessed if EQUED = 'N'.  S is   
            an input argument if FACT = 'F'; otherwise, S is an output   
            argument.  If FACT = 'F' and EQUED = 'Y', each element of S   
            must be positive.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y',   
            B is overwritten by diag(S) * B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X to the original 
  
            system of equations.  Note that if EQUED = 'Y', A and B are   
            modified on exit, and the solution to the equilibrated system 
  
            is inv(diag(S))*X.   

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

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= N: the leading minor of order i of A   
                        is not positive definite, so the factorization   
                        could not be completed, and the solution and   
                        error bounds could not be computed.   
                  = N+1: RCOND is less than machine precision.  The   
                        factorization has been completed, but the matrix 
  
                        is singular to working precision, and the   
                        solution and error bounds have not been computed. 
  

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = conjg(aji))   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static LONG DOUBLE amax, smin, smax;
    static int i, j;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE scond, anorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static long int equil, rcequ;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static long int nofact;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlansp_(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlansp(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlansp_(char *, char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dppcon_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qppcon(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qppcon_(char *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), dlaqsp_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlaqsp(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlaqsp_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, char *);
    static int infequ;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dppequ_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qppequ(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qppequ_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dpprfs_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qpprfs(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qpprfs_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dpptrf_(char *, int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qpptrf(char *, int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qpptrf_(char *, int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE smlnum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dpptrs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qpptrs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qpptrs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);



#define AP(I) ap[(I)-1]
#define AFP(I) afp[(I)-1]
#define S(I) s[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    if (nofact || equil) {
	*(unsigned char *)equed = 'N';
	rcequ = 0;
    } else {
	rcequ = lsame_(equed, "Y");

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
    } else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (lsame_(fact, "F") && ! (rcequ || lsame_(equed, "N"))) {
	*info = -7;
    } else {
	if (rcequ) {
	    smin = bignum;
	    smax = 0.;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		d__1 = smin, d__2 = S(j);
		smin = MIN(d__1,d__2);
/* Computing MAX */
		d__1 = smax, d__2 = S(j);
		smax = MAX(d__1,d__2);
/* L10: */
	    }
	    if (smin <= 0.) {
		*info = -8;
	    } else if (*n > 0) {
		scond = MAX(smin,smlnum) / MIN(smax,bignum);
	    } else {
		scond = 1.;
	    }
	}
	if (*info == 0) {
	    if (*ldb < MAX(1,*n)) {
		*info = -10;
	    } else if (*ldx < MAX(1,*n)) {
		*info = -12;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPPSVX", &i__1);
	return;
    }

    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A.
 */


#ifdef PETSC_PREFIX_SUFFIX
	dppequ_(uplo, n, &AP(1), &S(1), &scond, &amax, &infequ);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qppequ(uplo, n, &AP(1), &S(1), &scond, &amax, &infequ);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qppequ_(uplo, n, &AP(1), &S(1), &scond, &amax, &infequ);
#endif

	if (infequ == 0) {

/*           Equilibrate the matrix. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlaqsp_(uplo, n, &AP(1), &S(1), &scond, &amax, equed);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaqsp(uplo, n, &AP(1), &S(1), &scond, &amax, equed);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaqsp_(uplo, n, &AP(1), &S(1), &scond, &amax, equed);
#endif

	    rcequ = lsame_(equed, "Y");
	}
    }

/*     Scale the right-hand side. */

    if (rcequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		B(i,j) = S(i) * B(i,j);
/* L20: */
	    }
/* L30: */
	}
    }

    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U'*U or A = L*L'. */

	i__1 = *n * (*n + 1) / 2;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dpptrf_(uplo, n, &AFP(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpptrf(uplo, n, &AFP(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpptrf_(uplo, n, &AFP(1), info);
#endif


/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {
		*rcond = 0.;
	    }
	    return;
	}
    }

/*     Compute the norm of the matrix A. */


#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlansp_("I", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlansp("I", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlansp_("I", uplo, n, &AP(1), &WORK(1));
#endif


/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dppcon_(uplo, n, &AFP(1), &anorm, rcond, &WORK(1), &IWORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qppcon(uplo, n, &AFP(1), &anorm, rcond, &WORK(1), &IWORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qppcon_(uplo, n, &AFP(1), &anorm, rcond, &WORK(1), &IWORK(1), info);
#endif


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
    dpptrs_(uplo, n, nrhs, &AFP(1), &X(1,1), ldx, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpptrs(uplo, n, nrhs, &AFP(1), &X(1,1), ldx, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpptrs_(uplo, n, nrhs, &AFP(1), &X(1,1), ldx, info);
#endif


/*     Use iterative refinement to improve the computed solution and   
       compute error bounds and backward error estimates for it. */


#ifdef PETSC_PREFIX_SUFFIX
    dpprfs_(uplo, n, nrhs, &AP(1), &AFP(1), &B(1,1), ldb, &X(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpprfs(uplo, n, nrhs, &AP(1), &AFP(1), &B(1,1), ldb, &X(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpprfs_(uplo, n, nrhs, &AP(1), &AFP(1), &B(1,1), ldb, &X(1,1), 
#endif

	    ldx, &FERR(1), &BERR(1), &WORK(1), &IWORK(1), info);

/*     Transform the solution matrix X to a solution of the original   
       system. */

    if (rcequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__2 = *n;
	    for (i = 1; i <= *n; ++i) {
		X(i,j) = S(i) * X(i,j);
/* L40: */
	    }
/* L50: */
	}
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) /= scond;
/* L60: */
	}
    }

    return;

/*     End of DPPSVX */

} /* dppsvx_ */

