#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsygv_(int *itype, char *jobz, char *uplo, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsygv(int *itype, char *jobz, char *uplo, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsygv_(int *itype, char *jobz, char *uplo, int *
#endif

	n, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *w, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYGV computes all the eigenvalues, and optionally, the eigenvectors 
  
    of a real generalized symmetric-definite eigenproblem, of the form   
    A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.   
    Here A and B are assumed to be symmetric and B is also   
    positive definite.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            Specifies the problem type to be solved:   
            = 1:  A*x = (lambda)*B*x   
            = 2:  A*B*x = (lambda)*x   
            = 3:  B*A*x = (lambda)*x   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangles of A and B are stored;   
            = 'L':  Lower triangles of A and B are stored.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   

            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            matrix Z of eigenvectors.  The eigenvectors are normalized   
            as follows:   
            if ITYPE = 1 or 2, Z**T*B*Z = I;   
            if ITYPE = 3, Z**T*inv(B)*Z = I.   
            If JOBZ = 'N', then on exit the upper triangle (if UPLO='U') 
  
            or the lower triangle (if UPLO='L') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the symmetric matrix B.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of B contains the   
            upper triangular part of the matrix B.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of B contains   
            the lower triangular part of the matrix B.   

            On exit, if INFO <= N, the part of B containing the matrix is 
  
            overwritten by the triangular factor U or L from the Cholesky 
  
            factorization B = U**T*U or B = L*L**T.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    W       (output) LONG DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= MAX(1,3*N-1).   
            For optimal efficiency, LWORK >= (NB+2)*N,   
            where NB is the blocksize for DSYTRD returned by ILAENV.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  DPOTRF or DSYEV returned an error code:   
               <= N:  if INFO = i, DSYEV failed to converge;   
                      i off-diagonal elements of an intermediate   
                      tridiagonal form did not converge to zero;   
               > N:   if INFO = N + i, for 1 <= i <= N, then the leading 
  
                      minor of order i of B is not positive definite.   
                      The factorization of B could not be completed and   
                      no eigenvalues or eigenvectors were computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b11 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int neig;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrmm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static char trans[1];

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrsm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyev_(char *, char *, int *, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyev(char *, char *, int *, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyev_(char *, char *, int *, LONG DOUBLE *
#endif

	    , int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static long int wantz;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dpotrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpotrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpotrf_(
#endif

	    char *, int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dsygst_(int *, char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsygst(int *, char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsygst_(int *, char *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, int *);



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    wantz = lsame_(jobz, "V");
    upper = lsame_(uplo, "U");

    *info = 0;
    if (*itype < 0 || *itype > 3) {
	*info = -1;
    } else if (! (wantz || lsame_(jobz, "N"))) {
	*info = -2;
    } else if (! (upper || lsame_(uplo, "L"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < MAX(1,*n)) {
	*info = -6;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYGV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Form a Cholesky factorization of B. */


#ifdef PETSC_PREFIX_SUFFIX
    dpotrf_(uplo, n, &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpotrf(uplo, n, &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpotrf_(uplo, n, &B(1,1), ldb, info);
#endif

    if (*info != 0) {
	*info = *n + *info;
	return;
    }

/*     Transform problem to standard eigenvalue problem and solve. */


#ifdef PETSC_PREFIX_SUFFIX
    dsygst_(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsygst(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsygst_(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dsyev_(jobz, uplo, n, &A(1,1), lda, &W(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsyev(jobz, uplo, n, &A(1,1), lda, &W(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsyev_(jobz, uplo, n, &A(1,1), lda, &W(1), &WORK(1), lwork, info);
#endif


    if (wantz) {

/*        Backtransform eigenvectors to the original problem. */

	neig = *n;
	if (*info > 0) {
	    neig = *info - 1;
	}
	if (*itype == 1 || *itype == 2) {

/*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;   
             backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
 */

	    if (upper) {
		*(unsigned char *)trans = 'N';
	    } else {
		*(unsigned char *)trans = 'T';
	    }


#ifdef PETSC_PREFIX_SUFFIX
	    dtrsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrsm("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrsm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif


	} else if (*itype == 3) {

/*           For B*A*x=(lambda)*x;   
             backtransform eigenvectors: x = L*y or U'*y */

	    if (upper) {
		*(unsigned char *)trans = 'T';
	    } else {
		*(unsigned char *)trans = 'N';
	    }


#ifdef PETSC_PREFIX_SUFFIX
	    dtrmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrmm("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrmm_("Left", uplo, trans, "Non-unit", n, &neig, &c_b11, &B(1,1), ldb, &A(1,1), lda);
#endif

	}
    }
    return;

/*     End of DSYGV */

} /* dsygv_ */

