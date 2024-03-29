#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbsvx_(char *fact, char *trans, int *n, int *kl,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbsvx(char *fact, char *trans, int *n, int *kl,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbsvx_(char *fact, char *trans, int *n, int *kl,
#endif

	 int *ku, int *nrhs, LONG DOUBLE *ab, int *ldab, 
	LONG DOUBLE *afb, int *ldafb, int *ipiv, char *equed, 
	LONG DOUBLE *r, LONG DOUBLE *c, LONG DOUBLE *b, int *ldb, LONG DOUBLE 
	*x, int *ldx, LONG DOUBLE *rcond, LONG DOUBLE *ferr, LONG DOUBLE *
	berr, LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGBSVX uses the LU factorization to compute the solution to a real   
    system of linear equations A * X = B, A**T * X = B, or A**H * X = B, 
  
    where A is a band matrix of order N with KL subdiagonals and KU   
    superdiagonals, and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed by this subroutine:   

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
          A = L * U,   
       where L is a product of permutation and unit lower triangular   
       matrices with KL subdiagonals, and U is upper triangular with   
       KL+KU superdiagonals.   

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
            = 'F':  On entry, AFB and IPIV contain the factored form of   
                    A.  If EQUED is not 'N', the matrix A has been   
                    equilibrated with scaling factors given by R and C.   
                    AB, AFB, and IPIV are not modified.   
            = 'N':  The matrix A will be copied to AFB and factored.   
            = 'E':  The matrix A will be equilibrated if necessary, then 
  
                    copied to AFB and factored.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations.   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Transpose)   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the matrix A in band storage, in rows 1 to KL+KU+1. 
  
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(KU+1+i-j,j) = A(i,j) for MAX(1,j-KU)<=i<=MIN(N,j+kl)   

            If FACT = 'F' and EQUED is not 'N', then A must have been   
            equilibrated by the scaling factors in R and/or C.  AB is not 
  
            modified if FACT = 'F' or 'N', or if FACT = 'E' and   
            EQUED = 'N' on exit.   

            On exit, if EQUED .ne. 'N', A is scaled as follows:   
            EQUED = 'R':  A := diag(R) * A   
            EQUED = 'C':  A := A * diag(C)   
            EQUED = 'B':  A := diag(R) * A * diag(C).   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KL+KU+1.   

    AFB     (input or output) LONG DOUBLE PRECISION array, dimension (LDAFB,N) 
  
            If FACT = 'F', then AFB is an input argument and on entry   
            contains details of the LU factorization of the band matrix   
            A, as computed by DGBTRF.  U is stored as an upper triangular 
  
            band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,   
            and the multipliers used during the factorization are stored 
  
            in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is 
  
            the factored form of the equilibrated matrix A.   

            If FACT = 'N', then AFB is an output argument and on exit   
            returns details of the LU factorization of A.   

            If FACT = 'E', then AFB is an output argument and on exit   
            returns details of the LU factorization of the equilibrated   
            matrix A (see the description of AB for the form of the   
            equilibrated matrix).   

    LDAFB   (input) INTEGER   
            The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. 
  

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains the pivot indices from the factorization A = L*U   
            as computed by DGBTRF; row i of the matrix was interchanged   
            with row IPIV(i).   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = L*U   
            of the original matrix A.   

            If FACT = 'E', then IPIV is an output argument and on exit   
            contains the pivot indices from the factorization A = L*U   
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
            On entry, the right hand side matrix B.   
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

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (3*N)   
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
                  <= N:  U(i,i) is exactly zero.  The factorization   
                         has been completed, but the factor U is exactly 
  
                         singular, so the solution and error bounds   
                         could not be computed.   
                  = N+1: RCOND is less than machine precision.  The   
                         factorization has been completed, but the   
                         matrix A is singular to working precision, and   
                         the solution and error bounds have not been   
                         computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static LONG DOUBLE amax;
    static char norm[1];
    static int i, j;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE rcmin, rcmax, anorm;

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
    static long int equil;
    static int P_j1, j2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *), dlangb_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlangb_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlangb_(char *, int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaqgb_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqgb(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaqgb_(int *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, char *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dgbcon_(char *, int *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgbcon(char *, int *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgbcon_(char *, int *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, int *);
    static LONG DOUBLE colcnd;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlantb_(char *, char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlantb(char *, char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlantb_(char *, char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgbequ_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbequ(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbequ_(int *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), dgbrfs_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qgbrfs(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qgbrfs_(
#endif

	    char *, int *, int *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), dgbtrf_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qgbtrf(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qgbtrf_(int *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, int 
	    *, int *);
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
    extern /* Subroutine */ void dgbtrs_(char *, int *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtrs(char *, int *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtrs_(char *, int *, int *, int 
#endif

	    *, int *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, int *);
    static int infequ;
    static long int colequ;
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

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define AFB(I,J) afb[(I)-1 + ((J)-1)* ( *ldafb)]
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
    } else if (*kl < 0) {
	*info = -4;
    } else if (*ku < 0) {
	*info = -5;
    } else if (*nrhs < 0) {
	*info = -6;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -8;
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
	*info = -10;
    } else if (lsame_(fact, "F") && ! (rowequ || colequ || lsame_(
	    equed, "N"))) {
	*info = -12;
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
		*info = -13;
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
		*info = -14;
	    } else if (*n > 0) {
		colcnd = MAX(rcmin,smlnum) / MIN(rcmax,bignum);
	    } else {
		colcnd = 1.;
	    }
	}
	if (*info == 0) {
	    if (*ldb < MAX(1,*n)) {
		*info = -16;
	    } else if (*ldx < MAX(1,*n)) {
		*info = -18;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBSVX", &i__1);
	return;
    }

    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A.
 */


#ifdef PETSC_PREFIX_SUFFIX
	dgbequ_(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgbequ(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgbequ_(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd, &
#endif

		colcnd, &amax, &infequ);
	if (infequ == 0) {

/*           Equilibrate the matrix. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlaqgb_(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaqgb(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaqgb_(n, n, kl, ku, &AB(1,1), ldab, &R(1), &C(1), &rowcnd,
#endif

		     &colcnd, &amax, equed);
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

/*        Compute the LU factorization of the band matrix A. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__2 = j - *ku;
	    P_j1 = MAX(i__2,1);
/* Computing MIN */
	    i__2 = j + *kl;
	    j2 = MIN(i__2,*n);
	    i__2 = j2 - P_j1 + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__2, &AB(*ku+1-j+P_j1,j), &c__1, &AFB(*kl+*ku+1-j+P_j1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__2, &AB(*ku+1-j+P_j1,j), &c__1, &AFB(*kl+*ku+1-j+P_j1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__2, &AB(*ku+1-j+P_j1,j), &c__1, &AFB(*kl+*ku+1-j+P_j1,j), &c__1);
#endif

/* L70: */
	}


#ifdef PETSC_PREFIX_SUFFIX
	dgbtrf_(n, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgbtrf(n, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgbtrf_(n, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), info);
#endif


/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {

/*              Compute the reciprocal pivot growth factor of 
the   
                leading rank-deficient INFO columns of A. */

		anorm = 0.;
		i__1 = *info;
		for (j = 1; j <= *info; ++j) {
/* Computing MAX */
		    i__2 = *ku + 2 - j;
/* Computing MIN */
		    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
		    i__3 = MIN(i__4,i__5);
		    for (i = MAX(*ku+2-j,1); i <= MIN(*n+*ku+1-j,*kl+*ku+1); ++i) {
/* Computing MAX */
			d__2 = anorm, d__3 = (d__1 = AB(i,j), ABS(
				d__1));
			anorm = MAX(d__2,d__3);
/* L80: */
		    }
/* L90: */
		}
/* Computing MIN */
		i__3 = *info - 1, i__2 = *kl + *ku;
		i__1 = MIN(i__3,i__2);
/* Computing MAX */
		i__4 = 1, i__5 = *kl + *ku + 2 - *info;

#ifdef PETSC_PREFIX_SUFFIX
		rpvgrw = dlantb_("M", "U", "N", info, &i__1, &AFB(MAX(1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		rpvgrw = qlantb("M", "U", "N", info, &i__1, &AFB(MAX(1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		rpvgrw = qlantb_("M", "U", "N", info, &i__1, &AFB(MAX(1,
#endif

			*kl+*ku+2-*info),1), ldafb, &WORK(1));
		if (rpvgrw == 0.) {
		    rpvgrw = 1.;
		} else {
		    rpvgrw = anorm / rpvgrw;
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
    anorm = dlangb_(norm, n, kl, ku, &AB(1,1), ldab, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlangb(norm, n, kl, ku, &AB(1,1), ldab, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlangb_(norm, n, kl, ku, &AB(1,1), ldab, &WORK(1));
#endif

    i__1 = *kl + *ku;

#ifdef PETSC_PREFIX_SUFFIX
    rpvgrw = dlantb_("M", "U", "N", n, &i__1, &AFB(1,1), ldafb, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    rpvgrw = qlantb("M", "U", "N", n, &i__1, &AFB(1,1), ldafb, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    rpvgrw = qlantb_("M", "U", "N", n, &i__1, &AFB(1,1), ldafb, &WORK(
#endif

	    1));
    if (rpvgrw == 0.) {
	rpvgrw = 1.;
    } else {

#ifdef PETSC_PREFIX_SUFFIX
	rpvgrw = dlangb_("M", n, kl, ku, &AB(1,1), ldab, &WORK(1)) /
#endif
#ifdef Q_C_PREFIX_SUFFIX
	rpvgrw = qlangb("M", n, kl, ku, &AB(1,1), ldab, &WORK(1)) /
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	rpvgrw = qlangb_("M", n, kl, ku, &AB(1,1), ldab, &WORK(1)) /
#endif

		 rpvgrw;
    }

/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dgbcon_(norm, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), &anorm, rcond,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgbcon(norm, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), &anorm, rcond,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgbcon_(norm, n, kl, ku, &AFB(1,1), ldafb, &IPIV(1), &anorm, rcond,
#endif

	     &WORK(1), &IWORK(1), info);

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
    dgbtrs_(trans, n, kl, ku, nrhs, &AFB(1,1), ldafb, &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgbtrs(trans, n, kl, ku, nrhs, &AFB(1,1), ldafb, &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgbtrs_(trans, n, kl, ku, nrhs, &AFB(1,1), ldafb, &IPIV(1), &X(1,1), ldx, info);
#endif


/*     Use iterative refinement to improve the computed solution and   
       compute error bounds and backward error estimates for it. */


#ifdef PETSC_PREFIX_SUFFIX
    dgbrfs_(trans, n, kl, ku, nrhs, &AB(1,1), ldab, &AFB(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgbrfs(trans, n, kl, ku, nrhs, &AB(1,1), ldab, &AFB(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgbrfs_(trans, n, kl, ku, nrhs, &AB(1,1), ldab, &AFB(1,1), 
#endif

	    ldafb, &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &
	    BERR(1), &WORK(1), &IWORK(1), info);

/*     Transform the solution matrix X to a solution of the original   
       system. */

    if (notran) {
	if (colequ) {
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__3 = *n;
		for (i = 1; i <= *n; ++i) {
		    X(i,j) = C(i) * X(i,j);
/* L100: */
		}
/* L110: */
	    }
	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		FERR(j) /= colcnd;
/* L120: */
	    }
	}
    } else if (rowequ) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    i__3 = *n;
	    for (i = 1; i <= *n; ++i) {
		X(i,j) = R(i) * X(i,j);
/* L130: */
	    }
/* L140: */
	}
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) /= rowcnd;
/* L150: */
	}
    }

    WORK(1) = rpvgrw;
    return;

/*     End of DGBSVX */

} /* dgbsvx_ */

