#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgelsx_(int *m, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgelsx(int *m, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgelsx_(int *m, int *n, int *nrhs, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, int *
	jpvt, LONG DOUBLE *rcond, int *rank, LONG DOUBLE *work, int *
	info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGELSX computes the minimum-norm solution to a real linear least   
    squares problem:   
        minimize || A * X - B ||   
    using a complete orthogonal factorization of A.  A is an M-by-N   
    matrix which may be rank-deficient.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution   
    matrix X.   

    The routine first computes a QR factorization with column pivoting:   
        A * P = Q * [ R11 R12 ]   
                    [  0  R22 ]   
    with R11 defined as the largest leading submatrix whose estimated   
    condition number is less than 1/RCOND.  The order of R11, RANK,   
    is the effective rank of A.   

    Then, R22 is considered to be negligible, and R12 is annihilated   
    by orthogonal transformations from the right, arriving at the   
    complete orthogonal factorization:   
       A * P = Q * [ T11 0 ] * Z   
                   [  0  0 ]   
    The minimum-norm solution is then   
       X = P * Z' [ inv(T11)*Q1'*B ]   
                  [        0       ]   
    where Q1 consists of the first RANK columns of Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of   
            columns of matrices B and X. NRHS >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A has been overwritten by details of its   
            complete orthogonal factorization.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the M-by-NRHS right hand side matrix B.   
            On exit, the N-by-NRHS solution matrix X.   
            If m >= n and RANK = n, the residual sum-of-squares for   
            the solution in the i-th column is given by the sum of   
            squares of elements N+1:M in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,M,N).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(i) .ne. 0, the i-th column of A is an   
            initial column, otherwise it is a free column.  Before   
            the QR factorization of A, all initial columns are   
            permuted to the leading positions; only the remaining   
            free columns are moved as a result of column pivoting   
            during the factorization.   
            On exit, if JPVT(i) = k, then the i-th column of A*P   
            was the k-th column of A.   

    RCOND   (input) LONG DOUBLE PRECISION   
            RCOND is used to determine the effective rank of A, which   
            is defined as the order of the largest leading triangular   
            submatrix R11 in the QR factorization with pivoting of A,   
            whose estimated condition number < 1/RCOND.   

    RANK    (output) INTEGER   
            The effective rank of A, i.e., the order of the submatrix   
            R11.  This is the same as the order of the submatrix T11   
            in the complete orthogonal factorization of A.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
                        (MAX( MIN(M,N)+3*N, 2*MIN(M,N)+NRHS )),   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__0 = 0;
    static LONG DOUBLE c_b13 = 0.;
    static int c__2 = 2;
    static int c__1 = 1;
    static LONG DOUBLE c_b36 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1;
    /* Local variables */
    static LONG DOUBLE anrm, bnrm, smin, smax;
    static int i, j, k, iascl, ibscl, ismin, ismax;
    static LONG DOUBLE c1, c2;

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

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dlaic1_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaic1(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaic1_(
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE s1, s2, t1, t2;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorm2r_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorm2r(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorm2r_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), dlabad_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qlabad(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qlabad_(
#endif

	    LONG DOUBLE *, LONG DOUBLE *);

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
    static int mn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlascl_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dgeqpf_(int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qgeqpf(int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qgeqpf_(int *, int *, 
#endif

	    LONG DOUBLE *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlaset_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlaset(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlaset_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, 
	    int *);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatzm_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatzm(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatzm_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *, LONG DOUBLE *);
    static LONG DOUBLE sminpr, smaxpr, smlnum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtzrqf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtzrqf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtzrqf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);



#define JPVT(I) jpvt[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    mn = MIN(*m,*n);
    ismin = mn + 1;
    ismax = (mn << 1) + 1;

/*     Test the input arguments. */

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*m);
	if (*ldb < MAX(i__1,*n)) {
	    *info = -7;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGELSX", &i__1);
	return;
    }

/*     Quick return if possible   

   Computing MIN */
    i__1 = MIN(*m,*n);
    if (MIN(i__1,*nrhs) == 0) {
	*rank = 0;
	return;
    }

/*     Get machine parameters */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S") / dlamch_("P");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S") / dlamch_("P");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S") / dlamch_("P");
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


/*     Scale A, B if max elements outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", m, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", m, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", m, n, &A(1,1), lda, &WORK(1));
#endif

    iascl = 0;
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
#endif

		info);
	iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
#endif

		info);
	iascl = 2;
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

	i__1 = MAX(*m,*n);

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif

	*rank = 0;
	goto L100;
    }


#ifdef PETSC_PREFIX_SUFFIX
    bnrm = dlange_("M", m, nrhs, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnrm = qlange("M", m, nrhs, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnrm = qlange_("M", m, nrhs, &B(1,1), ldb, &WORK(1));
#endif

    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &B(1,1), ldb,
#endif

		 info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &B(1,1), ldb,
#endif

		 info);
	ibscl = 2;
    }

/*     Compute QR factorization with column pivoting of A:   
          A * P = Q * R */


#ifdef PETSC_PREFIX_SUFFIX
    dgeqpf_(m, n, &A(1,1), lda, &JPVT(1), &WORK(1), &WORK(mn + 1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqpf(m, n, &A(1,1), lda, &JPVT(1), &WORK(1), &WORK(mn + 1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqpf_(m, n, &A(1,1), lda, &JPVT(1), &WORK(1), &WORK(mn + 1), info);
#endif


/*     workspace 3*N. Details of Householder rotations stored   
       in WORK(1:MN).   

       Determine RANK using incremental condition estimation */

    WORK(ismin) = 1.;
    WORK(ismax) = 1.;
    smax = (d__1 = A(1,1), ABS(d__1));
    smin = smax;
    if ((d__1 = A(1,1), ABS(d__1)) == 0.) {
	*rank = 0;
	i__1 = MAX(*m,*n);

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &B(1,1), ldb);
#endif

	goto L100;
    } else {
	*rank = 1;
    }

L10:
    if (*rank < mn) {
	i = *rank + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlaic1_(&c__2, rank, &WORK(ismin), &smin, &A(1,i), &A(i,i), &sminpr, &s1, &c1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaic1(&c__2, rank, &WORK(ismin), &smin, &A(1,i), &A(i,i), &sminpr, &s1, &c1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaic1_(&c__2, rank, &WORK(ismin), &smin, &A(1,i), &A(i,i), &sminpr, &s1, &c1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlaic1_(&c__1, rank, &WORK(ismax), &smax, &A(1,i), &A(i,i), &smaxpr, &s2, &c2);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaic1(&c__1, rank, &WORK(ismax), &smax, &A(1,i), &A(i,i), &smaxpr, &s2, &c2);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaic1_(&c__1, rank, &WORK(ismax), &smax, &A(1,i), &A(i,i), &smaxpr, &s2, &c2);
#endif


	if (smaxpr * *rcond <= sminpr) {
	    i__1 = *rank;
	    for (i = 1; i <= *rank; ++i) {
		WORK(ismin + i - 1) = s1 * WORK(ismin + i - 1);
		WORK(ismax + i - 1) = s2 * WORK(ismax + i - 1);
/* L20: */
	    }
	    WORK(ismin + *rank) = c1;
	    WORK(ismax + *rank) = c2;
	    smin = sminpr;
	    smax = smaxpr;
	    ++(*rank);
	    goto L10;
	}
    }

/*     Logically partition R = [ R11 R12 ]   
                               [  0  R22 ]   
       where R11 = R(1:RANK,1:RANK)   

       [R11,R12] = [ T11, 0 ] * Y */

    if (*rank < *n) {

#ifdef PETSC_PREFIX_SUFFIX
	dtzrqf_(rank, n, &A(1,1), lda, &WORK(mn + 1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtzrqf(rank, n, &A(1,1), lda, &WORK(mn + 1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtzrqf_(rank, n, &A(1,1), lda, &WORK(mn + 1), info);
#endif

    }

/*     Details of Householder rotations stored in WORK(MN+1:2*MN)   

       B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
    dorm2r_("Left", "Transpose", m, nrhs, &mn, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qorm2r("Left", "Transpose", m, nrhs, &mn, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qorm2r_("Left", "Transpose", m, nrhs, &mn, &A(1,1), lda, &WORK(1), &
#endif

	    B(1,1), ldb, &WORK((mn << 1) + 1), info);

/*     workspace NRHS   

       B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
    dtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrsm("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &
#endif

	    A(1,1), lda, &B(1,1), ldb);

    i__1 = *n;
    for (i = *rank + 1; i <= *n; ++i) {
	i__2 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    B(i,j) = 0.;
/* L30: */
	}
/* L40: */
    }

/*     B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS) */

    if (*rank < *n) {
	i__1 = *rank;
	for (i = 1; i <= *rank; ++i) {
	    i__2 = *n - *rank + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlatzm_("Left", &i__2, nrhs, &A(i,*rank+1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatzm("Left", &i__2, nrhs, &A(i,*rank+1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatzm_("Left", &i__2, nrhs, &A(i,*rank+1), lda, &
#endif

		    WORK(mn + i), &B(i,1), &B(*rank+1,1), ldb,
		     &WORK((mn << 1) + 1));
/* L50: */
	}
    }

/*     workspace NRHS   

       B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK((mn << 1) + i) = 1.;
/* L60: */
	}
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK((mn << 1) + i) == 1.) {
		if (JPVT(i) != i) {
		    k = i;
		    t1 = B(k,j);
		    t2 = B(JPVT(k),j);
L70:
		    B(JPVT(k),j) = t1;
		    WORK((mn << 1) + k) = 0.;
		    t1 = t2;
		    k = JPVT(k);
		    t2 = B(JPVT(k),j);
		    if (JPVT(k) != i) {
			goto L70;
		    }
		    B(i,j) = t1;
		    WORK((mn << 1) + k) = 0.;
		}
	    }
/* L80: */
	}
/* L90: */
    }

/*     Undo scaling */

    if (iascl == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &B(1,1), ldb,
#endif

		 info);

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &A(1,1), 
#endif

		lda, info);
    } else if (iascl == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &B(1,1), ldb,
#endif

		 info);

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &A(1,1), 
#endif

		lda, info);
    }
    if (ibscl == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif

		 info);
    } else if (ibscl == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &B(1,1), ldb,
#endif

		 info);
    }

L100:

    return;

/*     End of DGELSX */

} /* dgelsx_ */

