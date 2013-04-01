#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgels_(char *trans, int *m, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgels(char *trans, int *m, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgels_(char *trans, int *m, int *n, int *
#endif

	nrhs, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGELS solves overdetermined or underdetermined real linear systems   
    involving an M-by-N matrix A, or its transpose, using a QR or LQ   
    factorization of A.  It is assumed that A has full rank.   

    The following options are provided:   

    1. If TRANS = 'N' and m >= n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A*X ||.   

    2. If TRANS = 'N' and m < n:  find the minimum norm solution of   
       an underdetermined system A * X = B.   

    3. If TRANS = 'T' and m >= n:  find the minimum norm solution of   
       an undetermined system A**T * X = B.   

    4. If TRANS = 'T' and m < n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A**T * X ||.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution   
    matrix X.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
            = 'N': the linear system involves A;   
            = 'T': the linear system involves A**T.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of   
            columns of the matrices B and X. NRHS >=0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit,   
              if M >= N, A is overwritten by details of its QR   
                         factorization as returned by DGEQRF;   
              if M <  N, A is overwritten by details of its LQ   
                         factorization as returned by DGELQF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the matrix B of right hand side vectors, stored   
            columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS   
            if TRANS = 'T'.   
            On exit, B is overwritten by the solution vectors, stored   
            columnwise:   
            if TRANS = 'N' and m >= n, rows 1 to n of B contain the least 
  
            squares solution vectors; the residual sum of squares for the 
  
            solution in each column is given by the sum of squares of   
            elements N+1 to M in that column;   
            if TRANS = 'N' and m < n, rows 1 to N of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'T' and m >= n, rows 1 to M of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'T' and m < n, rows 1 to M of B contain the   
            least squares solution vectors; the residual sum of squares   
            for the solution in each column is given by the sum of   
            squares of elements M+1 to N in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,M,N).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            LWORK >= MIN(M,N) + MAX(1,M,N,NRHS).   
            For optimal performance,   
            LWORK >= MIN(M,N) + MAX(1,M,N,NRHS) * NB   
            where NB is the optimum block size.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static LONG DOUBLE c_b33 = 0.;
    static int c__0 = 0;
    static LONG DOUBLE c_b61 = 1.;
    
    /* System generated locals */
    int   i__1, i__2, i__3;
    /* Local variables */
    static LONG DOUBLE anrm, bnrm;
    static int brow;
    static long int tpsd;
    static int i, j, iascl, ibscl;
    extern long int lsame_(char *, char *);

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
    static int wsize;
    static LONG DOUBLE rwork[1];

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static int nb;

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
    extern /* Subroutine */ void dgelqf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelqf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelqf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, int *),

#ifdef PETSC_PREFIX_SUFFIX
	     dgeqrf_(int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     qgeqrf(int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     qgeqrf_(int *, int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), dlaset_(char *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlaset(char *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlaset_(char *,
#endif

	     int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static int scllen;
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormlq_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormlq(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormlq_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_(char *, char *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, int *);
    static LONG DOUBLE smlnum;



#define RWORK(I) rwork[(I)]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    mn = MIN(*m,*n);
    if (! (lsame_(trans, "N") || lsame_(trans, "T"))) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < MAX(1,*m)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*m);
	if (*ldb < MAX(i__1,*n)) {
	    *info = -8;
	} else /* if(complicated condition) */ {
/* Computing MAX   
   Computing MAX */
	    i__3 = MAX(*m,*n);
	    i__1 = 1, i__2 = mn + MAX(i__3,*nrhs);
	    if (*lwork < MAX(i__1,i__2)) {
		*info = -10;
	    }
	}
    }

/*     Figure out optimal block size */

    if (*info == 0 || *info == -10) {

	tpsd = 1;
	if (lsame_(trans, "N")) {
	    tpsd = 0;
	}

	if (*m >= *n) {
	    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "DORMQR", "LN", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = MAX(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "DORMQR", "LT", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = MAX(i__1,i__2);
	    }
	} else {
	    nb = ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "DORMLQ", "LT", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = MAX(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_(&c__1, "DORMLQ", "LN", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = MAX(i__1,i__2);
	    }
	}

/* Computing MAX */
	i__1 = MAX(*m,*n);
	wsize = mn + MAX(i__1,*nrhs) * nb;
	WORK(1) = (LONG DOUBLE) wsize;

    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGELS ", &i__1);
	return;
    }

/*     Quick return if possible   

   Computing MIN */
    i__1 = MIN(*m,*n);
    if (MIN(i__1,*nrhs) == 0) {
	i__1 = MAX(*m,*n);

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif

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


/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", m, n, &A(1,1), lda, rwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", m, n, &A(1,1), lda, rwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", m, n, &A(1,1), lda, rwork);
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
	dlaset_("F", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("F", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("F", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
#endif

	goto L50;
    }

    brow = *m;
    if (tpsd) {
	brow = *n;
    }

#ifdef PETSC_PREFIX_SUFFIX
    bnrm = dlange_("M", &brow, nrhs, &B(1,1), ldb, rwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnrm = qlange("M", &brow, nrhs, &B(1,1), ldb, rwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnrm = qlange_("M", &brow, nrhs, &B(1,1), ldb, rwork);
#endif

    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &B(1,1), 
#endif

		ldb, info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */


#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &B(1,1), 
#endif

		ldb, info);
	ibscl = 2;
    }

    if (*m >= *n) {

/*        compute QR factorization of A */

	i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	dgeqrf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgeqrf(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgeqrf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif

		;

/*        workspace at least N, optimally N*NB */

	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B ||   

             B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_("Left", "Transpose", m, nrhs, n, &A(1,1), lda, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr("Left", "Transpose", m, nrhs, n, &A(1,1), lda, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_("Left", "Transpose", m, nrhs, n, &A(1,1), lda, &WORK(
#endif

		    1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info)
		    ;

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
	    dtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &
#endif

		    c_b61, &A(1,1), lda, &B(1,1), ldb);

	    scllen = *n;

	} else {

/*           Overdetermined system of equations A' * X = B   

             B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
	    dtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b61, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrsm("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b61, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b61, 
#endif

		    &A(1,1), lda, &B(1,1), ldb);

/*           B(N+1:M,1:NRHS) = ZERO */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *m;
		for (i = *n + 1; i <= *m; ++i) {
		    B(i,j) = 0.;
/* L10: */
		}
/* L20: */
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_("Left", "No transpose", m, nrhs, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr("Left", "No transpose", m, nrhs, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_("Left", "No transpose", m, nrhs, n, &A(1,1), lda, &
#endif

		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *m;

	}

    } else {

/*        Compute LQ factorization of A */

	i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	dgelqf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgelqf(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgelqf_(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
#endif

		;

/*        workspace at least M, optimally M*NB. */

	if (! tpsd) {

/*           underdetermined system of equations A * X = B   

             B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
	    dtrsm_("Left", "Lower", "No transpose", "Non-unit", m, nrhs, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrsm("Left", "Lower", "No transpose", "Non-unit", m, nrhs, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrsm_("Left", "Lower", "No transpose", "Non-unit", m, nrhs, &
#endif

		    c_b61, &A(1,1), lda, &B(1,1), ldb);

/*           B(M+1:N,1:NRHS) = 0 */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *n;
		for (i = *m + 1; i <= *n; ++i) {
		    B(i,j) = 0.;
/* L30: */
		}
/* L40: */
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	    dormlq_("Left", "Transpose", n, nrhs, m, &A(1,1), lda, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormlq("Left", "Transpose", n, nrhs, m, &A(1,1), lda, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormlq_("Left", "Transpose", n, nrhs, m, &A(1,1), lda, &WORK(
#endif

		    1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info)
		    ;

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *n;

	} else {

/*           overdetermined system min || A' * X - B ||   

             B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;

#ifdef PETSC_PREFIX_SUFFIX
	    dormlq_("Left", "No transpose", n, nrhs, m, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormlq("Left", "No transpose", n, nrhs, m, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormlq_("Left", "No transpose", n, nrhs, m, &A(1,1), lda, &
#endif

		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS) */


#ifdef PETSC_PREFIX_SUFFIX
	    dtrsm_("Left", "Lower", "Transpose", "Non-unit", m, nrhs, &c_b61, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrsm("Left", "Lower", "Transpose", "Non-unit", m, nrhs, &c_b61, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrsm_("Left", "Lower", "Transpose", "Non-unit", m, nrhs, &c_b61, 
#endif

		    &A(1,1), lda, &B(1,1), ldb);

	    scllen = *m;

	}

    }

/*     Undo scaling */

    if (iascl == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &B(1,1)
#endif

		, ldb, info);
    } else if (iascl == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &B(1,1)
#endif

		, ldb, info);
    }
    if (ibscl == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &B(1,1)
#endif

		, ldb, info);
    } else if (ibscl == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &B(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &B(1,1)
#endif

		, ldb, info);
    }

L50:
    WORK(1) = (LONG DOUBLE) wsize;

    return;

/*     End of DGELS */

} /* dgels_ */

