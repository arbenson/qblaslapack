#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgelss_(int *m, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgelss(int *m, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgelss_(int *m, int *n, int *nrhs, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *
	s, LONG DOUBLE *rcond, int *rank, LONG DOUBLE *work, int *lwork,
	 int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGELSS computes the minimum norm solution to a real linear least   
    squares problem:   

    Minimize 2-norm(| b - A*x |).   

    using the singular value decomposition (SVD) of A. A is an M-by-N   
    matrix which may be rank-deficient.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix 
  
    X.   

    The effective rank of A is determined by treating as zero those   
    singular values which are less than RCOND times the largest singular 
  
    value.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X. NRHS >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the first MIN(m,n) rows of A are overwritten with   
            its right singular vectors, stored rowwise.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the M-by-NRHS right hand side matrix B.   
            On exit, B is overwritten by the N-by-NRHS solution   
            matrix X.  If m >= n and RANK = n, the residual   
            sum-of-squares for the solution in the i-th column is given   
            by the sum of squares of elements n+1:m in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,MAX(M,N)). 
  

    S       (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The singular values of A in decreasing order.   
            The condition number of A in the 2-norm = S(1)/S(MIN(m,n)).   

    RCOND   (input) LONG DOUBLE PRECISION   
            RCOND is used to determine the effective rank of A.   
            Singular values S(i) <= RCOND*S(1) are treated as zero.   
            If RCOND < 0, machine precision is used instead.   

    RANK    (output) INTEGER   
            The effective rank of A, i.e., the number of singular values 
  
            which are greater than RCOND*S(1).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= 1, and also:   
            LWORK >= 3*MIN(M,N) + MAX( 2*MIN(M,N), MAX(M,N), NRHS )   
            For good performance, LWORK should generally be larger.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  the algorithm for computing the SVD failed to converge; 
  
                  if INFO = i, i off-diagonal elements of an intermediate 
  
                  bidiagonal form did not converge to zero.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__6 = 6;
    static int c_n1 = -1;
    static int c__1 = 1;
    static int c__0 = 0;
    static LONG DOUBLE c_b74 = 0.;
    static LONG DOUBLE c_b108 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1;
    /* Local variables */
    static LONG DOUBLE anrm, bnrm;
    static int itau;
    static LONG DOUBLE vdum[1];
    static int i;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemm_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int iascl, ibscl;

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
	    LONG DOUBLE *, LONG DOUBLE *, int *), drscl_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qrscl(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qrscl_(int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static int chunk;
    static LONG DOUBLE sfmin;
    static int minmn;

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
    static int maxmn, itaup, itauq, mnthr, iwork;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabad_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static int bl, ie, il;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgebrd_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebrd(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebrd_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, int *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int mm;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlange_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlange(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlange_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *);
    static int bdspac;

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
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), dlacpy_(char *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlacpy(char *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qlacpy_(char *,
#endif

	     int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    xerbla_(char *, int *), dbdsqr_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    xerbla_(char *, int *), qbdsqr(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    xerbla_(char *, int *), qbdsqr_(char *, int *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), dorgbr_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qorgbr(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qorgbr_(char *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static LONG DOUBLE bignum;
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormbr_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormbr(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormbr_(char *, char *, char *, int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), dormlq_(char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qormlq(char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qormlq_(char *, char *, int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *);
    static int ldwork;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormqr_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormqr(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormqr_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *);
    static int minwrk, maxwrk;
    static LONG DOUBLE smlnum, eps, thr;



#define VDUM(I) vdum[(I)]
#define S(I) s[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    minmn = MIN(*m,*n);
    maxmn = MAX(*m,*n);
    mnthr = ilaenv_(&c__6, "DGELSS", " ", m, n, nrhs, &c_n1, 6L, 1L);
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*ldb < MAX(1,maxmn)) {
	*info = -7;
    }

/*     Compute workspace   
        (Note: Comments in the code beginning "Workspace:" describe the   
         minimal amount of workspace needed at that point in the code,   
         as well as the preferred amount for good performance.   
         NB refers to the optimal block size for the immediately   
         following subroutine, as returned by ILAENV.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	maxwrk = 0;
	mm = *m;
	if (*m >= *n && *m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than co
lumns */

	    mm = *n;
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, "DGEQRF", " ", m, 
		    n, &c_n1, &c_n1, 6L, 1L);
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, "DORMQR", "LT", 
		    m, nrhs, n, &c_n1, 6L, 2L);
	    maxwrk = MAX(i__1,i__2);
	}
	if (*m >= *n) {

/*           Path 1 - overdetermined or exactly determined   

             Compute workspace neede for DBDSQR   

   Computing MAX */
	    i__1 = 1, i__2 = *n * 5 - 4;
	    bdspac = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, "DGEBRD"
		    , " ", &mm, n, &c_n1, &c_n1, 6L, 1L);
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, "DORMBR", 
		    "QLT", &mm, nrhs, n, &c_n1, 6L, 3L);
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, "DORGBR",
		     "P", n, n, n, &c_n1, 6L, 1L);
	    maxwrk = MAX(i__1,i__2);
	    maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * *nrhs;
	    maxwrk = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = MAX(i__1,i__2);
	    minwrk = MAX(i__1,bdspac);
	    maxwrk = MAX(minwrk,maxwrk);
	}
	if (*n > *m) {

/*           Compute workspace neede for DBDSQR   

   Computing MAX */
	    i__1 = 1, i__2 = *m * 5 - 4;
	    bdspac = MAX(i__1,i__2);
/* Computing MAX */
	    i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = MAX(i__1,i__2);
	    minwrk = MAX(i__1,bdspac);
	    if (*n >= mnthr) {

/*              Path 2a - underdetermined, with many more colu
mns   
                than rows */

		maxwrk = *m + *m * ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, 
			&c_n1, 6L, 1L);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m << 1) * 
			ilaenv_(&c__1, "DGEBRD", " ", m, m, &c_n1, &c_n1, 6L, 
			1L);
		maxwrk = MAX(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + *nrhs * ilaenv_(&
			c__1, "DORMBR", "QLT", m, nrhs, m, &c_n1, 6L, 3L);
		maxwrk = MAX(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + (*m - 1) * 
			ilaenv_(&c__1, "DORGBR", "P", m, m, m, &c_n1, 6L, 1L);
		maxwrk = MAX(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
		maxwrk = MAX(i__1,i__2);
		if (*nrhs > 1) {
/* Computing MAX */
		    i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
		    maxwrk = MAX(i__1,i__2);
		} else {
/* Computing MAX */
		    i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
		    maxwrk = MAX(i__1,i__2);
		}
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m + *nrhs * ilaenv_(&c__1, "DORMLQ", 
			"LT", n, nrhs, m, &c_n1, 6L, 2L);
		maxwrk = MAX(i__1,i__2);
	    } else {

/*              Path 2 - underdetermined */

		maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, "DGEBRD", " ", m,
			 n, &c_n1, &c_n1, 6L, 1L);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, "DORMBR"
			, "QLT", m, nrhs, m, &c_n1, 6L, 3L);
		maxwrk = MAX(i__1,i__2);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *m * 3 + *m * ilaenv_(&c__1, "DORGBR", 
			"P", m, n, m, &c_n1, 6L, 1L);
		maxwrk = MAX(i__1,i__2);
		maxwrk = MAX(maxwrk,bdspac);
/* Computing MAX */
		i__1 = maxwrk, i__2 = *n * *nrhs;
		maxwrk = MAX(i__1,i__2);
	    }
	}
	maxwrk = MAX(minwrk,maxwrk);
	WORK(1) = (LONG DOUBLE) maxwrk;
    }

    minwrk = MAX(minwrk,1);
    if (*lwork < minwrk) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGELSS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	*rank = 0;
	return;
    }

/*     Get machine parameters */


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
    sfmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sfmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sfmin = qlamch_("S");
#endif

    smlnum = sfmin / eps;
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


/*     Scale A if max element outside range [SMLNUM,BIGNUM] */


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
	dlaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("F", &i__1, nrhs, &c_b74, &c_b74, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(1,1), ldb);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("F", &minmn, &c__1, &c_b74, &c_b74, &S(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("F", &minmn, &c__1, &c_b74, &c_b74, &S(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("F", &minmn, &c__1, &c_b74, &c_b74, &S(1), &c__1);
#endif

	*rank = 0;
	goto L70;
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */


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

/*     Overdetermined case */

    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

	mm = *m;
	if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than co
lumns */

	    mm = *n;
	    itau = 1;
	    iwork = itau + *n;

/*           Compute A=Q*R   
             (Workspace: need 2*N, prefer N+N*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgeqrf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgeqrf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__1,
#endif

		     info);

/*           Multiply B by transpose(Q)   
             (Workspace: need N+NRHS, prefer N+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_("L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr("L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_("L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif


/*           Zero out below R */

	    if (*n > 1) {
		i__1 = *n - 1;
		i__2 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlaset_("L", &i__1, &i__2, &c_b74, &c_b74, &A(2,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaset("L", &i__1, &i__2, &c_b74, &c_b74, &A(2,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaset_("L", &i__1, &i__2, &c_b74, &c_b74, &A(2,1), 
#endif

			lda);
	    }
	}

	ie = 1;
	itauq = ie + *n;
	itaup = itauq + *n;
	iwork = itaup + *n;

/*        Bidiagonalize R in A   
          (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

	i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgebrd_(&mm, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgebrd(&mm, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgebrd_(&mm, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif

		WORK(itaup), &WORK(iwork), &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R
   
          (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

	i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dormbr_("Q", "L", "T", &mm, nrhs, n, &A(1,1), lda, &WORK(itauq), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qormbr("Q", "L", "T", &mm, nrhs, n, &A(1,1), lda, &WORK(itauq), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qormbr_("Q", "L", "T", &mm, nrhs, n, &A(1,1), lda, &WORK(itauq), 
#endif

		&B(1,1), ldb, &WORK(iwork), &i__1, info);

/*        Generate right bidiagonalizing vectors of R in A   
          (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

	i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(iwork), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgbr("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(iwork), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgbr_("P", n, n, n, &A(1,1), lda, &WORK(itaup), &WORK(iwork), &
#endif

		i__1, info);
	iwork = ie + *n;

/*        Perform bidiagonal QR iteration   
            multiply B by transpose of left singular vectors   
            compute right singular vectors in A   
          (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
	dbdsqr_("U", n, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qbdsqr("U", n, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qbdsqr_("U", n, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), lda, 
#endif

		vdum, &c__1, &B(1,1), ldb, &WORK(iwork), info);
	if (*info != 0) {
	    goto L70;
	}

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	d__1 = *rcond * S(1);
	thr = MAX(d__1,sfmin);
	if (*rcond < 0.) {
/* Computing MAX */
	    d__1 = eps * S(1);
	    thr = MAX(d__1,sfmin);
	}
	*rank = 0;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (S(i) > thr) {

#ifdef PETSC_PREFIX_SUFFIX
		drscl_(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrscl(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrscl_(nrhs, &S(i), &B(i,1), ldb);
#endif

		++(*rank);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaset("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), ldb);
#endif

	    }
/* L10: */
	}

/*        Multiply B by right singular vectors   
          (Workspace: need N, prefer N*NRHS) */

	if (*lwork >= *ldb * *nrhs && *nrhs > 1) {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemm_("T", "N", n, nrhs, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemm("T", "N", n, nrhs, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemm_("T", "N", n, nrhs, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("G", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("G", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("G", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif

	} else if (*nrhs > 1) {
	    chunk = *lwork / *n;
	    i__1 = *nrhs;
	    i__2 = chunk;
	    for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		i__3 = *nrhs - i + 1;
		bl = MIN(i__3,chunk);

#ifdef PETSC_PREFIX_SUFFIX
		dgemm_("T", "N", n, &bl, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), n);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemm("T", "N", n, &bl, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), n);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemm_("T", "N", n, &bl, n, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), n);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("G", n, &bl, &WORK(1), n, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("G", n, &bl, &WORK(1), n, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("G", n, &bl, &WORK(1), n, &B(1,1), ldb);
#endif

/* L20: */
	    }
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("T", n, n, &c_b108, &A(1,1), lda, &B(1,1), &c__1,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("T", n, n, &c_b108, &A(1,1), lda, &B(1,1), &c__1,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("T", n, n, &c_b108, &A(1,1), lda, &B(1,1), &c__1,
#endif

		     &c_b74, &WORK(1), &c__1);

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif

	}

    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__2 = *m, i__1 = (*m << 1) - 4, i__2 = MAX(i__2,i__1), i__2 = MAX(
		i__2,*nrhs), i__1 = *n - *m * 3;
	if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + MAX(i__2,i__1)) {

/*        Path 2a - underdetermined, with many more columns than r
ows   
          and sufficient workspace for an efficient algorithm */

	    ldwork = *m;
/* Computing MAX   
   Computing MAX */
	    i__3 = *m, i__4 = (*m << 1) - 4, i__3 = MAX(i__3,i__4), i__3 = 
		    MAX(i__3,*nrhs), i__4 = *n - *m * 3;
	    i__2 = (*m << 2) + *m * *lda + MAX(i__3,i__4), i__1 = *m * *lda + 
		    *m + *m * *nrhs;
	    if (*lwork >= MAX(i__2,i__1)) {
		ldwork = *lda;
	    }
	    itau = 1;
	    iwork = *m + 1;

/*        Compute A=L*Q   
          (Workspace: need 2*M, prefer M+M*NB) */

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__2,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgelqf(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__2,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgelqf_(m, n, &A(1,1), lda, &WORK(itau), &WORK(iwork), &i__2,
#endif

		     info);
	    il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */


#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("L", m, m, &A(1,1), lda, &WORK(il), &ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("L", m, m, &A(1,1), lda, &WORK(il), &ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("L", m, m, &A(1,1), lda, &WORK(il), &ldwork);
#endif

	    i__2 = *m - 1;
	    i__1 = *m - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlaset_("U", &i__2, &i__1, &c_b74, &c_b74, &WORK(il + ldwork), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaset("U", &i__2, &i__1, &c_b74, &c_b74, &WORK(il + ldwork), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaset_("U", &i__2, &i__1, &c_b74, &c_b74, &WORK(il + ldwork), &
#endif

		    ldwork);
	    ie = il + ldwork * *m;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL)   
          (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgebrd_(m, m, &WORK(il), &ldwork, &S(1), &WORK(ie), &WORK(itauq), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgebrd(m, m, &WORK(il), &ldwork, &S(1), &WORK(ie), &WORK(itauq), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgebrd_(m, m, &WORK(il), &ldwork, &S(1), &WORK(ie), &WORK(itauq), 
#endif

		    &WORK(itaup), &WORK(iwork), &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors 
of L   
          (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) 
*/

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormbr_("Q", "L", "T", m, nrhs, m, &WORK(il), &ldwork, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormbr("Q", "L", "T", m, nrhs, m, &WORK(il), &ldwork, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormbr_("Q", "L", "T", m, nrhs, m, &WORK(il), &ldwork, &WORK(
#endif

		    itauq), &B(1,1), ldb, &WORK(iwork), &i__2, info);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) 
  
          (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */

	    i__2 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorgbr_("P", m, m, m, &WORK(il), &ldwork, &WORK(itaup), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorgbr("P", m, m, m, &WORK(il), &ldwork, &WORK(itaup), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorgbr_("P", m, m, m, &WORK(il), &ldwork, &WORK(itaup), &WORK(
#endif

		    iwork), &i__2, info);
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration,   
             computing right singular vectors of L in WORK(IL) and
   
             multiplying B by transpose of left singular vectors 
  
          (Workspace: need M*M+M+BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
	    dbdsqr_("U", m, m, &c__0, nrhs, &S(1), &WORK(ie), &WORK(il), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qbdsqr("U", m, m, &c__0, nrhs, &S(1), &WORK(ie), &WORK(il), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qbdsqr_("U", m, m, &c__0, nrhs, &S(1), &WORK(ie), &WORK(il), &
#endif

		    ldwork, &A(1,1), lda, &B(1,1), ldb, &WORK(iwork)
		    , info);
	    if (*info != 0) {
		goto L70;
	    }

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	    d__1 = *rcond * S(1);
	    thr = MAX(d__1,sfmin);
	    if (*rcond < 0.) {
/* Computing MAX */
		d__1 = eps * S(1);
		thr = MAX(d__1,sfmin);
	    }
	    *rank = 0;
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		if (S(i) > thr) {

#ifdef PETSC_PREFIX_SUFFIX
		    drscl_(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrscl(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrscl_(nrhs, &S(i), &B(i,1), ldb);
#endif

		    ++(*rank);
		} else {

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif

			    ldb);
		}
/* L30: */
	    }
	    iwork = ie;

/*        Multiply B by right singular vectors of L in WORK(IL)   
          (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */

	    if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {

#ifdef PETSC_PREFIX_SUFFIX
		dgemm_("T", "N", m, nrhs, m, &c_b108, &WORK(il), &ldwork, &B(1,1), ldb, &c_b74, &WORK(iwork), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemm("T", "N", m, nrhs, m, &c_b108, &WORK(il), &ldwork, &B(1,1), ldb, &c_b74, &WORK(iwork), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemm_("T", "N", m, nrhs, m, &c_b108, &WORK(il), &ldwork, &B(1,1), ldb, &c_b74, &WORK(iwork), ldb);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("G", m, nrhs, &WORK(iwork), ldb, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("G", m, nrhs, &WORK(iwork), ldb, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("G", m, nrhs, &WORK(iwork), ldb, &B(1,1), ldb);
#endif

	    } else if (*nrhs > 1) {
		chunk = (*lwork - iwork + 1) / *m;
		i__2 = *nrhs;
		i__1 = chunk;
		for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		    i__3 = *nrhs - i + 1;
		    bl = MIN(i__3,chunk);

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("T", "N", m, &bl, m, &c_b108, &WORK(il), &ldwork, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("T", "N", m, &bl, m, &c_b108, &WORK(il), &ldwork, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("T", "N", m, &bl, m, &c_b108, &WORK(il), &ldwork, &
#endif

			    B(1,i), ldb, &c_b74, &WORK(iwork), n);

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("G", m, &bl, &WORK(iwork), n, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("G", m, &bl, &WORK(iwork), n, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("G", m, &bl, &WORK(iwork), n, &B(1,1), ldb);
#endif

/* L40: */
		}
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("T", m, m, &c_b108, &WORK(il), &ldwork, &B(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("T", m, m, &c_b108, &WORK(il), &ldwork, &B(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("T", m, m, &c_b108, &WORK(il), &ldwork, &B(1,1),
#endif

			 &c__1, &c_b74, &WORK(iwork), &c__1);

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(m, &WORK(iwork), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(m, &WORK(iwork), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(m, &WORK(iwork), &c__1, &B(1,1), &c__1);
#endif

	    }

/*        Zero out below first M rows of B */

	    i__1 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
	    dlaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(*m+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaset("F", &i__1, nrhs, &c_b74, &c_b74, &B(*m+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaset_("F", &i__1, nrhs, &c_b74, &c_b74, &B(*m+1,1), 
#endif

		    ldb);
	    iwork = itau + *m;

/*        Multiply transpose(Q) by B   
          (Workspace: need M+NRHS, prefer M+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormlq_("L", "T", n, nrhs, m, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormlq("L", "T", n, nrhs, m, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormlq_("L", "T", n, nrhs, m, &A(1,1), lda, &WORK(itau), &B(1,1), ldb, &WORK(iwork), &i__1, info);
#endif


	} else {

/*        Path 2 - remaining underdetermined cases */

	    ie = 1;
	    itauq = ie + *m;
	    itaup = itauq + *m;
	    iwork = itaup + *m;

/*        Bidiagonalize A   
          (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgebrd(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgebrd_(m, n, &A(1,1), lda, &S(1), &WORK(ie), &WORK(itauq), &
#endif

		    WORK(itaup), &WORK(iwork), &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors 
  
          (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormbr_("Q", "L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itauq)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormbr("Q", "L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itauq)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormbr_("Q", "L", "T", m, nrhs, n, &A(1,1), lda, &WORK(itauq)
#endif

		    , &B(1,1), ldb, &WORK(iwork), &i__1, info);

/*        Generate right bidiagonalizing vectors in A   
          (Workspace: need 4*M, prefer 3*M+M*NB) */

	    i__1 = *lwork - iwork + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorgbr("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorgbr_("P", m, n, m, &A(1,1), lda, &WORK(itaup), &WORK(
#endif

		    iwork), &i__1, info);
	    iwork = ie + *m;

/*        Perform bidiagonal QR iteration,   
             computing right singular vectors of A in A and   
             multiplying B by transpose of left singular vectors 
  
          (Workspace: need BDSPAC) */


#ifdef PETSC_PREFIX_SUFFIX
	    dbdsqr_("L", m, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qbdsqr("L", m, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qbdsqr_("L", m, n, &c__0, nrhs, &S(1), &WORK(ie), &A(1,1), 
#endif

		    lda, vdum, &c__1, &B(1,1), ldb, &WORK(iwork), info);
	    if (*info != 0) {
		goto L70;
	    }

/*        Multiply B by reciprocals of singular values   

   Computing MAX */
	    d__1 = *rcond * S(1);
	    thr = MAX(d__1,sfmin);
	    if (*rcond < 0.) {
/* Computing MAX */
		d__1 = eps * S(1);
		thr = MAX(d__1,sfmin);
	    }
	    *rank = 0;
	    i__1 = *m;
	    for (i = 1; i <= *m; ++i) {
		if (S(i) > thr) {

#ifdef PETSC_PREFIX_SUFFIX
		    drscl_(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrscl(nrhs, &S(i), &B(i,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrscl_(nrhs, &S(i), &B(i,1), ldb);
#endif

		    ++(*rank);
		} else {

#ifdef PETSC_PREFIX_SUFFIX
		    dlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaset("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaset_("F", &c__1, nrhs, &c_b74, &c_b74, &B(i,1), 
#endif

			    ldb);
		}
/* L50: */
	    }

/*        Multiply B by right singular vectors of A   
          (Workspace: need N, prefer N*NRHS) */

	    if (*lwork >= *ldb * *nrhs && *nrhs > 1) {

#ifdef PETSC_PREFIX_SUFFIX
		dgemm_("T", "N", n, nrhs, m, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemm("T", "N", n, nrhs, m, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemm_("T", "N", n, nrhs, m, &c_b108, &A(1,1), lda, &B(1,1), ldb, &c_b74, &WORK(1), ldb);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("F", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("F", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("F", n, nrhs, &WORK(1), ldb, &B(1,1), ldb);
#endif

	    } else if (*nrhs > 1) {
		chunk = *lwork / *n;
		i__1 = *nrhs;
		i__2 = chunk;
		for (i = 1; chunk < 0 ? i >= *nrhs : i <= *nrhs; i += chunk) {
/* Computing MIN */
		    i__3 = *nrhs - i + 1;
		    bl = MIN(i__3,chunk);

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("T", "N", n, &bl, m, &c_b108, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("T", "N", n, &bl, m, &c_b108, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("T", "N", n, &bl, m, &c_b108, &A(1,1), lda, &
#endif

			    B(1,i), ldb, &c_b74, &WORK(1), n);

#ifdef PETSC_PREFIX_SUFFIX
		    dlacpy_("F", n, &bl, &WORK(1), n, &B(1,i), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlacpy("F", n, &bl, &WORK(1), n, &B(1,i), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlacpy_("F", n, &bl, &WORK(1), n, &B(1,i), ldb);
#endif

/* L60: */
		}
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("T", m, n, &c_b108, &A(1,1), lda, &B(1,1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("T", m, n, &c_b108, &A(1,1), lda, &B(1,1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("T", m, n, &c_b108, &A(1,1), lda, &B(1,1), &
#endif

			c__1, &c_b74, &WORK(1), &c__1);

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(n, &WORK(1), &c__1, &B(1,1), &c__1);
#endif

	    }
	}
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
	dlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &S(1), &
#endif

		minmn, info);
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
	dlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &S(1), &
#endif

		minmn, info);
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

L70:
    WORK(1) = (LONG DOUBLE) maxwrk;
    return;

/*     End of DGELSS */

} /* dgelss_ */

