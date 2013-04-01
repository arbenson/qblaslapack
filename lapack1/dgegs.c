#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgegs_(char *jobvsl, char *jobvsr, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgegs(char *jobvsl, char *jobvsr, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgegs_(char *jobvsl, char *jobvsr, int *n, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *
	alphar, LONG DOUBLE *alphai, LONG DOUBLE *beta, LONG DOUBLE *vsl, 
	int *ldvsl, LONG DOUBLE *vsr, int *ldvsr, LONG DOUBLE *work, 
	int *lwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEGS computes for a pair of N-by-N real nonsymmetric matrices A, B: 
  
    the generalized eigenvalues (alphar +/- alphai*i, beta), the real   
    Schur form (A, B), and optionally left and/or right Schur vectors   
    (VSL and VSR).   

    (If only the generalized eigenvalues are needed, use the driver DGEGV 
  
    instead.)   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "Matrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    The (generalized) Schur form of a pair of matrices is the result of   
    multiplying both matrices on the left by one orthogonal matrix and   
    both on the right by another orthogonal matrix, these two orthogonal 
  
    matrices being chosen so as to bring the pair of matrices into   
    (real) Schur form.   

    A pair of matrices A, B is in generalized real Schur form if B is   
    upper triangular with non-negative diagonal and A is block upper   
    triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond   
    to real generalized eigenvalues, while 2-by-2 blocks of A will be   
    "standardized" by making the corresponding elements of B have the   
    form:   
            [  a  0  ]   
            [  0  b  ]   

    and the pair of corresponding 2-by-2 blocks in A and B will   
    have a complex conjugate pair of generalized eigenvalues.   

    The left and right Schur vectors are the columns of VSL and VSR,   
    respectively, where VSL and VSR are the orthogonal matrices   
    which reduce A and B to Schur form:   

    Schur form of (A,B) = ( (VSL)**T A (VSR), (VSL)**T B (VSR) )   

    Arguments   
    =========   

    JOBVSL  (input) CHARACTER*1   
            = 'N':  do not compute the left Schur vectors;   
            = 'V':  compute the left Schur vectors.   

    JOBVSR  (input) CHARACTER*1   
            = 'N':  do not compute the right Schur vectors;   
            = 'V':  compute the right Schur vectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VSL, and VSR.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose generalized 
  
            eigenvalues and (optionally) Schur vectors are to be   
            computed.   
            On exit, the generalized Schur form of A.   
            Note: to avoid overflow, the Frobenius norm of the matrix   
            A should be less than the overflow threshold.   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) Schur vectors are   
            to be computed.   
            On exit, the generalized Schur form of B.   
            Note: to avoid overflow, the Frobenius norm of the matrix   
            B should be less than the overflow threshold.   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= MAX(1,N).   

    ALPHAR  (output) LONG DOUBLE PRECISION array, dimension (N)   
    ALPHAI  (output) LONG DOUBLE PRECISION array, dimension (N)   
    BETA    (output) LONG DOUBLE PRECISION array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,   
            j=1,...,N  and  BETA(j),j=1,...,N  are the diagonals of the   
            complex Schur form (A,B) that would result if the 2-by-2   
            diagonal blocks of the real Schur form of (A,B) were further 
  
            reduced to triangular form using 2-by-2 complex unitary   
            transformations.  If ALPHAI(j) is zero, then the j-th   
            eigenvalue is real; if positive, then the j-th and (j+1)-st   
            eigenvalues are a complex conjugate pair, with ALPHAI(j+1)   
            negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero. 
  
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VSL     (output) LONG DOUBLE PRECISION array, dimension (LDVSL,N)   
            If JOBVSL = 'V', VSL will contain the left Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSL = 'N'.   

    LDVSL   (input) INTEGER   
            The leading dimension of the matrix VSL. LDVSL >=1, and   
            if JOBVSL = 'V', LDVSL >= N.   

    VSR     (output) LONG DOUBLE PRECISION array, dimension (LDVSR,N)   
            If JOBVSR = 'V', VSR will contain the right Schur vectors.   
            (See "Purpose", above.)   
            Not referenced if JOBVSR = 'N'.   

    LDVSR   (input) INTEGER   
            The leading dimension of the matrix VSR. LDVSR >= 1, and   
            if JOBVSR = 'V', LDVSR >= N.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= MAX(1,4*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR   
            The optimal LWORK is  2*N + N*(NB+1).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  (A,B) are not in Schur   
                  form, but ALPHAR(j), ALPHAI(j), and BETA(j) should   
                  be correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from DGGBAL   
                  =N+2: error return from DGEQRF   
                  =N+3: error return from DORMQR   
                  =N+4: error return from DORGQR   
                  =N+5: error return from DGGHRD   
                  =N+6: error return from DHGEQZ (other than failed   
                                                  iteration)   
                  =N+7: error return from DGGBAK (computing VSL)   
                  =N+8: error return from DGGBAK (computing VSR)   
                  =N+9: error return from DLASCL (various places)   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c_n1 = -1;
    static LONG DOUBLE c_b23 = 0.;
    static LONG DOUBLE c_b24 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static LONG DOUBLE anrm, bnrm;
    static int itau;
    extern long int lsame_(char *, char *);
    static int ileft, iinfo, icols;
    static long int ilvsl;
    static int iwork;
    static long int ilvsr;
    static int irows;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dggbak_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qggbak(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qggbak_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dggbal_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qggbal(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qggbal_(char *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *);

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
    extern /* Subroutine */ void dgghrd_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgghrd(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgghrd_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), dlascl_(char *, int *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qlascl(char *, int *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *), qlascl_(char *, int *, int *, LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    int *);
    static long int ilascl, ilbscl;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgeqrf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqrf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqrf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static LONG DOUBLE safmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dhgeqz_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qhgeqz(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qhgeqz_(char *, char *, char *, int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    int *);
    static int ijobvl, iright, ijobvr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorgqr_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgqr(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgqr_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    int *);
    static LONG DOUBLE anrmto;
    static int lwkmin;
    static LONG DOUBLE bnrmto;

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
    static LONG DOUBLE smlnum;
    static int lwkopt, ihi, ilo;
    static LONG DOUBLE eps;



#define ALPHAR(I) alphar[(I)-1]
#define ALPHAI(I) alphai[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define VSL(I,J) vsl[(I)-1 + ((J)-1)* ( *ldvsl)]
#define VSR(I,J) vsr[(I)-1 + ((J)-1)* ( *ldvsr)]

    if (lsame_(jobvsl, "N")) {
	ijobvl = 1;
	ilvsl = 0;
    } else if (lsame_(jobvsl, "V")) {
	ijobvl = 2;
	ilvsl = 1;
    } else {
	ijobvl = -1;
	ilvsl = 0;
    }

    if (lsame_(jobvsr, "N")) {
	ijobvr = 1;
	ilvsr = 0;
    } else if (lsame_(jobvsr, "V")) {
	ijobvr = 2;
	ilvsr = 1;
    } else {
	ijobvr = -1;
	ilvsr = 0;
    }

/*     Test the input arguments   

   Computing MAX */
    i__1 = *n << 2;
    lwkmin = MAX(i__1,1);
    lwkopt = lwkmin;
    *info = 0;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    } else if (*ldvsl < 1 || (ilvsl && *ldvsl < *n)) {
	*info = -12;
    } else if (*ldvsr < 1 || (ilvsr && *ldvsr < *n)) {
	*info = -14;
    } else if (*lwork < lwkmin) {
	*info = -16;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEGS ", &i__1);
	return;
    }

/*     Quick return if possible */

    WORK(1) = (LONG DOUBLE) lwkopt;
    if (*n == 0) {
	return;
    }

/*     Get machine constants */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("E") * dlamch_("B");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("E") * dlamch_("B");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("E") * dlamch_("B");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    safmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    safmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    safmin = qlamch_("S");
#endif

    smlnum = *n * safmin / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", n, n, &A(1,1), lda, &WORK(1));
#endif

    ilascl = 0;
    if (anrm > 0. && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = 1;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = 1;
    }

    if (ilascl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &A(1,1), lda, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */


#ifdef PETSC_PREFIX_SUFFIX
    bnrm = dlange_("M", n, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnrm = qlange("M", n, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnrm = qlange_("M", n, n, &B(1,1), ldb, &WORK(1));
#endif

    ilbscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = 1;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = 1;
    }

    if (ilbscl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &B(1,1), ldb, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Workspace layout:  (2*N words -- "work..." not actually used)   
          left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dggbal_("P", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qggbal("P", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qggbal_("P", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif

	    ileft), &WORK(iright), &WORK(iwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L10;
    }

/*     Reduce B to triangular form, and initialize VSL and/or VSR   
       Workspace layout:  ("work..." must have at least N words)   
          left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
    dgeqrf_(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqrf(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqrf_(&irows, &icols, &B(ilo,ilo), ldb, &WORK(itau), &WORK(
#endif

	    iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) WORK(iwork) + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L10;
    }

    i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
    dormqr_("L", "T", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormqr("L", "T", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormqr_("L", "T", &irows, &icols, &irows, &B(ilo,ilo), ldb, &
#endif

	    WORK(itau), &A(ilo,ilo), lda, &WORK(iwork), &i__1, &
	    iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) WORK(iwork) + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L10;
    }

    if (ilvsl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b23, &c_b24, &VSL(1,1), ldvsl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b23, &c_b24, &VSL(1,1), ldvsl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b23, &c_b24, &VSL(1,1), ldvsl);
#endif

	i__1 = irows - 1;
	i__2 = irows - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VSL(ilo+1,ilo), ldvsl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VSL(ilo+1,ilo), ldvsl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VSL(ilo+1,ilo), ldvsl);
#endif

	i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
	dorgqr_(&irows, &irows, &irows, &VSL(ilo,ilo), ldvsl, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgqr(&irows, &irows, &irows, &VSL(ilo,ilo), ldvsl, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgqr_(&irows, &irows, &irows, &VSL(ilo,ilo), ldvsl, &
#endif

		WORK(itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (int) WORK(iwork) + iwork - 1;
	    lwkopt = MAX(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L10;
	}
    }

    if (ilvsr) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b23, &c_b24, &VSR(1,1), ldvsr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b23, &c_b24, &VSR(1,1), ldvsr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b23, &c_b24, &VSR(1,1), ldvsr);
#endif

    }

/*     Reduce to generalized Hessenberg form */


#ifdef PETSC_PREFIX_SUFFIX
    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgghrd(jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif

	    ldb, &VSL(1,1), ldvsl, &VSR(1,1), ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = *n + 5;
	goto L10;
    }

/*     Perform QZ algorithm, computing Schur vectors if desired   
       Workspace layout:  ("work..." must have at least 1 word)   
          left_permutation, right_permutation, work... */

    iwork = itau;
    i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VSL(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qhgeqz("S", jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VSL(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VSL(1,1)
#endif

	    , ldvsl, &VSR(1,1), ldvsr, &WORK(iwork), &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (int) WORK(iwork) + iwork - 1;
	lwkopt = MAX(i__1,i__2);
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= *n) {
	    *info = iinfo;
	} else if (iinfo > *n && iinfo <= *n << 1) {
	    *info = iinfo - *n;
	} else {
	    *info = *n + 6;
	}
	goto L10;
    }

/*     Apply permutation to VSL and VSR */

    if (ilvsl) {

#ifdef PETSC_PREFIX_SUFFIX
	dggbak_("P", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSL(1,1), ldvsl, &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qggbak("P", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSL(1,1), ldvsl, &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qggbak_("P", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSL(1,1), ldvsl, &iinfo);
#endif

	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {

#ifdef PETSC_PREFIX_SUFFIX
	dggbak_("P", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSR(1,1), ldvsr, &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qggbak("P", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSR(1,1), ldvsr, &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qggbak_("P", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &VSR(1,1), ldvsr, &iinfo);
#endif

	if (iinfo != 0) {
	    *info = *n + 8;
	    goto L10;
	}
    }

/*     Undo scaling */

    if (ilascl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &A(1,1), lda, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAR(1), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAR(1), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAR(1), n, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAI(1), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAI(1), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &ALPHAI(1), n, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}
    }

    if (ilbscl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &B(1,1), ldb, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &BETA(1), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &BETA(1), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &BETA(1), n, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 9;
	    return;
	}
    }

L10:
    WORK(1) = (LONG DOUBLE) lwkopt;

    return;

/*     End of DGEGS */

} /* dgegs_ */

