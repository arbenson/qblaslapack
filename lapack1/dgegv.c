#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgegv_(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgegv(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgegv_(char *jobvl, char *jobvr, int *n, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *alphar, 
	LONG DOUBLE *alphai, LONG DOUBLE *beta, LONG DOUBLE *vl, int *ldvl, 
	LONG DOUBLE *vr, int *ldvr, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEGV computes for a pair of n-by-n real nonsymmetric matrices A and 
  
    B, the generalized eigenvalues (alphar +/- alphai*i, beta), and   
    optionally, the left and/or right generalized eigenvectors (VL and   
    VR).   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "Matrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    A right generalized eigenvector corresponding to a generalized   
    eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such   
    that  (A - w B) r = 0 .  A left generalized eigenvector is a vector   
    l such that l**H * (A - w B) = 0, where l**H is the   
    conjugate-transpose of l.   

    Note: this routine performs "full balancing" on A and B -- see   
    "Further Details", below.   

    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N':  do not compute the left generalized eigenvectors;   
            = 'V':  compute the left generalized eigenvectors.   

    JOBVR   (input) CHARACTER*1   
            = 'N':  do not compute the right generalized eigenvectors;   
            = 'V':  compute the right generalized eigenvectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VL, and VR.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of A on exit, see "Further   
            Details", below.)   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of B on exit, see "Further   
            Details", below.)   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= MAX(1,N).   

    ALPHAR  (output) LONG DOUBLE PRECISION array, dimension (N)   
    ALPHAI  (output) LONG DOUBLE PRECISION array, dimension (N)   
    BETA    (output) LONG DOUBLE PRECISION array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  If ALPHAI(j) is zero, then   
            the j-th eigenvalue is real; if positive, then the j-th and   
            (j+1)-st eigenvalues are a complex conjugate pair, with   
            ALPHAI(j+1) negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero. 
  
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VL      (output) LONG DOUBLE PRECISION array, dimension (LDVL,N)   
            If JOBVL = 'V', the left generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have ABS(real part) + ABS(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVL = 'N'.   

    LDVL    (input) INTEGER   
            The leading dimension of the matrix VL. LDVL >= 1, and   
            if JOBVL = 'V', LDVL >= N.   

    VR      (output) LONG DOUBLE PRECISION array, dimension (LDVR,N)   
            If JOBVL = 'V', the right generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have ABS(real part) + ABS(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVR = 'N'.   

    LDVR    (input) INTEGER   
            The leading dimension of the matrix VR. LDVR >= 1, and   
            if JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= MAX(1,8*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR; 
  
            The optimal LWORK is:   
                2*N + MAX( 6*N, N*(NB+1) ).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  No eigenvectors have been   
                  calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)   
                  should be correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from DGGBAL   
                  =N+2: error return from DGEQRF   
                  =N+3: error return from DORMQR   
                  =N+4: error return from DORGQR   
                  =N+5: error return from DGGHRD   
                  =N+6: error return from DHGEQZ (other than failed   
                                                  iteration)   
                  =N+7: error return from DTGEVC   
                  =N+8: error return from DGGBAK (computing VL)   
                  =N+9: error return from DGGBAK (computing VR)   
                  =N+10: error return from DLASCL (various calls)   

    Further Details   
    ===============   

    Balancing   
    ---------   

    This driver calls DGGBAL to both permute and scale rows and columns   
    of A and B.  The permutations PL and PR are chosen so that PL*A*PR   
    and PL*B*R will be upper triangular except for the diagonal blocks   
    A(i:j,i:j) and B(i:j,i:j), with i and j as close together as   
    possible.  The diagonal scaling matrices DL and DR are chosen so   
    that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to   
    one (except for the elements that start out zero.)   

    After the eigenvalues and eigenvectors of the balanced matrices   
    have been computed, DGGBAK transforms the eigenvectors back to what   
    they would have been (in perfect arithmetic) if they had not been   
    balanced.   

    Contents of A and B on Exit   
    -------- -- - --- - -- ----   

    If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or   
    both), then on exit the arrays A and B will contain the real Schur   
    form[*] of the "balanced" versions of A and B.  If no eigenvectors   
    are computed, then only the diagonal blocks will be correct.   

    [*] See DHGEQZ, DGEGS, or read the book "Matrix Computations",   
        by Golub & van Loan, pub. by Johns Hopkins U. Press.   

    ===================================================================== 
  


       Decode the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c_n1 = -1;
    static LONG DOUBLE c_b14 = 1.;
    static LONG DOUBLE c_b25 = 0.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Local variables */
    static LONG DOUBLE absb, anrm, bnrm;
    static int itau;
    static LONG DOUBLE temp;
    static long int ilvl, ilvr;
    static LONG DOUBLE anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern long int lsame_(char *, char *);
    static int ileft, iinfo, icols, iwork, irows, jc;

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
    static int in;

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
    static int jr;
    static LONG DOUBLE salfai;

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
    static LONG DOUBLE salfar;

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

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE safmax;
    static char chtemp[1];
    static long int ldumma[1];

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

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dtgevc_(char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qtgevc(char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qtgevc_(char *, char *, 
#endif

	    long int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    int *, int *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);
    static int ijobvl, iright;
    static long int ilimit;
    static int ijobvr;

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
    static LONG DOUBLE onepls;
    static int lwkmin;

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
    static int lwkopt, ihi, ilo;
    static LONG DOUBLE eps;
    static long int ilv;



#define ALPHAR(I) alphar[(I)-1]
#define ALPHAI(I) alphai[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define VL(I,J) vl[(I)-1 + ((J)-1)* ( *ldvl)]
#define VR(I,J) vr[(I)-1 + ((J)-1)* ( *ldvr)]

    if (lsame_(jobvl, "N")) {
	ijobvl = 1;
	ilvl = 0;
    } else if (lsame_(jobvl, "V")) {
	ijobvl = 2;
	ilvl = 1;
    } else {
	ijobvl = -1;
	ilvl = 0;
    }

    if (lsame_(jobvr, "N")) {
	ijobvr = 1;
	ilvr = 0;
    } else if (lsame_(jobvr, "V")) {
	ijobvr = 2;
	ilvr = 1;
    } else {
	ijobvr = -1;
	ilvr = 0;
    }
    ilv = ilvl || ilvr;

/*     Test the input arguments   

   Computing MAX */
    i__1 = *n << 3;
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
    } else if (*ldvl < 1 || (ilvl && *ldvl < *n)) {
	*info = -12;
    } else if (*ldvr < 1 || (ilvr && *ldvr < *n)) {
	*info = -14;
    } else if (*lwork < lwkmin) {
	*info = -16;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEGV ", &i__1);
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

    safmin += safmin;
    safmax = 1. / safmin;
    onepls = eps * 4 + 1.;

/*     Scale A */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlange_("M", n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlange("M", n, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlange_("M", n, n, &A(1,1), lda, &WORK(1));
#endif

    anrm1 = anrm;
    anrm2 = 1.;
    if (anrm < 1.) {
	if (safmax * anrm < 1.) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }

    if (anrm > 0.) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &anrm, &c_b14, n, n, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &anrm, &c_b14, n, n, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &anrm, &c_b14, n, n, &A(1,1), lda, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return;
	}
    }

/*     Scale B */


#ifdef PETSC_PREFIX_SUFFIX
    bnrm = dlange_("M", n, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnrm = qlange("M", n, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnrm = qlange_("M", n, n, &B(1,1), ldb, &WORK(1));
#endif

    bnrm1 = bnrm;
    bnrm2 = 1.;
    if (bnrm < 1.) {
	if (safmax * bnrm < 1.) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }

    if (bnrm > 0.) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_("G", &c_n1, &c_n1, &bnrm, &c_b14, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl("G", &c_n1, &c_n1, &bnrm, &c_b14, n, n, &B(1,1), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_("G", &c_n1, &c_n1, &bnrm, &c_b14, n, n, &B(1,1), ldb, &
#endif

		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Workspace layout:  (8*N words -- "work" requires 6*N words)   
          left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dggbal_("B", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qggbal("B", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qggbal_("B", n, &A(1,1), lda, &B(1,1), ldb, &ilo, &ihi, &WORK(
#endif

	    ileft), &WORK(iright), &WORK(iwork), &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L120;
    }

/*     Reduce B to triangular form, and initialize VL and/or VR   
       Workspace layout:  ("work..." must have at least N words)   
          left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = *n + 1 - ilo;
    } else {
	icols = irows;
    }
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
	goto L120;
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
	goto L120;
    }

    if (ilvl) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b25, &c_b14, &VL(1,1), ldvl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b25, &c_b14, &VL(1,1), ldvl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b25, &c_b14, &VL(1,1), ldvl);
#endif

	i__1 = irows - 1;
	i__2 = irows - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VL(ilo+1,ilo), ldvl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VL(ilo+1,ilo), ldvl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("L", &i__1, &i__2, &B(ilo+1,ilo), ldb, &VL(ilo+1,ilo), ldvl);
#endif

	i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
	dorgqr_(&irows, &irows, &irows, &VL(ilo,ilo), ldvl, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgqr(&irows, &irows, &irows, &VL(ilo,ilo), ldvl, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgqr_(&irows, &irows, &irows, &VL(ilo,ilo), ldvl, &WORK(
#endif

		itau), &WORK(iwork), &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (int) WORK(iwork) + iwork - 1;
	    lwkopt = MAX(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L120;
	}
    }

    if (ilvr) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b25, &c_b14, &VR(1,1), ldvr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b25, &c_b14, &VR(1,1), ldvr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b25, &c_b14, &VR(1,1), ldvr);
#endif

    }

/*     Reduce to generalized Hessenberg form */

    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */


#ifdef PETSC_PREFIX_SUFFIX
	dgghrd_(jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgghrd(jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgghrd_(jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), 
#endif

		ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    } else {

#ifdef PETSC_PREFIX_SUFFIX
	dgghrd_("N", "N", &irows, &c__1, &irows, &A(ilo,ilo), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgghrd("N", "N", &irows, &c__1, &irows, &A(ilo,ilo), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgghrd_("N", "N", &irows, &c__1, &irows, &A(ilo,ilo), lda, 
#endif

		&B(ilo,ilo), ldb, &VL(1,1), ldvl, &VR(1,1), ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = *n + 5;
	goto L120;
    }

/*     Perform QZ algorithm   
       Workspace layout:  ("work..." must have at least 1 word)   
          left_permutation, right_permutation, work... */

    iwork = itau;
    if (ilv) {
	*(unsigned char *)chtemp = 'S';
    } else {
	*(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;

#ifdef PETSC_PREFIX_SUFFIX
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VL(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qhgeqz(chtemp, jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VL(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &A(1,1), lda, &B(1,1), ldb, &ALPHAR(1), &ALPHAI(1), &BETA(1), &VL(1,1), 
#endif

	    ldvl, &VR(1,1), ldvr, &WORK(iwork), &i__1, &iinfo);
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
	goto L120;
    }

    if (ilv) {

/*        Compute Eigenvectors  (DTGEVC requires 6*N words of workspac
e) */

	if (ilvl) {
	    if (ilvr) {
		*(unsigned char *)chtemp = 'B';
	    } else {
		*(unsigned char *)chtemp = 'L';
	    }
	} else {
	    *(unsigned char *)chtemp = 'R';
	}


#ifdef PETSC_PREFIX_SUFFIX
	dtgevc_(chtemp, "B", ldumma, n, &A(1,1), lda, &B(1,1), ldb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtgevc(chtemp, "B", ldumma, n, &A(1,1), lda, &B(1,1), ldb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtgevc_(chtemp, "B", ldumma, n, &A(1,1), lda, &B(1,1), ldb, 
#endif

		&VL(1,1), ldvl, &VR(1,1), ldvr, n, &in, &WORK(
		iwork), &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L120;
	}

/*        Undo balancing on VL and VR, rescale */

	if (ilvl) {

#ifdef PETSC_PREFIX_SUFFIX
	    dggbak_("B", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qggbak("B", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qggbak_("B", "L", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif

		    VL(1,1), ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 8;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		if (ALPHAI(jc) < 0.) {
		    goto L50;
		}
		temp = 0.;
		if (ALPHAI(jc) == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = VL(jr,jc), 
				ABS(d__1));
			temp = MAX(d__2,d__3);
/* L10: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = VL(jr,jc), 
				ABS(d__1)) + (d__2 = VL(jr,jc+1), ABS(d__2));
			temp = MAX(d__3,d__4);
/* L20: */
		    }
		}
		if (temp < safmin) {
		    goto L50;
		}
		temp = 1. / temp;
		if (ALPHAI(jc) == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VL(jr,jc) *= temp;
/* L30: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VL(jr,jc) *= temp;
			VL(jr,jc+1) *= temp;
/* L40: */
		    }
		}
L50:
		;
	    }
	}
	if (ilvr) {

#ifdef PETSC_PREFIX_SUFFIX
	    dggbak_("B", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qggbak("B", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qggbak_("B", "R", n, &ilo, &ihi, &WORK(ileft), &WORK(iright), n, &
#endif

		    VR(1,1), ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 9;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= *n; ++jc) {
		if (ALPHAI(jc) < 0.) {
		    goto L100;
		}
		temp = 0.;
		if (ALPHAI(jc) == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = VR(jr,jc), 
				ABS(d__1));
			temp = MAX(d__2,d__3);
/* L60: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = VR(jr,jc), 
				ABS(d__1)) + (d__2 = VR(jr,jc+1), ABS(d__2));
			temp = MAX(d__3,d__4);
/* L70: */
		    }
		}
		if (temp < safmin) {
		    goto L100;
		}
		temp = 1. / temp;
		if (ALPHAI(jc) == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,jc) *= temp;
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			VR(jr,jc) *= temp;
			VR(jr,jc+1) *= temp;
/* L90: */
		    }
		}
L100:
		;
	    }
	}

/*        End of eigenvector calculation */

    }

/*     Undo scaling in alpha, beta   

       Note: this does not give the alpha and beta for the unscaled   
       problem.   

       Un-scaling is limited to avoid underflow in alpha and beta   
       if they are significant. */

    i__1 = *n;
    for (jc = 1; jc <= *n; ++jc) {
	absar = (d__1 = ALPHAR(jc), ABS(d__1));
	absai = (d__1 = ALPHAI(jc), ABS(d__1));
	absb = (d__1 = BETA(jc), ABS(d__1));
	salfar = anrm * ALPHAR(jc);
	salfai = anrm * ALPHAI(jc);
	sbeta = bnrm * BETA(jc);
	ilimit = 0;
	scale = 1.;

/*        Check for significant underflow in ALPHAI   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = MAX(d__1,d__2), d__2 = eps *
		 absb;
	if (ABS(salfai) < safmin && absai >= MAX(d__1,d__2)) {
	    ilimit = 1;
/* Computing MAX */
	    d__1 = onepls * safmin, d__2 = anrm2 * absai;
	    scale = onepls * safmin / anrm1 / MAX(d__1,d__2);

	} else if (salfai == 0.) {

/*           If insignificant underflow in ALPHAI, then make the 
  
             conjugate eigenvalue real. */

	    if (ALPHAI(jc) < 0. && jc > 1) {
		ALPHAI(jc - 1) = 0.;
	    } else if (ALPHAI(jc) > 0. && jc < *n) {
		ALPHAI(jc + 1) = 0.;
	    }
	}

/*        Check for significant underflow in ALPHAR   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absai, d__1 = MAX(d__1,d__2), d__2 = eps *
		 absb;
	if (ABS(salfar) < safmin && absar >= MAX(d__1,d__2)) {
	    ilimit = 1;
/* Computing MAX   
   Computing MAX */
	    d__3 = onepls * safmin, d__4 = anrm2 * absar;
	    d__1 = scale, d__2 = onepls * safmin / anrm1 / MAX(d__3,d__4);
	    scale = MAX(d__1,d__2);
	}

/*        Check for significant underflow in BETA   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = MAX(d__1,d__2), d__2 = eps *
		 absai;
	if (ABS(sbeta) < safmin && absb >= MAX(d__1,d__2)) {
	    ilimit = 1;
/* Computing MAX   
   Computing MAX */
	    d__3 = onepls * safmin, d__4 = bnrm2 * absb;
	    d__1 = scale, d__2 = onepls * safmin / bnrm1 / MAX(d__3,d__4);
	    scale = MAX(d__1,d__2);
	}

/*        Check for possible overflow when limiting scaling */

	if (ilimit) {
/* Computing MAX */
	    d__1 = ABS(salfar), d__2 = ABS(salfai), d__1 = MAX(d__1,d__2), 
		    d__2 = ABS(sbeta);
	    temp = scale * safmin * MAX(d__1,d__2);
	    if (temp > 1.) {
		scale /= temp;
	    }
	    if (scale < 1.) {
		ilimit = 0;
	    }
	}

/*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */

	if (ilimit) {
	    salfar = scale * ALPHAR(jc) * anrm;
	    salfai = scale * ALPHAI(jc) * anrm;
	    sbeta = scale * BETA(jc) * bnrm;
	}
	ALPHAR(jc) = salfar;
	ALPHAI(jc) = salfai;
	BETA(jc) = sbeta;
/* L110: */
    }

L120:
    WORK(1) = (LONG DOUBLE) lwkopt;

    return;

/*     End of DGEGV */

} /* dgegv_ */

