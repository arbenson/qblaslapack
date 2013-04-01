#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsyevd_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsyevd(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsyevd_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *w, LONG DOUBLE *work, int *lwork, 
	int *iwork, int *liwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A. If eigenvectors are desired, it uses a   
    divide and conquer algorithm.   

    The divide and conquer algorithm makes very mild assumptions about   
    floating point arithmetic. It will work on machines with a guard   
    digit in add/subtract, or on those binary machines without guard   
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or   
    Cray-2. It could conceivably fail on hexadecimal or decimal machines 
  
    without guard digits, but we know of none.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') 
  
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    W       (output) LONG DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If N <= 1,               LWORK must be at least 1.   
            If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.   
            If JOBZ = 'V' and N > 1, LWORK must be at least   
                           1 + 5*N + 2*N*lg N + 3*N**2,   
                           where lg( N ) = smallest int k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 2 + 5*N.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__2 = 2;
    static int c__0 = 0;
    static LONG DOUBLE c_b15 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int inde;
    static LONG DOUBLE anrm, rmin, rmax;
    static int lopt;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE sigma;
    extern long int lsame_(char *, char *);
    static int iinfo, lwmin, liopt;
    static long int lower, wantz;
    static int indwk2, llwrk2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int iscale;

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
	    int *, int *), dstedc_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qstedc(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qstedc_(char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *, int *, int *, int *), dlacpy_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *, int *, int *, int *), qlacpy(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *, int *, int *, int *), qlacpy_(
#endif

	    char *, int *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
	    *, int *);
    static LONG DOUBLE safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum;
    static int indtau;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlansy_(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlansy(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlansy_(char *, char *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *);
    static int indwrk, liwmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormtr_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormtr(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormtr_(char *, char *, char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, int *), dsytrd_(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, int *), qsytrd(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, int *), qsytrd_(char *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *,
	     int *);
    static int llwork;
    static LONG DOUBLE smlnum;
    static int lgn;
    static LONG DOUBLE eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	liwmin = 1;
	lwmin = 1;
	lopt = lwmin;
	liopt = liwmin;
    } else {
	lgn = (int) (log((LONG DOUBLE) (*n)) / log(2.));
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (wantz) {
	    liwmin = *n * 5 + 2;
/* Computing 2nd power */
	    i__1 = *n;
	    lwmin = *n * 5 + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	} else {
	    liwmin = 1;
	    lwmin = (*n << 1) + 1;
	}
	lopt = lwmin;
	liopt = liwmin;
    }
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*lwork < lwmin) {
	*info = -8;
    } else if (*liwork < liwmin) {
	*info = -10;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	W(1) = A(1,1);
	if (wantz) {
	    A(1,1) = 1.;
	}
	goto L10;
    }

/*     Get machine constants. */


#ifdef PETSC_PREFIX_SUFFIX
    safmin = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    safmin = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    safmin = qlamch_("Safe minimum");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Precision");
#endif

    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */


#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlansy_("M", uplo, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlansy("M", uplo, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlansy_("M", uplo, n, &A(1,1), lda, &WORK(1));
#endif

    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm;
    }
    if (iscale == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dlascl_(uplo, &c__0, &c__0, &c_b15, &sigma, n, n, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl(uplo, &c__0, &c__0, &c_b15, &sigma, n, n, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_(uplo, &c__0, &c__0, &c_b15, &sigma, n, n, &A(1,1), lda, 
#endif

		info);
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;


#ifdef PETSC_PREFIX_SUFFIX
    dsytrd_(uplo, n, &A(1,1), lda, &W(1), &WORK(inde), &WORK(indtau), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsytrd(uplo, n, &A(1,1), lda, &W(1), &WORK(inde), &WORK(indtau), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsytrd_(uplo, n, &A(1,1), lda, &W(1), &WORK(inde), &WORK(indtau), &
#endif

	    WORK(indwrk), &llwork, &iinfo);
    lopt = (int) ((*n << 1) + WORK(indwrk));

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call DORMTR to multiply it by the   
       Householder transformations stored in A. */

    if (! wantz) {

#ifdef PETSC_PREFIX_SUFFIX
	dsterf_(n, &W(1), &WORK(inde), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsterf(n, &W(1), &WORK(inde), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsterf_(n, &W(1), &WORK(inde), info);
#endif

    } else {

#ifdef PETSC_PREFIX_SUFFIX
	dstedc_("I", n, &W(1), &WORK(inde), &WORK(indwrk), n, &WORK(indwk2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qstedc("I", n, &W(1), &WORK(inde), &WORK(indwrk), n, &WORK(indwk2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qstedc_("I", n, &W(1), &WORK(inde), &WORK(indwrk), n, &WORK(indwk2), &
#endif

		llwrk2, &IWORK(1), liwork, info);

#ifdef PETSC_PREFIX_SUFFIX
	dormtr_("L", uplo, "N", n, n, &A(1,1), lda, &WORK(indtau), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qormtr("L", uplo, "N", n, n, &A(1,1), lda, &WORK(indtau), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qormtr_("L", uplo, "N", n, n, &A(1,1), lda, &WORK(indtau), &WORK(
#endif

		indwrk), n, &WORK(indwk2), &llwrk2, &iinfo);

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("A", n, n, &WORK(indwrk), n, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("A", n, n, &WORK(indwrk), n, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("A", n, n, &WORK(indwrk), n, &A(1,1), lda);
#endif

/* Computing MAX   
   Computing 2nd power */
	i__3 = *n;
	i__1 = lopt, i__2 = *n * 5 + 1 + (*n << 1) * lgn + i__3 * i__3 * 3;
	lopt = MAX(i__1,i__2);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	d__1 = 1. / sigma;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__1, &W(1), &c__1);
#endif

    }
L10:
    if (*lwork > 0) {
	WORK(1) = (LONG DOUBLE) lopt;
    }
    if (*liwork > 0) {
	IWORK(1) = liopt;
    }
    return;

/*     End of DSYEVD */

} /* dsyevd_ */

