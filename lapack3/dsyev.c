#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsyev_(char *jobz, char *uplo, int *n, LONG DOUBLE *a,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsyev(char *jobz, char *uplo, int *n, LONG DOUBLE *a,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsyev_(char *jobz, char *uplo, int *n, LONG DOUBLE *a,
#endif

	 int *lda, LONG DOUBLE *w, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYEV computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A.   

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

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= MAX(1,3*N-1).   
            For optimal efficiency, LWORK >= (NB+2)*N,   
            where NB is the blocksize for DSYTRD returned by ILAENV.   

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
    static int c__0 = 0;
    static LONG DOUBLE c_b12 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int inde;
    static LONG DOUBLE anrm;
    static int imax;
    static LONG DOUBLE rmin, rmax;
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
    static int iinfo;
    static long int lower, wantz;

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
	    int *, int *);
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
    static int indwrk;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorgtr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgtr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgtr_(char *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), dsteqr_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qsteqr(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qsteqr_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dsytrd_(char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsytrd(char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsytrd_(char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static int llwork;
    static LONG DOUBLE smlnum, eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 3 - 1;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -8;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYEV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return;
    }

    if (*n == 1) {
	W(1) = A(1,1);
	WORK(1) = 3.;
	if (wantz) {
	    A(1,1) = 1.;
	}
	return;
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
	dlascl_(uplo, &c__0, &c__0, &c_b12, &sigma, n, n, &A(1,1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlascl(uplo, &c__0, &c__0, &c_b12, &sigma, n, n, &A(1,1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlascl_(uplo, &c__0, &c__0, &c_b12, &sigma, n, n, &A(1,1), lda, 
#endif

		info);
    }

/*     Call DSYTRD to reduce symmetric matrix to tridiagonal form. */

    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;

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
       DORGTR to generate the orthogonal matrix, then call DSTEQR. */

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
	dorgtr_(uplo, n, &A(1,1), lda, &WORK(indtau), &WORK(indwrk), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgtr(uplo, n, &A(1,1), lda, &WORK(indtau), &WORK(indwrk), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgtr_(uplo, n, &A(1,1), lda, &WORK(indtau), &WORK(indwrk), &
#endif

		llwork, &iinfo);

#ifdef PETSC_PREFIX_SUFFIX
	dsteqr_(jobz, n, &W(1), &WORK(inde), &A(1,1), lda, &WORK(indtau),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsteqr(jobz, n, &W(1), &WORK(inde), &A(1,1), lda, &WORK(indtau),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsteqr_(jobz, n, &W(1), &WORK(inde), &A(1,1), lda, &WORK(indtau),
#endif

		 info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *n;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&imax, &d__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&imax, &d__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&imax, &d__1, &W(1), &c__1);
#endif

    }

/*     Set WORK(1) to optimal workspace size.   

   Computing MAX */
    i__1 = *n * 3 - 1;
    WORK(1) = (LONG DOUBLE) MAX(i__1,lopt);

    return;

/*     End of DSYEV */

} /* dsyev_ */

