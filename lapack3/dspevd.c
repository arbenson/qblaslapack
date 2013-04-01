#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspevd_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspevd(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspevd_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif

	ap, LONG DOUBLE *w, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, 
	int *lwork, int *iwork, int *liwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSPEVD computes all the eigenvalues and, optionally, eigenvectors   
    of a real symmetric matrix A in packed storage. If eigenvectors are   
    desired, it uses a divide and conquer algorithm.   

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

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  

            On exit, AP is overwritten by values generated during the   
            reduction to tridiagonal form.  If UPLO = 'U', the diagonal   
            and first superdiagonal of the tridiagonal matrix T overwrite 
  
            the corresponding elements of A, and if UPLO = 'L', the   
            diagonal and first subdiagonal of T overwrite the   
            corresponding elements of A.   

    W       (output) LONG DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    Z       (output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with W(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= MAX(1,N).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If N <= 1,               LWORK must be at least 1.   
            If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.   
            If JOBZ = 'V' and N > 1, LWORK must be at least   
                           ( 1 + 5*N + 2*N*lg N + 2*N**2 ),   
                           where lg( N ) = smallest int k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 2 + 5*N.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__2 = 2;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int inde;
    static LONG DOUBLE anrm, rmin, rmax;

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
    static int iinfo, lwmin;
    static long int wantz;

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
    extern /* Subroutine */ void dstedc_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qstedc(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qstedc_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    int *, int *, int *);
    static LONG DOUBLE safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
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
    static int indwrk, liwmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsptrd_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsptrd(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsptrd_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dopmtr_(char *, char *, char *, int *, int *, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qopmtr(char *, char *, char *, int *, int *, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qopmtr_(char *, char *, char *, int *, int *, LONG DOUBLE *
#endif

	    , LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static int llwork;
    static LONG DOUBLE smlnum;
    static int lgn;
    static LONG DOUBLE eps;



#define AP(I) ap[(I)-1]
#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	liwmin = 1;
	lwmin = 1;
    } else {
	lgn = (int) (log((LONG DOUBLE) (*n)) / log(2.));
	if (pow(c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (pow(c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (wantz) {
	    liwmin = *n * 5 + 2;
/* Computing 2nd power */
	    i__1 = *n;
	    lwmin = *n * 5 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 1);
	} else {
	    liwmin = 1;
	    lwmin = *n << 1;
	}
    }
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lsame_(uplo, "U") || lsame_(uplo, "L"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -7;
    } else if (*lwork < lwmin) {
	*info = -9;
    } else if (*liwork < liwmin) {
	*info = -11;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	W(1) = AP(1);
	if (wantz) {
	    Z(1,1) = 1.;
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
    anrm = dlansp_("M", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlansp("M", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlansp_("M", uplo, n, &AP(1), &WORK(1));
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
	i__1 = *n * (*n + 1) / 2;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__1, &sigma, &AP(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__1, &sigma, &AP(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__1, &sigma, &AP(1), &c__1);
#endif

    }

/*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. 
*/

    inde = 1;
    indtau = inde + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dsptrd_(uplo, n, &AP(1), &W(1), &WORK(inde), &WORK(indtau), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsptrd(uplo, n, &AP(1), &W(1), &WORK(inde), &WORK(indtau), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsptrd_(uplo, n, &AP(1), &W(1), &WORK(inde), &WORK(indtau), &iinfo);
#endif


/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call   
       DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the   
       tridiagonal matrix, then call DOPMTR to multiply it by the   
       Householder transformations represented in AP. */

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
	indwrk = indtau + *n;
	llwork = *lwork - indwrk + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dstedc_("I", n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qstedc("I", n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qstedc_("I", n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), 
#endif

		&llwork, &IWORK(1), liwork, info);

#ifdef PETSC_PREFIX_SUFFIX
	dopmtr_("L", uplo, "N", n, n, &AP(1), &WORK(indtau), &Z(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qopmtr("L", uplo, "N", n, n, &AP(1), &WORK(indtau), &Z(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qopmtr_("L", uplo, "N", n, n, &AP(1), &WORK(indtau), &Z(1,1), 
#endif

		ldz, &WORK(indwrk), &iinfo);
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
	WORK(1) = (LONG DOUBLE) lwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return;

/*     End of DSPEVD */

} /* dspevd_ */

