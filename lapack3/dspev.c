#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspev_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspev(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspev_(char *jobz, char *uplo, int *n, LONG DOUBLE *
#endif

	ap, LONG DOUBLE *w, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPEV computes all the eigenvalues and, optionally, eigenvectors of a 
  
    real symmetric matrix A in packed storage.   

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

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int inde;
    static LONG DOUBLE anrm;
    static int imax;
    static LONG DOUBLE rmin, rmax;

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
    static int indwrk;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dopgtr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qopgtr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qopgtr_(char *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), dsptrd_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qsptrd(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qsptrd_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dsteqr_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qsteqr(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qsteqr_(char *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static LONG DOUBLE smlnum, eps;



#define AP(I) ap[(I)-1]
#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lsame_(uplo, "U") || lsame_(uplo, "L"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -7;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPEV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	W(1) = AP(1);
	if (wantz) {
	    Z(1,1) = 1.;
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
       DOPGTR to generate the orthogonal matrix, then call DSTEQR. */

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

#ifdef PETSC_PREFIX_SUFFIX
	dopgtr_(uplo, n, &AP(1), &WORK(indtau), &Z(1,1), ldz, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qopgtr(uplo, n, &AP(1), &WORK(indtau), &Z(1,1), ldz, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qopgtr_(uplo, n, &AP(1), &WORK(indtau), &Z(1,1), ldz, &WORK(
#endif

		indwrk), &iinfo);

#ifdef PETSC_PREFIX_SUFFIX
	dsteqr_(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indtau),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsteqr(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indtau),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsteqr_(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indtau),
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

    return;

/*     End of DSPEV */

} /* dspev_ */

