#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsbev_(char *jobz, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsbev(char *jobz, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsbev_(char *jobz, char *uplo, int *n, int *kd, 
#endif

	LONG DOUBLE *ab, int *ldab, LONG DOUBLE *w, LONG DOUBLE *z, int *
	ldz, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSBEV computes all the eigenvalues and, optionally, eigenvectors of   
    a real symmetric band matrix A.   

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

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB, N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  

            On exit, AB is overwritten by values generated during the   
            reduction to tridiagonal form.  If UPLO = 'U', the first   
            superdiagonal and the diagonal of the tridiagonal matrix T   
            are returned in rows KD and KD+1 of AB, and if UPLO = 'L',   
            the diagonal and first subdiagonal of T are returned in the   
            first two rows of AB.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD + 1.   

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

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (MAX(1,3*N-2)) 
  

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
    static LONG DOUBLE c_b11 = 1.;
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

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlansb_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlansb(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlansb_(char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *);
    static LONG DOUBLE safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsbtrd_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsbtrd(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsbtrd_(char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *), dsterf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *), qsterf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *, LONG DOUBLE *, int *), qsterf_(
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int indwrk;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsteqr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsteqr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsteqr_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static LONG DOUBLE smlnum, eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*kd < 0) {
	*info = -4;
    } else if (*ldab < *kd + 1) {
	*info = -6;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -9;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSBEV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	W(1) = AB(1,1);
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
    anrm = dlansb_("M", uplo, n, kd, &AB(1,1), ldab, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlansb("M", uplo, n, kd, &AB(1,1), ldab, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlansb_("M", uplo, n, kd, &AB(1,1), ldab, &WORK(1));
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
	if (lower) {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("B", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("B", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("B", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif

		    info);
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("Q", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("Q", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("Q", kd, kd, &c_b11, &sigma, n, n, &AB(1,1), ldab, 
#endif

		    info);
	}
    }

/*     Call DSBTRD to reduce symmetric band matrix to tridiagonal form. */

    inde = 1;
    indwrk = inde + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dsbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsbtrd(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif


/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR. 
*/

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
	dsteqr_(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsteqr(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsteqr_(jobz, n, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk),
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

/*     End of DSBEV */

} /* dsbev_ */

