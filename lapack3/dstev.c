#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dstev_(char *jobz, int *n, LONG DOUBLE *d, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qstev(char *jobz, int *n, LONG DOUBLE *d, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qstev_(char *jobz, int *n, LONG DOUBLE *d, LONG DOUBLE 
#endif

	*e, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEV computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric tridiagonal matrix A.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A, stored in elements 1 to N-1 of E; E(N) need not   
            be set, but is used by the routine.   
            On exit, the contents of E are destroyed.   

    Z       (output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal   
            eigenvectors of the matrix A, with the i-th column of Z   
            holding the eigenvector associated with D(i).   
            If JOBZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (MAX(1,2*N-2)) 
  
            If JOBZ = 'N', WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of E did not converge to zero.   

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
    static int imax;
    static LONG DOUBLE rmin, rmax, tnrm;

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
    extern LONG DOUBLE dlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif


#ifdef PETSC_PREFIX_SUFFIX
	     int *), dsteqr_(char *, int *, LONG DOUBLE *, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *), qsteqr(char *, int *, LONG DOUBLE *, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *), qsteqr_(char *, int *, LONG DOUBLE *, LONG DOUBLE *
#endif

	    , LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static LONG DOUBLE smlnum, eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -6;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEV ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
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

    iscale = 0;

#ifdef PETSC_PREFIX_SUFFIX
    tnrm = dlanst_("M", n, &D(1), &E(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    tnrm = qlanst("M", n, &D(1), &E(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    tnrm = qlanst_("M", n, &D(1), &E(1));
#endif

    if (tnrm > 0. && tnrm < rmin) {
	iscale = 1;
	sigma = rmin / tnrm;
    } else if (tnrm > rmax) {
	iscale = 1;
	sigma = rmax / tnrm;
    }
    if (iscale == 1) {

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &sigma, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &sigma, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &sigma, &D(1), &c__1);
#endif

	i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__1, &sigma, &E(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__1, &sigma, &E(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__1, &sigma, &E(1), &c__1);
#endif

    }

/*     For eigenvalues only, call DSTERF.  For eigenvalues and   
       eigenvectors, call DSTEQR. */

    if (! wantz) {

#ifdef PETSC_PREFIX_SUFFIX
	dsterf_(n, &D(1), &E(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsterf(n, &D(1), &E(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsterf_(n, &D(1), &E(1), info);
#endif

    } else {

#ifdef PETSC_PREFIX_SUFFIX
	dsteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsteqr("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif

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
	dscal_(&imax, &d__1, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&imax, &d__1, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&imax, &d__1, &D(1), &c__1);
#endif

    }

    return;

/*     End of DSTEV */

} /* dstev_ */

