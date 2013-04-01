#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dstevd_(char *jobz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qstevd(char *jobz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qstevd_(char *jobz, int *n, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int 
	*lwork, int *iwork, int *liwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric tridiagonal matrix. If eigenvectors are desired, it   
    uses a divide and conquer algorithm.   

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

    WORK    (workspace/output) LONG DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1 then LWORK must be at least   
                           ( 1 + 3*N + 2*N*lg N + 2*N**2 ),   
                           where lg( N ) = smallest int k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1 then LIWORK must be at least 2+5*N. 
  

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
    static int c__2 = 2;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
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
    static int lwmin;
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

	     int *);
    static int liwmin;
    static LONG DOUBLE smlnum;
    static int lgn;
    static LONG DOUBLE eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");

    *info = 0;
    liwmin = 1;
    lwmin = 1;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -6;
    } else if (*n > 1 && wantz) {
	lgn = (int) (log((LONG DOUBLE) (*n)) / log(2.));
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
/* Computing 2nd power */
	i__1 = *n;
	lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 1);
	liwmin = *n * 5 + 2;
	if (*lwork < lwmin) {
	    *info = -8;
	} else if (*liwork < liwmin) {
	    *info = -10;
	}
    } else if (*lwork < 1) {
	*info = -8;
    } else if (*liwork < 1) {
	*info = -10;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEVD", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
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
       eigenvectors, call DSTEDC. */

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
	dstedc_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), lwork, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qstedc("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), lwork, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qstedc_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), lwork, &
#endif

		IWORK(1), liwork, info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
	d__1 = 1. / sigma;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &d__1, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &d__1, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &d__1, &D(1), &c__1);
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

/*     End of DSTEVD */

} /* dstevd_ */

