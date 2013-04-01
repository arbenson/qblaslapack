#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dstevx_(char *jobz, char *range, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qstevx(char *jobz, char *range, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qstevx_(char *jobz, char *range, int *n, LONG DOUBLE *
#endif

	d, LONG DOUBLE *e, LONG DOUBLE *vl, LONG DOUBLE *vu, int *il, 
	int *iu, LONG DOUBLE *abstol, int *m, LONG DOUBLE *w, 
	LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int *iwork, 
	int *ifail, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEVX computes selected eigenvalues and, optionally, eigenvectors   
    of a real symmetric tridiagonal matrix A.  Eigenvalues and   
    eigenvectors can be selected by specifying either a range of values   
    or a range of indices for the desired eigenvalues.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    RANGE   (input) CHARACTER*1   
            = 'A': all eigenvalues will be found.   
            = 'V': all eigenvalues in the half-open interval (VL,VU]   
                   will be found.   
            = 'I': the IL-th through IU-th eigenvalues will be found.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.   
            On exit, D may be multiplied by a constant factor chosen   
            to avoid over/underflow in computing the eigenvalues.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A in elements 1 to N-1 of E; E(N) need not be set.   
            On exit, E may be multiplied by a constant factor chosen   
            to avoid over/underflow in computing the eigenvalues.   

    VL      (input) LONG DOUBLE PRECISION   
    VU      (input) LONG DOUBLE PRECISION   
            If RANGE='V', the lower and upper bounds of the interval to   
            be searched for eigenvalues. VL < VU.   
            Not referenced if RANGE = 'A' or 'I'.   

    IL      (input) INTEGER   
    IU      (input) INTEGER   
            If RANGE='I', the indices (in ascending order) of the   
            smallest and largest eigenvalues to be returned.   
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.   
            Not referenced if RANGE = 'A' or 'V'.   

    ABSTOL  (input) LONG DOUBLE PRECISION   
            The absolute error tolerance for the eigenvalues.   
            An approximate eigenvalue is accepted as converged   
            when it is determined to lie in an interval [a,b]   
            of width less than or equal to   

                    ABSTOL + EPS *   MAX( |a|,|b| ) ,   

            where EPS is the machine precision.  If ABSTOL is less   
            than or equal to zero, then  EPS*|T|  will be used in   
            its place, where |T| is the 1-norm of the tridiagonal   
            matrix.   

            Eigenvalues will be computed most accurately when ABSTOL is   
            set to twice the underflow threshold 2*DLAMCH('S'), not zero. 
  
            If this routine returns with INFO>0, indicating that some   
            eigenvectors did not converge, try setting ABSTOL to   
            2*DLAMCH('S').   

            See "Computing Small Singular Values of Bidiagonal Matrices   
            with Guaranteed High Relative Accuracy," by Demmel and   
            Kahan, LAPACK Working Note #3.   

    M       (output) INTEGER   
            The total number of eigenvalues found.  0 <= M <= N.   
            If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.   

    W       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The first M elements contain the selected eigenvalues in   
            ascending order.   

    Z       (output) LONG DOUBLE PRECISION array, dimension (LDZ, MAX(1,M) )   
            If JOBZ = 'V', then if INFO = 0, the first M columns of Z   
            contain the orthonormal eigenvectors of the matrix A   
            corresponding to the selected eigenvalues, with the i-th   
            column of Z holding the eigenvector associated with W(i).   
            If an eigenvector fails to converge (INFO > 0), then that   
            column of Z contains the latest approximation to the   
            eigenvector, and the index of the eigenvector is returned   
            in IFAIL.  If JOBZ = 'N', then Z is not referenced.   
            Note: the user must ensure that at least MAX(1,M) columns are 
  
            supplied in the array Z; if RANGE = 'V', the exact value of M 
  
            is not known in advance and an upper bound must be used.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (5*N)   

    IWORK   (workspace) INTEGER array, dimension (5*N)   

    IFAIL   (output) INTEGER array, dimension (N)   
            If JOBZ = 'V', then if INFO = 0, the first M elements of   
            IFAIL are zero.  If INFO > 0, then IFAIL contains the   
            indices of the eigenvectors that failed to converge.   
            If JOBZ = 'N', then IFAIL is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, then i eigenvectors failed to converge.   
                  Their indices are stored in array IFAIL.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static int imax;
    static LONG DOUBLE rmin, rmax, tnrm;
    static int itmp1, i, j;

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
    static char order[1];

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dswap_(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap_(int *, LONG DOUBLE *, int 
#endif

	    *, LONG DOUBLE *, int *);
    static long int wantz;
    static int jj;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static long int alleig, indeig;
    static int iscale, indibl;
    static long int valeig;
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

    static int indisp;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dstein_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qstein(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qstein_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dsterf_(int *, LONG DOUBLE *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsterf(int *, LONG DOUBLE *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsterf_(int *, LONG DOUBLE *, LONG DOUBLE *, int *);
#endif

    static int indiwo;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dstebz_(char *, char *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qstebz(char *, char *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qstebz_(char *, char *, int *, LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    int *, LONG DOUBLE *, int *, int *);
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
    static int nsplit;
    static LONG DOUBLE smlnum, eps, vll, vuu, tmp1;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]
#define IFAIL(I) ifail[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && *il < 1) {
	*info = -8;
    } else if (indeig && (*iu < MIN(*n,*il) || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -14;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEVX", &i__1);
	return;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    W(1) = D(1);
	} else {
	    if (*vl < D(1) && *vu >= D(1)) {
		*m = 1;
		W(1) = D(1);
	    }
	}
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
/* Computing MIN */
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = MIN(d__1,d__2);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    if (valeig) {
	vll = *vl;
	vuu = *vu;
    }

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

	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     If all eigenvalues are desired and ABSTOL is less than zero, then 
  
       call DSTERF or SSTEQR.  If this fails for some eigenvalue, then   
       try DSTEBZ. */

    if ((alleig || (indeig && *il == 1 && *iu == *n)) && *abstol <= 0.) {

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &D(1), &c__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &D(1), &c__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &D(1), &c__1, &W(1), &c__1);
#endif

	i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(&i__1, &E(1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(&i__1, &E(1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(&i__1, &E(1), &c__1, &WORK(1), &c__1);
#endif

	indwrk = *n + 1;
	if (! wantz) {

#ifdef PETSC_PREFIX_SUFFIX
	    dsterf_(n, &W(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsterf(n, &W(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsterf_(n, &W(1), &WORK(1), info);
#endif

	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_("I", n, &W(1), &WORK(1), &Z(1,1), ldz, &WORK(indwrk),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr("I", n, &W(1), &WORK(1), &Z(1,1), ldz, &WORK(indwrk),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_("I", n, &W(1), &WORK(1), &Z(1,1), ldz, &WORK(indwrk),
#endif

		     info);
	    if (*info == 0) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    IFAIL(i) = 0;
/* L10: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L20;
	}
	*info = 0;
    }

/*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indwrk = 1;
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &D(1), &E(1), m, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qstebz(range, order, n, &vll, &vuu, il, iu, abstol, &D(1), &E(1), m, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qstebz_(range, order, n, &vll, &vuu, il, iu, abstol, &D(1), &E(1), m, &
#endif

	    nsplit, &W(1), &IWORK(indibl), &IWORK(indisp), &WORK(indwrk), &
	    IWORK(indiwo), info);

    if (wantz) {

#ifdef PETSC_PREFIX_SUFFIX
	dstein_(n, &D(1), &E(1), m, &W(1), &IWORK(indibl), &IWORK(indisp), &Z(1,1), ldz, &WORK(indwrk), &IWORK(indiwo), &IFAIL(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qstein(n, &D(1), &E(1), m, &W(1), &IWORK(indibl), &IWORK(indisp), &Z(1,1), ldz, &WORK(indwrk), &IWORK(indiwo), &IFAIL(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qstein_(n, &D(1), &E(1), m, &W(1), &IWORK(indibl), &IWORK(indisp), &Z(1,1), ldz, &WORK(indwrk), &IWORK(indiwo), &IFAIL(1), 
#endif

		info);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L20:
    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
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

/*     If eigenvalues are not in order, then sort them, along with   
       eigenvectors. */

    if (wantz) {
	i__1 = *m - 1;
	for (j = 1; j <= *m-1; ++j) {
	    i = 0;
	    tmp1 = W(j);
	    i__2 = *m;
	    for (jj = j + 1; jj <= *m; ++jj) {
		if (W(jj) < tmp1) {
		    i = jj;
		    tmp1 = W(jj);
		}
/* L30: */
	    }

	    if (i != 0) {
		itmp1 = IWORK(indibl + i - 1);
		W(i) = W(j);
		IWORK(indibl + i - 1) = IWORK(indibl + j - 1);
		W(j) = tmp1;
		IWORK(indibl + j - 1) = itmp1;

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(n, &Z(1,i), &c__1, &Z(1,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(n, &Z(1,i), &c__1, &Z(1,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(n, &Z(1,i), &c__1, &Z(1,j), &
#endif

			c__1);
		if (*info != 0) {
		    itmp1 = IFAIL(i);
		    IFAIL(i) = IFAIL(j);
		    IFAIL(j) = itmp1;
		}
	    }
/* L40: */
	}
    }

    return;

/*     End of DSTEVX */

} /* dstevx_ */

