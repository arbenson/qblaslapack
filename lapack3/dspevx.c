#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspevx_(char *jobz, char *range, char *uplo, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspevx(char *jobz, char *range, char *uplo, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspevx_(char *jobz, char *range, char *uplo, int *n, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *vl, LONG DOUBLE *vu, int *il, int *
	iu, LONG DOUBLE *abstol, int *m, LONG DOUBLE *w, LONG DOUBLE *z, 
	int *ldz, LONG DOUBLE *work, int *iwork, int *ifail, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSPEVX computes selected eigenvalues and, optionally, eigenvectors   
    of a real symmetric matrix A in packed storage.  Eigenvalues/vectors 
  
    can be selected by specifying either a range of values or a range of 
  
    indices for the desired eigenvalues.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    RANGE   (input) CHARACTER*1   
            = 'A': all eigenvalues will be found;   
            = 'V': all eigenvalues in the half-open interval (VL,VU]   
                   will be found;   
            = 'I': the IL-th through IU-th eigenvalues will be found.   

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

            where EPS is the machine precision.  If ABSTOL is less than   
            or equal to zero, then  EPS*|T|  will be used in its place,   
            where |T| is the 1-norm of the tridiagonal matrix obtained   
            by reducing AP to tridiagonal form.   

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
            If INFO = 0, the selected eigenvalues in ascending order.   

    Z       (output) LONG DOUBLE PRECISION array, dimension (LDZ, MAX(1,M))   
            If JOBZ = 'V', then if INFO = 0, the first M columns of Z   
            contain the orthonormal eigenvectors of the matrix A   
            corresponding to the selected eigenvalues, with the i-th   
            column of Z holding the eigenvector associated with W(i).   
            If an eigenvector fails to converge, then that column of Z   
            contains the latest approximation to the eigenvector, and the 
  
            index of the eigenvector is returned in IFAIL.   
            If JOBZ = 'N', then Z is not referenced.   
            Note: the user must ensure that at least MAX(1,M) columns are 
  
            supplied in the array Z; if RANGE = 'V', the exact value of M 
  
            is not known in advance and an upper bound must be used.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            JOBZ = 'V', LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (8*N)   

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
    static int indd, inde;
    static LONG DOUBLE anrm;
    static int imax;
    static LONG DOUBLE rmin, rmax;
    static int itmp1, i, j, indee;

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
    static LONG DOUBLE abstll, bignum;

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
    static int indtau, indisp;

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

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dopmtr_(char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qopmtr(char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qopmtr_(char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *);
    static int nsplit;
    static LONG DOUBLE smlnum, eps, vll, vuu, tmp1;



#define AP(I) ap[(I)-1]
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
    } else if (! (lsame_(uplo, "L") || lsame_(uplo, "U"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
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
	xerbla_("DSPEVX", &i__1);
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
	    W(1) = AP(1);
	} else {
	    if (*vl < AP(1) && *vu >= AP(1)) {
		*m = 1;
		W(1) = AP(1);
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
    abstll = *abstol;
    if (valeig) {
	vll = *vl;
	vuu = *vu;
    }

#ifdef PETSC_PREFIX_SUFFIX
    anrm = dlansp_("M", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anrm = qlansp("M", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anrm = qlansp_("M", uplo, n, &AP(1), &WORK(1));
#endif

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

	if (*abstol > 0.) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form. 
*/

    indtau = 1;
    inde = indtau + *n;
    indd = inde + *n;
    indwrk = indd + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dsptrd_(uplo, n, &AP(1), &WORK(indd), &WORK(inde), &WORK(indtau), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsptrd(uplo, n, &AP(1), &WORK(indd), &WORK(inde), &WORK(indtau), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsptrd_(uplo, n, &AP(1), &WORK(indd), &WORK(inde), &WORK(indtau), &iinfo);
#endif


/*     If all eigenvalues are desired and ABSTOL is less than or equal   
       to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails   
       for some eigenvalue, then try DSTEBZ. */

    if ((alleig || (indeig && *il == 1 && *iu == *n)) && *abstol <= 0.) {

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &WORK(indd), &c__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &WORK(indd), &c__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &WORK(indd), &c__1, &W(1), &c__1);
#endif

	indee = indwrk + (*n << 1);
	if (! wantz) {
	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dsterf_(n, &W(1), &WORK(indee), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsterf(n, &W(1), &WORK(indee), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsterf_(n, &W(1), &WORK(indee), info);
#endif

	} else {

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
	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &WORK(inde), &c__1, &WORK(indee), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_(jobz, n, &W(1), &WORK(indee), &Z(1,1), ldz, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr(jobz, n, &W(1), &WORK(indee), &Z(1,1), ldz, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_(jobz, n, &W(1), &WORK(indee), &Z(1,1), ldz, &WORK(
#endif

		    indwrk), info);
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
    indibl = 1;
    indisp = indibl + *n;
    indiwo = indisp + *n;

#ifdef PETSC_PREFIX_SUFFIX
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &WORK(indd), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qstebz(range, order, n, &vll, &vuu, il, iu, &abstll, &WORK(indd), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &WORK(indd), &WORK(
#endif

	    inde), m, &nsplit, &W(1), &IWORK(indibl), &IWORK(indisp), &WORK(
	    indwrk), &IWORK(indiwo), info);

    if (wantz) {

#ifdef PETSC_PREFIX_SUFFIX
	dstein_(n, &WORK(indd), &WORK(inde), m, &W(1), &IWORK(indibl), &IWORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qstein(n, &WORK(indd), &WORK(inde), m, &W(1), &IWORK(indibl), &IWORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qstein_(n, &WORK(indd), &WORK(inde), m, &W(1), &IWORK(indibl), &IWORK(
#endif

		indisp), &Z(1,1), ldz, &WORK(indwrk), &IWORK(indiwo), &
		IFAIL(1), info);

/*        Apply orthogonal matrix used in reduction to tridiagonal   
          form to eigenvectors returned by DSTEIN. */


#ifdef PETSC_PREFIX_SUFFIX
	dopmtr_("L", uplo, "N", n, m, &AP(1), &WORK(indtau), &Z(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qopmtr("L", uplo, "N", n, m, &AP(1), &WORK(indtau), &Z(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qopmtr_("L", uplo, "N", n, m, &AP(1), &WORK(indtau), &Z(1,1), 
#endif

		ldz, &WORK(indwrk), info);
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

/*     End of DSPEVX */

} /* dspevx_ */

