#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsbevd_(char *jobz, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsbevd(char *jobz, char *uplo, int *n, int *kd, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsbevd_(char *jobz, char *uplo, int *n, int *kd, 
#endif

	LONG DOUBLE *ab, int *ldab, LONG DOUBLE *w, LONG DOUBLE *z, int *
	ldz, LONG DOUBLE *work, int *lwork, int *iwork, int *
	liwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSBEVD computes all the eigenvalues and, optionally, eigenvectors of 
  
    a real symmetric band matrix A. If eigenvectors are desired, it uses 
  
    a divide and conquer algorithm.   

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

    WORK    (workspace/output) LONG DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            IF N <= 1,                LWORK must be at least 1.   
            If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.   
            If JOBZ  = 'V' and N > 2, LWORK must be at least   
                           ( 1 + 4*N + 2*N*lg N + 3*N**2 ),   
                           where lg( N ) = smallest int k such   
                                           that 2**k >= N.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array LIWORK.   
            If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 2, LIWORK must be at least 2 + 5*N.   

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
    static LONG DOUBLE c_b14 = 1.;
    static LONG DOUBLE c_b21 = 0.;
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


#ifdef PETSC_PREFIX_SUFFIX
	    int *), dgemm_(char *, char *, int *, int *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qgemm(char *, char *, int *, int *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qgemm_(char *, char *, int *, int *, int *
#endif

	    , LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE sigma;
    extern long int lsame_(char *, char *);
    static int iinfo, lwmin;
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

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *, int *), dlacpy_(char *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *, int *), qlacpy(char *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *, int *), qlacpy_(char *, int 
#endif

	    *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
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
    static int indwrk, liwmin;
    static LONG DOUBLE smlnum;
    static int lgn;
    static LONG DOUBLE eps;



#define W(I) w[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");

    *info = 0;
    if (*n <= 1) {
	lgn = 0;
	liwmin = 1;
	lwmin = 1;
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
	    lwmin = (*n << 2) + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	} else {
	    liwmin = 1;
	    lwmin = *n << 1;
	}
    }
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
    } else if (*lwork < lwmin) {
	*info = -11;
    } else if (*liwork < liwmin) {
	*info = -13;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSBEVD ", &i__1);
	goto L10;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L10;
    }

    if (*n == 1) {
	W(1) = AB(1,1);
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
	    dlascl_("B", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("B", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("B", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif

		    info);
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("Q", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("Q", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("Q", kd, kd, &c_b14, &sigma, n, n, &AB(1,1), ldab, 
#endif

		    info);
	}
    }

/*     Call DSBTRD to reduce symmetric band matrix to tridiagonal form. */

    inde = 1;
    indwrk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwrk2 = *lwork - indwk2 + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dsbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsbtrd(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsbtrd_(jobz, uplo, n, kd, &AB(1,1), ldab, &W(1), &WORK(inde), &Z(1,1), ldz, &WORK(indwrk), &iinfo);
#endif


/*     For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEDC. 
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
	dgemm_("N", "N", n, n, n, &c_b14, &Z(1,1), ldz, &WORK(indwrk), n,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemm("N", "N", n, n, n, &c_b14, &Z(1,1), ldz, &WORK(indwrk), n,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemm_("N", "N", n, n, n, &c_b14, &Z(1,1), ldz, &WORK(indwrk), n,
#endif

		 &c_b21, &WORK(indwk2), n);

#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("A", n, n, &WORK(indwk2), n, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("A", n, n, &WORK(indwk2), n, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("A", n, n, &WORK(indwk2), n, &Z(1,1), ldz);
#endif

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

/*     End of DSBEVD */

} /* dsbevd_ */

