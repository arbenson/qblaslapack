#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtpcon_(char *norm, char *uplo, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtpcon(char *norm, char *uplo, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtpcon_(char *norm, char *uplo, char *diag, int *n, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *rcond, LONG DOUBLE *work, int *iwork, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTPCON estimates the reciprocal of the condition number of a packed   
    triangular matrix A, in either the 1-norm or the infinity-norm.   

    The norm of A is computed and an estimate is obtained for   
    norm(inv(A)), then the reciprocal of the condition number is   
    computed as   
       RCOND = 1 / ( norm(A) * norm(inv(A)) ).   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            If DIAG = 'U', the diagonal elements of A are not referenced 
  
            and are assumed to be 1.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Local variables */
    static int kase, kase1;
    static LONG DOUBLE scale;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drscl_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrscl(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrscl_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE anorm;
    static long int upper;
    static LONG DOUBLE xnorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, int *);
    static int ix;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlantp_(char *, char *, char *, int *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlantp(char *, char *, char *, int *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlantp_(char *, char *, char *, int *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *);
    static LONG DOUBLE ainvnm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatps_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatps(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatps_(char *, char *, char *, char *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *);
    static long int onenrm;
    static char normin[1];
    static LONG DOUBLE smlnum;
    static long int nounit;



#define IWORK(I) iwork[(I)-1]
#define WORK(I) work[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    nounit = lsame_(diag, "N");

    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPCON", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	*rcond = 1.;
	return;
    }

    *rcond = 0.;

#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("Safe minimum") * (LONG DOUBLE) MAX(1,*n);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("Safe minimum") * (LONG DOUBLE) MAX(1,*n);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("Safe minimum") * (LONG DOUBLE) MAX(1,*n);
#endif


/*     Compute the norm of the triangular matrix A. */


#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlantp_(norm, uplo, diag, n, &AP(1), &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlantp(norm, uplo, diag, n, &AP(1), &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlantp_(norm, uplo, diag, n, &AP(1), &WORK(1));
#endif


/*     Continue only if ANORM > 0. */

    if (anorm > 0.) {

/*        Estimate the norm of the inverse of A. */

	ainvnm = 0.;
	*(unsigned char *)normin = 'N';
	if (onenrm) {
	    kase1 = 1;
	} else {
	    kase1 = 2;
	}
	kase = 0;
L10:

#ifdef PETSC_PREFIX_SUFFIX
	dlacon_(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacon(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacon_(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif

	if (kase != 0) {
	    if (kase == kase1) {

/*              Multiply by inv(A). */


#ifdef PETSC_PREFIX_SUFFIX
		dlatps_(uplo, "No transpose", diag, normin, n, &AP(1), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlatps(uplo, "No transpose", diag, normin, n, &AP(1), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlatps_(uplo, "No transpose", diag, normin, n, &AP(1), &WORK(
#endif

			1), &scale, &WORK((*n << 1) + 1), info);
	    } else {

/*              Multiply by inv(A'). */


#ifdef PETSC_PREFIX_SUFFIX
		dlatps_(uplo, "Transpose", diag, normin, n, &AP(1), &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlatps(uplo, "Transpose", diag, normin, n, &AP(1), &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlatps_(uplo, "Transpose", diag, normin, n, &AP(1), &WORK(1), 
#endif

			&scale, &WORK((*n << 1) + 1), info);
	    }
	    *(unsigned char *)normin = 'Y';

/*           Multiply by 1/SCALE if doing so will not cause overfl
ow. */

	    if (scale != 1.) {

#ifdef PETSC_PREFIX_SUFFIX
		ix = idamax_(n, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		ix = iqamax(n, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		ix = iqamax_(n, &WORK(1), &c__1);
#endif

		xnorm = (d__1 = WORK(ix), ABS(d__1));
		if (scale < xnorm * smlnum || scale == 0.) {
		    goto L20;
		}

#ifdef PETSC_PREFIX_SUFFIX
		drscl_(n, &scale, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrscl(n, &scale, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrscl_(n, &scale, &WORK(1), &c__1);
#endif

	    }
	    goto L10;
	}

/*        Compute the estimate of the reciprocal condition number. */

	if (ainvnm != 0.) {
	    *rcond = 1. / anorm / ainvnm;
	}
    }

L20:
    return;

/*     End of DTPCON */

} /* dtpcon_ */

