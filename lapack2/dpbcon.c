#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbcon_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbcon(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbcon_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif

	ab, int *ldab, LONG DOUBLE *anorm, LONG DOUBLE *rcond, LONG DOUBLE *
	work, int *iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPBCON estimates the reciprocal of the condition number (in the   
    1-norm) of a real symmetric positive definite band matrix using the   
    Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangular factor stored in AB;   
            = 'L':  Lower triangular factor stored in AB.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T of the band matrix A, stored in the 
  
            first KD+1 rows of the array.  The j-th column of U or L is   
            stored in the j-th column of the array AB as follows:   
            if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=MIN(n,j+kd). 
  

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    ANORM   (input) LONG DOUBLE PRECISION   
            The 1-norm (or infinity-norm) of the symmetric band matrix A. 
  

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an   
            estimate of the 1-norm of inv(A) computed in this routine.   

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
    int  i__1;
    LONG DOUBLE d__1;
    /* Local variables */
    static int kase;
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
    static long int upper;

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
    static LONG DOUBLE scalel;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatbs_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatbs(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatbs_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE scaleu;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ainvnm;
    static char normin[1];
    static LONG DOUBLE smlnum;



#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*ldab < *kd + 1) {
	*info = -5;
    } else if (*anorm < 0.) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBCON", &i__1);
	return;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return;
    } else if (*anorm == 0.) {
	return;
    }


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("Safe minimum");
#endif


/*     Estimate the 1-norm of the inverse. */

    kase = 0;
    *(unsigned char *)normin = 'N';
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
	if (upper) {

/*           Multiply by inv(U'). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Upper", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Upper", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Upper", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif

		     info);
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(U). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Upper", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Upper", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif

		     info);
	} else {

/*           Multiply by inv(L). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Lower", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Lower", "No transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scalel, &WORK((*n << 1) + 1),
#endif

		     info);
	    *(unsigned char *)normin = 'Y';

/*           Multiply by inv(L'). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Lower", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Lower", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Lower", "Transpose", "Non-unit", normin, n, kd, &AB(1,1), ldab, &WORK(1), &scaleu, &WORK((*n << 1) + 1),
#endif

		     info);
	}

/*        Multiply by 1/SCALE if doing so will not cause overflow. */

	scale = scalel * scaleu;
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

	    if (scale < (d__1 = WORK(ix), ABS(d__1)) * smlnum || scale == 0.) 
		    {
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

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

L20:

    return;

/*     End of DPBCON */

} /* dpbcon_ */

