#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgecon_(char *norm, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgecon(char *norm, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgecon_(char *norm, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *anorm, LONG DOUBLE *rcond, LONG DOUBLE *work, int *
	iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGECON estimates the reciprocal of the condition number of a general 
  
    real matrix A, in either the 1-norm or the infinity-norm, using   
    the LU factorization computed by DGETRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as   
       RCOND = 1 / ( norm(A) * norm(inv(A)) ).   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    ANORM   (input) LONG DOUBLE PRECISION   
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (4*N)   

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

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE sl;
    static int ix;

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

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE su;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ainvnm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatrs_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrs(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrs_(char *, char *, char *, char *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *);
    static long int onenrm;
    static char normin[1];
    static LONG DOUBLE smlnum;



#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    } else if (*anorm < 0.) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGECON", &i__1);
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


/*     Estimate the norm of inv(A). */

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

/*           Multiply by inv(L). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatrs_("Lower", "No transpose", "Unit", normin, n, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrs("Lower", "No transpose", "Unit", normin, n, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrs_("Lower", "No transpose", "Unit", normin, n, &A(1,1), 
#endif

		    lda, &WORK(1), &sl, &WORK((*n << 1) + 1), info);

/*           Multiply by inv(U). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatrs_("Upper", "No transpose", "Non-unit", normin, n, &A(1,1), lda, &WORK(1), &su, &WORK(*n * 3 + 1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrs("Upper", "No transpose", "Non-unit", normin, n, &A(1,1), lda, &WORK(1), &su, &WORK(*n * 3 + 1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrs_("Upper", "No transpose", "Non-unit", normin, n, &A(1,1), lda, &WORK(1), &su, &WORK(*n * 3 + 1), info);
#endif

	} else {

/*           Multiply by inv(U'). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatrs_("Upper", "Transpose", "Non-unit", normin, n, &A(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrs("Upper", "Transpose", "Non-unit", normin, n, &A(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrs_("Upper", "Transpose", "Non-unit", normin, n, &A(1,1),
#endif

		     lda, &WORK(1), &su, &WORK(*n * 3 + 1), info);

/*           Multiply by inv(L'). */


#ifdef PETSC_PREFIX_SUFFIX
	    dlatrs_("Lower", "Transpose", "Unit", normin, n, &A(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrs("Lower", "Transpose", "Unit", normin, n, &A(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrs_("Lower", "Transpose", "Unit", normin, n, &A(1,1), 
#endif

		    lda, &WORK(1), &sl, &WORK((*n << 1) + 1), info);
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. 
*/

	scale = sl * su;
	*(unsigned char *)normin = 'Y';
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

/*     End of DGECON */

} /* dgecon_ */

