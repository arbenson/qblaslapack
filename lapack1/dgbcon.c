#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbcon_(char *norm, int *n, int *kl, int *ku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbcon(char *norm, int *n, int *kl, int *ku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbcon_(char *norm, int *n, int *kl, int *ku,
#endif

	 LONG DOUBLE *ab, int *ldab, int *ipiv, LONG DOUBLE *anorm, 
	LONG DOUBLE *rcond, LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGBCON estimates the reciprocal of the condition number of a real   
    general band matrix A, in either the 1-norm or the infinity-norm,   
    using the LU factorization computed by DGBTRF.   

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

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            Details of the LU factorization of the band matrix A, as   
            computed by DGBTRF.  U is stored as an upper triangular band 
  
            matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and   
            the multipliers used during the factorization are stored in   
            rows KL+KU+2 to 2*KL+KU+1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= N, row i of the matrix was   
            interchanged with row IPIV(i).   

    ANORM   (input) LONG DOUBLE PRECISION   
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Local variables */
    static int kase;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static int kase1, j;
    static LONG DOUBLE t, scale;
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
    static long int lnoti;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static int kd;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int lm, jp, ix;

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
	    LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    static LONG DOUBLE ainvnm;
    static long int onenrm;
    static char normin[1];
    static LONG DOUBLE smlnum;



#define IPIV(I) ipiv[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -6;
    } else if (*anorm < 0.) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBCON", &i__1);
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
    kd = *kl + *ku + 1;
    lnoti = *kl > 0;
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

	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = MIN(i__2,i__3);
		    jp = IPIV(j);
		    t = WORK(jp);
		    if (jp != j) {
			WORK(jp) = WORK(j);
			WORK(j) = t;
		    }
		    d__1 = -t;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&lm, &d__1, &AB(kd+1,j), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&lm, &d__1, &AB(kd+1,j), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&lm, &d__1, &AB(kd+1,j), &c__1, &
#endif

			    WORK(j + 1), &c__1);
/* L20: */
		}
	    }

/*           Multiply by inv(U). */

	    i__1 = *kl + *ku;

#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
#endif

		    AB(1,1), ldab, &WORK(1), &scale, &WORK((*n << 1) + 
		    1), info);
	} else {

/*           Multiply by inv(U'). */

	    i__1 = *kl + *ku;

#ifdef PETSC_PREFIX_SUFFIX
	    dlatbs_("Upper", "Transpose", "Non-unit", normin, n, &i__1, &AB(1,1), ldab, &WORK(1), &scale, &WORK((*n << 1) + 1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatbs("Upper", "Transpose", "Non-unit", normin, n, &i__1, &AB(1,1), ldab, &WORK(1), &scale, &WORK((*n << 1) + 1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatbs_("Upper", "Transpose", "Non-unit", normin, n, &i__1, &AB(1,1), ldab, &WORK(1), &scale, &WORK((*n << 1) + 1), 
#endif

		    info);

/*           Multiply by inv(L'). */

	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = MIN(i__1,i__2);

#ifdef PETSC_PREFIX_SUFFIX
		    WORK(j) -= ddot_(&lm, &AB(kd+1,j), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    WORK(j) -= qdot(&lm, &AB(kd+1,j), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    WORK(j) -= qdot_(&lm, &AB(kd+1,j), &c__1, &
#endif

			    WORK(j + 1), &c__1);
		    jp = IPIV(j);
		    if (jp != j) {
			t = WORK(jp);
			WORK(jp) = WORK(j);
			WORK(j) = t;
		    }
/* L30: */
		}
	    }
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

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
		goto L40;
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

L40:
    return;

/*     End of DGBCON */

} /* dgbcon_ */

