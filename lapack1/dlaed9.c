#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed9_(int *k, int *kstart, int *kstop, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed9(int *k, int *kstart, int *kstop, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed9_(int *k, int *kstart, int *kstop, 
#endif

	int *n, LONG DOUBLE *d, LONG DOUBLE *q, int *ldq, LONG DOUBLE *
	rho, LONG DOUBLE *dlamda, LONG DOUBLE *w, LONG DOUBLE *s, int *lds, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED9 finds the roots of the secular equation, as defined by the   
    values in D, Z, and RHO, between KSTART and KSTOP.  It makes the   
    appropriate calls to DLAED4 and then stores the new matrix of   
    eigenvectors for use in calculating the next level of Z vectors.   

    Arguments   
    =========   

    K       (input) INTEGER   
            The number of terms in the rational function to be solved by 
  
            DLAED4.  K >= 0.   

    KSTART  (input) INTEGER   
    KSTOP   (input) INTEGER   
            The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP   
            are to be computed.  1 <= KSTART <= KSTOP <= K.   

    N       (input) INTEGER   
            The number of rows and columns in the Q matrix.   
            N >= K (delation may result in N > K).   

    D       (output) LONG DOUBLE PRECISION array, dimension (N)   
            D(I) contains the updated eigenvalues   
            for KSTART <= I <= KSTOP.   

    Q       (workspace) LONG DOUBLE PRECISION array, dimension (LDQ,N)   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= MAX( 1, N ).   

    RHO     (input) LONG DOUBLE PRECISION   
            The value of the parameter in the rank one update equation.   
            RHO >= 0 required.   

    DLAMDA  (input) LONG DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the old roots   
            of the deflated updating problem.  These are the poles   
            of the secular equation.   

    W       (input) LONG DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the components   
            of the deflation-adjusted updating vector.   

    S       (output) LONG DOUBLE PRECISION array, dimension (LDS, K)   
            Will contain the eigenvectors of the repaired matrix which   
            will be stored for subsequent Z vector calculation and   
            multiplied by the previously accumulated eigenvectors   
            to update the system.   

    LDS     (input) INTEGER   
            The leading dimension of S.  LDS >= MAX( 1, K ).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, an eigenvalue did not converge   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

    static int i, j;

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
	    LONG DOUBLE *, int *), dlaed4_(int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaed4(int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaed4_(int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);



#define D(I) d[(I)-1]
#define DLAMDA(I) dlamda[(I)-1]
#define W(I) w[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define S(I,J) s[(I)-1 + ((J)-1)* ( *lds)]

    *info = 0;

    if (*k < 0) {
	*info = -1;
    } else if (*kstart < 1 || *kstart > MAX(1,*k)) {
	*info = -2;
    } else if (MAX(1,*kstop) < *kstart || *kstop > MAX(1,*k)) {
	*info = -3;
    } else if (*n < *k) {
	*info = -4;
    } else if (*ldq < MAX(1,*k)) {
	*info = -7;
    } else if (*lds < MAX(1,*k)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED9", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*k == 0) {
	return;
    }

/*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can   
       be computed with high relative accuracy (barring over/underflow). 
  
       This is a problem on machines without a guard digit in   
       add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).   
       The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),   
       which on any of these machines zeros out the bottommost   
       bit of DLAMDA(I) if it is 1; this makes the subsequent   
       subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation   
       occurs. On binary machines with a guard digit (almost all   
       machines) it does not change DLAMDA(I) at all. On hexadecimal   
       and decimal machines with a guard digit, it slightly   
       changes the bottommost bits of DLAMDA(I). It does not account   
       for hexadecimal or decimal machines without guard digits   
       (we know of none). We use a subroutine call to compute   
       2*DLAMBDA(I) to prevent optimizing compilers from eliminating   
       this code. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
	DLAMDA(i) = dlamc3_(&DLAMDA(i), &DLAMDA(i)) - DLAMDA(i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	DLAMDA(i) = qlamc3(&DLAMDA(i), &DLAMDA(i)) - DLAMDA(i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	DLAMDA(i) = qlamc3_(&DLAMDA(i), &DLAMDA(i)) - DLAMDA(i);
#endif

/* L10: */
    }

    i__1 = *kstop;
    for (j = *kstart; j <= *kstop; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaed4_(k, &j, &DLAMDA(1), &W(1), &Q(1,j), rho, &D(j), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaed4(k, &j, &DLAMDA(1), &W(1), &Q(1,j), rho, &D(j), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaed4_(k, &j, &DLAMDA(1), &W(1), &Q(1,j), rho, &D(j), 
#endif

		info);

/*        If the zero finder fails, the computation is terminated. */

	if (*info != 0) {
	    goto L120;
	}
/* L20: */
    }

    if (*k == 1 || *k == 2) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    i__2 = *k;
	    for (j = 1; j <= *k; ++j) {
		S(j,i) = Q(j,i);
/* L30: */
	    }
/* L40: */
	}
	goto L120;
    }

/*     Compute updated W. */


#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(k, &W(1), &c__1, &S(1,1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(k, &W(1), &c__1, &S(1,1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(k, &W(1), &c__1, &S(1,1), &c__1);
#endif


/*     Initialize W(I) = Q(I,I) */

    i__1 = *ldq + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dcopy_(k, &Q(1,1), &i__1, &W(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qcopy(k, &Q(1,1), &i__1, &W(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qcopy_(k, &Q(1,1), &i__1, &W(1), &c__1);
#endif

    i__1 = *k;
    for (j = 1; j <= *k; ++j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    W(i) *= Q(i,j) / (DLAMDA(i) - DLAMDA(j));
/* L50: */
	}
	i__2 = *k;
	for (i = j + 1; i <= *k; ++i) {
	    W(i) *= Q(i,j) / (DLAMDA(i) - DLAMDA(j));
/* L60: */
	}
/* L70: */
    }
    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	d__1 = sqrt(-W(i));
	W(i) = SIGN(d__1, S(i,1));
/* L80: */
    }

/*     Compute eigenvectors of the modified rank-1 modification. */

    i__1 = *k;
    for (j = 1; j <= *k; ++j) {
	i__2 = *k;
	for (i = 1; i <= *k; ++i) {
	    Q(i,j) = W(i) / Q(i,j);
/* L90: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	temp = dnrm2_(k, &Q(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	temp = qnrm2(k, &Q(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	temp = qnrm2_(k, &Q(1,j), &c__1);
#endif

	i__2 = *k;
	for (i = 1; i <= *k; ++i) {
	    S(i,j) = Q(i,j) / temp;
/* L100: */
	}
/* L110: */
    }

L120:
    return;

/*     End of DLAED9 */

} /* dlaed9_ */

