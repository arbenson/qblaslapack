#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed3_(int *k, int *kstart, int *kstop, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed3(int *k, int *kstart, int *kstop, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed3_(int *k, int *kstart, int *kstop, 
#endif

	int *n, LONG DOUBLE *d, LONG DOUBLE *q, int *ldq, LONG DOUBLE *
	rho, int *cutpnt, LONG DOUBLE *dlamda, LONG DOUBLE *q2, int *
	ldq2, int *indxc, int *ctot, LONG DOUBLE *w, LONG DOUBLE *s, 
	int *lds, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED3 finds the roots of the secular equation, as defined by the   
    values in D, W, and RHO, between KSTART and KSTOP.  It makes the   
    appropriate calls to DLAED4 and then updates the eigenvectors by   
    multiplying the matrix of eigenvectors of the pair of eigensystems   
    being combined by the matrix of eigenvectors of the K-by-K system   
    which is solved here.   

    This code makes very mild assumptions about floating point   
    arithmetic. It will work on machines with a guard digit in   
    add/subtract, or on those binary machines without guard digits   
    which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.   
    It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.   

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
            N >= K (deflation may result in N>K).   

    D       (output) LONG DOUBLE PRECISION array, dimension (N)   
            D(I) contains the updated eigenvalues for   
            KSTART <= I <= KSTOP.   

    Q       (output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            Initially the first K columns are used as workspace.   
            On output the columns KSTART to KSTOP contain   
            the updated eigenvectors.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= MAX(1,N).   

    RHO     (input) LONG DOUBLE PRECISION   
            The value of the parameter in the rank one update equation.   
            RHO >= 0 required.   

    CUTPNT  (input) INTEGER   
            The location of the last eigenvalue in the leading submatrix. 
  
            MIN(1,N) <= CUTPNT <= N.   

    DLAMDA  (input/output) LONG DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the old roots   
            of the deflated updating problem.  These are the poles   
            of the secular equation. May be changed on output by   
            having lowest order bit set to zero on Cray X-MP, Cray Y-MP, 
  
            Cray-2, or Cray C-90, as described above.   

    Q2      (input) LONG DOUBLE PRECISION array, dimension (LDQ2, N)   
            The first K columns of this matrix contain the non-deflated   
            eigenvectors for the split problem.   

    LDQ2    (input) INTEGER   
            The leading dimension of the array Q2.  LDQ2 >= MAX(1,N).   

    INDXC   (input) INTEGER array, dimension (N)   
            The permutation used to arrange the columns of the deflated   
            Q matrix into three groups:  the first group contains   
            non-zero elements only at and above CUTPNT, the second   
            contains non-zero elements only below CUTPNT, and the third   
            is dense.  The rows of the eigenvectors found by DLAED4   
            must be likewise permuted before the matrix multiply can take 
  
            place.   

    CTOT    (input) INTEGER array, dimension (4)   
            A count of the total number of the various types of columns   
            in Q, as described in INDXC.  The fourth column type is any   
            column which has been deflated.   

    W       (input/output) LONG DOUBLE PRECISION array, dimension (K)   
            The first K elements of this array contain the components   
            of the deflation-adjusted updating vector. Destroyed on   
            output.   

    S       (workspace) LONG DOUBLE PRECISION array, dimension (LDS, K)   
            Will contain the eigenvectors of the repaired matrix which   
            will be multiplied by the previously accumulated eigenvectors 
  
            to update the system.   

    LDS     (input) INTEGER   
            The leading dimension of S.  LDS >= MAX(1,K).   

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
    static LONG DOUBLE c_b22 = 1.;
    static LONG DOUBLE c_b23 = 0.;
    
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
    extern /* Subroutine */ void dgemm_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *),

#ifdef PETSC_PREFIX_SUFFIX
	     dcopy_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     qcopy(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     qcopy_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int 
#endif

	    *);
    static int ktemp, parts;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaed4_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed4(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed4_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *)
	    ;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3(LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamc3_(LONG DOUBLE *, LONG DOUBLE *);
#endif

    static int jc;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);



#define D(I) d[(I)-1]
#define DLAMDA(I) dlamda[(I)-1]
#define INDXC(I) indxc[(I)-1]
#define CTOT(I) ctot[(I)-1]
#define W(I) w[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define Q2(I,J) q2[(I)-1 + ((J)-1)* ( *ldq2)]
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
    } else if (*ldq < MAX(1,*n)) {
	*info = -7;
    } else if (*ldq2 < MAX(1,*n)) {
	*info = -12;
    } else if (*lds < MAX(1,*k)) {
	*info = -17;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED3", &i__1);
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

    ktemp = *kstop - *kstart + 1;
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
	    goto L130;
	}
/* L20: */
    }

    if (*k == 1 || *k == 2) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    i__2 = *k;
	    for (j = 1; j <= *k; ++j) {
		jc = INDXC(j);
		S(j,i) = Q(jc,i);
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
	    jc = INDXC(i);
	    S(i,j) = Q(jc,j) / temp;
/* L100: */
	}
/* L110: */
    }

/*     Compute the updated eigenvectors. */

L120:

    parts = 0;
    if (CTOT(1) > 0) {
	++parts;

#ifdef PETSC_PREFIX_SUFFIX
	dgemm_("N", "N", cutpnt, &ktemp, &CTOT(1), &c_b22, &Q2(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemm("N", "N", cutpnt, &ktemp, &CTOT(1), &c_b22, &Q2(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemm_("N", "N", cutpnt, &ktemp, &CTOT(1), &c_b22, &Q2(1,1), 
#endif

		ldq2, &S(1,*kstart), lds, &c_b23, &Q(1,*kstart), ldq);
    }
    if (CTOT(2) > 0) {
	parts += 2;
	i__1 = *n - *cutpnt;

#ifdef PETSC_PREFIX_SUFFIX
	dgemm_("N", "N", &i__1, &ktemp, &CTOT(2), &c_b22, &Q2(*cutpnt+1,CTOT(1)+1), ldq2, &S(CTOT(1)+1,*kstart), lds, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemm("N", "N", &i__1, &ktemp, &CTOT(2), &c_b22, &Q2(*cutpnt+1,CTOT(1)+1), ldq2, &S(CTOT(1)+1,*kstart), lds, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemm_("N", "N", &i__1, &ktemp, &CTOT(2), &c_b22, &Q2(*cutpnt+1,CTOT(1)+1), ldq2, &S(CTOT(1)+1,*kstart), lds, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif

    }
    if (parts == 1) {
	i__1 = *n - *cutpnt;

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("A", &i__1, &ktemp, &c_b23, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("A", &i__1, &ktemp, &c_b23, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("A", &i__1, &ktemp, &c_b23, &c_b23, &Q(*cutpnt+1,*kstart), ldq);
#endif

    }
    if (parts == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("A", cutpnt, &ktemp, &c_b23, &c_b23, &Q(1,*kstart),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("A", cutpnt, &ktemp, &c_b23, &c_b23, &Q(1,*kstart),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("A", cutpnt, &ktemp, &c_b23, &c_b23, &Q(1,*kstart),
#endif

		 ldq);
    }
    if (CTOT(3) > 0) {
	if (parts > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemm_("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b22, &Q(1,*kstart), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemm("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b22, &Q(1,*kstart), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemm_("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b22, &Q(1,*kstart), 
#endif

		    ldq);
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dgemm_("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b23, &Q(1,*kstart), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemm("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b23, &Q(1,*kstart), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemm_("N", "N", n, &ktemp, &CTOT(3), &c_b22, &Q2(1,CTOT(1)+1+CTOT(2)), ldq2, &S(CTOT(1)+1+CTOT(2),*kstart), lds, &c_b23, &Q(1,*kstart), 
#endif

		    ldq);
	}
    }

L130:
    return;

/*     End of DLAED3 */

} /* dlaed3_ */

