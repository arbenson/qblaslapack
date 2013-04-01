#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed8_(int *icompq, int *k, int *n, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed8(int *icompq, int *k, int *n, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed8_(int *icompq, int *k, int *n, int 
#endif

	*qsiz, LONG DOUBLE *d, LONG DOUBLE *q, int *ldq, int *indxq, 
	LONG DOUBLE *rho, int *cutpnt, LONG DOUBLE *z, LONG DOUBLE *dlamda, 
	LONG DOUBLE *q2, int *ldq2, LONG DOUBLE *w, int *perm, int *
	givptr, int *givcol, LONG DOUBLE *givnum, int *indxp, int *
	indx, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED8 merges the two sets of eigenvalues together into a single   
    sorted set.  Then it tries to deflate the size of the problem.   
    There are two ways in which deflation can occur:  when two or more   
    eigenvalues are close together or if there is a tiny element in the   
    Z vector.  For each such occurrence the order of the related secular 
  
    equation problem is reduced by one.   

    Arguments   
    =========   

    ICOMPQ  (input) INTEGER   
            = 0:  Compute eigenvalues only.   
            = 1:  Compute eigenvectors of original dense symmetric matrix 
  
                  also.  On entry, Q contains the orthogonal matrix used 
  
                  to reduce the original matrix to tridiagonal form.   

    K      (output) INTEGER   
           The number of non-deflated eigenvalues, and the order of the   
           related secular equation.   

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    QSIZ   (input) INTEGER   
           The dimension of the orthogonal matrix used to reduce   
           the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. 
  

    D      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
           On entry, the eigenvalues of the two submatrices to be   
           combined.  On exit, the trailing (N-K) updated eigenvalues   
           (those which were deflated) sorted into increasing order.   

    Q      (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
           If ICOMPQ = 0, Q is not referenced.  Otherwise,   
           on entry, Q contains the eigenvectors of the partially solved 
  
           system which has been previously updated in matrix   
           multiplies with other partially solved eigensystems.   
           On exit, Q contains the trailing (N-K) updated eigenvectors   
           (those which were deflated) in its last N-K columns.   

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  LDQ >= MAX(1,N).   

    INDXQ  (input) INTEGER array, dimension (N)   
           The permutation which separately sorts the two sub-problems   
           in D into ascending order.  Note that elements in the second   
           half of this permutation must first have CUTPNT added to   
           their values in order to be accurate.   

    RHO    (input/output) LONG DOUBLE PRECISION   
           On entry, the off-diagonal element associated with the rank-1 
  
           cut which originally split the two submatrices which are now   
           being recombined.   
           On exit, RHO has been modified to the value required by   
           DLAED3.   

    CUTPNT (input) INTEGER   
           The location of the last eigenvalue in the leading   
           sub-matrix.  MIN(1,N) <= CUTPNT <= N.   

    Z      (input) LONG DOUBLE PRECISION array, dimension (N)   
           On entry, Z contains the updating vector (the last row of   
           the first sub-eigenvector matrix and the first row of the   
           second sub-eigenvector matrix).   
           On exit, the contents of Z are destroyed by the updating   
           process.   

    DLAMDA (output) LONG DOUBLE PRECISION array, dimension (N)   
           A copy of the first K eigenvalues which will be used by   
           DLAED3 to form the secular equation.   

    Q2     (output) LONG DOUBLE PRECISION array, dimension (LDQ2,N)   
           If ICOMPQ = 0, Q2 is not referenced.  Otherwise,   
           a copy of the first K eigenvectors which will be used by   
           DLAED7 in a matrix multiply (DGEMM) to update the new   
           eigenvectors.   

    LDQ2   (input) INTEGER   
           The leading dimension of the array Q2.  LDQ2 >= MAX(1,N).   

    W      (output) LONG DOUBLE PRECISION array, dimension (N)   
           The first k values of the final deflation-altered z-vector and 
  
           will be passed to DLAED3.   

    PERM   (output) INTEGER array, dimension (N)   
           The permutations (from deflation and sorting) to be applied   
           to each eigenblock.   

    GIVPTR (output) INTEGER   
           The number of Givens rotations which took place in this   
           subproblem.   

    GIVCOL (output) INTEGER array, dimension (2, N)   
           Each pair of numbers indicates a pair of columns to take place 
  
           in a Givens rotation.   

    GIVNUM (output) LONG DOUBLE PRECISION array, dimension (2, N)   
           Each number indicates the S value to be used in the   
           corresponding Givens rotation.   

    INDXP  (workspace) INTEGER array, dimension (N)   
           The permutation used to place deflated values of D at the end 
  
           of the array.  INDXP(1:K) points to the nondeflated D-values   
           and INDXP(K+1:N) points to the deflated eigenvalues.   

    INDX   (workspace) INTEGER array, dimension (N)   
           The permutation used to sort the contents of D into ascending 
  
           order.   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  



       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static LONG DOUBLE c_b3 = -1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static int jlam, imax, jmax;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE c;
    static int i, j;
    static LONG DOUBLE s, t;

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
	    int *), dcopy_(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qcopy(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qcopy_(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif

	    *, int *);
    static int k2, n1, n2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlamch_(char *);
#endif

    static int jp;

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
    extern /* Subroutine */ void dlamrg_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamrg(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlamrg_(int *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *, int *), dlacpy_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *, int *), qlacpy(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *, int *), qlacpy_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    static int n1p1;
    static LONG DOUBLE eps, tau, tol;


    --d;
    q_dim1 = *ldq;
    q_offset = q_dim1 + 1;
    q -= q_offset;
    --indxq;
    --z;
    --dlamda;
    q2_dim1 = *ldq2;
    q2_offset = q2_dim1 + 1;
    q2 -= q2_offset;
    --w;
    --perm;
    givcol -= 3;
    givnum -= 3;
    --indxp;
    --indx;

    /* Function Body */
    *info = 0;

    if (*icompq < 0 || *icompq > 1) {
	*info = -1;
    } else if (*n < 0) {
	*info = -3;
    } else if (*icompq == 1 && *qsiz < *n) {
	*info = -4;
    } else if (*ldq < MAX(1,*n)) {
	*info = -7;
    } else if (*cutpnt < MIN(1,*n) || *cutpnt > *n) {
	*info = -10;
    } else if (*ldq2 < MAX(1,*n)) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED8", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    n1 = *cutpnt;
    n2 = *n - n1;
    n1p1 = n1 + 1;

    if (*rho < 0.) {

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&n2, &c_b3, &z[n1p1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&n2, &c_b3, &z[n1p1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&n2, &c_b3, &z[n1p1], &c__1);
#endif

    }

/*     Normalize z so that norm(z) = 1 */

    t = 1. / sqrt(2.);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	indx[j] = j;
/* L10: */
    }

#ifdef PETSC_PREFIX_SUFFIX
    dscal_(n, &t, &z[1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qscal(n, &t, &z[1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qscal_(n, &t, &z[1], &c__1);
#endif

    *rho = (d__1 = *rho * 2., ABS(d__1));

/*     Sort the eigenvalues into increasing order */

    i__1 = *n;
    for (i = *cutpnt + 1; i <= i__1; ++i) {
	indxq[i] += *cutpnt;
/* L20: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dlamda[i] = d[indxq[i]];
	w[i] = z[indxq[i]];
/* L30: */
    }
    i = 1;
    j = *cutpnt + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dlamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlamrg(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
#endif

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	d[i] = dlamda[indx[i]];
	z[i] = w[indx[i]];
/* L40: */
    }

/*     Calculate the allowable deflation tolerence */


#ifdef PETSC_PREFIX_SUFFIX
    imax = idamax_(n, &z[1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    imax = iqamax(n, &z[1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    imax = iqamax_(n, &z[1], &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    jmax = idamax_(n, &d[1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    jmax = iqamax(n, &d[1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    jmax = iqamax_(n, &d[1], &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif

    tol = eps * 8. * (d__1 = d[jmax], ABS(d__1));

/*     If the rank-1 modifier is small enough, no more needs to be done   
       except to reorganize Q so that its columns correspond with the   
       elements in D. */

    if (*rho * (d__1 = z[imax], ABS(d__1)) <= tol) {
	*k = 0;
	if (*icompq == 0) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		perm[j] = indxq[indx[j]];
/* L50: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		perm[j] = indxq[indx[j]];

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
#endif

			+ 1], &c__1);
/* L60: */
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq);
#endif

	}
	return;
    }

/*     If there are multiple eigenvalues then the problem deflates.  Here 
  
       the number of equal eigenvalues are found.  As each equal   
       eigenvalue is found, an elementary reflector is computed to rotate 
  
       the corresponding eigensubspace so that the corresponding   
       components of Z are zero in this new basis. */

    *k = 0;
    *givptr = 0;
    k2 = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (*rho * (d__1 = z[j], ABS(d__1)) <= tol) {

/*           Deflate due to small z component. */

	    --k2;
	    indxp[k2] = j;
	    if (j == *n) {
		goto L110;
	    }
	} else {
	    jlam = j;
	    goto L80;
	}
/* L70: */
    }
L80:
    ++j;
    if (j > *n) {
	goto L100;
    }
    if (*rho * (d__1 = z[j], ABS(d__1)) <= tol) {

/*        Deflate due to small z component. */

	--k2;
	indxp[k2] = j;
    } else {

/*        Check if eigenvalues are close enough to allow deflation. */

	s = z[jlam];
	c = z[j];

/*        Find sqrt(a**2+b**2) without overflow or   
          destructive underflow. */


#ifdef PETSC_PREFIX_SUFFIX
	tau = dlapy2_(&c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	tau = qlapy2(&c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	tau = qlapy2_(&c, &s);
#endif

	t = d[j] - d[jlam];
	c /= tau;
	s = -s / tau;
	if ((d__1 = t * c * s, ABS(d__1)) <= tol) {

/*           Deflation is possible. */

	    z[j] = tau;
	    z[jlam] = 0.;

/*           Record the appropriate Givens rotation */

	    ++(*givptr);
	    givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
	    givcol[(*givptr << 1) + 2] = indxq[indx[j]];
	    givnum[(*givptr << 1) + 1] = c;
	    givnum[(*givptr << 1) + 2] = s;
	    if (*icompq == 1) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[
#endif

			indxq[indx[j]] * q_dim1 + 1], &c__1, &c, &s);
	    }
	    t = d[jlam] * c * c + d[j] * s * s;
	    d[j] = d[jlam] * s * s + d[j] * c * c;
	    d[jlam] = t;
	    --k2;
	    i = 1;
L90:
	    if (k2 + i <= *n) {
		if (d[jlam] < d[indxp[k2 + i]]) {
		    indxp[k2 + i - 1] = indxp[k2 + i];
		    indxp[k2 + i] = jlam;
		    ++i;
		    goto L90;
		} else {
		    indxp[k2 + i - 1] = jlam;
		}
	    } else {
		indxp[k2 + i - 1] = jlam;
	    }
	    jlam = j;
	} else {
	    ++(*k);
	    w[*k] = z[jlam];
	    dlamda[*k] = d[jlam];
	    indxp[*k] = jlam;
	    jlam = j;
	}
    }
    goto L80;
L100:

/*     Record the last eigenvalue. */

    ++(*k);
    w[*k] = z[jlam];
    dlamda[*k] = d[jlam];
    indxp[*k] = jlam;

L110:

/*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA   
       and Q2 respectively.  The eigenvalues/vectors which were not   
       deflated go into the first K slots of DLAMDA and Q2 respectively, 
  
       while those which were deflated go into the last N - K slots. */

    if (*icompq == 0) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jp = indxp[j];
	    dlamda[j] = d[jp];
	    perm[j] = indxq[indx[jp]];
/* L120: */
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jp = indxp[j];
	    dlamda[j] = d[jp];
	    perm[j] = indxq[indx[jp]];

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1]
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1]
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1]
#endif

		    , &c__1);
/* L130: */
	}
    }

/*     The deflated eigenvalues and their corresponding vectors go back   
       into the last N - K slots of D and Q respectively. */

    if (*k < *n) {
	if (*icompq == 0) {
	    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif

	} else {
	    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &dlamda[*k + 1], &c__1, &d[*k + 1], &c__1);
#endif

	    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*
#endif

		    k + 1) * q_dim1 + 1], ldq);
	}
    }

    return;

/*     End of DLAED8 */

} /* dlaed8_ */

