#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed0_(int *icompq, int *qsiz, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed0(int *icompq, int *qsiz, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed0_(int *icompq, int *qsiz, int *n, 
#endif

	LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *q, int *ldq, LONG DOUBLE 
	*qstore, int *ldqs, LONG DOUBLE *work, int *iwork, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAED0 computes all eigenvalues and corresponding eigenvectors of a   
    symmetric tridiagonal matrix using the divide and conquer method.   

    Arguments   
    =========   

    ICOMPQ  (input) INTEGER   
            = 0:  Compute eigenvalues only.   
            = 1:  Compute eigenvectors of original dense symmetric matrix 
  
                  also.  On entry, Q contains the orthogonal matrix used 
  
                  to reduce the original matrix to tridiagonal form.   
            = 2:  Compute eigenvalues and eigenvectors of tridiagonal   
                  matrix.   

    QSIZ   (input) INTEGER   
           The dimension of the orthogonal matrix used to reduce   
           the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. 
  

    N      (input) INTEGER   
           The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D      (input/output) LONG DOUBLE PRECISION array, dimension (N)   
           On entry, the main diagonal of the tridiagonal matrix.   
           On exit, its eigenvalues.   

    E      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
           The off-diagonal elements of the tridiagonal matrix.   
           On exit, E has been destroyed.   

    Q      (input/output) LONG DOUBLE PRECISION array, dimension (LDQ, N)   
           On entry, Q must contain an N-by-N orthogonal matrix.   
           If ICOMPQ = 0    Q is not referenced.   
           If ICOMPQ = 1    On entry, Q is a subset of the columns of the 
  
                            orthogonal matrix used to reduce the full   
                            matrix to tridiagonal form corresponding to   
                            the subset of the full matrix which is being 
  
                            decomposed at this time.   
           If ICOMPQ = 2    On entry, Q will be the identity matrix.   
                            On exit, Q contains the eigenvectors of the   
                            tridiagonal matrix.   

    LDQ    (input) INTEGER   
           The leading dimension of the array Q.  If eigenvectors are   
           desired, then  LDQ >= MAX(1,N).  In any case,  LDQ >= 1.   

    QSTORE (workspace) LONG DOUBLE PRECISION array, dimension (LDQS, N)   
           Referenced only when ICOMPQ = 1.  Used to store parts of   
           the eigenvector matrix when the updating matrix multiplies   
           take place.   

    LDQS   (input) INTEGER   
           The leading dimension of the array QSTORE.  If ICOMPQ = 1,   
           then  LDQS >= MAX(1,N).  In any case,  LDQS >= 1.   

    WORK   (workspace) LONG DOUBLE PRECISION array,   
                                  dimension (1 + 3*N + 2*N*lg N + 2*N**2) 
  
                          ( lg( N ) = smallest int k   
                                      such that 2^k >= N )   

    IWORK  (workspace) INTEGER array,   
           If ICOMPQ = 0 or 1, the dimension of IWORK must be at least   
                          6 + 6*N + 5*N*lg N.   
                          ( lg( N ) = smallest int k   
                                      such that 2^k >= N )   
           If ICOMPQ = 2, the dimension of IWORK must be at least   
                          2 + 5*N.   

    INFO   (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  The algorithm failed to compute an eigenvalue while   
                  working on the submatrix lying in rows and columns   
                  INFO/(N+1) through mod(INFO,N+1).   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__2 = 2;
    static LONG DOUBLE c_b16 = 1.;
    static LONG DOUBLE c_b17 = 0.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp;
    static int curr, i, j, k;

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
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int iperm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static int indxq, iwrem;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaed1_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed1(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed1_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *);
    static int iqptr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaed7_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed7(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed7_(int *, int *, int *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *, int *, int *, int *, LONG DOUBLE 
	    *, LONG DOUBLE *, int *, int *);
    static int tlvls, iq;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static int igivcl;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static int igivnm, submat, curprb, subpbs, igivpt;

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
    static int curlvl, matsiz, iprmpt, lgn, msd2, smm1, spm1, spm2;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define QSTORE(I,J) qstore[(I)-1 + ((J)-1)* ( *ldqs)]

    *info = 0;

    if (*icompq < 0 || *icompq > 2) {
	*info = -1;
    } else if (*icompq == 1 && *qsiz < MAX(0,*n)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldq < MAX(1,*n)) {
	*info = -7;
    } else if (*ldqs < MAX(1,*n)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED0", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Determine the size and placement of the submatrices, and save in   
       the leading elements of IWORK. */

    IWORK(1) = *n;
    subpbs = 1;
    tlvls = 0;
L10:
    if (IWORK(subpbs) > 25) {
	for (j = subpbs; j >= 1; --j) {
	    IWORK(j * 2) = (IWORK(j) + 1) / 2;
	    IWORK((j << 1) - 1) = IWORK(j) / 2;
/* L20: */
	}
	++tlvls;
	subpbs <<= 1;
	goto L10;
    }
    i__1 = subpbs;
    for (j = 2; j <= subpbs; ++j) {
	IWORK(j) += IWORK(j - 1);
/* L30: */
    }

/*     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1 
  
       using rank-1 modifications (cuts). */

    spm1 = subpbs - 1;
    i__1 = spm1;
    for (i = 1; i <= spm1; ++i) {
	submat = IWORK(i) + 1;
	smm1 = submat - 1;
	D(smm1) -= (d__1 = E(smm1), ABS(d__1));
	D(submat) -= (d__1 = E(smm1), ABS(d__1));
/* L40: */
    }

    indxq = (*n << 2) + 3;
    if (*icompq != 2) {

/*        Set up workspaces for eigenvalues only/accumulate new vector
s   
          routine */

	temp = log((LONG DOUBLE) (*n)) / log(2.);
	lgn = (int) temp;
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	if (pow((LONG DOUBLE)c__2, (LONG DOUBLE)lgn) < *n) {
	    ++lgn;
	}
	iprmpt = indxq + *n + 1;
	iperm = iprmpt + *n * lgn;
	iqptr = iperm + *n * lgn;
	igivpt = iqptr + *n + 2;
	igivcl = igivpt + *n * lgn;

	igivnm = 1;
	iq = igivnm + (*n << 1) * lgn;
/* Computing 2nd power */
	i__1 = *n;
	iwrem = iq + i__1 * i__1 + 1;

/*        Initialize pointers */

	i__1 = subpbs;
	for (i = 0; i <= subpbs; ++i) {
	    IWORK(iprmpt + i) = 1;
	    IWORK(igivpt + i) = 1;
/* L50: */
	}
	IWORK(iqptr) = 1;
    }

/*     Solve each submatrix eigenproblem at the bottom of the divide and 
  
       conquer tree. */

    curr = 0;
    i__1 = spm1;
    for (i = 0; i <= spm1; ++i) {
	if (i == 0) {
	    submat = 1;
	    matsiz = IWORK(1);
	} else {
	    submat = IWORK(i) + 1;
	    matsiz = IWORK(i + 1) - IWORK(i);
	}
	if (*icompq == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_("I", &matsiz, &D(submat), &E(submat), &Q(submat,submat), ldq, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr("I", &matsiz, &D(submat), &E(submat), &Q(submat,submat), ldq, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_("I", &matsiz, &D(submat), &E(submat), &Q(submat,submat), ldq, &WORK(1), info);
#endif

	    if (*info != 0) {
		goto L130;
	    }
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_("I", &matsiz, &D(submat), &E(submat), &WORK(iq - 1 + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr("I", &matsiz, &D(submat), &E(submat), &WORK(iq - 1 + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_("I", &matsiz, &D(submat), &E(submat), &WORK(iq - 1 + 
#endif

		    IWORK(iqptr + curr)), &matsiz, &WORK(1), info);
	    if (*info != 0) {
		goto L130;
	    }
	    if (*icompq == 1) {

#ifdef PETSC_PREFIX_SUFFIX
		dgemm_("N", "N", qsiz, &matsiz, &matsiz, &c_b16, &Q(1,submat), ldq, &WORK(iq - 1 + IWORK(iqptr + curr)),
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemm("N", "N", qsiz, &matsiz, &matsiz, &c_b16, &Q(1,submat), ldq, &WORK(iq - 1 + IWORK(iqptr + curr)),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemm_("N", "N", qsiz, &matsiz, &matsiz, &c_b16, &Q(1,submat), ldq, &WORK(iq - 1 + IWORK(iqptr + curr)),
#endif

			 &matsiz, &c_b17, &QSTORE(1,submat), 
			ldqs);
	    }
/* Computing 2nd power */
	    i__2 = matsiz;
	    IWORK(iqptr + curr + 1) = IWORK(iqptr + curr) + i__2 * i__2;
	    ++curr;
	}
	k = 1;
	i__2 = IWORK(i + 1);
	for (j = submat; j <= IWORK(i+1); ++j) {
	    IWORK(indxq + j) = k;
	    ++k;
/* L60: */
	}
/* L70: */
    }

/*     Successively merge eigensystems of adjacent submatrices   
       into eigensystem for the corresponding larger matrix.   

       while ( SUBPBS > 1 ) */

    curlvl = 1;
L80:
    if (subpbs > 1) {
	spm2 = subpbs - 2;
	i__1 = spm2;
	for (i = 0; i <= spm2; i += 2) {
	    if (i == 0) {
		submat = 1;
		matsiz = IWORK(2);
		msd2 = IWORK(1);
		curprb = 0;
	    } else {
		submat = IWORK(i) + 1;
		matsiz = IWORK(i + 2) - IWORK(i);
		msd2 = matsiz / 2;
		++curprb;
	    }

/*     Merge lower order eigensystems (of size MSD2 and MATSIZ - M
SD2)   
       into an eigensystem of size MATSIZ.   
       DLAED1 is used only for the full eigensystem of a tridiagon
al   
       matrix.   
       DLAED7 handles the cases in which eigenvalues only or eigen
values   
       and eigenvectors of a full symmetric matrix (which was redu
ced to   
       tridiagonal form) are desired. */

	    if (*icompq == 2) {

#ifdef PETSC_PREFIX_SUFFIX
		dlaed1_(&matsiz, &D(submat), &Q(submat,submat), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaed1(&matsiz, &D(submat), &Q(submat,submat), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaed1_(&matsiz, &D(submat), &Q(submat,submat), 
#endif

			ldq, &IWORK(indxq + submat), &E(submat + msd2 - 1), &
			msd2, &WORK(1), &IWORK(subpbs + 1), info);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dlaed7_(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &D(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaed7(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &D(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaed7_(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &D(
#endif

			submat), &QSTORE(1,submat), ldqs, &
			IWORK(indxq + submat), &E(submat + msd2 - 1), &msd2, &
			WORK(iq), &IWORK(iqptr), &IWORK(iprmpt), &IWORK(iperm)
			, &IWORK(igivpt), &IWORK(igivcl), &WORK(igivnm), &
			WORK(iwrem), &IWORK(subpbs + 1), info);
	    }
	    if (*info != 0) {
		goto L130;
	    }
	    IWORK(i / 2 + 1) = IWORK(i + 2);
/* L90: */
	}
	subpbs /= 2;
	++curlvl;
	goto L80;
    }

/*     end while   

       Re-merge the eigenvalues/vectors which were deflated at the final 
  
       merge step. */

    if (*icompq == 1) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    j = IWORK(indxq + i);
	    WORK(i) = D(j);

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(qsiz, &QSTORE(1,j), &c__1, &Q(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(qsiz, &QSTORE(1,j), &c__1, &Q(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(qsiz, &QSTORE(1,j), &c__1, &Q(1,i), &c__1);
#endif

/* L100: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif

    } else if (*icompq == 2) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    j = IWORK(indxq + i);
	    WORK(i) = D(j);

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(n, &Q(1,j), &c__1, &WORK(*n * i + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(n, &Q(1,j), &c__1, &WORK(*n * i + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(n, &Q(1,j), &c__1, &WORK(*n * i + 1), &c__1);
#endif

/* L110: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_("A", n, n, &WORK(*n + 1), n, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy("A", n, n, &WORK(*n + 1), n, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_("A", n, n, &WORK(*n + 1), n, &Q(1,1), ldq);
#endif

    } else {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    j = IWORK(indxq + i);
	    WORK(i) = D(j);
/* L120: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &WORK(1), &c__1, &D(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &WORK(1), &c__1, &D(1), &c__1);
#endif

    }
    goto L140;

L130:
    *info = submat * (*n + 1) + submat + matsiz - 1;

L140:
    return;

/*     End of DLAED0 */

} /* dlaed0_ */

