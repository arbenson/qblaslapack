#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dstedc_(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qstedc(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qstedc_(char *compz, int *n, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int 
	*lwork, int *iwork, int *liwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEDC computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the divide and conquer method.   
    The eigenvectors of a full or band real symmetric matrix can also be 
  
    found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this   
    matrix to tridiagonal form.   

    This code makes very mild assumptions about floating point   
    arithmetic. It will work on machines with a guard digit in   
    add/subtract, or on those binary machines without guard digits   
    which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.   
    It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.  See DLAED3 for details.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'I':  Compute eigenvectors of tridiagonal matrix also.   
            = 'V':  Compute eigenvectors of original dense symmetric   
                    matrix also.  On entry, Z contains the orthogonal   
                    matrix used to reduce the original matrix to   
                    tridiagonal form.   

    N       (input) INTEGER   
            The dimension of the symmetric tridiagonal matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, the subdiagonal elements of the tridiagonal matrix. 
  
            On exit, E has been destroyed.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ,N)   
            On entry, if COMPZ = 'V', then Z contains the orthogonal   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original symmetric matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If  COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If eigenvectors are desired, then LDZ >= MAX(1,N).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if LWORK > 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.   
            If COMPZ = 'V' and N > 1 then LWORK must be at least   
                           ( 1 + 3*N + 2*N*lg N + 3*N**2 ),   
                           where lg( N ) = smallest int k such   
                           that 2**k >= N.   
            If COMPZ = 'I' and N > 1 then LWORK must be at least   
                           ( 1 + 3*N + 2*N*lg N + 2*N**2 ).   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK. 
  

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.   
            If COMPZ = 'V' and N > 1 then LIWORK must be at least   
                           ( 6 + 6*N + 5*N*lg N ).   
            If COMPZ = 'I' and N > 1 then LIWORK must be at least   
                           ( 2 + 5*N ).   

    INFO    (output) INTEGER   
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
    static LONG DOUBLE c_b12 = 0.;
    static LONG DOUBLE c_b13 = 1.;
    static int c__0 = 0;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE tiny;
    static int i, j, k, m;
    static LONG DOUBLE p;

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
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dswap_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static int lwmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaed0_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed0(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaed0_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     int *, LONG DOUBLE *, int *, int *);
    static int start, ii;

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
    extern /* Subroutine */ void dlascl_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlascl_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dlacpy_(char *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qlacpy(char *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qlacpy_(char *, int *, int 
#endif

	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlaset_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaset(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaset_(char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *), xerbla_(char *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsterf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif


#ifdef PETSC_PREFIX_SUFFIX
	     int *), dlasrt_(char *, int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *), qlasrt(char *, int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *), qlasrt_(char *, int *, LONG DOUBLE *, int *);
#endif

    static int liwmin, icompz;

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
    static LONG DOUBLE orgnrm;
    static int dtrtrw, storez, end, lgn;
    static LONG DOUBLE eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (*n <= 1 || icompz <= 0) {
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
	if (icompz == 1) {
/* Computing 2nd power */
	    i__1 = *n;
	    lwmin = *n * 3 + 1 + (*n << 1) * lgn + i__1 * i__1 * 3;
	    liwmin = *n * 6 + 6 + *n * 5 * lgn;
	} else if (icompz == 2) {
/* Computing 2nd power */
	    i__1 = *n;
	    lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 1);
	    liwmin = *n * 5 + 2;
	}
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < MAX(1,*n))) {
	*info = -6;
    } else if (*lwork < lwmin) {
	*info = -8;
    } else if (*liwork < liwmin) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEDC", &i__1);
	goto L50;
    }

/*     Quick return if possible */

    if (*n == 0) {
	goto L50;
    }
    if (*n == 1) {
	if (icompz != 0) {
	    Z(1,1) = 1.;
	}
	goto L50;
    }

/*     If the following conditional clause is removed, then the routine   
       will use the Divide and Conquer routine to compute only the   
       eigenvalues, which requires (3N + 3N**2) real workspace and   
       (2 + 5N + 2N lg(N)) int workspace.   
       Since on many architectures DSTERF is much faster than any other   
       algorithm for finding eigenvalues only, it is used here   
       as the default.   

       If COMPZ = 'N', use DSTERF to compute the eigenvalues. */

    if (icompz == 0) {

#ifdef PETSC_PREFIX_SUFFIX
	dsterf_(n, &D(1), &E(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsterf(n, &D(1), &E(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsterf_(n, &D(1), &E(1), info);
#endif

	goto L50;
    }

/*     If N is smaller than the minimum divide size (SMLSIZ+1), then   
       solve the problem with another solver. */

    if (*n <= 25) {
	if (icompz == 0) {

#ifdef PETSC_PREFIX_SUFFIX
	    dsterf_(n, &D(1), &E(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsterf(n, &D(1), &E(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsterf_(n, &D(1), &E(1), info);
#endif

	    goto L50;
	} else if (icompz == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_("I", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif

	    goto L50;
	} else {

#ifdef PETSC_PREFIX_SUFFIX
	    dsteqr_("V", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsteqr("V", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsteqr_("V", n, &D(1), &E(1), &Z(1,1), ldz, &WORK(1), info);
#endif

	    goto L50;
	}
    }

/*     If COMPZ = 'V', the Z matrix must be stored elsewhere for later   
       use. */

    if (icompz == 1) {
	storez = *n * *n + 1;
    } else {
	storez = 1;
    }

    if (icompz == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif

    }

/*     Scale. */


#ifdef PETSC_PREFIX_SUFFIX
    orgnrm = dlanst_("M", n, &D(1), &E(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    orgnrm = qlanst("M", n, &D(1), &E(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    orgnrm = qlanst_("M", n, &D(1), &E(1));
#endif

    if (orgnrm == 0.) {
	goto L50;
    }


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif


    start = 1;

/*     while ( START <= N ) */

L10:
    if (start <= *n) {

/*     Let END be the position of the next subdiagonal entry such that
   
       E( END ) <= TINY or END = N if no such subdiagonal exists.  The
   
       matrix identified by the elements between START and END   
       constitutes an independent sub-problem. */

	end = start;
L20:
	if (end < *n) {
	    d__1 = D(end);d__2 = D(end + 1);
            tiny = eps * sqrt(( ABS(d__1))) * sqrt((ABS(d__2)));
	    if ((d__1 = E(end), ABS(d__1)) > tiny) {
		++end;
		goto L20;
	    }
	}

/*        (Sub) Problem determined.  Compute its size and solve it. */

	m = end - start + 1;
	if (m == 1) {
	    start = end + 1;
	    goto L10;
	}
	if (m > 25) {
	    *info = 25;

/*           Scale. */


#ifdef PETSC_PREFIX_SUFFIX
	    orgnrm = dlanst_("M", &m, &D(start), &E(start));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    orgnrm = qlanst("M", &m, &D(start), &E(start));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    orgnrm = qlanst_("M", &m, &D(start), &E(start));
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &m, &c__1, &D(start), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &orgnrm, &c_b13, &m, &c__1, &D(start), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &m, &c__1, &D(start), 
#endif

		    &m, info);
	    i__1 = m - 1;
	    i__2 = m - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &i__1, &c__1, &E(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &orgnrm, &c_b13, &i__1, &c__1, &E(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &orgnrm, &c_b13, &i__1, &c__1, &E(
#endif

		    start), &i__2, info);

	    if (icompz == 1) {
		dtrtrw = 1;
	    } else {
		dtrtrw = start;
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    dlaed0_(&icompz, n, &m, &D(start), &E(start), &Z(dtrtrw,start), ldz, &WORK(1), n, &WORK(storez), &IWORK(1), info)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlaed0(&icompz, n, &m, &D(start), &E(start), &Z(dtrtrw,start), ldz, &WORK(1), n, &WORK(storez), &IWORK(1), info)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlaed0_(&icompz, n, &m, &D(start), &E(start), &Z(dtrtrw,start), ldz, &WORK(1), n, &WORK(storez), &IWORK(1), info)
#endif

		    ;
	    if (*info != 0) {
		*info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m 
			+ 1) + start - 1;
		goto L50;
	    }

/*           Scale back. */


#ifdef PETSC_PREFIX_SUFFIX
	    dlascl_("G", &c__0, &c__0, &c_b13, &orgnrm, &m, &c__1, &D(start), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlascl("G", &c__0, &c__0, &c_b13, &orgnrm, &m, &c__1, &D(start), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlascl_("G", &c__0, &c__0, &c_b13, &orgnrm, &m, &c__1, &D(start), 
#endif

		    &m, info);

	} else {
	    if (icompz == 1) {

/*     Since QR won't update a Z matrix which is larger than t
he   
       length of D, we must solve the sub-problem in a workspa
ce and   
       then multiply back into Z. */


#ifdef PETSC_PREFIX_SUFFIX
		dsteqr_("I", &m, &D(start), &E(start), &WORK(1), &m, &WORK(m *
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsteqr("I", &m, &D(start), &E(start), &WORK(1), &m, &WORK(m *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsteqr_("I", &m, &D(start), &E(start), &WORK(1), &m, &WORK(m *
#endif

			 m + 1), info);

#ifdef PETSC_PREFIX_SUFFIX
		dlacpy_("A", n, &m, &Z(1,start), ldz, &WORK(storez)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlacpy("A", n, &m, &Z(1,start), ldz, &WORK(storez)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlacpy_("A", n, &m, &Z(1,start), ldz, &WORK(storez)
#endif

			, n);

#ifdef PETSC_PREFIX_SUFFIX
		dgemm_("N", "N", n, &m, &m, &c_b13, &WORK(storez), ldz, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemm("N", "N", n, &m, &m, &c_b13, &WORK(storez), ldz, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemm_("N", "N", n, &m, &m, &c_b13, &WORK(storez), ldz, &WORK(
#endif

			1), &m, &c_b12, &Z(1,start), ldz);
	    } else if (icompz == 2) {

#ifdef PETSC_PREFIX_SUFFIX
		dsteqr_("I", &m, &D(start), &E(start), &Z(start,start), ldz, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsteqr("I", &m, &D(start), &E(start), &Z(start,start), ldz, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsteqr_("I", &m, &D(start), &E(start), &Z(start,start), ldz, &WORK(1), info);
#endif

	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dsterf_(&m, &D(start), &E(start), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsterf(&m, &D(start), &E(start), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsterf_(&m, &D(start), &E(start), info);
#endif

	    }
	    if (*info != 0) {
		*info = start * (*n + 1) + end;
		goto L50;
	    }
	}

	start = end + 1;
	goto L10;
    }

/*     endwhile   

       If the problem split any number of times, then the eigenvalues   
       will not be properly ordered.  Here we permute the eigenvalues   
       (and the associated eigenvectors) into ascending order. */

    if (m != *n) {
	if (icompz == 0) {

/*        Use Quick Sort */


#ifdef PETSC_PREFIX_SUFFIX
	    dlasrt_("I", n, &D(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasrt("I", n, &D(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasrt_("I", n, &D(1), info);
#endif


	} else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	    i__1 = *n;
	    for (ii = 2; ii <= *n; ++ii) {
		i = ii - 1;
		k = i;
		p = D(i);
		i__2 = *n;
		for (j = ii; j <= *n; ++j) {
		    if (D(j) < p) {
			k = j;
			p = D(j);
		    }
/* L30: */
		}
		if (k != i) {
		    D(k) = D(i);
		    D(i) = p;

#ifdef PETSC_PREFIX_SUFFIX
		    dswap_(n, &Z(1,i), &c__1, &Z(1,k), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qswap(n, &Z(1,i), &c__1, &Z(1,k), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qswap_(n, &Z(1,i), &c__1, &Z(1,k), &
#endif

			    c__1);
		}
/* L40: */
	    }
	}
    }

L50:
    if (*lwork > 0) {
	WORK(1) = (LONG DOUBLE) lwmin;
    }
    if (*liwork > 0) {
	IWORK(1) = liwmin;
    }
    return;

/*     End of DSTEDC */

} /* dstedc_ */

