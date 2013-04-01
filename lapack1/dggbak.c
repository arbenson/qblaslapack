#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggbak_(char *job, char *side, int *n, int *ilo, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggbak(char *job, char *side, int *n, int *ilo, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggbak_(char *job, char *side, int *n, int *ilo, 
#endif

	int *ihi, LONG DOUBLE *lscale, LONG DOUBLE *rscale, int *m, 
	LONG DOUBLE *v, int *ldv, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAK forms the right or left eigenvectors of a real generalized   
    eigenvalue problem A*x = lambda*B*x, by backward transformation on   
    the computed eigenvectors of the balanced pair of matrices output by 
  
    DGGBAL.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the type of backward transformation required:   
            = 'N':  do nothing, return immediately;   
            = 'P':  do backward transformation for permutation only;   
            = 'S':  do backward transformation for scaling only;   
            = 'B':  do backward transformations for both permutation and 
  
                    scaling.   
            JOB must be the same as the argument JOB supplied to DGGBAL. 
  

    SIDE    (input) CHARACTER*1   
            = 'R':  V contains right eigenvectors;   
            = 'L':  V contains left eigenvectors.   

    N       (input) INTEGER   
            The number of rows of the matrix V.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            The ints ILO and IHI determined by DGGBAL.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    LSCALE  (input) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the left side of A and B, as returned by DGGBAL.   

    RSCALE  (input) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the right side of A and B, as returned by DGGBAL.   

    M       (input) INTEGER   
            The number of columns of the matrix V.  M >= 0.   

    V       (input/output) LONG DOUBLE PRECISION array, dimension (LDV,M)   
            On entry, the matrix of right or left eigenvectors to be   
            transformed, as returned by DTGEVC.   
            On exit, V is overwritten by the transformed eigenvectors.   

    LDV     (input) INTEGER   
            The leading dimension of the matrix V. LDV >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. Ward, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int i, k;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
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
    static long int leftv;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int rightv;


#define LSCALE(I) lscale[(I)-1]
#define RSCALE(I) rscale[(I)-1]

#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]

    rightv = lsame_(side, "R");
    leftv = lsame_(side, "L");

    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(
	    job, "S") && ! lsame_(job, "B")) {
	*info = -1;
    } else if (! rightv && ! leftv) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi < *ilo || *ihi > MAX(1,*n)) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*ldv < MAX(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGBAK", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }
    if (*m == 0) {
	return;
    }
    if (lsame_(job, "N")) {
	return;
    }

    if (*ilo == *ihi) {
	goto L30;
    }

/*     Backward balance */

    if (lsame_(job, "S") || lsame_(job, "B")) {

/*        Backward transformation on right eigenvectors */

	if (rightv) {
	    i__1 = *ihi;
	    for (i = *ilo; i <= *ihi; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(m, &RSCALE(i), &V(i,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(m, &RSCALE(i), &V(i,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(m, &RSCALE(i), &V(i,1), ldv);
#endif

/* L10: */
	    }
	}

/*        Backward transformation on left eigenvectors */

	if (leftv) {
	    i__1 = *ihi;
	    for (i = *ilo; i <= *ihi; ++i) {

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(m, &LSCALE(i), &V(i,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(m, &LSCALE(i), &V(i,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(m, &LSCALE(i), &V(i,1), ldv);
#endif

/* L20: */
	    }
	}
    }

/*     Backward permutation */

L30:
    if (lsame_(job, "P") || lsame_(job, "B")) {

/*        Backward permutation on right eigenvectors */

	if (rightv) {
	    if (*ilo == 1) {
		goto L50;
	    }

	    for (i = *ilo - 1; i >= 1; --i) {
		k = (int) RSCALE(i);
		if (k == i) {
		    goto L40;
		}

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif

L40:
		;
	    }

L50:
	    if (*ihi == *n) {
		goto L70;
	    }
	    i__1 = *n;
	    for (i = *ihi + 1; i <= *n; ++i) {
		k = (int) RSCALE(i);
		if (k == i) {
		    goto L60;
		}

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif

L60:
		;
	    }
	}

/*        Backward permutation on left eigenvectors */

L70:
	if (leftv) {
	    if (*ilo == 1) {
		goto L90;
	    }
	    for (i = *ilo - 1; i >= 1; --i) {
		k = (int) LSCALE(i);
		if (k == i) {
		    goto L80;
		}

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif

L80:
		;
	    }

L90:
	    if (*ihi == *n) {
		goto L110;
	    }
	    i__1 = *n;
	    for (i = *ihi + 1; i <= *n; ++i) {
		k = (int) LSCALE(i);
		if (k == i) {
		    goto L100;
		}

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &V(i,1), ldv, &V(k,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &V(i,1), ldv, &V(k,1), ldv);
#endif

L100:
		;
	    }
	}
    }

L110:

    return;

/*     End of DGGBAK */

} /* dggbak_ */

