#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgebak_(char *job, char *side, int *n, int *ilo, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgebak(char *job, char *side, int *n, int *ilo, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgebak_(char *job, char *side, int *n, int *ilo, 
#endif

	int *ihi, LONG DOUBLE *scale, int *m, LONG DOUBLE *v, int *
	ldv, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEBAK forms the right or left eigenvectors of a real general matrix 
  
    by backward transformation on the computed eigenvectors of the   
    balanced matrix output by DGEBAL.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the type of backward transformation required:   
            = 'N', do nothing, return immediately;   
            = 'P', do backward transformation for permutation only;   
            = 'S', do backward transformation for scaling only;   
            = 'B', do backward transformations for both permutation and   
                   scaling.   
            JOB must be the same as the argument JOB supplied to DGEBAL. 
  

    SIDE    (input) CHARACTER*1   
            = 'R':  V contains right eigenvectors;   
            = 'L':  V contains left eigenvectors.   

    N       (input) INTEGER   
            The number of rows of the matrix V.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            The ints ILO and IHI determined by DGEBAL.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    SCALE   (input) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutation and scaling factors, as returned   
            by DGEBAL.   

    M       (input) INTEGER   
            The number of columns of the matrix V.  M >= 0.   

    V       (input/output) LONG DOUBLE PRECISION array, dimension (LDV,M)   
            On entry, the matrix of right or left eigenvectors to be   
            transformed, as returned by DHSEIN or DTREVC.   
            On exit, V is overwritten by the transformed eigenvectors.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Decode and Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int i, k;
    static LONG DOUBLE s;

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
    static int ii;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int rightv;


#define SCALE(I) scale[(I)-1]

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
    } else if (*ilo < 1 || *ilo > MAX(1,*n)) {
	*info = -4;
    } else if (*ihi < MIN(*ilo,*n) || *ihi > *n) {
	*info = -5;
    } else if (*m < 0) {
	*info = -7;
    } else if (*ldv < MAX(1,*n)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEBAK", &i__1);
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

	if (rightv) {
	    i__1 = *ihi;
	    for (i = *ilo; i <= *ihi; ++i) {
		s = SCALE(i);

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(m, &s, &V(i,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(m, &s, &V(i,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(m, &s, &V(i,1), ldv);
#endif

/* L10: */
	    }
	}

	if (leftv) {
	    i__1 = *ihi;
	    for (i = *ilo; i <= *ihi; ++i) {
		s = 1. / SCALE(i);

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(m, &s, &V(i,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(m, &s, &V(i,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(m, &s, &V(i,1), ldv);
#endif

/* L20: */
	    }
	}

    }

/*     Backward permutation   

       For  I = ILO-1 step -1 until 1,   
                IHI+1 step 1 until N do -- */

L30:
    if (lsame_(job, "P") || lsame_(job, "B")) {
	if (rightv) {
	    i__1 = *n;
	    for (ii = 1; ii <= *n; ++ii) {
		i = ii;
		if (i >= *ilo && i <= *ihi) {
		    goto L40;
		}
		if (i < *ilo) {
		    i = *ilo - ii;
		}
		k = (int) SCALE(i);
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
	}

	if (leftv) {
	    i__1 = *n;
	    for (ii = 1; ii <= *n; ++ii) {
		i = ii;
		if (i >= *ilo && i <= *ihi) {
		    goto L50;
		}
		if (i < *ilo) {
		    i = *ilo - ii;
		}
		k = (int) SCALE(i);
		if (k == i) {
		    goto L50;
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

L50:
		;
	    }
	}
    }

    return;

/*     End of DGEBAK */

} /* dgebak_ */

