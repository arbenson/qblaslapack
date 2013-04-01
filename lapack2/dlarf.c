#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlarf_(char *side, int *m, int *n, LONG DOUBLE *v,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlarf(char *side, int *m, int *n, LONG DOUBLE *v,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlarf_(char *side, int *m, int *n, LONG DOUBLE *v,
#endif

	 int *incv, LONG DOUBLE *tau, LONG DOUBLE *c, int *ldc, 
	LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARF applies a real elementary reflector H to a real m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) LONG DOUBLE PRECISION array, dimension   
                       (1 + (M-1)*ABS(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*ABS(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) LONG DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
  
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= MAX(1,M).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b4 = 1.;
    static LONG DOUBLE c_b5 = 0.;
    static int c__1 = 1;
    
    /* System generated locals */
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dger_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qger(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qger_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);



#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(side, "L")) {

/*        Form  H * C */

	if (*tau != 0.) {

/*           w := C' * v */


#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("Transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), incv, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("Transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), incv, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("Transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), incv, &
#endif

		    c_b5, &WORK(1), &c__1);

/*           C := C - v * w' */

	    d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	    dger_(m, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qger(m, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qger_(m, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
#endif

		    ldc);
	}
    } else {

/*        Form  C * H */

	if (*tau != 0.) {

/*           w := C * v */


#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("No transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("No transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("No transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), 
#endif

		    incv, &c_b5, &WORK(1), &c__1);

/*           C := C - w * v' */

	    d__1 = -(*tau);

#ifdef PETSC_PREFIX_SUFFIX
	    dger_(m, n, &d__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qger(m, n, &d__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qger_(m, n, &d__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
#endif

		    ldc);
	}
    }
    return;

/*     End of DLARF */

} /* dlarf_ */

