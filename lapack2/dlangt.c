#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlangt_(char *norm, int *n, LONG DOUBLE *dl, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlangt(char *norm, int *n, LONG DOUBLE *dl, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlangt_(char *norm, int *n, LONG DOUBLE *dl, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *du)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLANGT  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    real tridiagonal matrix A.   

    Description   
    ===========   

    DLANGT returns the value   

       DLANGT = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  MAX(ABS(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANGT as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANGT is   
            set to zero.   

    DL      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) sub-diagonal elements of A.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) super-diagonal elements of A.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE ret_val, d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    /* Local variables */
    static int i;
    static LONG DOUBLE scale;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE anorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlassq_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlassq(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlassq_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE sum;



#define DU(I) du[(I)-1]
#define D(I) d[(I)-1]
#define DL(I) dl[(I)-1]


    if (*n <= 0) {
	anorm = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	anorm = (d__1 = D(*n), ABS(d__1));
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = DL(i), ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = D(i), ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = DU(i), ABS(d__1));
	    anorm = MAX(d__2,d__3);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = ABS(D(1));
	} else {
/* Computing MAX */
	    d__3 = ABS(D(1)) + ABS(DL(1)), d__4 = (d__1 = D(*n), ABS(d__1)) + 
		    (d__2 = DU(*n - 1), ABS(d__2));
	    anorm = MAX(d__3,d__4);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = D(i), ABS(d__1)) + (d__2 = DL(i),
			 ABS(d__2)) + (d__3 = DU(i - 1), ABS(d__3));
		anorm = MAX(d__4,d__5);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (*n == 1) {
	    anorm = ABS(D(1));
	} else {
/* Computing MAX */
	    d__3 = ABS(D(1)) + ABS(DU(1)), d__4 = (d__1 = D(*n), ABS(d__1)) + 
		    (d__2 = DL(*n - 1), ABS(d__2));
	    anorm = MAX(d__3,d__4);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		d__4 = anorm, d__5 = (d__1 = D(i), ABS(d__1)) + (d__2 = DU(i),
			 ABS(d__2)) + (d__3 = DL(i - 1), ABS(d__3));
		anorm = MAX(d__4,d__5);
/* L30: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;

#ifdef PETSC_PREFIX_SUFFIX
	dlassq_(n, &D(1), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlassq(n, &D(1), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlassq_(n, &D(1), &c__1, &scale, &sum);
#endif

	if (*n > 1) {
	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlassq_(&i__1, &DL(1), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlassq(&i__1, &DL(1), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlassq_(&i__1, &DL(1), &c__1, &scale, &sum);
#endif

	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlassq_(&i__1, &DU(1), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlassq(&i__1, &DU(1), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlassq_(&i__1, &DU(1), &c__1, &scale, &sum);
#endif

	}
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of DLANGT */

} /* dlangt_ */

