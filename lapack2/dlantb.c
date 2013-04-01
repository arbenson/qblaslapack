#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlantb_(char *norm, char *uplo, char *diag, int *n, int *k,
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlantb(char *norm, char *uplo, char *diag, int *n, int *k,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlantb_(char *norm, char *uplo, char *diag, int *n, int *k,
#endif

	 LONG DOUBLE *ab, int *ldab, LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANTB  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the element of  largest absolute value  of an 
  
    n by n triangular band matrix A,  with ( k + 1 ) diagonals.   

    Description   
    ===========   

    DLANTB returns the value   

       DLANTB = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANTB as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANTB is   
            set to zero.   

    K       (input) INTEGER   
            The number of super-diagonals of the matrix A if UPLO = 'U', 
  
            or the number of sub-diagonals of the matrix A if UPLO = 'L'. 
  
            K >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangular band matrix A, stored in the   
            first k+1 rows of AB.  The j-th column of A is stored   
            in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for MAX(1,j-k)<=i<=j;   
            if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=MIN(n,j+k).   
            Note that when DIAG = 'U', the elements of the array AB   
            corresponding to the diagonal elements of the matrix A are   
            not referenced, but are assumed to be one.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= K+1.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static int i, j, l;
    static LONG DOUBLE scale;
    static long int udiag;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE value;

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



#define WORK(I) work[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	if (lsame_(diag, "U")) {
	    value = 1.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MAX */
		    i__2 = *k + 2 - j;
		    i__3 = *k;
		    for (i = MAX(*k+2-j,1); i <= *k; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
				d__1));
			value = MAX(d__2,d__3);
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__2 = *n + 1 - j, i__4 = *k + 1;
		    i__3 = MIN(i__2,i__4);
		    for (i = 2; i <= MIN(*n+1-j,*k+1); ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
				d__1));
			value = MAX(d__2,d__3);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    value = 0.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MAX */
		    i__3 = *k + 2 - j;
		    i__2 = *k + 1;
		    for (i = MAX(*k+2-j,1); i <= *k+1; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
				d__1));
			value = MAX(d__2,d__3);
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = MIN(i__3,i__4);
		    for (i = 1; i <= MIN(*n+1-j,*k+1); ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
				d__1));
			value = MAX(d__2,d__3);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
/* Computing MAX */
		    i__2 = *k + 2 - j;
		    i__3 = *k;
		    for (i = MAX(*k+2-j,1); i <= *k; ++i) {
			sum += (d__1 = AB(i,j), ABS(d__1));
/* L90: */
		    }
		} else {
		    sum = 0.;
/* Computing MAX */
		    i__3 = *k + 2 - j;
		    i__2 = *k + 1;
		    for (i = MAX(*k+2-j,1); i <= *k+1; ++i) {
			sum += (d__1 = AB(i,j), ABS(d__1));
/* L100: */
		    }
		}
		value = MAX(value,sum);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = MIN(i__3,i__4);
		    for (i = 2; i <= MIN(*n+1-j,*k+1); ++i) {
			sum += (d__1 = AB(i,j), ABS(d__1));
/* L120: */
		    }
		} else {
		    sum = 0.;
/* Computing MIN */
		    i__3 = *n + 1 - j, i__4 = *k + 1;
		    i__2 = MIN(i__3,i__4);
		    for (i = 1; i <= MIN(*n+1-j,*k+1); ++i) {
			sum += (d__1 = AB(i,j), ABS(d__1));
/* L130: */
		    }
		}
		value = MAX(value,sum);
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = *k + 1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *k;
		    i__4 = j - 1;
		    for (i = MAX(1,j-*k); i <= j-1; ++i) {
			WORK(i) += (d__1 = AB(l+i,j), ABS(d__1))
				;
/* L160: */
		    }
/* L170: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.;
/* L180: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = *k + 1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *k;
		    i__3 = j;
		    for (i = MAX(1,j-*k); i <= j; ++i) {
			WORK(i) += (d__1 = AB(l+i,j), ABS(d__1))
				;
/* L190: */
		    }
/* L200: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L210: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__2 = j + *k;
		    i__3 = MIN(i__4,i__2);
		    for (i = j + 1; i <= MIN(*n,j+*k); ++i) {
			WORK(i) += (d__1 = AB(l+i,j), ABS(d__1))
				;
/* L220: */
		    }
/* L230: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 0.;
/* L240: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__2 = j + *k;
		    i__3 = MIN(i__4,i__2);
		    for (i = j; i <= MIN(*n,j+*k); ++i) {
			WORK(i) += (d__1 = AB(l+i,j), ABS(d__1))
				;
/* L250: */
		    }
/* L260: */
		}
	    }
	}
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = MAX(d__1,d__2);
/* L270: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (LONG DOUBLE) (*n);
		if (*k > 0) {
		    i__1 = *n;
		    for (j = 2; j <= *n; ++j) {
/* Computing MIN */
			i__4 = j - 1;
			i__3 = MIN(i__4,*k);
/* Computing MAX */
			i__2 = *k + 2 - j;

#ifdef PETSC_PREFIX_SUFFIX
			dlassq_(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlassq(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlassq_(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, 
#endif

				&scale, &sum);
/* L280: */
		    }
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__4 = j, i__2 = *k + 1;
		    i__3 = MIN(i__4,i__2);
/* Computing MAX */
		    i__5 = *k + 2 - j;

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__3, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif

			    scale, &sum);
/* L290: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (LONG DOUBLE) (*n);
		if (*k > 0) {
		    i__1 = *n - 1;
		    for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
			i__4 = *n - j;
			i__3 = MIN(i__4,*k);

#ifdef PETSC_PREFIX_SUFFIX
			dlassq_(&i__3, &AB(2,j), &c__1, &scale, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlassq(&i__3, &AB(2,j), &c__1, &scale, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlassq_(&i__3, &AB(2,j), &c__1, &scale, &
#endif

				sum);
/* L300: */
		    }
		}
	    } else {
		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		    i__4 = *n - j + 1, i__2 = *k + 1;
		    i__3 = MIN(i__4,i__2);

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__3, &AB(1,j), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__3, &AB(1,j), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__3, &AB(1,j), &c__1, &scale, &sum);
#endif

/* L310: */
		}
	    }
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANTB */

} /* dlantb_ */

