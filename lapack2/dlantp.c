#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlantp_(char *norm, char *uplo, char *diag, int *n, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlantp(char *norm, char *uplo, char *diag, int *n, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlantp_(char *norm, char *uplo, char *diag, int *n, LONG DOUBLE 
#endif

	*ap, LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANTP  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    triangular matrix A, supplied in packed form.   

    Description   
    ===========   

    DLANTP returns the value   

       DLANTP = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANTP as described   
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
            The order of the matrix A.  N >= 0.  When N = 0, DLANTP is   
            set to zero.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            Note that when DIAG = 'U', the elements of the array AP   
            corresponding to the diagonal elements of the matrix A are   
            not referenced, but are assumed to be one.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1, i__2;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static int i, j, k;
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
#define AP(I) ap[(I)-1]


    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	k = 1;
	if (lsame_(diag, "U")) {
	    value = 1.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + j - 2;
		    for (i = k; i <= k+j-2; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
			value = MAX(d__2,d__3);
/* L10: */
		    }
		    k += j;
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + *n - j;
		    for (i = k + 1; i <= k+*n-j; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
			value = MAX(d__2,d__3);
/* L30: */
		    }
		    k = k + *n - j + 1;
/* L40: */
		}
	    }
	} else {
	    value = 0.;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + j - 1;
		    for (i = k; i <= k+j-1; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
			value = MAX(d__2,d__3);
/* L50: */
		    }
		    k += j;
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = k + *n - j;
		    for (i = k; i <= k+*n-j; ++i) {
/* Computing MAX */
			d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
			value = MAX(d__2,d__3);
/* L70: */
		    }
		    k = k + *n - j + 1;
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	k = 1;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
		    i__2 = k + j - 2;
		    for (i = k; i <= k+j-2; ++i) {
			sum += (d__1 = AP(i), ABS(d__1));
/* L90: */
		    }
		} else {
		    sum = 0.;
		    i__2 = k + j - 1;
		    for (i = k; i <= k+j-1; ++i) {
			sum += (d__1 = AP(i), ABS(d__1));
/* L100: */
		    }
		}
		k += j;
		value = MAX(value,sum);
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (udiag) {
		    sum = 1.;
		    i__2 = k + *n - j;
		    for (i = k + 1; i <= k+*n-j; ++i) {
			sum += (d__1 = AP(i), ABS(d__1));
/* L120: */
		    }
		} else {
		    sum = 0.;
		    i__2 = k + *n - j;
		    for (i = k; i <= k+*n-j; ++i) {
			sum += (d__1 = AP(i), ABS(d__1));
/* L130: */
		    }
		}
		k = k + *n - j + 1;
		value = MAX(value,sum);
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	k = 1;
	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(i) = 1.;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			WORK(i) += (d__1 = AP(k), ABS(d__1));
			++k;
/* L160: */
		    }
		    ++k;
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
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			WORK(i) += (d__1 = AP(k), ABS(d__1));
			++k;
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
		    ++k;
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			WORK(i) += (d__1 = AP(k), ABS(d__1));
			++k;
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
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			WORK(i) += (d__1 = AP(k), ABS(d__1));
			++k;
/* L250: */
		    }
/* L260: */
		}
	    }
	}
	value = 0.;
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
		k = 2;
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
		    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif

		    k += j;
/* L280: */
		}
	    } else {
		scale = 0.;
		sum = 1.;
		k = 1;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&j, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&j, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&j, &AP(k), &c__1, &scale, &sum);
#endif

		    k += j;
/* L290: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		scale = 1.;
		sum = (LONG DOUBLE) (*n);
		k = 2;
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
		    i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif

		    k = k + *n - j + 1;
/* L300: */
		}
	    } else {
		scale = 0.;
		sum = 1.;
		k = 1;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n - j + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif

		    k = k + *n - j + 1;
/* L310: */
		}
	    }
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANTP */

} /* dlantp_ */

