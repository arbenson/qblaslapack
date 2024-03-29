#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlangb_(char *norm, int *n, int *kl, int *ku, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlangb(char *norm, int *n, int *kl, int *ku, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlangb_(char *norm, int *n, int *kl, int *ku, 
#endif

	LONG DOUBLE *ab, int *ldab, LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANGB  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the element of  largest absolute value  of an 
  
    n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals. 
  

    Description   
    ===========   

    DLANGB returns the value   

       DLANGB = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANGB as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANGB is   
            set to zero.   

    KL      (input) INTEGER   
            The number of sub-diagonals of the matrix A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of super-diagonals of the matrix A.  KU >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The band matrix A, stored in rows 1 to KL+KU+1.  The j-th   
            column of A is stored in the j-th column of the array AB as   
            follows:   
            AB(ku+1+i-j,j) = A(i,j) for MAX(1,j-ku)<=i<=MIN(n,j+kl).   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KL+KU+1.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  



    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5, i__6;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static int i, j, k, l;
    static LONG DOUBLE scale;
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

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__2 = *ku + 2 - j;
/* Computing MIN */
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
	    i__3 = MIN(i__4,i__5);
	    for (i = MAX(*ku+2-j,1); i <= MIN(*n+*ku+1-j,*kl+*ku+1); ++i) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = AB(i,j), ABS(d__1));
		value = MAX(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = 0.;
/* Computing MAX */
	    i__3 = *ku + 2 - j;
/* Computing MIN */
	    i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
	    i__2 = MIN(i__4,i__5);
	    for (i = MAX(*ku+2-j,1); i <= MIN(*n+*ku+1-j,*kl+*ku+1); ++i) {
		sum += (d__1 = AB(i,j), ABS(d__1));
/* L30: */
	    }
	    value = MAX(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK(i) = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    k = *ku + 1 - j;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
	    i__5 = *n, i__6 = j + *kl;
	    i__4 = MIN(i__5,i__6);
	    for (i = MAX(1,j-*ku); i <= MIN(*n,j+*kl); ++i) {
		WORK(i) += (d__1 = AB(k+i,j), ABS(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = MAX(d__1,d__2);
/* L80: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__4 = 1, i__2 = j - *ku;
	    l = MAX(i__4,i__2);
	    k = *ku + 1 - j + l;
/* Computing MIN */
	    i__2 = *n, i__3 = j + *kl;
	    i__4 = MIN(i__2,i__3) - l + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlassq_(&i__4, &AB(k,j), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlassq(&i__4, &AB(k,j), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlassq_(&i__4, &AB(k,j), &c__1, &scale, &sum);
#endif

/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANGB */

} /* dlangb_ */

