#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlamrg_(int *n1, int *n2, LONG DOUBLE *a, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlamrg(int *n1, int *n2, LONG DOUBLE *a, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlamrg_(int *n1, int *n2, LONG DOUBLE *a, int 
#endif

	*dtrd1, int *dtrd2, int *index)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLAMRG will create a permutation list which will merge the elements   
    of A (which is composed of two independently sorted sets) into a   
    single set which is sorted in ascending order.   

    Arguments   
    =========   

    N1     (input) INTEGER   
    N2     (input) INTEGER   
           These arguements contain the respective lengths of the two   
           sorted lists to be merged.   

    A      (input) LONG DOUBLE PRECISION array, dimension (N1+N2)   
           The first N1 elements of A contain a list of numbers which   
           are sorted in either ascending or descending order.  Likewise 
  
           for the final N2 elements.   

    DTRD1  (input) INTEGER   
    DTRD2  (input) INTEGER   
           These are the strides to be taken through the array A.   
           Allowable strides are 1 and -1.  They indicate whether a   
           subset of A is sorted in ascending (DTRDx = 1) or descending   
           (DTRDx = -1) order.   

    INDEX  (output) INTEGER array, dimension (N1+N2)   
           On exit this array will contain a permutation such that   
           if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be   
           sorted in ascending order.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    /* Local variables */
    static int i, ind1, ind2, n1sv, n2sv;


#define INDEX(I) index[(I)-1]
#define A(I) a[(I)-1]


    n1sv = *n1;
    n2sv = *n2;
    if (*dtrd1 > 0) {
	ind1 = 1;
    } else {
	ind1 = *n1;
    }
    if (*dtrd2 > 0) {
	ind2 = *n1 + 1;
    } else {
	ind2 = *n1 + *n2;
    }
    i = 1;
/*     while ( (N1SV > 0) & (N2SV > 0) ) */
L10:
    if (n1sv > 0 && n2sv > 0) {
	if (A(ind1) <= A(ind2)) {
	    INDEX(i) = ind1;
	    ++i;
	    ind1 += *dtrd1;
	    --n1sv;
	} else {
	    INDEX(i) = ind2;
	    ++i;
	    ind2 += *dtrd2;
	    --n2sv;
	}
	goto L10;
    }
/*     end while */
    if (n1sv == 0) {
	i__1 = n2sv;
	for (n1sv = 1; n1sv <= n2sv; ++n1sv) {
	    INDEX(i) = ind2;
	    ++i;
	    ind2 += *dtrd2;
/* L20: */
	}
    } else {
/*     N2SV .EQ. 0 */
	i__1 = n1sv;
	for (n2sv = 1; n2sv <= n1sv; ++n2sv) {
	    INDEX(i) = ind1;
	    ++i;
	    ind1 += *dtrd1;
/* L30: */
	}
    }

    return;

/*     End of DLAMRG */

} /* dlamrg_ */

