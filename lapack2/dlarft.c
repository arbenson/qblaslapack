#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlarft_(char *direct, char *storev, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlarft(char *direct, char *storev, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlarft_(char *direct, char *storev, int *n, int *
#endif

	k, LONG DOUBLE *v, int *ldv, LONG DOUBLE *tau, LONG DOUBLE *t, 
	int *ldt)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFT forms the triangular factor T of a real block reflector H   
    of order n, which is defined as a product of k elementary reflectors. 
  

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; 
  

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. 
  

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) LONG DOUBLE PRECISION array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= MAX(1,N); if STOREV = 'R', LDV >= K. 
  

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) LONG DOUBLE PRECISION array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is 
  
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define 
  
    the H(i) is best illustrated by the following example with n = 5 and 
  
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': 
  

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 ) 
  
                     ( v1  1    )                     (     1 v2 v2 v2 ) 
  
                     ( v1 v2  1 )                     (        1 v3 v3 ) 
  
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': 
  

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) 
  
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    ) 
  
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) 
  
                     (     1 v3 )   
                     (        1 )   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b8 = 0.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Local variables */
    static int i, j;
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

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dtrmv_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qtrmv(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qtrmv_(char *, 
#endif

	    char *, char *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static LONG DOUBLE vii;



#define TAU(I) tau[(I)-1]

#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    if (*n == 0) {
	return;
    }

    if (lsame_(direct, "F")) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__2 = i;
		for (j = 1; j <= i; ++j) {
		    T(j,i) = 0.;
/* L10: */
		}
	    } else {

/*              general case */

		vii = V(i,i);
		V(i,i) = 1.;
		if (lsame_(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' 
* V(i:n,i) */

		    i__2 = *n - i + 1;
		    i__3 = i - 1;
		    d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("Transpose", &i__2, &i__3, &d__1, &V(i,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("Transpose", &i__2, &i__3, &d__1, &V(i,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("Transpose", &i__2, &i__3, &d__1, &V(i,1), 
#endif

			    ldv, &V(i,i), &c__1, &c_b8, &T(1,i), &c__1);
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) *
 V(i,i:n)' */

		    i__2 = i - 1;
		    i__3 = *n - i + 1;
		    d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("No transpose", &i__2, &i__3, &d__1, &V(1,i), ldv, &V(i,i), ldv, &c_b8, &T(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("No transpose", &i__2, &i__3, &d__1, &V(1,i), ldv, &V(i,i), ldv, &c_b8, &T(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("No transpose", &i__2, &i__3, &d__1, &V(1,i), ldv, &V(i,i), ldv, &c_b8, &T(1,i), &c__1);
#endif

		}
		V(i,i) = vii;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmv("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
#endif

		T(i,i) = TAU(i);
	    }
/* L20: */
	}
    } else {
	for (i = *k; i >= 1; --i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i; j <= *k; ++j) {
		    T(j,i) = 0.;
/* L30: */
		}
	    } else {

/*              general case */

		if (i < *k) {
		    if (lsame_(storev, "C")) {
			vii = V(*n-*k+i,i);
			V(*n-*k+i,i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1
:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i;
			i__2 = *k - i;
			d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
			dgemv_("Transpose", &i__1, &i__2, &d__1, &V(1,i+1), ldv, &V(1,i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemv("Transpose", &i__1, &i__2, &d__1, &V(1,i+1), ldv, &V(1,i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemv_("Transpose", &i__1, &i__2, &d__1, &V(1,i+1), ldv, &V(1,i), &c__1, &
#endif

				c_b8, &T(i+1,i), &c__1);
			V(*n-*k+i,i) = vii;
		    } else {
			vii = V(i,*n-*k+i);
			V(i,*n-*k+i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k
+i) * V(i,1:n-k+i)' */

			i__1 = *k - i;
			i__2 = *n - *k + i;
			d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
			dgemv_("No transpose", &i__1, &i__2, &d__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b8, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemv("No transpose", &i__1, &i__2, &d__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b8, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemv_("No transpose", &i__1, &i__2, &d__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b8, &
#endif

				T(i+1,i), &c__1);
			V(i,*n-*k+i) = vii;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,
i) */

		    i__1 = *k - i;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmv_("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmv("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmv_("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
#endif

			    , &c__1);
		}
		T(i,i) = TAU(i);
	    }
/* L40: */
	}
    }
    return;

/*     End of DLARFT */

} /* dlarft_ */

