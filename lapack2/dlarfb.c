#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlarfb_(char *side, char *trans, char *direct, char *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlarfb(char *side, char *trans, char *direct, char *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlarfb_(char *side, char *trans, char *direct, char *
#endif

	storev, int *m, int *n, int *k, LONG DOUBLE *v, int *
	ldv, LONG DOUBLE *t, int *ldt, LONG DOUBLE *c, int *ldc, 
	LONG DOUBLE *work, int *ldwork)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFB applies a real block reflector H or its transpose H' to a   
    real m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H' from the Left   
            = 'R': apply H or H' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H' (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) LONG DOUBLE PRECISION array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L' 
  
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R' 
  
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= MAX(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= MAX(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) LONG DOUBLE PRECISION array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= MAX(1,M).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= MAX(1,N);   
            if SIDE = 'R', LDWORK >= MAX(1,M).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b14 = 1.;
    static LONG DOUBLE c_b25 = -1.;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int i, j;

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
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dtrmm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qtrmm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qtrmm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static char transt[1];




#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]
#define WORK(I,J) work[(I)-1 + ((J)-1)* ( *ldwork)]

    if (*m <= 0 || *n <= 0) {
	return;
    }

    if (lsame_(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, "C")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(n, &C(j,1), ldc, &WORK(1,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(n, &C(j,1), ldc, &WORK(1,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(n, &C(j,1), ldc, &WORK(1,j), &
#endif

			    c__1);
/* L10: */
		}

/*              W := W * V1 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif

			 &V(1,1), ldv, &WORK(1,1), ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif

			    C(*k+1,1), ldc, &V(*k+1,1), ldv,
			     &c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif

			    V(*k+1,1), ldv, &WORK(1,1), 
			    ldwork, &c_b14, &C(*k+1,1), ldc)
			    ;
		}

/*              W := W * V1' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif

			V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif

/* L40: */
		}

/*              W := W * V1 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif

			 &V(1,1), ldv, &WORK(1,1), ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "No transpose", m, k, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "No transpose", m, k, &i__1, &
#endif

			    c_b14, &C(1,*k+1), ldc, &V(*k+1,1), ldv, &c_b14, &WORK(1,1), 
			    ldwork);
		}

/*              W := W * T  or  W * T' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif

			    WORK(1,1), ldwork, &V(*k+1,1), 
			    ldv, &c_b14, &C(1,*k+1), ldc);
		}

/*              W := W * V1' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif

			V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif

/* L70: */
		}

/*              W := W * V2 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif

			 &V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
#endif

			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
#endif

			    V(1,1), ldv, &WORK(1,1), ldwork, &
			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif

			V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif

/* L100: */
		}

/*              W := W * V2 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif

			 &V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "No transpose", m, k, &i__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "No transpose", m, k, &i__1, &
#endif

			    c_b14, &C(1,1), ldc, &V(1,1), ldv, &
			    c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T  or  W * T' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
#endif

			    WORK(1,1), ldwork, &V(1,1), ldv, &
			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif

			V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, "R")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(n, &C(j,1), ldc, &WORK(1,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(n, &C(j,1), ldc, &WORK(1,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(n, &C(j,1), ldc, &WORK(1,j), &
#endif

			    c__1);
/* L130: */
		}

/*              W := W * V1' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
#endif

			V(1,1), ldv, &WORK(1,1), ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(*k+1,1), ldc, &V(1,*k+1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(*k+1,1), ldc, &V(1,*k+1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(*k+1,1), ldc, &V(1,*k+1), 
#endif

			    ldv, &c_b14, &WORK(1,1), ldwork);
		}

/*              W := W * T'  or  W * T */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,*k+1), ldv, &WORK(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,*k+1), ldv, &WORK(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,*k+1), ldv, &WORK(1,1), 
#endif

			    ldwork, &c_b14, &C(*k+1,1), ldc);
		}

/*              W := W * V1 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
#endif

			 &V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
#endif

/* L160: */
		}

/*              W := W * V1' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
#endif

			V(1,1), ldv, &WORK(1,1), ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif

			    C(1,*k+1), ldc, &V(1,*k+1), ldv, &c_b14, &WORK(1,1), 
			    ldwork);
		}

/*              W := W * T  or  W * T' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "No transpose", m, &i__1, k, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "No transpose", m, &i__1, k, &
#endif

			    c_b25, &WORK(1,1), ldwork, &V(1,*k+1), ldv, &c_b14, &C(1,*k+1), ldc);
		}

/*              W := W * V1 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
#endif

			 &V(1,1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
#endif

/* L190: */
		}

/*              W := W * V2' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
#endif

			V(1,*m-*k+1), ldv, &WORK(1,1)
			, ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(1,1), ldc, &V(1,1), ldv, &c_b14, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(1,1), ldc, &V(1,1), ldv, &c_b14, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(1,1), ldc, &V(1,1), ldv, &c_b14, &WORK(1,1), ldwork);
#endif

		}

/*              W := W * T'  or  W * T */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,1), ldv, &WORK(1,1), ldwork, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,1), ldv, &WORK(1,1), ldwork, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,1), ldv, &WORK(1,1), ldwork, &
#endif

			    c_b14, &C(1,1), ldc);
		}

/*              W := W * V2 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
#endif

			 &V(1,*m-*k+1), ldv, &WORK(1,1), ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *n;
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
#endif

/* L220: */
		}

/*              W := W * V2' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
#endif

			V(1,*n-*k+1), ldv, &WORK(1,1)
			, ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
#endif

			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork);
		}

/*              W := W * T  or  W * T' */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork);
#endif


/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemm("No transpose", "No transpose", m, &i__1, k, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemm_("No transpose", "No transpose", m, &i__1, k, &
#endif

			    c_b25, &WORK(1,1), ldwork, &V(1,1), 
			    ldv, &c_b14, &C(1,1), ldc);
		}

/*              W := W * V2 */


#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
#endif

			 &V(1,*n-*k+1), ldv, &WORK(1,1), ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return;

/*     End of DLARFB */

} /* dlarfb_ */

