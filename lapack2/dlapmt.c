#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlapmt_(long int *forwrd, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlapmt(long int *forwrd, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlapmt_(long int *forwrd, int *m, int *n, 
#endif

	LONG DOUBLE *x, int *ldx, int *k)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAPMT rearranges the columns of the M by N matrix X as specified   
    by the permutation K(1),K(2),...,K(N) of the ints 1,...,N.   
    If FORWRD = .TRUE.,  forward permutation:   

         X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.   

    If FORWRD = .FALSE., backward permutation:   

         X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.   

    Arguments   
    =========   

    FORWRD  (input) LOGICAL   
            = .TRUE., forward permutation   
            = .FALSE., backward permutation   

    M       (input) INTEGER   
            The number of rows of the matrix X. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix X. N >= 0.   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (LDX,N)   
            On entry, the M by N matrix X.   
            On exit, X contains the permuted matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X, LDX >= MAX(1,M).   

    K       (input) INTEGER array, dimension (N)   
            On entry, K contains the permutation vector.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static LONG DOUBLE temp;
    static int i, j, ii, in;


#define K(I) k[(I)-1]

#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    if (*n <= 1) {
	return;
    }

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	K(i) = -K(i);
/* L10: */
    }

    if (*forwrd) {

/*        Forward permutation */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {

	    if (K(i) > 0) {
		goto L40;
	    }

	    j = i;
	    K(j) = -K(j);
	    in = K(j);

L20:
	    if (K(in) > 0) {
		goto L40;
	    }

	    i__2 = *m;
	    for (ii = 1; ii <= *m; ++ii) {
		temp = X(ii,j);
		X(ii,j) = X(ii,in);
		X(ii,in) = temp;
/* L30: */
	    }

	    K(in) = -K(in);
	    j = in;
	    in = K(in);
	    goto L20;

L40:

/* L50: */
	    ;
	}

    } else {

/*        Backward permutation */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {

	    if (K(i) > 0) {
		goto L80;
	    }

	    K(i) = -K(i);
	    j = K(i);
L60:
	    if (j == i) {
		goto L80;
	    }

	    i__2 = *m;
	    for (ii = 1; ii <= *m; ++ii) {
		temp = X(ii,i);
		X(ii,i) = X(ii,j);
		X(ii,j) = temp;
/* L70: */
	    }

	    K(j) = -K(j);
	    j = K(j);
	    goto L60;

L80:

/* L90: */
	    ;
	}

    }

    return;

/*     End of DLAPMT */

} /* dlapmt_ */

