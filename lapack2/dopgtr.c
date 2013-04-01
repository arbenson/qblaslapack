#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dopgtr_(char *uplo, int *n, LONG DOUBLE *ap, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qopgtr(char *uplo, int *n, LONG DOUBLE *ap, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qopgtr_(char *uplo, int *n, LONG DOUBLE *ap, 
#endif

	LONG DOUBLE *tau, LONG DOUBLE *q, int *ldq, LONG DOUBLE *work, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DOPGTR generates a real orthogonal matrix Q which is defined as the   
    product of n-1 elementary reflectors H(i) of order n, as returned by 
  
    DSPTRD using packed storage:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangular packed storage used in previous   
                   call to DSPTRD;   
            = 'L': Lower triangular packed storage used in previous   
                   call to DSPTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The vectors which define the elementary reflectors, as   
            returned by DSPTRD.   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSPTRD.   

    Q       (output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            The N-by-N orthogonal matrix Q.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N-1)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i, j;
    extern long int lsame_(char *, char *);
    static int iinfo;
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorg2l_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorg2l(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorg2l_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dorg2r_(int *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorg2r(int *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorg2r_(int *, int *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static int ij;
    extern /* Subroutine */ void xerbla_(char *, int *);


#define AP(I) ap[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldq < MAX(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DOPGTR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (upper) {

/*        Q was determined by a call to DSPTRD with UPLO = 'U'   

          Unpack the vectors which define the elementary reflectors an
d   
          set the last row and column of Q equal to those of the unit 
  
          matrix */

	ij = 2;
	i__1 = *n - 1;
	for (j = 1; j <= *n-1; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		Q(i,j) = AP(ij);
		++ij;
/* L10: */
	    }
	    ij += 2;
	    Q(*n,j) = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    Q(i,*n) = 0.;
/* L30: */
	}
	Q(*n,*n) = 1.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorg2l_(&i__1, &i__2, &i__3, &Q(1,1), ldq, &TAU(1), &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorg2l(&i__1, &i__2, &i__3, &Q(1,1), ldq, &TAU(1), &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorg2l_(&i__1, &i__2, &i__3, &Q(1,1), ldq, &TAU(1), &WORK(1), &
#endif

		iinfo);

    } else {

/*        Q was determined by a call to DSPTRD with UPLO = 'L'.   

          Unpack the vectors which define the elementary reflectors an
d   
          set the first row and column of Q equal to those of the unit
   
          matrix */

	Q(1,1) = 1.;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    Q(i,1) = 0.;
/* L40: */
	}
	ij = 3;
	i__1 = *n;
	for (j = 2; j <= *n; ++j) {
	    Q(1,j) = 0.;
	    i__2 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		Q(i,j) = AP(ij);
		++ij;
/* L50: */
	    }
	    ij += 2;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorg2r_(&i__1, &i__2, &i__3, &Q(2,2), ldq, &TAU(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorg2r(&i__1, &i__2, &i__3, &Q(2,2), ldq, &TAU(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorg2r_(&i__1, &i__2, &i__3, &Q(2,2), ldq, &TAU(1), 
#endif

		    &WORK(1), &iinfo);
	}
    }
    return;

/*     End of DOPGTR */

} /* dopgtr_ */

