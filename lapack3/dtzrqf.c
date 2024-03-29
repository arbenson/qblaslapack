#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtzrqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtzrqf(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtzrqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A   
    to upper triangular form by means of orthogonal transformations.   

    The upper trapezoidal matrix A is factored as   

       A = ( R  0 ) * Z,   

    where Z is an N-by-N orthogonal matrix and R is an M-by-M upper   
    triangular matrix.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= M.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the leading M-by-N upper trapezoidal part of the   
            array A must contain the matrix to be factorized.   
            On exit, the leading M-by-M upper triangular part of A   
            contains the upper triangular matrix R, and elements M+1 to   
            N of the first M rows of A, with the array TAU, represent the 
  
            orthogonal matrix Z as a product of M elementary reflectors. 
  

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (M)   
            The scalar factors of the elementary reflectors.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The factorization is obtained by Householder's method.  The kth   
    transformation matrix, Z( k ), which is used to introduce zeros into 
  
    the ( m - k + 1 )th row of A, is given in the form   

       Z( k ) = ( I     0   ),   
                ( 0  T( k ) )   

    where   

       T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),   
                                                   (   0    )   
                                                   ( z( k ) )   

    tau is a scalar and z( k ) is an ( n - m ) element vector.   
    tau and z( k ) are chosen to annihilate the elements of the kth row   
    of X.   

    The scalar tau is returned in the kth element of TAU and the vector   
    u( k ) in the kth row of A, such that the elements of z( k ) are   
    in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in 
  
    the upper triangular part of A.   

    Z is given by   

       Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b8 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
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
    static int i, k;

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
	    LONG DOUBLE *, LONG DOUBLE *, int *), dcopy_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy_(int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), daxpy_(int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy(int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy_(int 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *)
	    ;
    static int m1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *), xerbla_(char *, int *);



#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < *m) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTZRQF", &i__1);
	return;
    }

/*     Perform the factorization. */

    if (*m == 0) {
	return;
    }
    if (*m == *n) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    TAU(i) = 0.;
/* L10: */
	}
    } else {
/* Computing MIN */
	i__1 = *m + 1;
	m1 = MIN(i__1,*n);
	for (k = *m; k >= 1; --k) {

/*           Use a Householder reflection to zero the kth row of A
.   
             First set up the reflection. */

	    i__1 = *n - *m + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarfg_(&i__1, &A(k,k), &A(k,m1), lda, &TAU(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfg(&i__1, &A(k,k), &A(k,m1), lda, &TAU(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfg_(&i__1, &A(k,k), &A(k,m1), lda, &TAU(
#endif

		    k));

	    if (TAU(k) != 0. && k > 1) {

/*              We now perform the operation  A := A*P( k ). 
  

                Use the first ( k - 1 ) elements of TAU to sto
re  a( k ),   
                where  a( k ) consists of the first ( k - 1 ) 
elements of   
                the  kth column  of  A.  Also  let  B  denote 
 the  first   
                ( k - 1 ) rows of the last ( n - m ) columns o
f A. */

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(1,k), &c__1, &TAU(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(1,k), &c__1, &TAU(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(1,k), &c__1, &TAU(1), &c__1);
#endif


/*              Form   w = a( k ) + B*z( k )  in TAU. */

		i__1 = k - 1;
		i__2 = *n - *m;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i__1, &i__2, &c_b8, &A(1,m1), lda, &A(k,m1), lda, &c_b8, &TAU(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i__1, &i__2, &c_b8, &A(1,m1), lda, &A(k,m1), lda, &c_b8, &TAU(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i__1, &i__2, &c_b8, &A(1,m1), lda, &A(k,m1), lda, &c_b8, &TAU(1), &
#endif

			c__1);

/*              Now form  a( k ) := a( k ) - tau*w   
                and       B      := B      - tau*w*z( k )'. */

		i__1 = k - 1;
		d__1 = -TAU(k);

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__1, &d__1, &TAU(1), &c__1, &A(1,k), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__1, &d__1, &TAU(1), &c__1, &A(1,k), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__1, &d__1, &TAU(1), &c__1, &A(1,k), &
#endif

			c__1);
		i__1 = k - 1;
		i__2 = *n - *m;
		d__1 = -TAU(k);

#ifdef PETSC_PREFIX_SUFFIX
		dger_(&i__1, &i__2, &d__1, &TAU(1), &c__1, &A(k,m1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qger(&i__1, &i__2, &d__1, &TAU(1), &c__1, &A(k,m1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qger_(&i__1, &i__2, &d__1, &TAU(1), &c__1, &A(k,m1)
#endif

			, lda, &A(1,m1), lda);
	    }
/* L20: */
	}
    }

    return;

/*     End of DTZRQF */

} /* dtzrqf_ */

