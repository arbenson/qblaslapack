#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsyr2k_(char *uplo, char *trans, int *n, int *k, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsyr2k(char *uplo, char *trans, int *n, int *k, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsyr2k_(char *uplo, char *trans, int *n, int *k, 
#endif

	LONG DOUBLE *alpha, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
	int *ldb, LONG DOUBLE *beta, LONG DOUBLE *c, int *ldc)
{


    /* System generated locals */
    int i__1, i__2, 
	    i__3;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp1, temp2;
    static int i, j, l;
    extern long int lsame_(char *, char *);
    static int nrowa;
    static long int upper;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSYR2K  performs one of the symmetric rank 2k operations   

       C := alpha*A*B' + alpha*B*A' + beta*C,   

    or   

       C := alpha*A'*B + alpha*B'*A + beta*C,   

    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix 
  
    and  A and B  are  n by k  matrices  in the  first  case  and  k by n 
  
    matrices in the second case.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of the  array  C  is to be  referenced  as 
  
             follows:   

                UPLO = 'U' or 'u'   Only the  upper triangular part of  C 
  
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the  lower triangular part of  C 
  
                                    is to be referenced.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as 
  
             follows:   

                TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +   
                                          beta*C.   

                TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   

                TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +   
                                          beta*C.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number 
  
             of  columns  of the  matrices  A and B,  and on  entry  with 
  
             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number 
  
             of rows of the matrices  A and B.  K must be at least  zero. 
  
             Unchanged on exit.   

    ALPHA  - LONG DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  A  must contain the matrix  A,  otherwise 
  
             the leading  k by n  part of the array  A  must contain  the 
  
             matrix A.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDA must be at least  MAX( 1, n ), otherwise  LDA must 
  
             be at least  MAX( 1, k ).   
             Unchanged on exit.   

    B      - LONG DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  B  must contain the matrix  B,  otherwise 
  
             the leading  k by n  part of the array  B  must contain  the 
  
             matrix B.   
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDB must be at least  MAX( 1, n ), otherwise  LDB must 
  
             be at least  MAX( 1, k ).   
             Unchanged on exit.   

    BETA   - LONG DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   

    C      - LONG DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n 
  
             upper triangular part of the array C must contain the upper 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             lower triangular part of C is not referenced.  On exit, the 
  
             upper triangular part of the array  C is overwritten by the 
  
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n 
  
             lower triangular part of the array C must contain the lower 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             upper triangular part of C is not referenced.  On exit, the 
  
             lower triangular part of the array  C is overwritten by the 
  
             lower triangular part of the updated matrix.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             MAX( 1, n ).   
             Unchanged on exit.   


    Level 3 Blas routine.   


    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame_(uplo, "U");

    info = 0;
    if (! upper && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < MAX(1,nrowa)) {
	info = 7;
    } else if (*ldb < MAX(1,nrowa)) {
	info = 9;
    } else if (*ldc < MAX(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("DSYR2K", &info);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0 || ((*alpha == 0. || *k == 0) && *beta == 1.)) {
	return;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return;
    }

/*     Start the operations. */

    if (lsame_(trans, "N")) {

/*        Form  C := alpha*A*B' + alpha*B*A' + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0. || B(j,l) != 0.) {
			temp1 = *alpha * B(j,l);
			temp2 = *alpha * A(j,l);
			i__3 = j;
			for (i = 1; i <= j; ++i) {
			    C(i,j) = C(i,j) + A(i,l) * temp1 + B(i,l) * 
				    temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0. || B(j,l) != 0.) {
			temp1 = *alpha * B(j,l);
			temp2 = *alpha * A(j,l);
			i__3 = *n;
			for (i = j; i <= *n; ++i) {
			    C(i,j) = C(i,j) + A(i,l) * temp1 + B(i,l) * 
				    temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*B + alpha*B'*A + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = j;
		for (i = 1; i <= j; ++i) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			temp1 += A(l,i) * B(l,j);
			temp2 += B(l,i) * A(l,j);
/* L190: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + *
				alpha * temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *n;
		for (i = j; i <= *n; ++i) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			temp1 += A(l,i) * B(l,j);
			temp2 += B(l,i) * A(l,j);
/* L220: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp1 + *alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + *
				alpha * temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return;

/*     End of DSYR2K. */

} /* dsyr2k_ */

