#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsymm_(char *side, char *uplo, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsymm(char *side, char *uplo, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsymm_(char *side, char *uplo, int *m, int *n, 
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
    static int i, j, k;
    extern long int lsame_(char *, char *);
    static int nrowa;
    static long int upper;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSYMM  performs one of the matrix-matrix operations   

       C := alpha*A*B + beta*C,   

    or   

       C := alpha*B*A + beta*C,   

    where alpha and beta are scalars,  A is a symmetric matrix and  B and 
  
    C are  m by n matrices.   

    Parameters   
    ==========   

    SIDE   - CHARACTER*1.   
             On entry,  SIDE  specifies whether  the  symmetric matrix  A 
  
             appears on the  left or right  in the  operation as follows: 
  

                SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,   

                SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,   

             Unchanged on exit.   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of  the  symmetric  matrix   A  is  to  be 
  
             referenced as follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of the 
  
                                    symmetric matrix is to be referenced. 
  

                UPLO = 'L' or 'l'   Only the lower triangular part of the 
  
                                    symmetric matrix is to be referenced. 
  

             Unchanged on exit.   

    M      - INTEGER.   
             On entry,  M  specifies the number of rows of the matrix  C. 
  
             M  must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix C. 
  
             N  must be at least zero.   
             Unchanged on exit.   

    ALPHA  - LONG DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
  
             m  when  SIDE = 'L' or 'l'  and is  n otherwise.   
             Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of 
  
             the array  A  must contain the  symmetric matrix,  such that 
  
             when  UPLO = 'U' or 'u', the leading m by m upper triangular 
  
             part of the array  A  must contain the upper triangular part 
  
             of the  symmetric matrix and the  strictly  lower triangular 
  
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l', 
  
             the leading  m by m  lower triangular part  of the  array  A 
  
             must  contain  the  lower triangular part  of the  symmetric 
  
             matrix and the  strictly upper triangular part of  A  is not 
  
             referenced.   
             Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of 
  
             the array  A  must contain the  symmetric matrix,  such that 
  
             when  UPLO = 'U' or 'u', the leading n by n upper triangular 
  
             part of the array  A  must contain the upper triangular part 
  
             of the  symmetric matrix and the  strictly  lower triangular 
  
             part of  A  is not referenced,  and when  UPLO = 'L' or 'l', 
  
             the leading  n by n  lower triangular part  of the  array  A 
  
             must  contain  the  lower triangular part  of the  symmetric 
  
             matrix and the  strictly upper triangular part of  A  is not 
  
             referenced.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
  
             LDA must be at least  MAX( 1, m ), otherwise  LDA must be at 
  
             least  MAX( 1, n ).   
             Unchanged on exit.   

    B      - LONG DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry, the leading  m by n part of the array  B  must 
  
             contain the matrix B.   
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   LDB  must  be  at  least 
  
             MAX( 1, m ).   
             Unchanged on exit.   

    BETA   - LONG DOUBLE PRECISION.   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is 
  
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   

    C      - LONG DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must 
  
             contain the matrix  C,  except when  beta  is zero, in which 
  
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n updated 
  
             matrix.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             MAX( 1, m ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Set NROWA as the number of rows of A.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (lsame_(side, "L")) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    upper = lsame_(uplo, "U");

/*     Test the input parameters. */

    info = 0;
    if (! lsame_(side, "L") && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < MAX(1,nrowa)) {
	info = 7;
    } else if (*ldb < MAX(1,*m)) {
	info = 9;
    } else if (*ldc < MAX(1,*m)) {
	info = 12;
    }
    if (info != 0) {
	xerbla_("DSYMM ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
	return;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = *beta * C(i,j);
/* L30: */
		}
/* L40: */
	    }
	}
	return;
    }

/*     Start the operations. */

    if (lsame_(side, "L")) {

/*        Form  C := alpha*A*B + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp1 = *alpha * B(i,j);
		    temp2 = 0.;
		    i__3 = i - 1;
		    for (k = 1; k <= i-1; ++k) {
			C(k,j) += temp1 * A(k,i);
			temp2 += B(k,j) * A(k,i);
/* L50: */
		    }
		    if (*beta == 0.) {
			C(i,j) = temp1 * A(i,i) + *
				alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + temp1 
				* A(i,i) + *alpha * temp2;
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		for (i = *m; i >= 1; --i) {
		    temp1 = *alpha * B(i,j);
		    temp2 = 0.;
		    i__2 = *m;
		    for (k = i + 1; k <= *m; ++k) {
			C(k,j) += temp1 * A(k,i);
			temp2 += B(k,j) * A(k,i);
/* L80: */
		    }
		    if (*beta == 0.) {
			C(i,j) = temp1 * A(i,i) + *
				alpha * temp2;
		    } else {
			C(i,j) = *beta * C(i,j) + temp1 
				* A(i,i) + *alpha * temp2;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
    } else {

/*        Form  C := alpha*B*A + beta*C. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    temp1 = *alpha * A(j,j);
	    if (*beta == 0.) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = temp1 * B(i,j);
/* L110: */
		}
	    } else {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = *beta * C(i,j) + temp1 * B(i,j);
/* L120: */
		}
	    }
	    i__2 = j - 1;
	    for (k = 1; k <= j-1; ++k) {
		if (upper) {
		    temp1 = *alpha * A(k,j);
		} else {
		    temp1 = *alpha * A(j,k);
		}
		i__3 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) += temp1 * B(i,k);
/* L130: */
		}
/* L140: */
	    }
	    i__2 = *n;
	    for (k = j + 1; k <= *n; ++k) {
		if (upper) {
		    temp1 = *alpha * A(j,k);
		} else {
		    temp1 = *alpha * A(k,j);
		}
		i__3 = *m;
		for (i = 1; i <= *m; ++i) {
		    C(i,j) += temp1 * B(i,k);
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
    }

    return;

/*     End of DSYMM . */

} /* dsymm_ */

