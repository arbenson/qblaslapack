#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrtri_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrtri(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrtri_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif

	a, int *lda, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTRTRI computes the inverse of a real upper or lower triangular   
    matrix A.   

    This is the Level 3 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of the array A contains 
  
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of the array A contains 
  
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   
            On exit, the (triangular) inverse of the original matrix, in 
  
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, A(i,i) is exactly zero.  The triangular   
                 matrix is singular and its inverse can not be computed. 
  

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static LONG DOUBLE c_b18 = 1.;
    static LONG DOUBLE c_b22 = -1.;
    
    /* System generated locals */
    int  i__1, i__3, i__4, i__5;
    char ch__1[3];
    /* Builtin functions   */
    /* Local variables */
    static int j;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrmm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dtrsm_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qtrsm(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qtrsm_(
#endif

	    char *, char *, char *, char *, int *, int *, LONG DOUBLE *
	    , LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrti2_(char *, char *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrti2(char *, char *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrti2_(char *, char *, int *, LONG DOUBLE 
#endif

	    *, int *, int *);
    static int jb, nb, nn;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static long int nounit;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTRTRI", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (A(*info,*info) == 0.) {
		return;
	    }
/* L10: */
	}
	*info = 0;
    }

/*     Determine the block size for this environment.   

   Writing concatenation */
    /* i__2[0] = 1, a__1[0] = uplo;
    i__2[1] = 1, a__1[1] = diag;
    s_cat(ch__1, a__1, i__2, &c__2, 2L); */
    ch__1[0] = *uplo; ch__1[1] = *diag; ch__1[2] = 0;
    nb = ilaenv_(&c__1, "DTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, 6L, 2L);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */


#ifdef PETSC_PREFIX_SUFFIX
	dtrti2_(uplo, diag, n, &A(1,1), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrti2(uplo, diag, n, &A(1,1), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrti2_(uplo, diag, n, &A(1,1), lda, info);
#endif

    } else {

/*        Use blocked code */

	if (upper) {

/*           Compute inverse of upper triangular matrix */

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; nb < 0 ? j >= *n : j <= *n; j += nb) {
/* Computing MIN */
		i__4 = nb, i__5 = *n - j + 1;
		jb = MIN(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

		i__4 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtrmm_("Left", "Upper", "No transpose", diag, &i__4, &jb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmm("Left", "Upper", "No transpose", diag, &i__4, &jb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmm_("Left", "Upper", "No transpose", diag, &i__4, &jb, &
#endif

			c_b18, &A(1,1), lda, &A(1,j), lda);
		i__4 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtrsm_("Right", "Upper", "No transpose", diag, &i__4, &jb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrsm("Right", "Upper", "No transpose", diag, &i__4, &jb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrsm_("Right", "Upper", "No transpose", diag, &i__4, &jb, &
#endif

			c_b22, &A(j,j), lda, &A(1,j), 
			lda);

/*              Compute inverse of current diagonal block */


#ifdef PETSC_PREFIX_SUFFIX
		dtrti2_("Upper", diag, &jb, &A(j,j), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrti2("Upper", diag, &jb, &A(j,j), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrti2_("Upper", diag, &jb, &A(j,j), lda, info);
#endif

/* L20: */
	    }
	} else {

/*           Compute inverse of lower triangular matrix */

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; -nb < 0 ? j >= 1 : j <= 1; j += -nb) {
/* Computing MIN */
		i__1 = nb, i__4 = *n - j + 1;
		jb = MIN(i__1,i__4);
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block co
lumn */

		    i__1 = *n - j - jb + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmm_("Left", "Lower", "No transpose", diag, &i__1, &jb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmm("Left", "Lower", "No transpose", diag, &i__1, &jb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmm_("Left", "Lower", "No transpose", diag, &i__1, &jb, 
#endif

			    &c_b18, &A(j+jb,j+jb), lda, &A(j+jb,j), lda);
		    i__1 = *n - j - jb + 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrsm_("Right", "Lower", "No transpose", diag, &i__1, &jb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrsm("Right", "Lower", "No transpose", diag, &i__1, &jb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrsm_("Right", "Lower", "No transpose", diag, &i__1, &jb,
#endif

			     &c_b22, &A(j,j), lda, &A(j+jb,j), lda);
		}

/*              Compute inverse of current diagonal block */


#ifdef PETSC_PREFIX_SUFFIX
		dtrti2_("Lower", diag, &jb, &A(j,j), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrti2("Lower", diag, &jb, &A(j,j), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrti2_("Lower", diag, &jb, &A(j,j), lda, info);
#endif

/* L30: */
	    }
	}
    }

    return;

/*     End of DTRTRI */

} /* dtrtri_ */

