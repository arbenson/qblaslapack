#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbtrf_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbtrf(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbtrf_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif

	ab, int *ldab, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPBTRF computes the Cholesky factorization of a real symmetric   
    positive definite band matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U**T*U or A = L*L**T of the band   
            matrix A, in the same storage format as A.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    Further Details   
    ===============   

    The band storage scheme is illustrated by the following example, when 
  
    N = 6, KD = 2, and UPLO = 'U':   

    On entry:                       On exit:   

        *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46   
        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56   
       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66   

    Similarly, if UPLO = 'L' the format of A is as follows:   

    On entry:                       On exit:   

       a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66   
       a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *   
       a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *   

    Array elements marked * are not used by the routine.   

    Contributed by   
    Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static LONG DOUBLE c_b18 = 1.;
    static LONG DOUBLE c_b21 = -1.;
    static int c__33 = 33;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static LONG DOUBLE work[1056]	/* was [33][32] */;
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
    extern /* Subroutine */ void dtrsm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static int i2, i3;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyrk_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyrk(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyrk_(char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *), dpbtf2_(char *, int *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *), qpbtf2(char *, int *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *), qpbtf2_(char *, int *, int *,
#endif


#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, int *, int *), dpotf2_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, int *, int *), qpotf2(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, int *, int *), qpotf2_(char *, 
#endif

	    int *, LONG DOUBLE *, int *, int *);
    static int ib, nb, ii, jj;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);



#define WORK(I) work[(I)]
#define WAS(I) was[(I)]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*ldab < *kd + 1) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBTRF", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Determine the block size for this environment */

    nb = ilaenv_(&c__1, "DPBTRF", uplo, n, kd, &c_n1, &c_n1, 6L, 1L);

/*     The block size must not exceed the semi-bandwidth KD, and must not 
  
       exceed the limit set by the size of the local array WORK. */

    nb = MIN(nb,32);

    if (nb <= 1 || nb > *kd) {

/*        Use unblocked code */


#ifdef PETSC_PREFIX_SUFFIX
	dpbtf2_(uplo, n, kd, &AB(1,1), ldab, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpbtf2(uplo, n, kd, &AB(1,1), ldab, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpbtf2_(uplo, n, kd, &AB(1,1), ldab, info);
#endif

    } else {

/*        Use blocked code */

	if (lsame_(uplo, "U")) {

/*           Compute the Cholesky factorization of a symmetric ban
d   
             matrix, given the upper triangle of the matrix in ban
d   
             storage.   

             Zero the upper triangle of the work array. */

	    i__1 = nb;
	    for (j = 1; j <= nb; ++j) {
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    WORK(i + j * 33 - 34) = 0.;
/* L10: */
		}
/* L20: */
	    }

/*           Process the band matrix one diagonal block at a time.
 */

	    i__1 = *n;
	    i__2 = nb;
	    for (i = 1; nb < 0 ? i >= *n : i <= *n; i += nb) {
/* Computing MIN */
		i__3 = nb, i__4 = *n - i + 1;
		ib = MIN(i__3,i__4);

/*              Factorize the diagonal block */

		i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dpotf2_(uplo, &ib, &AB(*kd+1,i), &i__3, &ii)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qpotf2(uplo, &ib, &AB(*kd+1,i), &i__3, &ii)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qpotf2_(uplo, &ib, &AB(*kd+1,i), &i__3, &ii)
#endif

			;
		if (ii != 0) {
		    *info = i + ii - 1;
		    goto L150;
		}
		if (i + ib <= *n) {

/*                 Update the relevant part of the trailin
g submatrix.   
                   If A11 denotes the diagonal block which
 has just been   
                   factorized, then we need to update the 
remaining   
                   blocks in the diagram:   

                      A11   A12   A13   
                            A22   A23   
                                  A33   

                   The numbers of rows and columns in the 
partitioning   
                   are IB, I2, I3 respectively. The blocks
 A12, A22 and   
                   A23 are empty if IB = KD. The upper tri
angle of A13   
                   lies outside the band.   

   Computing MIN */
		    i__3 = *kd - ib, i__4 = *n - i - ib + 1;
		    i2 = MIN(i__3,i__4);
/* Computing MIN */
		    i__3 = ib, i__4 = *n - i - *kd + 1;
		    i3 = MIN(i__3,i__4);

		    if (i2 > 0) {

/*                    Update A12 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif

				&i2, &c_b18, &AB(*kd+1,i), &
				i__3, &AB(*kd+1-ib,i+ib), 
				&i__4);

/*                    Update A22 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyrk_("Upper", "Transpose", &i2, &ib, &c_b21, &AB(*kd+1-ib,i+ib), &i__3, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyrk("Upper", "Transpose", &i2, &ib, &c_b21, &AB(*kd+1-ib,i+ib), &i__3, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyrk_("Upper", "Transpose", &i2, &ib, &c_b21, &AB(*kd+1-ib,i+ib), &i__3, &
#endif

				c_b18, &AB(*kd+1,i+ib), &
				i__4);
		    }

		    if (i3 > 0) {

/*                    Copy the lower triangle of A13 i
nto the work array. */

			i__3 = i3;
			for (jj = 1; jj <= i3; ++jj) {
			    i__4 = ib;
			    for (ii = jj; ii <= ib; ++ii) {
				WORK(ii + jj * 33 - 34) = AB(ii-jj+1,jj+i+*kd-1);
/* L30: */
			    }
/* L40: */
			}

/*                    Update A13 (in the work array). 
*/

			i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Left", "Upper", "Transpose", "Non-unit", &ib, 
#endif

				&i3, &c_b18, &AB(*kd+1,i), &
				i__3, work, &c__33);

/*                    Update A23 */

			if (i2 > 0) {
			    i__3 = *ldab - 1;
			    i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dgemm_("Transpose", "No Transpose", &i2, &i3, &ib,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qgemm("Transpose", "No Transpose", &i2, &i3, &ib,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qgemm_("Transpose", "No Transpose", &i2, &i3, &ib,
#endif

				     &c_b21, &AB(*kd+1-ib,i+ib), &i__3, work, &c__33, &c_b18, &
				    AB(ib+1,i+*kd), &i__4);
			}

/*                    Update A33 */

			i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyrk_("Upper", "Transpose", &i3, &ib, &c_b21, work, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyrk("Upper", "Transpose", &i3, &ib, &c_b21, work, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyrk_("Upper", "Transpose", &i3, &ib, &c_b21, work, &
#endif

				c__33, &c_b18, &AB(*kd+1,i+*kd), &i__3);

/*                    Copy the lower triangle of A13 b
ack into place. */

			i__3 = i3;
			for (jj = 1; jj <= i3; ++jj) {
			    i__4 = ib;
			    for (ii = jj; ii <= ib; ++ii) {
				AB(ii-jj+1,jj+i+*kd-1)
					 = WORK(ii + jj * 33 - 34);
/* L50: */
			    }
/* L60: */
			}
		    }
		}
/* L70: */
	    }
	} else {

/*           Compute the Cholesky factorization of a symmetric ban
d   
             matrix, given the lower triangle of the matrix in ban
d   
             storage.   

             Zero the lower triangle of the work array. */

	    i__2 = nb;
	    for (j = 1; j <= nb; ++j) {
		i__1 = nb;
		for (i = j + 1; i <= nb; ++i) {
		    WORK(i + j * 33 - 34) = 0.;
/* L80: */
		}
/* L90: */
	    }

/*           Process the band matrix one diagonal block at a time.
 */

	    i__2 = *n;
	    i__1 = nb;
	    for (i = 1; nb < 0 ? i >= *n : i <= *n; i += nb) {
/* Computing MIN */
		i__3 = nb, i__4 = *n - i + 1;
		ib = MIN(i__3,i__4);

/*              Factorize the diagonal block */

		i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dpotf2_(uplo, &ib, &AB(1,i), &i__3, &ii);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qpotf2(uplo, &ib, &AB(1,i), &i__3, &ii);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qpotf2_(uplo, &ib, &AB(1,i), &i__3, &ii);
#endif

		if (ii != 0) {
		    *info = i + ii - 1;
		    goto L150;
		}
		if (i + ib <= *n) {

/*                 Update the relevant part of the trailin
g submatrix.   
                   If A11 denotes the diagonal block which
 has just been   
                   factorized, then we need to update the 
remaining   
                   blocks in the diagram:   

                      A11   
                      A21   A22   
                      A31   A32   A33   

                   The numbers of rows and columns in the 
partitioning   
                   are IB, I2, I3 respectively. The blocks
 A21, A22 and   
                   A32 are empty if IB = KD. The lower tri
angle of A31   
                   lies outside the band.   

   Computing MIN */
		    i__3 = *kd - ib, i__4 = *n - i - ib + 1;
		    i2 = MIN(i__3,i__4);
/* Computing MIN */
		    i__3 = ib, i__4 = *n - i - *kd + 1;
		    i3 = MIN(i__3,i__4);

		    if (i2 > 0) {

/*                    Update A21 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i2,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Right", "Lower", "Transpose", "Non-unit", &i2,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Right", "Lower", "Transpose", "Non-unit", &i2,
#endif

				 &ib, &c_b18, &AB(1,i), &i__3, &
				AB(ib+1,i), &i__4);

/*                    Update A22 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyrk_("Lower", "No Transpose", &i2, &ib, &c_b21, &AB(ib+1,i), &i__3, &c_b18, &AB(1,i+ib), &i__4);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyrk("Lower", "No Transpose", &i2, &ib, &c_b21, &AB(ib+1,i), &i__3, &c_b18, &AB(1,i+ib), &i__4);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyrk_("Lower", "No Transpose", &i2, &ib, &c_b21, &AB(ib+1,i), &i__3, &c_b18, &AB(1,i+ib), &i__4);
#endif

		    }

		    if (i3 > 0) {

/*                    Copy the upper triangle of A31 i
nto the work array. */

			i__3 = ib;
			for (jj = 1; jj <= ib; ++jj) {
			    i__4 = MIN(jj,i3);
			    for (ii = 1; ii <= MIN(jj,i3); ++ii) {
				WORK(ii + jj * 33 - 34) = AB(*kd+1-jj+ii,jj+i-1);
/* L100: */
			    }
/* L110: */
			}

/*                    Update A31 (in the work array). 
*/

			i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i3,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Right", "Lower", "Transpose", "Non-unit", &i3,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Right", "Lower", "Transpose", "Non-unit", &i3,
#endif

				 &ib, &c_b18, &AB(1,i), &i__3, 
				work, &c__33);

/*                    Update A32 */

			if (i2 > 0) {
			    i__3 = *ldab - 1;
			    i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dgemm_("No transpose", "Transpose", &i3, &i2, &ib,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qgemm("No transpose", "Transpose", &i3, &i2, &ib,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qgemm_("No transpose", "Transpose", &i3, &i2, &ib,
#endif

				     &c_b21, work, &c__33, &AB(ib+1,i), &i__3, &c_b18, &AB(*kd+1-ib,i+ib), &i__4);
			}

/*                    Update A33 */

			i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyrk_("Lower", "No Transpose", &i3, &ib, &c_b21, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyrk("Lower", "No Transpose", &i3, &ib, &c_b21, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyrk_("Lower", "No Transpose", &i3, &ib, &c_b21, 
#endif

				work, &c__33, &c_b18, &AB(1,i+*kd), &i__3);

/*                    Copy the upper triangle of A31 b
ack into place. */

			i__3 = ib;
			for (jj = 1; jj <= ib; ++jj) {
			    i__4 = MIN(jj,i3);
			    for (ii = 1; ii <= MIN(jj,i3); ++ii) {
				AB(*kd+1-jj+ii,jj+i-1)
					 = WORK(ii + jj * 33 - 34);
/* L120: */
			    }
/* L130: */
			}
		    }
		}
/* L140: */
	    }
	}
    }
    return;

L150:
    return;

/*     End of DPBTRF */

} /* dpbtrf_ */

