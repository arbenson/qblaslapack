#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbtrf_(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbtrf(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbtrf_(int *m, int *n, int *kl, int *ku,
#endif

	 LONG DOUBLE *ab, int *ldab, int *ipiv, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGBTRF computes an LU factorization of a real m-by-n band matrix A   
    using partial pivoting with row interchanges.   

    This is the blocked version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the matrix A in band storage, in rows KL+1 to   
            2*KL+KU+1; rows 1 to KL of the array need not be set.   
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(kl+ku+1+i-j,j) = A(i,j) for MAX(1,j-ku)<=i<=MIN(m,j+kl)   

            On exit, details of the factorization: U is stored as an   
            upper triangular band matrix with KL+KU superdiagonals in   
            rows 1 to KL+KU+1, and the multipliers used during the   
            factorization are stored in rows KL+KU+2 to 2*KL+KU+1.   
            See below for further details.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (output) INTEGER array, dimension (MIN(M,N))   
            The pivot indices; for 1 <= i <= MIN(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = +i, U(i,i) is exactly zero. The factorization 
  
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used 
  
                 to solve a system of equations.   

    Further Details   
    ===============   

    The band storage scheme is illustrated by the following example, when 
  
    M = N = 6, KL = 2, KU = 1:   

    On entry:                       On exit:   

        *    *    *    +    +    +       *    *    *   u14  u25  u36   
        *    *    +    +    +    +       *    *   u13  u24  u35  u46   
        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56   
       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66   
       a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *   
       a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *   

    Array elements marked * are not used by the routine; elements marked 
  
    + need not be set on entry, but are required by the routine to store 
  
    elements of U because of fill-in resulting from the row interchanges. 
  

    ===================================================================== 
  


       KV is the number of superdiagonals in the factor U, allowing for   
       fill-in   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c__65 = 65;
    static LONG DOUBLE c_b18 = -1.;
    static LONG DOUBLE c_b31 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5, i__6;
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
    static LONG DOUBLE temp;
    static int i, j;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *), dgemm_(char *, char *, int *, int *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qgemm(char *, char *, int *, int *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qgemm_(char *, char *, int *, int *, int *
#endif

	    , LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dcopy_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy_(
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dswap_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qswap(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qswap_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *
#endif

	    );
    static LONG DOUBLE work13[4160]	/* was [65][64] */, work31[4160]	
	    /* was [65][64] */;

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
    static int i2, i3, j2, j3, k2;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgbtf2_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtf2(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtf2_(int *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, int *, int *);
    static int jb, nb, ii, jj, jm, ip, jp, km, ju, kv;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    static int nw;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaswp_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaswp(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaswp_(int *, LONG DOUBLE *, int *, 
#endif

	    int *, int *, int *, int *);



#define WORK13(I) work13[(I)]
#define WAS(I) was[(I)]
#define IPIV(I) ipiv[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    kv = *ku + *kl;

/*     Test the input parameters. */

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < *kl + kv + 1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBTRF", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return;
    }

/*     Determine the block size for this environment */

    nb = ilaenv_(&c__1, "DGBTRF", " ", m, n, kl, ku, 6L, 1L);

/*     The block size must not exceed the limit set by the size of the   
       local arrays WORK13 and WORK31. */

    nb = MIN(nb,64);

    if (nb <= 1 || nb > *kl) {

/*        Use unblocked code */


#ifdef PETSC_PREFIX_SUFFIX
	dgbtf2_(m, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgbtf2(m, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgbtf2_(m, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif

    } else {

/*        Use blocked code   

          Zero the superdiagonal elements of the work array WORK13 */

	i__1 = nb;
	for (j = 1; j <= nb; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		WORK13(i + j * 65 - 66) = 0.;
/* L10: */
	    }
/* L20: */
	}

/*        Zero the subdiagonal elements of the work array WORK31 */

	i__1 = nb;
	for (j = 1; j <= nb; ++j) {
	    i__2 = nb;
	    for (i = j + 1; i <= nb; ++i) {
		work31[i + j * 65 - 66] = 0.;
/* L30: */
	    }
/* L40: */
	}

/*        Gaussian elimination with partial pivoting   

          Set fill-in elements in columns KU+2 to KV to zero */

	i__1 = MIN(kv,*n);
	for (j = *ku + 2; j <= MIN(kv,*n); ++j) {
	    i__2 = *kl;
	    for (i = kv - j + 2; i <= *kl; ++i) {
		AB(i,j) = 0.;
/* L50: */
	    }
/* L60: */
	}

/*        JU is the index of the last column affected by the current 
  
          stage of the factorization */

	ju = 1;

	i__1 = MIN(*m,*n);
	i__2 = nb;
	for (j = 1; nb < 0 ? j >= MIN(*m,*n) : j <= MIN(*m,*n); j += nb) {
/* Computing MIN */
	    i__3 = nb, i__4 = MIN(*m,*n) - j + 1;
	    jb = MIN(i__3,i__4);

/*           The active part of the matrix is partitioned   

                A11   A12   A13   
                A21   A22   A23   
                A31   A32   A33   

             Here A11, A21 and A31 denote the current block of JB 
columns   
             which is about to be factorized. The number of rows i
n the   
             partitioning are JB, I2, I3 respectively, and the num
bers   
             of columns are JB, J2, J3. The superdiagonal elements
 of A13   
             and the subdiagonal elements of A31 lie outside the b
and.   

   Computing MIN */
	    i__3 = *kl - jb, i__4 = *m - j - jb + 1;
	    i2 = MIN(i__3,i__4);
/* Computing MIN */
	    i__3 = jb, i__4 = *m - j - *kl + 1;
	    i3 = MIN(i__3,i__4);

/*           J2 and J3 are computed after JU has been updated.   

             Factorize the current block of JB columns */

	    i__3 = j + jb - 1;
	    for (jj = j; jj <= j+jb-1; ++jj) {

/*              Set fill-in elements in column JJ+KV to zero 
*/

		if (jj + kv <= *n) {
		    i__4 = *kl;
		    for (i = 1; i <= *kl; ++i) {
			AB(i,jj+kv) = 0.;
/* L70: */
		    }
		}

/*              Find pivot and test for singularity. KM is the
 number of   
                subdiagonal elements in the current column.   

   Computing MIN */
		i__4 = *kl, i__5 = *m - jj;
		km = MIN(i__4,i__5);
		i__4 = km + 1;

#ifdef PETSC_PREFIX_SUFFIX
		jp = idamax_(&i__4, &AB(kv+1,jj), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		jp = iqamax(&i__4, &AB(kv+1,jj), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		jp = iqamax_(&i__4, &AB(kv+1,jj), &c__1);
#endif

		IPIV(jj) = jp + jj - j;
		if (AB(kv+jp,jj) != 0.) {
/* Computing MAX   
   Computing MIN */
		    i__6 = jj + *ku + jp - 1;
		    i__4 = ju, i__5 = MIN(i__6,*n);
		    ju = MAX(i__4,i__5);
		    if (jp != 1) {

/*                    Apply interchange to columns J t
o J+JB-1 */

			if (jp + jj - 1 < j + *kl) {

			    i__4 = *ldab - 1;
			    i__5 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dswap_(&jb, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qswap(&jb, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qswap_(&jb, &AB(kv+1+jj-j,j), &
#endif

				    i__4, &AB(kv+jp+jj-j,j),
				     &i__5);
			} else {

/*                       The interchange affects c
olumns J to JJ-1 of A31   
                         which are stored in the w
ork array WORK31 */

			    i__4 = jj - j;
			    i__5 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dswap_(&i__4, &AB(kv+1+jj-j,j), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qswap(&i__4, &AB(kv+1+jj-j,j), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qswap_(&i__4, &AB(kv+1+jj-j,j), 
#endif

				    &i__5, &work31[jp + jj - j - *kl - 1], &
				    c__65);
			    i__4 = j + jb - jj;
			    i__5 = *ldab - 1;
			    i__6 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    dswap_(&i__4, &AB(kv+1,jj), &i__5, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qswap(&i__4, &AB(kv+1,jj), &i__5, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qswap_(&i__4, &AB(kv+1,jj), &i__5, &
#endif

				    AB(kv+jp,jj), &i__6);
			}
		    }

/*                 Compute multipliers */

		    d__1 = 1. / AB(kv+1,jj);

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(&km, &d__1, &AB(kv+2,jj), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(&km, &d__1, &AB(kv+2,jj), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(&km, &d__1, &AB(kv+2,jj), &c__1);
#endif


/*                 Update trailing submatrix within the ba
nd and within   
                   the current block. JM is the index of t
he last column   
                   which needs to be updated.   

   Computing MIN */
		    i__4 = ju, i__5 = j + jb - 1;
		    jm = MIN(i__4,i__5);
		    if (jm > jj) {
			i__4 = jm - jj;
			i__5 = *ldab - 1;
			i__6 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dger_(&km, &i__4, &c_b18, &AB(kv+2,jj), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qger(&km, &i__4, &c_b18, &AB(kv+2,jj), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qger_(&km, &i__4, &c_b18, &AB(kv+2,jj), 
#endif

				&c__1, &AB(kv,jj+1), &i__5, &
				AB(kv+1,jj+1), &i__6);
		    }
		} else {

/*                 If pivot is zero, set INFO to the index
 of the pivot   
                   unless a zero pivot has already been fo
und. */

		    if (*info == 0) {
			*info = jj;
		    }
		}

/*              Copy current column of A31 into the work array
 WORK31   

   Computing MIN */
		i__4 = jj - j + 1;
		nw = MIN(i__4,i3);
		if (nw > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(&nw, &AB(kv+*kl+1-jj+j,jj), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(&nw, &AB(kv+*kl+1-jj+j,jj), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(&nw, &AB(kv+*kl+1-jj+j,jj), &
#endif

			    c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
		}
/* L80: */
	    }
	    if (j + jb <= *n) {

/*              Apply the row interchanges to the other blocks
.   

   Computing MIN */
		i__3 = ju - j + 1;
		j2 = MIN(i__3,kv) - jb;
/* Computing MAX */
		i__3 = 0, i__4 = ju - j - kv + 1;
		j3 = MAX(i__3,i__4);

/*              Use DLASWP to apply the row interchanges to A1
2, A22, and   
                A32. */

		i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlaswp_(&j2, &AB(kv+1-jb,j+jb), &i__3, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlaswp(&j2, &AB(kv+1-jb,j+jb), &i__3, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlaswp_(&j2, &AB(kv+1-jb,j+jb), &i__3, &
#endif

			c__1, &jb, &IPIV(j), &c__1);

/*              Adjust the pivot indices. */

		i__3 = j + jb - 1;
		for (i = j; i <= j+jb-1; ++i) {
		    IPIV(i) = IPIV(i) + j - 1;
/* L90: */
		}

/*              Apply the row interchanges to A13, A23, and A3
3   
                columnwise. */

		k2 = j - 1 + jb + j2;
		i__3 = j3;
		for (i = 1; i <= j3; ++i) {
		    jj = k2 + i;
		    i__4 = j + jb - 1;
		    for (ii = j + i - 1; ii <= j+jb-1; ++ii) {
			ip = IPIV(ii);
			if (ip != ii) {
			    temp = AB(kv+1+ii-jj,jj);
			    AB(kv+1+ii-jj,jj) = AB(kv+1+ip-jj,jj);
			    AB(kv+1+ip-jj,jj) = temp;
			}
/* L100: */
		    }
/* L110: */
		}

/*              Update the relevant part of the trailing subma
trix */

		if (j2 > 0) {

/*                 Update A12 */

		    i__3 = *ldab - 1;
		    i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrsm("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2, 
#endif

			    &c_b31, &AB(kv+1,j), &i__3, &AB(kv+1-jb,j+jb), &i__4);

		    if (i2 > 0) {

/*                    Update A22 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;
			i__5 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("No transpose", "No transpose", &i2, &j2, &jb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("No transpose", "No transpose", &i2, &j2, &jb, 
#endif

				&c_b18, &AB(kv+1+jb,j), &i__3,
				 &AB(kv+1-jb,j+jb), &i__4,
				 &c_b31, &AB(kv+1,j+jb), &
				i__5);
		    }

		    if (i3 > 0) {

/*                    Update A32 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("No transpose", "No transpose", &i3, &j2, &jb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("No transpose", "No transpose", &i3, &j2, &jb, 
#endif

				&c_b18, work31, &c__65, &AB(kv+1-jb,j+jb), &i__3, &c_b31, &AB(kv+*kl+1-jb,j+jb), &i__4);
		    }
		}

		if (j3 > 0) {

/*                 Copy the lower triangle of A13 into the
 work array   
                   WORK13 */

		    i__3 = j3;
		    for (jj = 1; jj <= j3; ++jj) {
			i__4 = jb;
			for (ii = jj; ii <= jb; ++ii) {
			    WORK13(ii + jj * 65 - 66) = AB(ii-jj+1,jj+j+kv-1);
/* L120: */
			}
/* L130: */
		    }

/*                 Update A13 in the work array */

		    i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrsm("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3, 
#endif

			    &c_b31, &AB(kv+1,j), &i__3, work13, 
			    &c__65);

		    if (i2 > 0) {

/*                    Update A23 */

			i__3 = *ldab - 1;
			i__4 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("No transpose", "No transpose", &i2, &j3, &jb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("No transpose", "No transpose", &i2, &j3, &jb, 
#endif

				&c_b18, &AB(kv+1+jb,j), &i__3,
				 work13, &c__65, &c_b31, &AB(jb+1,j+kv), &i__4);
		    }

		    if (i3 > 0) {

/*                    Update A33 */

			i__3 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qgemm("No transpose", "No transpose", &i3, &j3, &jb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qgemm_("No transpose", "No transpose", &i3, &j3, &jb, 
#endif

				&c_b18, work31, &c__65, work13, &c__65, &
				c_b31, &AB(*kl+1,j+kv), &
				i__3);
		    }

/*                 Copy the lower triangle of A13 back int
o place */

		    i__3 = j3;
		    for (jj = 1; jj <= j3; ++jj) {
			i__4 = jb;
			for (ii = jj; ii <= jb; ++ii) {
			    AB(ii-jj+1,jj+j+kv-1) = 
				    WORK13(ii + jj * 65 - 66);
/* L140: */
			}
/* L150: */
		    }
		}
	    } else {

/*              Adjust the pivot indices. */

		i__3 = j + jb - 1;
		for (i = j; i <= j+jb-1; ++i) {
		    IPIV(i) = IPIV(i) + j - 1;
/* L160: */
		}
	    }

/*           Partially undo the interchanges in the current block 
to   
             restore the upper triangular form of A31 and copy the
 upper   
             triangle of A31 back into place */

	    i__3 = j;
	    for (jj = j + jb - 1; jj >= j; --jj) {
		jp = IPIV(jj) - jj + 1;
		if (jp != 1) {

/*                 Apply interchange to columns J to JJ-1 
*/

		    if (jp + jj - 1 < j + *kl) {

/*                    The interchange does not affect 
A31 */

			i__4 = jj - j;
			i__5 = *ldab - 1;
			i__6 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dswap_(&i__4, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qswap(&i__4, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qswap_(&i__4, &AB(kv+1+jj-j,j), &
#endif

				i__5, &AB(kv+jp+jj-j,j), &
				i__6);
		    } else {

/*                    The interchange does affect A31 
*/

			i__4 = jj - j;
			i__5 = *ldab - 1;

#ifdef PETSC_PREFIX_SUFFIX
			dswap_(&i__4, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qswap(&i__4, &AB(kv+1+jj-j,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qswap_(&i__4, &AB(kv+1+jj-j,j), &
#endif

				i__5, &work31[jp + jj - j - *kl - 1], &c__65);
		    }
		}

/*              Copy the current column of A31 back into place
   

   Computing MIN */
		i__4 = i3, i__5 = jj - j + 1;
		nw = MIN(i__4,i__5);
		if (nw > 0) {

#ifdef PETSC_PREFIX_SUFFIX
		    dcopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &AB(kv+*kl+1-jj+j,jj), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qcopy(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &AB(kv+*kl+1-jj+j,jj), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qcopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &AB(kv+*kl+1-jj+j,jj), &c__1);
#endif

		}
/* L170: */
	    }
/* L180: */
	}
    }

    return;

/*     End of DGBTRF */

} /* dgbtrf_ */

