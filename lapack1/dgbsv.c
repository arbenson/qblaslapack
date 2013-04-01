#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbsv_(int *n, int *kl, int *ku, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbsv(int *n, int *kl, int *ku, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbsv_(int *n, int *kl, int *ku, int *
#endif

	nrhs, LONG DOUBLE *ab, int *ldab, int *ipiv, LONG DOUBLE *b, 
	int *ldb, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGBSV computes the solution to a real system of linear equations   
    A * X = B, where A is a band matrix of order N with KL subdiagonals   
    and KU superdiagonals, and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as A = L * U, where L is a product of permutation   
    and unit lower triangular matrices with KL subdiagonals, and U is   
    upper triangular with KL+KU superdiagonals.  The factored form of A   
    is then used to solve the system of equations A * X = B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the matrix A in band storage, in rows KL+1 to   
            2*KL+KU+1; rows 1 to KL of the array need not be set.   
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(KL+KU+1+i-j,j) = A(i,j) for MAX(1,j-KU)<=i<=MIN(N,j+KL)   
            On exit, details of the factorization: U is stored as an   
            upper triangular band matrix with KL+KU superdiagonals in   
            rows 1 to KL+KU+1, and the multipliers used during the   
            factorization are stored in rows KL+KU+2 to 2*KL+KU+1.   
            See below for further details.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization 
  
                  has been completed, but the factor U is exactly   
                  singular, and the solution has not been computed.   

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
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgbtrf_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtrf(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbtrf_(int *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    xerbla_(char *, int *), dgbtrs_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    xerbla_(char *, int *), qgbtrs(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    xerbla_(char *, int *), qgbtrs_(char *, int *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, int 
	    *, LONG DOUBLE *, int *, int *);


#define IPIV(I) ipiv[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*kl < 0) {
	*info = -2;
    } else if (*ku < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -6;
    } else if (*ldb < MAX(*n,1)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBSV ", &i__1);
	return;
    }

/*     Compute the LU factorization of the band matrix A. */


#ifdef PETSC_PREFIX_SUFFIX
    dgbtrf_(n, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgbtrf(n, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgbtrf_(n, n, kl, ku, &AB(1,1), ldab, &IPIV(1), info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dgbtrs_("No transpose", n, kl, ku, nrhs, &AB(1,1), ldab, &IPIV(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgbtrs("No transpose", n, kl, ku, nrhs, &AB(1,1), ldab, &IPIV(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgbtrs_("No transpose", n, kl, ku, nrhs, &AB(1,1), ldab, &IPIV(
#endif

		1), &B(1,1), ldb, info);
    }
    return;

/*     End of DGBSV */

} /* dgbsv_ */

