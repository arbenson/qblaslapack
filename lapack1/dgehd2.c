#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgehd2_(int *n, int *ilo, int *ihi, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgehd2(int *n, int *ilo, int *ihi, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgehd2_(int *n, int *ilo, int *ihi, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *work, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DGEHD2 reduces a real general matrix A to upper Hessenberg form H by 
  
    an orthogonal similarity transformation:  Q' * A * Q = H .   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows   
            and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally   
            set by a previous call to DGEBAL; otherwise they should be   
            set to 1 and N respectively. See Further Details.   
            1 <= ILO <= IHI <= MAX(1,N).   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the n by n general matrix to be reduced.   
            On exit, the upper triangle and the first subdiagonal of A   
            are overwritten with the upper Hessenberg matrix H, and the   
            elements below the first subdiagonal, with the array TAU,   
            represent the orthogonal matrix Q as a product of elementary 
  
            reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of (ihi-ilo) elementary   
    reflectors   

       Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on   
    exit in A(i+2:ihi,i), and tau in TAU(i).   

    The contents of A are illustrated by the following example, with   
    n = 7, ilo = 2 and ihi = 6:   

    on entry,                        on exit,   

    ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      a   h   h   h   h   a )   
    (     a   a   a   a   a   a )    (      h   h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )   
    (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )   
    (                         a )    (                          a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an 
  
    element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarf_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *), dlarfg_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfg(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfg_(int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *), xerbla_(char *, int *);
    static LONG DOUBLE aii;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > MAX(1,*n)) {
	*info = -2;
    } else if (*ihi < MIN(*ilo,*n) || *ihi > *n) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEHD2", &i__1);
	return;
    }

    i__1 = *ihi - 1;
    for (i = *ilo; i <= *ihi-1; ++i) {

/*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
 */

	i__2 = *ihi - i;
/* Computing MIN */
	i__3 = i + 2;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), 
#endif

		&c__1, &TAU(i));
	aii = A(i+1,i);
	A(i+1,i) = 1.;

/*        Apply H(i) to A(1:ihi,i+1:ihi) from the right */

	i__2 = *ihi - i;

#ifdef PETSC_PREFIX_SUFFIX
	dlarf_("Right", ihi, &i__2, &A(i+1,i), &c__1, &TAU(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarf("Right", ihi, &i__2, &A(i+1,i), &c__1, &TAU(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarf_("Right", ihi, &i__2, &A(i+1,i), &c__1, &TAU(i), &
#endif

		A(1,i+1), lda, &WORK(1));

/*        Apply H(i) to A(i+1:ihi,i+1:n) from the left */

	i__2 = *ihi - i;
	i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
	dlarf_("Left", &i__2, &i__3, &A(i+1,i), &c__1, &TAU(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarf("Left", &i__2, &i__3, &A(i+1,i), &c__1, &TAU(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarf_("Left", &i__2, &i__3, &A(i+1,i), &c__1, &TAU(i), &
#endif

		A(i+1,i+1), lda, &WORK(1));

	A(i+1,i) = aii;
/* L10: */
    }

    return;

/*     End of DGEHD2 */

} /* dgehd2_ */

