#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgeqpf_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgeqpf(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgeqpf_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *jpvt, LONG DOUBLE *tau, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK test routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGEQPF computes a QR factorization with column pivoting of a   
    real M-by-N matrix A: A*P = Q*R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A. N >= 0   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the upper triangle of the array contains the   
            MIN(M,N)-by-N upper triangular matrix R; the elements   
            below the diagonal, together with the array TAU,   
            represent the orthogonal matrix Q as a product of   
            MIN(m,n) elementary reflectors.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    JPVT    (input/output) INTEGER array, dimension (N)   
            On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted 
  
            to the front of A*P (a leading column); if JPVT(i) = 0,   
            the i-th column of A is a free column.   
            On exit, if JPVT(i) = k, then the i-th column of A*P   
            was the k-th column of A.   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(n)   

    Each H(i) has the form   

       H = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i). 
  

    The matrix P is represented in jpvt as follows: If   
       jpvt(j) = i   
    then the jth column of P is the ith canonical unit vector.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dnrm2_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qnrm2_(int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE temp2;
    static int i, j;

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
	    LONG DOUBLE *);
    static int itemp;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dswap_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dgeqr2_(int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qgeqr2(int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qgeqr2_(int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dorm2r_(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorm2r(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorm2r_(char *, char *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static int ma, mn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE aii;
    static int pvt;



#define JPVT(I) jpvt[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQPF", &i__1);
	return;
    }

    mn = MIN(*m,*n);

/*     Move initial columns up front */

    itemp = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (JPVT(i) != 0) {
	    if (i != itemp) {

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &A(1,i), &c__1, &A(1,itemp), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &A(1,i), &c__1, &A(1,itemp), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &A(1,i), &c__1, &A(1,itemp), &
#endif

			c__1);
		JPVT(i) = JPVT(itemp);
		JPVT(itemp) = i;
	    } else {
		JPVT(i) = i;
	    }
	    ++itemp;
	} else {
	    JPVT(i) = i;
	}
/* L10: */
    }
    --itemp;

/*     Compute the QR factorization and update remaining columns */

    if (itemp > 0) {
	ma = MIN(itemp,*m);

#ifdef PETSC_PREFIX_SUFFIX
	dgeqr2_(m, &ma, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgeqr2(m, &ma, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgeqr2_(m, &ma, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif

	if (ma < *n) {
	    i__1 = *n - ma;

#ifdef PETSC_PREFIX_SUFFIX
	    dorm2r_("Left", "Transpose", m, &i__1, &ma, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorm2r("Left", "Transpose", m, &i__1, &ma, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorm2r_("Left", "Transpose", m, &i__1, &ma, &A(1,1), lda, &
#endif

		    TAU(1), &A(1,ma+1), lda, &WORK(1), info);
	}
    }

    if (itemp < mn) {

/*        Initialize partial column norms. The first n elements of   
          work store the exact column norms. */

	i__1 = *n;
	for (i = itemp + 1; i <= *n; ++i) {
	    i__2 = *m - itemp;

#ifdef PETSC_PREFIX_SUFFIX
	    WORK(i) = dnrm2_(&i__2, &A(itemp+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    WORK(i) = qnrm2(&i__2, &A(itemp+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    WORK(i) = qnrm2_(&i__2, &A(itemp+1,i), &c__1);
#endif

	    WORK(*n + i) = WORK(i);
/* L20: */
	}

/*        Compute factorization */

	i__1 = mn;
	for (i = itemp + 1; i <= mn; ++i) {

/*           Determine ith pivot column and swap if necessary */

	    i__2 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    pvt = i - 1 + idamax_(&i__2, &WORK(i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    pvt = i - 1 + iqamax(&i__2, &WORK(i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    pvt = i - 1 + iqamax_(&i__2, &WORK(i), &c__1);
#endif


	    if (pvt != i) {

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(m, &A(1,pvt), &c__1, &A(1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(m, &A(1,pvt), &c__1, &A(1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(m, &A(1,pvt), &c__1, &A(1,i), &
#endif

			c__1);
		itemp = JPVT(pvt);
		JPVT(pvt) = JPVT(i);
		JPVT(i) = itemp;
		WORK(pvt) = WORK(i);
		WORK(*n + pvt) = WORK(*n + i);
	    }

/*           Generate elementary reflector H(i) */

	    if (i < *m) {
		i__2 = *m - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfg_(&i__2, &A(i,i), &A(i+1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfg(&i__2, &A(i,i), &A(i+1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfg_(&i__2, &A(i,i), &A(i+1,i), &
#endif

			c__1, &TAU(i));
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dlarfg_(&c__1, &A(*m,*m), &A(*m,*m), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfg(&c__1, &A(*m,*m), &A(*m,*m), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfg_(&c__1, &A(*m,*m), &A(*m,*m), &
#endif

			c__1, &TAU(*m));
	    }

	    if (i < *n) {

/*              Apply H(i) to A(i:m,i+1:n) from the left */

		aii = A(i,i);
		A(i,i) = 1.;
		i__2 = *m - i + 1;
		i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dlarf_("LEFT", &i__2, &i__3, &A(i,i), &c__1, &TAU(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarf("LEFT", &i__2, &i__3, &A(i,i), &c__1, &TAU(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarf_("LEFT", &i__2, &i__3, &A(i,i), &c__1, &TAU(
#endif

			i), &A(i,i+1), lda, &WORK((*n << 1) + 
			1));
		A(i,i) = aii;
	    }

/*           Update partial column norms */

	    i__2 = *n;
	    for (j = i + 1; j <= *n; ++j) {
		if (WORK(j) != 0.) {
/* Computing 2nd power */
		    d__2 = (d__1 = A(i,j), ABS(d__1)) / WORK(j);
		    temp = 1. - d__2 * d__2;
		    temp = MAX(temp,0.);
/* Computing 2nd power */
		    d__1 = WORK(j) / WORK(*n + j);
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
		    if (temp2 == 1.) {
			if (*m - i > 0) {
			    i__3 = *m - i;

#ifdef PETSC_PREFIX_SUFFIX
			    WORK(j) = dnrm2_(&i__3, &A(i+1,j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    WORK(j) = qnrm2(&i__3, &A(i+1,j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    WORK(j) = qnrm2_(&i__3, &A(i+1,j), &
#endif

				    c__1);
			    WORK(*n + j) = WORK(j);
			} else {
			    WORK(j) = 0.;
			    WORK(*n + j) = 0.;
			}
		    } else {
			WORK(j) *= sqrt(temp);
		    }
		}
/* L30: */
	    }

/* L40: */
	}
    }
    return;

/*     End of DGEQPF */

} /* dgeqpf_ */

