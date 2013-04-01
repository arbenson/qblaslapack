#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlahrd_(int *n, int *k, int *nb, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlahrd(int *n, int *k, int *nb, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlahrd_(int *n, int *k, int *nb, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *t, int *ldt, 
	LONG DOUBLE *y, int *ldy)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAHRD reduces the first NB columns of a real general n-by-(n-k+1)   
    matrix A so that elements below the k-th subdiagonal are zero. The   
    reduction is performed by an orthogonal similarity transformation   
    Q' * A * Q. The routine returns the matrices V and T which determine 
  
    Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T. 
  

    This is an auxiliary routine called by DGEHRD.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.   

    K       (input) INTEGER   
            The offset for the reduction. Elements below the k-th   
            subdiagonal in the first NB columns are reduced to zero.   

    NB      (input) INTEGER   
            The number of columns to be reduced.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N-K+1) 
  
            On entry, the n-by-(n-k+1) general matrix A.   
            On exit, the elements on and above the k-th subdiagonal in   
            the first NB columns are overwritten with the corresponding   
            elements of the reduced matrix; the elements below the k-th   
            subdiagonal, with the array TAU, represent the matrix Q as a 
  
            product of elementary reflectors. The other columns of A are 
  
            unchanged. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (NB)   
            The scalar factors of the elementary reflectors. See Further 
  
            Details.   

    T       (output) LONG DOUBLE PRECISION array, dimension (NB,NB)   
            The upper triangular matrix T.   

    LDT     (input) INTEGER   
            The leading dimension of the array T.  LDT >= NB.   

    Y       (output) LONG DOUBLE PRECISION array, dimension (LDY,NB)   
            The n-by-nb matrix Y.   

    LDY     (input) INTEGER   
            The leading dimension of the array Y. LDY >= N.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of nb elementary reflectors 
  

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in   
    A(i+k+1:n,i), and tau in TAU(i).   

    The elements of the vectors v together form the (n-k+1)-by-nb matrix 
  
    V which is needed, with T and Y, to apply the transformation to the   
    unreduced part of the matrix, using an update of the form:   
    A := (I - V*T*V') * (A - Y*V').   

    The contents of A on exit are illustrated by the following example   
    with n = 7, k = 3 and nb = 2:   

       ( a   h   a   a   a )   
       ( a   h   a   a   a )   
       ( a   h   a   a   a )   
       ( h   h   a   a   a )   
       ( v1  h   a   a   a )   
       ( v1  v2  a   a   a )   
       ( v1  v2  a   a   a )   

    where a denotes an element of the original matrix A, h denotes a   
    modified element of the upper Hessenberg matrix H, and vi denotes an 
  
    element of the vector defining H(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b4 = -1.;
    static LONG DOUBLE c_b5 = 1.;
    static int c__1 = 1;
    static LONG DOUBLE c_b38 = 0.;
    
    /* System generated locals */
    int i__1, i__2, 
	    i__3;
    LONG DOUBLE d__1;
    /* Local variables */
    static int i;

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
	    int *), dgemv_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qgemv(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qgemv_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dcopy_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qcopy(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qcopy_(int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), daxpy_(int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qaxpy(int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qaxpy_(int *, LONG DOUBLE 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), dtrmv_(char 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrmv(char 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrmv_(char 
#endif

	    *, char *, char *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
	    *, int *);
    static LONG DOUBLE ei;

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



#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define Y(I,J) y[(I)-1 + ((J)-1)* ( *ldy)]

    if (*n <= 1) {
	return;
    }

    i__1 = *nb;
    for (i = 1; i <= *nb; ++i) {
	if (i > 1) {

/*           Update A(1:n,i)   

             Compute i-th column of A - Y * V' */

	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &A(*k+i-1,1), lda, &c_b5, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &A(*k+i-1,1), lda, &c_b5, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &A(*k+i-1,1), lda, &c_b5, &A(1,i), &c__1);
#endif


/*           Apply I - V * T' * V' to this column (call it b) from
 the   
             left, using the last column of T as workspace   

             Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
   
                      ( V2 )             ( b2 )   

             where V1 is unit lower triangular   

             w := V1' * b1 */

	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__2, &A(*k+1,i), &c__1, &T(1,*nb)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__2, &A(*k+1,i), &c__1, &T(1,*nb)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__2, &A(*k+1,i), &c__1, &T(1,*nb)
#endif

		    , &c__1);
	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtrmv_("Lower", "Transpose", "Unit", &i__2, &A(*k+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrmv("Lower", "Transpose", "Unit", &i__2, &A(*k+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrmv_("Lower", "Transpose", "Unit", &i__2, &A(*k+1,1), 
#endif

		    lda, &T(1,*nb), &c__1);

/*           w := w + V2'*b2 */

	    i__2 = *n - *k - i + 1;
	    i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda,
#endif

		     &A(*k+i,i), &c__1, &c_b5, &T(1,*nb), &c__1);

/*           w := T'*w */

	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtrmv_("Upper", "Transpose", "Non-unit", &i__2, &T(1,1), ldt,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrmv("Upper", "Transpose", "Non-unit", &i__2, &T(1,1), ldt,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrmv_("Upper", "Transpose", "Non-unit", &i__2, &T(1,1), ldt,
#endif

		     &T(1,*nb), &c__1);

/*           b2 := b2 - V2*w */

	    i__2 = *n - *k - i + 1;
	    i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("No transpose", &i__2, &i__3, &c_b4, &A(*k+i,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("No transpose", &i__2, &i__3, &c_b4, &A(*k+i,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("No transpose", &i__2, &i__3, &c_b4, &A(*k+i,1), 
#endif

		    lda, &T(1,*nb), &c__1, &c_b5, &A(*k+i,i), &c__1);

/*           b1 := b1 - V1*w */

	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtrmv_("Lower", "No transpose", "Unit", &i__2, &A(*k+1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrmv("Lower", "No transpose", "Unit", &i__2, &A(*k+1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrmv_("Lower", "No transpose", "Unit", &i__2, &A(*k+1,1)
#endif

		    , lda, &T(1,*nb), &c__1);
	    i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    daxpy_(&i__2, &c_b4, &T(1,*nb), &c__1, &A(*k+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qaxpy(&i__2, &c_b4, &T(1,*nb), &c__1, &A(*k+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qaxpy_(&i__2, &c_b4, &T(1,*nb), &c__1, &A(*k+1,i), &c__1);
#endif


	    A(*k+i-1,i-1) = ei;
	}

/*        Generate the elementary reflector H(i) to annihilate   
          A(k+i+1:n,i) */

	i__2 = *n - *k - i + 1;
/* Computing MIN */
	i__3 = *k + i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&i__2, &A(*k+i,i), &A(MIN(*k+i+1,*n),i),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&i__2, &A(*k+i,i), &A(MIN(*k+i+1,*n),i),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&i__2, &A(*k+i,i), &A(MIN(*k+i+1,*n),i),
#endif

		 &c__1, &TAU(i));
	ei = A(*k+i,i);
	A(*k+i,i) = 1.;

/*        Compute  Y(1:n,i) */

	i__2 = *n - *k - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("No transpose", n, &i__2, &c_b5, &A(1,i+1), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("No transpose", n, &i__2, &c_b5, &A(1,i+1), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("No transpose", n, &i__2, &c_b5, &A(1,i+1), lda,
#endif

		 &A(*k+i,i), &c__1, &c_b38, &Y(1,i), &
		c__1);
	i__2 = *n - *k - i + 1;
	i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda, &A(*k+i,i), &c__1, &c_b38, &T(1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda, &A(*k+i,i), &c__1, &c_b38, &T(1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("Transpose", &i__2, &i__3, &c_b5, &A(*k+i,1), lda, &A(*k+i,i), &c__1, &c_b38, &T(1,i), &
#endif

		c__1);
	i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemv_("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &T(1,i), &c__1, &c_b5, &Y(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemv("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &T(1,i), &c__1, &c_b5, &Y(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemv_("No transpose", n, &i__2, &c_b4, &Y(1,1), ldy, &T(1,i), &c__1, &c_b5, &Y(1,i), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dscal_(n, &TAU(i), &Y(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(n, &TAU(i), &Y(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(n, &TAU(i), &Y(1,i), &c__1);
#endif


/*        Compute T(1:i,i) */

	i__2 = i - 1;
	d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__2, &d__1, &T(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__2, &d__1, &T(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__2, &d__1, &T(1,i), &c__1);
#endif

	i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrmv("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrmv_("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, 
#endif

		&T(1,i), &c__1);
	T(i,i) = TAU(i);

/* L10: */
    }
    A(*k+*nb,*nb) = ei;

    return;

/*     End of DLAHRD */

} /* dlahrd_ */

