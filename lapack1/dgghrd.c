#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgghrd_(char *compq, char *compz, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgghrd(char *compq, char *compz, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgghrd_(char *compq, char *compz, int *n, int *
#endif

	ilo, int *ihi, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
	int *ldb, LONG DOUBLE *q, int *ldq, LONG DOUBLE *z, int *
	ldz, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGHRD reduces a pair of real matrices (A,B) to generalized upper   
    Hessenberg form using orthogonal transformations, where A is a   
    general matrix and B is upper triangular:  Q' * A * Z = H and   
    Q' * B * Z = T, where H is upper Hessenberg, T is upper triangular,   
    and Q and Z are orthogonal, and ' means transpose.   

    The orthogonal matrices Q and Z are determined as products of Givens 
  
    rotations.  They may either be formed explicitly, or they may be   
    postmultiplied into input matrices Q1 and Z1, so that   

         Q1 * A * Z1' = (Q1*Q) * H * (Z1*Z)'   
         Q1 * B * Z1' = (Q1*Q) * T * (Z1*Z)'   

    Arguments   
    =========   

    COMPQ   (input) CHARACTER*1   
            = 'N': do not compute Q;   
            = 'I': Q is initialized to the unit matrix, and the   
                   orthogonal matrix Q is returned;   
            = 'V': Q must contain an orthogonal matrix Q1 on entry,   
                   and the product Q1*Q is returned.   

    COMPZ   (input) CHARACTER*1   
            = 'N': do not compute Z;   
            = 'I': Z is initialized to the unit matrix, and the   
                   orthogonal matrix Z is returned;   
            = 'V': Z must contain an orthogonal matrix Z1 on entry,   
                   and the product Z1*Z is returned.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows and 
  
            columns 1:ILO-1 and IHI+1:N.  ILO and IHI are normally set   
            by a previous call to DGGBAL; otherwise they should be set   
            to 1 and N respectively.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the N-by-N general matrix to be reduced.   
            On exit, the upper triangle and the first subdiagonal of A   
            are overwritten with the upper Hessenberg matrix H, and the   
            rest is set to zero.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.   
            On exit, the upper triangular matrix T = Q' B Z.  The   
            elements below the diagonal are set to zero.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ, N)   
            If COMPQ='N':  Q is not referenced.   
            If COMPQ='I':  on entry, Q need not be set, and on exit it   
                           contains the orthogonal matrix Q, where Q'   
                           is the product of the Givens transformations   
                           which are applied to A and B on the left.   
            If COMPQ='V':  on entry, Q must contain an orthogonal matrix 
  
                           Q1, and on exit this is overwritten by Q1*Q.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            If COMPZ='N':  Z is not referenced.   
            If COMPZ='I':  on entry, Z need not be set, and on exit it   
                           contains the orthogonal matrix Z, which is   
                           the product of the Givens transformations   
                           which are applied to A and B on the right.   
            If COMPZ='V':  on entry, Z must contain an orthogonal matrix 
  
                           Z1, and on exit this is overwritten by Z1*Z.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.   
            LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    This routine reduces A to Hessenberg and B to triangular form by   
    an unblocked reduction, as described in _Matrix_Computations_,   
    by Golub and Van Loan (Johns Hopkins Press.)   

    ===================================================================== 
  


       Decode COMPQ   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b10 = 0.;
    static LONG DOUBLE c_b11 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int jcol;
    static LONG DOUBLE temp;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *);
    static int jrow;
    static LONG DOUBLE c, s;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *), xerbla_(char *, int *);
    static int icompq, icompz;
    static long int ilq, ilz;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    if (lsame_(compq, "N")) {
	ilq = 0;
	icompq = 1;
    } else if (lsame_(compq, "V")) {
	ilq = 1;
	icompq = 2;
    } else if (lsame_(compq, "I")) {
	ilq = 1;
	icompq = 3;
    } else {
	icompq = 0;
    }

/*     Decode COMPZ */

    if (lsame_(compz, "N")) {
	ilz = 0;
	icompz = 1;
    } else if (lsame_(compz, "V")) {
	ilz = 1;
	icompz = 2;
    } else if (lsame_(compz, "I")) {
	ilz = 1;
	icompz = 3;
    } else {
	icompz = 0;
    }

/*     Test the input parameters. */

    *info = 0;
    if (icompq <= 0) {
	*info = -1;
    } else if (icompz <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -5;
    } else if (*lda < MAX(1,*n)) {
	*info = -7;
    } else if (*ldb < MAX(1,*n)) {
	*info = -9;
    } else if ((ilq && *ldq < *n) || *ldq < 1) {
	*info = -11;
    } else if ((ilz && *ldz < *n) || *ldz < 1) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGHRD", &i__1);
	return;
    }

/*     Initialize Q and Z if desired. */

    if (icompq == 3) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b10, &c_b11, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b10, &c_b11, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b10, &c_b11, &Q(1,1), ldq);
#endif

    }
    if (icompz == 3) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b10, &c_b11, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b10, &c_b11, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b10, &c_b11, &Z(1,1), ldz);
#endif

    }

/*     Quick return if possible */

    if (*n <= 1) {
	return;
    }

/*     Zero out lower triangle of B */

    i__1 = *n - 1;
    for (jcol = 1; jcol <= *n-1; ++jcol) {
	i__2 = *n;
	for (jrow = jcol + 1; jrow <= *n; ++jrow) {
	    B(jrow,jcol) = 0.;
/* L10: */
	}
/* L20: */
    }

/*     Reduce A and B */

    i__1 = *ihi - 2;
    for (jcol = *ilo; jcol <= *ihi-2; ++jcol) {

	i__2 = jcol + 2;
	for (jrow = *ihi; jrow >= jcol+2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
 */

	    temp = A(jrow-1,jcol);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&temp, &A(jrow,jcol), &c, &s, &A(jrow-1,jcol));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&temp, &A(jrow,jcol), &c, &s, &A(jrow-1,jcol));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&temp, &A(jrow,jcol), &c, &s, &A(jrow-1,jcol));
#endif

	    A(jrow,jcol) = 0.;
	    i__3 = *n - jcol;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__3, &A(jrow-1,jcol+1), lda, &A(jrow,jcol+1), lda, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__3, &A(jrow-1,jcol+1), lda, &A(jrow,jcol+1), lda, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__3, &A(jrow-1,jcol+1), lda, &A(jrow,jcol+1), lda, &c, &s);
#endif

	    i__3 = *n + 2 - jrow;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__3, &B(jrow-1,jrow-1), ldb, &B(jrow,jrow-1), ldb, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__3, &B(jrow-1,jrow-1), ldb, &B(jrow,jrow-1), ldb, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__3, &B(jrow-1,jrow-1), ldb, &B(jrow,jrow-1), ldb, &c, &s);
#endif

	    if (ilq) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Q(1,jrow-1), &c__1, &Q(1,jrow), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Q(1,jrow-1), &c__1, &Q(1,jrow), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Q(1,jrow-1), &c__1, &Q(1,jrow), &c__1, &c, &s);
#endif

	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JR
OW-1) */

	    temp = B(jrow,jrow);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&temp, &B(jrow,jrow-1), &c, &s, &B(jrow,jrow));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&temp, &B(jrow,jrow-1), &c, &s, &B(jrow,jrow));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&temp, &B(jrow,jrow-1), &c, &s, &B(jrow,jrow));
#endif

	    B(jrow,jrow-1) = 0.;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(ihi, &A(1,jrow), &c__1, &A(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(ihi, &A(1,jrow), &c__1, &A(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(ihi, &A(1,jrow), &c__1, &A(1,jrow-1), &c__1, &c, &s);
#endif

	    i__3 = jrow - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__3, &B(1,jrow), &c__1, &B(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__3, &B(1,jrow), &c__1, &B(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__3, &B(1,jrow), &c__1, &B(1,jrow-1), &c__1, &c, &s);
#endif

	    if (ilz) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Z(1,jrow), &c__1, &Z(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Z(1,jrow), &c__1, &Z(1,jrow-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Z(1,jrow), &c__1, &Z(1,jrow-1), &c__1, &c, &s);
#endif

	    }
/* L30: */
	}
/* L40: */
    }

    return;

/*     End of DGGHRD */

} /* dgghrd_ */

