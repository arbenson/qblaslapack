#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dhgeqz_(char *job, char *compq, char *compz, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qhgeqz(char *job, char *compq, char *compz, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qhgeqz_(char *job, char *compq, char *compz, int *n, 
#endif

	int *ilo, int *ihi, LONG DOUBLE *a, int *lda, LONG DOUBLE *
	b, int *ldb, LONG DOUBLE *alphar, LONG DOUBLE *alphai, LONG DOUBLE *
	beta, LONG DOUBLE *q, int *ldq, LONG DOUBLE *z, int *ldz, 
	LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DHGEQZ implements a single-/LONG DOUBLE-shift version of the QZ method for 
  
    finding the generalized eigenvalues   

    w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation   

         det( A - w(i) B ) = 0   

    In addition, the pair A,B may be reduced to generalized Schur form:   
    B is upper triangular, and A is block upper triangular, where the   
    diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having 
  
    complex generalized eigenvalues (see the description of the argument 
  
    JOB.)   

    If JOB='S', then the pair (A,B) is simultaneously reduced to Schur   
    form by applying one orthogonal tranformation (usually called Q) on   
    the left and another (usually called Z) on the right.  The 2-by-2   
    upper-triangular diagonal blocks of B corresponding to 2-by-2 blocks 
  
    of A will be reduced to positive diagonal matrices.  (I.e.,   
    if A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) and   
    B(j+1,j+1) will be positive.)   

    If JOB='E', then at each iteration, the same transformations   
    are computed, but they are only applied to those parts of A and B   
    which are needed to compute ALPHAR, ALPHAI, and BETAR.   

    If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the orthogonal   
    transformations used to reduce (A,B) are accumulated into the arrays 
  
    Q and Z s.t.:   

         Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*   
         Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*   

    Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix 
  
         Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),   
         pp. 241--256.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will   
                   not necessarily be put into generalized Schur form.   
            = 'S': put A and B into generalized Schur form, as well   
                   as computing ALPHAR, ALPHAI, and BETA.   

    COMPQ   (input) CHARACTER*1   
            = 'N': do not modify Q.   
            = 'V': multiply the array Q on the right by the transpose of 
  
                   the orthogonal tranformation that is applied to the   
                   left side of A and B to reduce them to Schur form.   
            = 'I': like COMPQ='V', except that Q will be initialized to   
                   the identity first.   

    COMPZ   (input) CHARACTER*1   
            = 'N': do not modify Z.   
            = 'V': multiply the array Z on the right by the orthogonal   
                   tranformation that is applied to the right side of   
                   A and B to reduce them to Schur form.   
            = 'I': like COMPZ='V', except that Z will be initialized to   
                   the identity first.   

    N       (input) INTEGER   
            The order of the matrices A, B, Q, and Z.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows and 
  
            columns 1:ILO-1 and IHI+1:N.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the N-by-N upper Hessenberg matrix A.  Elements   
            below the subdiagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to generalized Schur form.   
            If JOB='E', then on exit A will have been destroyed.   
               The diagonal blocks will be correct, but the off-diagonal 
  
               portion will be meaningless.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX( 1, N ).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.  Elements   
            below the diagonal must be zero.  2-by-2 blocks in B   
            corresponding to 2-by-2 blocks in A will be reduced to   
            positive diagonal form.  (I.e., if A(j+1,j) is non-zero,   
            then B(j+1,j)=B(j,j+1)=0 and B(j,j) and B(j+1,j+1) will be   
            positive.)   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to Schur form.   
            If JOB='E', then on exit B will have been destroyed.   
               Elements corresponding to diagonal blocks of A will be   
               correct, but the off-diagonal portion will be meaningless. 
  

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX( 1, N ).   

    ALPHAR  (output) LONG DOUBLE PRECISION array, dimension (N)   
            ALPHAR(1:N) will be set to real parts of the diagonal   
            elements of A that would result from reducing A and B to   
            Schur form and then further reducing them both to triangular 
  
            form using unitary transformations s.t. the diagonal of B   
            was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block 
  
            (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=A(j,j).   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   

    ALPHAI  (output) LONG DOUBLE PRECISION array, dimension (N)   
            ALPHAI(1:N) will be set to imaginary parts of the diagonal   
            elements of A that would result from reducing A and B to   
            Schur form and then further reducing them both to triangular 
  
            form using unitary transformations s.t. the diagonal of B   
            was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block 
  
            (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=0.   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   

    BETA    (output) LONG DOUBLE PRECISION array, dimension (N)   
            BETA(1:N) will be set to the (real) diagonal elements of B   
            that would result from reducing A and B to Schur form and   
            then further reducing them both to triangular form using   
            unitary transformations s.t. the diagonal of B was   
            non-negative real.  Thus, if A(j,j) is in a 1-by-1 block   
            (i.e., A(j+1,j)=A(j,j+1)=0), then BETA(j)=B(j,j).   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   
            (Note that BETA(1:N) will always be non-negative, and no   
            BETAI is necessary.)   

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ, N)   
            If COMPQ='N', then Q will not be referenced.   
            If COMPQ='V' or 'I', then the transpose of the orthogonal   
               transformations which are applied to A and B on the left   
               will be applied to the array Q on the right.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= 1.   
            If COMPQ='V' or 'I', then LDQ >= N.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            If COMPZ='N', then Z will not be referenced.   
            If COMPZ='V' or 'I', then the orthogonal transformations   
               which are applied to A and B on the right will be applied 
  
               to the array Z on the right.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If COMPZ='V' or 'I', then LDZ >= N.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1,...,N: the QZ iteration did not converge.  (A,B) is not   
                       in Schur form, but ALPHAR(i), ALPHAI(i), and   
                       BETA(i), i=INFO+1,...,N should be correct.   
            = N+1,...,2*N: the shift calculation failed.  (A,B) is not   
                       in Schur form, but ALPHAR(i), ALPHAI(i), and   
                       BETA(i), i=INFO-N+1,...,N should be correct.   
            > 2*N:     various "impossible" errors.   

    Further Details   
    ===============   

    Iteration counters:   

    JITER  -- counts iterations.   
    IITER  -- counts iterations run since ILAST was last   
              changed.  This is therefore reset only when a 1-by-1 or   
              2-by-2 block deflates off the bottom.   

    ===================================================================== 
  

      $                     SAFETY = 1.0E+0 )   

       Decode JOB, COMPQ, COMPZ   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b12 = 0.;
    static LONG DOUBLE c_b13 = 1.;
    static int c__1 = 1;
    static int c__3 = 3;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE ad11l, ad12l, ad21l, ad22l, ad32l, wabs, atol, btol, 
	    temp;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *), dlag2_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *), qlag2(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *), qlag2_(
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *);
    static LONG DOUBLE temp2, s1inv, c;
    static int j;
    static LONG DOUBLE s, t, v[3], scale;
    extern long int lsame_(char *, char *);
    static int iiter, ilast, jiter;
    static LONG DOUBLE anorm, bnorm;
    static int maxit;
    static LONG DOUBLE tempi, tempr, s1, s2, u1, u2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlapy3_(LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2(LONG DOUBLE *, LONG DOUBLE *), dlapy3_(LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlapy2_(LONG DOUBLE *, LONG DOUBLE *), dlapy3_(LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasv2(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasv2_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static long int ilazr2;
    static LONG DOUBLE a11, a12, a21, a22, b11, b22, c12, c21;
    static int jc;
    static LONG DOUBLE an, bn, cl, cq, cr;
    static int in;
    static LONG DOUBLE ascale, bscale, u12, w11;
    static int jr;
    static LONG DOUBLE cz, sl, w12, w21, w22, wi;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE sr;

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
    static LONG DOUBLE vs, wr;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs(char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanhs_(char *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE safmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE safmax;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE eshift;
    static long int ilschr;
    static LONG DOUBLE b1a, b2a;
    static int icompq, ilastm;
    static LONG DOUBLE a1i;
    static int ischur;
    static LONG DOUBLE a2i, b1i;
    static long int ilazro;
    static int icompz, ifirst;
    static LONG DOUBLE b2i;
    static int ifrstm;
    static LONG DOUBLE a1r;
    static int istart;
    static long int ilpivt;
    static LONG DOUBLE a2r, b1r, b2r, wr2, ad11, ad12, ad21, ad22, c11i, c22i;
    static int jch;
    static LONG DOUBLE c11r, c22r, u12l;
    static long int ilq;
    static LONG DOUBLE tau, sqi;
    static long int ilz;
    static LONG DOUBLE ulp, sqr, szi, szr;



#define V(I) v[(I)]
#define ALPHAR(I) alphar[(I)-1]
#define ALPHAI(I) alphai[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    if (lsame_(job, "E")) {
	ilschr = 0;
	ischur = 1;
    } else if (lsame_(job, "S")) {
	ilschr = 1;
	ischur = 2;
    } else {
	ischur = 0;
    }

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

/*     Check Argument Values */

    *info = 0;
    if (ischur == 0) {
	*info = -1;
    } else if (icompq == 0) {
	*info = -2;
    } else if (icompz == 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1) {
	*info = -5;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -6;
    } else if (*lda < *n) {
	*info = -8;
    } else if (*ldb < *n) {
	*info = -10;
    } else if (*ldq < 1 || (ilq && *ldq < *n)) {
	*info = -15;
    } else if (*ldz < 1 || (ilz && *ldz < *n)) {
	*info = -17;
    } else if (*lwork < MAX(1,*n)) {
	*info = -19;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DHGEQZ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	WORK(1) = 1.;
	return;
    }

/*     Initialize Q and Z */

    if (icompq == 3) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b12, &c_b13, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b12, &c_b13, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b12, &c_b13, &Q(1,1), ldq);
#endif

    }
    if (icompz == 3) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b12, &c_b13, &Z(1,1), ldz);
#endif

    }

/*     Machine Constants */

    in = *ihi + 1 - *ilo;

#ifdef PETSC_PREFIX_SUFFIX
    safmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    safmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    safmin = qlamch_("S");
#endif

    safmax = 1. / safmin;

#ifdef PETSC_PREFIX_SUFFIX
    ulp = dlamch_("E") * dlamch_("B");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    ulp = qlamch("E") * dlamch_("B");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    ulp = qlamch_("E") * dlamch_("B");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlanhs_("F", &in, &A(*ilo,*ilo), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlanhs("F", &in, &A(*ilo,*ilo), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlanhs_("F", &in, &A(*ilo,*ilo), lda, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
    bnorm = dlanhs_("F", &in, &B(*ilo,*ilo), ldb, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnorm = qlanhs("F", &in, &B(*ilo,*ilo), ldb, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnorm = qlanhs_("F", &in, &B(*ilo,*ilo), ldb, &WORK(1));
#endif

/* Computing MAX */
    d__1 = safmin, d__2 = ulp * anorm;
    atol = MAX(d__1,d__2);
/* Computing MAX */
    d__1 = safmin, d__2 = ulp * bnorm;
    btol = MAX(d__1,d__2);
    ascale = 1. / MAX(safmin,anorm);
    bscale = 1. / MAX(safmin,bnorm);

/*     Set Eigenvalues IHI+1:N */

    i__1 = *n;
    for (j = *ihi + 1; j <= *n; ++j) {
	if (B(j,j) < 0.) {
	    if (ilschr) {
		i__2 = j;
		for (jr = 1; jr <= j; ++jr) {
		    A(jr,j) = -A(jr,j);
		    B(jr,j) = -B(jr,j);
/* L10: */
		}
	    } else {
		A(j,j) = -A(j,j);
		B(j,j) = -B(j,j);
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    Z(jr,j) = -Z(jr,j);
/* L20: */
		}
	    }
	}
	ALPHAR(j) = A(j,j);
	ALPHAI(j) = 0.;
	BETA(j) = B(j,j);
/* L30: */
    }

/*     If IHI < ILO, skip QZ steps */

    if (*ihi < *ilo) {
	goto L380;
    }

/*     MAIN QZ ITERATION LOOP   

       Initialize dynamic indices   

       Eigenvalues ILAST+1:N have been found.   
          Column operations modify rows IFRSTM:whatever.   
          Row operations modify columns whatever:ILASTM.   

       If only eigenvalues are being computed, then   
          IFRSTM is the row of the last splitting row above row ILAST;   
          this is always at least ILO.   
       IITER counts iterations since the last eigenvalue was found,   
          to tell when to use an extraordinary shift.   
       MAXIT is the maximum number of QZ sweeps allowed. */

    ilast = *ihi;
    if (ilschr) {
	ifrstm = 1;
	ilastm = *n;
    } else {
	ifrstm = *ilo;
	ilastm = *ihi;
    }
    iiter = 0;
    eshift = 0.;
    maxit = (*ihi - *ilo + 1) * 30;

    i__1 = maxit;
    for (jiter = 1; jiter <= maxit; ++jiter) {

/*        Split the matrix if possible.   

          Two tests:   
             1: A(j,j-1)=0  or  j=ILO   
             2: B(j,j)=0 */

	if (ilast == *ilo) {

/*           Special case: j=ILAST */

	    goto L80;
	} else {
	    if ((d__1 = A(ilast,ilast-1), ABS(d__1)) <= atol) {
		A(ilast,ilast-1) = 0.;
		goto L80;
	    }
	}

	if ((d__1 = B(ilast,ilast), ABS(d__1)) <= btol) {
	    B(ilast,ilast) = 0.;
	    goto L70;
	}

/*        General case: j<ILAST */

	i__2 = *ilo;
	for (j = ilast - 1; j >= *ilo; --j) {

/*           Test 1: for A(j,j-1)=0 or j=ILO */

	    if (j == *ilo) {
		ilazro = 1;
	    } else {
		if ((d__1 = A(j,j-1), ABS(d__1)) <= atol) {
		    A(j,j-1) = 0.;
		    ilazro = 1;
		} else {
		    ilazro = 0;
		}
	    }

/*           Test 2: for B(j,j)=0 */

	    if ((d__1 = B(j,j), ABS(d__1)) < btol) {
		B(j,j) = 0.;

/*              Test 1a: Check for 2 consecutive small subdiag
onals in A */

		ilazr2 = 0;
		if (! ilazro) {
		    temp = (d__1 = A(j,j-1), ABS(d__1));
		    temp2 = (d__1 = A(j,j), ABS(d__1));
		    tempr = MAX(temp,temp2);
		    if (tempr < 1. && tempr != 0.) {
			temp /= tempr;
			temp2 /= tempr;
		    }
		    if (temp * (ascale * (d__1 = A(j+1,j), ABS(
			    d__1))) <= temp2 * (ascale * atol)) {
			ilazr2 = 1;
		    }
		}

/*              If both tests pass (1 & 2), i.e., the leading 
diagonal   
                element of B in the block is zero, split a 1x1
 block off   
                at the top. (I.e., at the J-th row/column) The
 leading   
                diagonal element of the remainder can also be 
zero, so   
                this may have to be done repeatedly. */

		if (ilazro || ilazr2) {
		    i__3 = ilast - 1;
		    for (jch = j; jch <= ilast-1; ++jch) {
			temp = A(jch,jch);

#ifdef PETSC_PREFIX_SUFFIX
			dlartg_(&temp, &A(jch+1,jch), &c, &s, &A(jch,jch));
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlartg(&temp, &A(jch+1,jch), &c, &s, &A(jch,jch));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlartg_(&temp, &A(jch+1,jch), &c, &s, &A(jch,jch));
#endif

			A(jch+1,jch) = 0.;
			i__4 = ilastm - jch;

#ifdef PETSC_PREFIX_SUFFIX
			drot_(&i__4, &A(jch,jch+1), lda, &A(jch+1,jch+1), lda, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qrot(&i__4, &A(jch,jch+1), lda, &A(jch+1,jch+1), lda, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qrot_(&i__4, &A(jch,jch+1), lda, &A(jch+1,jch+1), lda, &c, &s);
#endif

			i__4 = ilastm - jch;

#ifdef PETSC_PREFIX_SUFFIX
			drot_(&i__4, &B(jch,jch+1), ldb, &B(jch+1,jch+1), ldb, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qrot(&i__4, &B(jch,jch+1), ldb, &B(jch+1,jch+1), ldb, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qrot_(&i__4, &B(jch,jch+1), ldb, &B(jch+1,jch+1), ldb, &c, &s);
#endif

			if (ilq) {

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif

			}
			if (ilazr2) {
			    A(jch,jch-1) *= c;
			}
			ilazr2 = 0;
			if ((d__1 = B(jch+1,jch+1), ABS(d__1)
				) >= btol) {
			    if (jch + 1 >= ilast) {
				goto L80;
			    } else {
				ifirst = jch + 1;
				goto L110;
			    }
			}
			B(jch+1,jch+1) = 0.;
/* L40: */
		    }
		    goto L70;
		} else {

/*                 Only test 2 passed -- chase the zero to
 B(ILAST,ILAST)   
                   Then process as in the case B(ILAST,ILA
ST)=0 */

		    i__3 = ilast - 1;
		    for (jch = j; jch <= ilast-1; ++jch) {
			temp = B(jch,jch+1);

#ifdef PETSC_PREFIX_SUFFIX
			dlartg_(&temp, &B(jch+1,jch+1), &c, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlartg(&temp, &B(jch+1,jch+1), &c, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlartg_(&temp, &B(jch+1,jch+1), &c, &
#endif

				s, &B(jch,jch+1));
			B(jch+1,jch+1) = 0.;
			if (jch < ilastm - 1) {
			    i__4 = ilastm - jch - 1;

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(&i__4, &B(jch,jch+2), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(&i__4, &B(jch,jch+2), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(&i__4, &B(jch,jch+2), ldb, &
#endif

				    B(jch+1,jch+2), ldb, &c, 
				    &s);
			}
			i__4 = ilastm - jch + 2;

#ifdef PETSC_PREFIX_SUFFIX
			drot_(&i__4, &A(jch,jch-1), lda, &A(jch+1,jch-1), lda, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qrot(&i__4, &A(jch,jch-1), lda, &A(jch+1,jch-1), lda, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qrot_(&i__4, &A(jch,jch-1), lda, &A(jch+1,jch-1), lda, &c, &s);
#endif

			if (ilq) {

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(n, &Q(1,jch), &c__1, &Q(1,jch+1), &c__1, &c, &s);
#endif

			}
			temp = A(jch+1,jch);

#ifdef PETSC_PREFIX_SUFFIX
			dlartg_(&temp, &A(jch+1,jch-1), &c, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qlartg(&temp, &A(jch+1,jch-1), &c, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qlartg_(&temp, &A(jch+1,jch-1), &c, &
#endif

				s, &A(jch+1,jch));
			A(jch+1,jch-1) = 0.;
			i__4 = jch + 1 - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
			drot_(&i__4, &A(ifrstm,jch), &c__1, &A(ifrstm,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qrot(&i__4, &A(ifrstm,jch), &c__1, &A(ifrstm,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qrot_(&i__4, &A(ifrstm,jch), &c__1, &A(ifrstm,jch-1), &c__1, &c, &s);
#endif

			i__4 = jch - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
			drot_(&i__4, &B(ifrstm,jch), &c__1, &B(ifrstm,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qrot(&i__4, &B(ifrstm,jch), &c__1, &B(ifrstm,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qrot_(&i__4, &B(ifrstm,jch), &c__1, &B(ifrstm,jch-1), &c__1, &c, &s);
#endif

			if (ilz) {

#ifdef PETSC_PREFIX_SUFFIX
			    drot_(n, &Z(1,jch), &c__1, &Z(1,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qrot(n, &Z(1,jch), &c__1, &Z(1,jch-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qrot_(n, &Z(1,jch), &c__1, &Z(1,jch-1), &c__1, &c, &s);
#endif

			}
/* L50: */
		    }
		    goto L70;
		}
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

		ifirst = j;
		goto L110;
	    }

/*           Neither test passed -- try next J   

   L60: */
	}

/*        (Drop-through is "impossible") */

	*info = *n + 1;
	goto L420;

/*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a   
          1x1 block. */

L70:
	temp = A(ilast,ilast);

#ifdef PETSC_PREFIX_SUFFIX
	dlartg_(&temp, &A(ilast,ilast-1), &c, &s, &A(ilast,ilast));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlartg(&temp, &A(ilast,ilast-1), &c, &s, &A(ilast,ilast));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlartg_(&temp, &A(ilast,ilast-1), &c, &s, &A(ilast,ilast));
#endif

	A(ilast,ilast-1) = 0.;
	i__2 = ilast - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
	drot_(&i__2, &A(ifrstm,ilast), &c__1, &A(ifrstm,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qrot(&i__2, &A(ifrstm,ilast), &c__1, &A(ifrstm,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qrot_(&i__2, &A(ifrstm,ilast), &c__1, &A(ifrstm,ilast-1), &c__1, &c, &s);
#endif

	i__2 = ilast - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
	drot_(&i__2, &B(ifrstm,ilast), &c__1, &B(ifrstm,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qrot(&i__2, &B(ifrstm,ilast), &c__1, &B(ifrstm,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qrot_(&i__2, &B(ifrstm,ilast), &c__1, &B(ifrstm,ilast-1), &c__1, &c, &s);
#endif

	if (ilz) {

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(n, &Z(1,ilast), &c__1, &Z(1,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(n, &Z(1,ilast), &c__1, &Z(1,ilast-1), &c__1, &c, &s);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(n, &Z(1,ilast), &c__1, &Z(1,ilast-1), &c__1, &c, &s);
#endif

	}

/*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,   
                                and BETA */

L80:
	if (B(ilast,ilast) < 0.) {
	    if (ilschr) {
		i__2 = ilast;
		for (j = ifrstm; j <= ilast; ++j) {
		    A(j,ilast) = -A(j,ilast);
		    B(j,ilast) = -B(j,ilast);
/* L90: */
		}
	    } else {
		A(ilast,ilast) = -A(ilast,ilast);
		B(ilast,ilast) = -B(ilast,ilast);
	    }
	    if (ilz) {
		i__2 = *n;
		for (j = 1; j <= *n; ++j) {
		    Z(j,ilast) = -Z(j,ilast);
/* L100: */
		}
	    }
	}
	ALPHAR(ilast) = A(ilast,ilast);
	ALPHAI(ilast) = 0.;
	BETA(ilast) = B(ilast,ilast);

/*        Go to next block -- exit if finished. */

	--ilast;
	if (ilast < *ilo) {
	    goto L380;
	}

/*        Reset counters */

	iiter = 0;
	eshift = 0.;
	if (! ilschr) {
	    ilastm = ilast;
	    if (ifrstm > ilast) {
		ifrstm = *ilo;
	    }
	}
	goto L350;

/*        QZ step   

          This iteration only involves rows/columns IFIRST:ILAST. We 
  
          assume IFIRST < ILAST, and that the diagonal of B is non-zer
o. */

L110:
	++iiter;
	if (! ilschr) {
	    ifrstm = ifirst;
	}

/*        Compute single shifts.   

          At this point, IFIRST < ILAST, and the diagonal elements of 
  
          B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in   
          magnitude) */

	if (iiter / 10 * 10 == iiter) {

/*           Exceptional shift.  Chosen for no particularly good r
eason.   
             (Single shift only.) */

	    if ((LONG DOUBLE) maxit * safmin * (d__1 = A(ilast-1,ilast), ABS(d__1)) < (d__2 = B(ilast-1,ilast-1), ABS(d__2))) {
		eshift += A(ilast-1,ilast) / B(ilast-1,ilast-1);
	    } else {
		eshift += 1. / (safmin * (LONG DOUBLE) maxit);
	    }
	    s1 = 1.;
	    wr = eshift;

	} else {

/*           Shifts based on the generalized eigenvalues of the   
             bottom-right 2x2 block of A and B. The first eigenval
ue   
             returned by DLAG2 is the Wilkinson shift (AEP p.512),
 */

	    d__1 = safmin * 100.;

#ifdef PETSC_PREFIX_SUFFIX
	    dlag2_(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &s2, &wr, &wr2, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlag2(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &s2, &wr, &wr2, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlag2_(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &s2, &wr, &wr2, &
#endif

		    wi);

/* Computing MAX   
   Computing MAX */
	    d__3 = 1., d__4 = ABS(wr), d__3 = MAX(d__3,d__4), d__4 = ABS(wi);
	    d__1 = s1, d__2 = safmin * MAX(d__3,d__4);
	    temp = MAX(d__1,d__2);
	    if (wi != 0.) {
		goto L200;
	    }
	}

/*        Fiddle with shift to avoid overflow */

	temp = MIN(ascale,1.) * (safmax * .5);
	if (s1 > temp) {
	    scale = temp / s1;
	} else {
	    scale = 1.;
	}

	temp = MIN(bscale,1.) * (safmax * .5);
	if (ABS(wr) > temp) {
/* Computing MIN */
	    d__1 = scale, d__2 = temp / ABS(wr);
	    scale = MIN(d__1,d__2);
	}
	s1 = scale * s1;
	wr = scale * wr;

/*        Now check for two consecutive small subdiagonals. */

	i__2 = ifirst + 1;
	for (j = ilast - 1; j >= ifirst+1; --j) {
	    istart = j;
	    temp = (d__1 = s1 * A(j,j-1), ABS(d__1));
	    temp2 = (d__1 = s1 * A(j,j) - wr * B(j,j), 
		    ABS(d__1));
	    tempr = MAX(temp,temp2);
	    if (tempr < 1. && tempr != 0.) {
		temp /= tempr;
		temp2 /= tempr;
	    }
	    if ((d__1 = ascale * A(j+1,j) * temp, ABS(d__1)) <= 
		    ascale * atol * temp2) {
		goto L130;
	    }
/* L120: */
	}

	istart = ifirst;
L130:

/*        Do an implicit single-shift QZ sweep.   

          Initial Q */

	temp = s1 * A(istart,istart) - wr * B(istart,istart);
	temp2 = s1 * A(istart+1,istart);

#ifdef PETSC_PREFIX_SUFFIX
	dlartg_(&temp, &temp2, &c, &s, &tempr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlartg(&temp, &temp2, &c, &s, &tempr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlartg_(&temp, &temp2, &c, &s, &tempr);
#endif


/*        Sweep */

	i__2 = ilast - 1;
	for (j = istart; j <= ilast-1; ++j) {
	    if (j > istart) {
		temp = A(j,j-1);

#ifdef PETSC_PREFIX_SUFFIX
		dlartg_(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlartg(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlartg_(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif

		A(j+1,j-1) = 0.;
	    }

	    i__3 = ilastm;
	    for (jc = j; jc <= ilastm; ++jc) {
		temp = c * A(j,jc) + s * A(j+1,jc);
		A(j+1,jc) = -s * A(j,jc) + c * A(j+1,jc);
		A(j,jc) = temp;
		temp2 = c * B(j,jc) + s * B(j+1,jc);
		B(j+1,jc) = -s * B(j,jc) + c * B(j+1,jc);
		B(j,jc) = temp2;
/* L140: */
	    }
	    if (ilq) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    temp = c * Q(jr,j) + s * Q(jr,j+1);
		    Q(jr,j+1) = -s * Q(jr,j) + c * 
			    Q(jr,j+1);
		    Q(jr,j) = temp;
/* L150: */
		}
	    }

	    temp = B(j+1,j+1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif

	    B(j+1,j) = 0.;

/* Computing MIN */
	    i__4 = j + 2;
	    i__3 = MIN(i__4,ilast);
	    for (jr = ifrstm; jr <= MIN(j+2,ilast); ++jr) {
		temp = c * A(jr,j+1) + s * A(jr,j);
		A(jr,j) = -s * A(jr,j+1) + c * A(jr,j);
		A(jr,j+1) = temp;
/* L160: */
	    }
	    i__3 = j;
	    for (jr = ifrstm; jr <= j; ++jr) {
		temp = c * B(jr,j+1) + s * B(jr,j);
		B(jr,j) = -s * B(jr,j+1) + c * B(jr,j);
		B(jr,j+1) = temp;
/* L170: */
	    }
	    if (ilz) {
		i__3 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    temp = c * Z(jr,j+1) + s * Z(jr,j);
		    Z(jr,j) = -s * Z(jr,j+1) + c * 
			    Z(jr,j);
		    Z(jr,j+1) = temp;
/* L180: */
		}
	    }
/* L190: */
	}

	goto L350;

/*        Use Francis LONG DOUBLE-shift   

          Note: the Francis LONG DOUBLE-shift should work with real shifts,
   
                but only if the block is at least 3x3.   
                This code may break if this point is reached with   
                a 2x2 block with real eigenvalues. */

L200:
	if (ifirst + 1 == ilast) {

/*           Special case -- 2x2 block with complex eigenvectors 
  

             Step 1: Standardize, that is, rotate so that   

                         ( B11  0  )   
                     B = (         )  with B11 non-negative.   
                         (  0  B22 ) */


#ifdef PETSC_PREFIX_SUFFIX
	    dlasv2_(&B(ilast-1,ilast-1), &B(ilast-1,ilast), &B(ilast,ilast), &b22, &b11, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasv2(&B(ilast-1,ilast-1), &B(ilast-1,ilast), &B(ilast,ilast), &b22, &b11, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasv2_(&B(ilast-1,ilast-1), &B(ilast-1,ilast), &B(ilast,ilast), &b22, &b11, &
#endif

		    sr, &cr, &sl, &cl);

	    if (b11 < 0.) {
		cr = -cr;
		sr = -sr;
		b11 = -b11;
		b22 = -b22;
	    }

	    i__2 = ilastm + 1 - ifirst;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__2, &A(ilast-1,ilast-1), lda, &A(ilast,ilast-1), lda, &cl, &sl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__2, &A(ilast-1,ilast-1), lda, &A(ilast,ilast-1), lda, &cl, &sl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__2, &A(ilast-1,ilast-1), lda, &A(ilast,ilast-1), lda, &cl, &sl);
#endif

	    i__2 = ilast + 1 - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
	    drot_(&i__2, &A(ifrstm,ilast-1), &c__1, &A(ifrstm,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qrot(&i__2, &A(ifrstm,ilast-1), &c__1, &A(ifrstm,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qrot_(&i__2, &A(ifrstm,ilast-1), &c__1, &A(ifrstm,ilast), &c__1, &cr, &sr);
#endif


	    if (ilast < ilastm) {
		i__2 = ilastm - ilast;

#ifdef PETSC_PREFIX_SUFFIX
		drot_(&i__2, &B(ilast-1,ilast+1), ldb, &B(ilast,ilast+1), lda, &cl, &sl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(&i__2, &B(ilast-1,ilast+1), ldb, &B(ilast,ilast+1), lda, &cl, &sl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(&i__2, &B(ilast-1,ilast+1), ldb, &B(ilast,ilast+1), lda, &cl, &sl);
#endif

	    }
	    if (ifrstm < ilast - 1) {
		i__2 = ifirst - ifrstm;

#ifdef PETSC_PREFIX_SUFFIX
		drot_(&i__2, &B(ifrstm,ilast-1), &c__1, &B(ifrstm,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(&i__2, &B(ifrstm,ilast-1), &c__1, &B(ifrstm,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(&i__2, &B(ifrstm,ilast-1), &c__1, &B(ifrstm,ilast), &c__1, &cr, &sr);
#endif

	    }

	    if (ilq) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Q(1,ilast-1), &c__1, &Q(1,ilast), &c__1, &cl, &sl);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Q(1,ilast-1), &c__1, &Q(1,ilast), &c__1, &cl, &sl);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Q(1,ilast-1), &c__1, &Q(1,ilast), &c__1, &cl, &sl);
#endif

	    }
	    if (ilz) {

#ifdef PETSC_PREFIX_SUFFIX
		drot_(n, &Z(1,ilast-1), &c__1, &Z(1,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(n, &Z(1,ilast-1), &c__1, &Z(1,ilast), &c__1, &cr, &sr);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(n, &Z(1,ilast-1), &c__1, &Z(1,ilast), &c__1, &cr, &sr);
#endif

	    }

	    B(ilast-1,ilast-1) = b11;
	    B(ilast-1,ilast) = 0.;
	    B(ilast,ilast-1) = 0.;
	    B(ilast,ilast) = b22;

/*           If B22 is negative, negate column ILAST */

	    if (b22 < 0.) {
		i__2 = ilast;
		for (j = ifrstm; j <= ilast; ++j) {
		    A(j,ilast) = -A(j,ilast);
		    B(j,ilast) = -B(j,ilast);
/* L210: */
		}

		if (ilz) {
		    i__2 = *n;
		    for (j = 1; j <= *n; ++j) {
			Z(j,ilast) = -Z(j,ilast);
/* L220: */
		    }
		}
	    }

/*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.) 
  

             Recompute shift */

	    d__1 = safmin * 100.;

#ifdef PETSC_PREFIX_SUFFIX
	    dlag2_(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &temp, &wr, &temp2,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlag2(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &temp, &wr, &temp2,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlag2_(&A(ilast-1,ilast-1), lda, &B(ilast-1,ilast-1), ldb, &d__1, &s1, &temp, &wr, &temp2,
#endif

		     &wi);

/*           If standardization has perturbed the shift onto real 
line,   
             do another (real single-shift) QR step. */

	    if (wi == 0.) {
		goto L350;
	    }
	    s1inv = 1. / s1;

/*           Do EISPACK (QZVAL) computation of alpha and beta */

	    a11 = A(ilast-1,ilast-1);
	    a21 = A(ilast,ilast-1);
	    a12 = A(ilast-1,ilast);
	    a22 = A(ilast,ilast);

/*           Compute complex Givens rotation on right   
             (Assume some element of C = (sA - wB) > unfl )   
                              __   
             (sA - wB) ( CZ   -SZ )   
                       ( SZ    CZ ) */

	    c11r = s1 * a11 - wr * b11;
	    c11i = -wi * b11;
	    c12 = s1 * a12;
	    c21 = s1 * a21;
	    c22r = s1 * a22 - wr * b22;
	    c22i = -wi * b22;

	    if (ABS(c11r) + ABS(c11i) + ABS(c12) > ABS(c21) + ABS(c22r) + ABS(
		    c22i)) {

#ifdef PETSC_PREFIX_SUFFIX
		t = dlapy3_(&c12, &c11r, &c11i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		t = qlapy3(&c12, &c11r, &c11i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		t = qlapy3_(&c12, &c11r, &c11i);
#endif

		cz = c12 / t;
		szr = -c11r / t;
		szi = -c11i / t;
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		cz = dlapy2_(&c22r, &c22i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		cz = qlapy2(&c22r, &c22i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		cz = qlapy2_(&c22r, &c22i);
#endif

		if (cz <= safmin) {
		    cz = 0.;
		    szr = 1.;
		    szi = 0.;
		} else {
		    tempr = c22r / cz;
		    tempi = c22i / cz;

#ifdef PETSC_PREFIX_SUFFIX
		    t = dlapy2_(&cz, &c21);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    t = qlapy2(&cz, &c21);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    t = qlapy2_(&cz, &c21);
#endif

		    cz /= t;
		    szr = -c21 * tempr / t;
		    szi = c21 * tempi / t;
		}
	    }

/*           Compute Givens rotation on left   

             (  CQ   SQ )   
             (  __      )  A or B   
             ( -SQ   CQ ) */

	    an = ABS(a11) + ABS(a12) + ABS(a21) + ABS(a22);
	    bn = ABS(b11) + ABS(b22);
	    wabs = ABS(wr) + ABS(wi);
	    if (s1 * an > wabs * bn) {
		cq = cz * b11;
		sqr = szr * b22;
		sqi = -szi * b22;
	    } else {
		a1r = cz * a11 + szr * a12;
		a1i = szi * a12;
		a2r = cz * a21 + szr * a22;
		a2i = szi * a22;

#ifdef PETSC_PREFIX_SUFFIX
		cq = dlapy2_(&a1r, &a1i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		cq = qlapy2(&a1r, &a1i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		cq = qlapy2_(&a1r, &a1i);
#endif

		if (cq <= safmin) {
		    cq = 0.;
		    sqr = 1.;
		    sqi = 0.;
		} else {
		    tempr = a1r / cq;
		    tempi = a1i / cq;
		    sqr = tempr * a2r + tempi * a2i;
		    sqi = tempi * a2r - tempr * a2i;
		}
	    }

#ifdef PETSC_PREFIX_SUFFIX
	    t = dlapy3_(&cq, &sqr, &sqi);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    t = qlapy3(&cq, &sqr, &sqi);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    t = qlapy3_(&cq, &sqr, &sqi);
#endif

	    cq /= t;
	    sqr /= t;
	    sqi /= t;

/*           Compute diagonal elements of QBZ */

	    tempr = sqr * szr - sqi * szi;
	    tempi = sqr * szi + sqi * szr;
	    b1r = cq * cz * b11 + tempr * b22;
	    b1i = tempi * b22;

#ifdef PETSC_PREFIX_SUFFIX
	    b1a = dlapy2_(&b1r, &b1i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    b1a = qlapy2(&b1r, &b1i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    b1a = qlapy2_(&b1r, &b1i);
#endif

	    b2r = cq * cz * b22 + tempr * b11;
	    b2i = -tempi * b11;

#ifdef PETSC_PREFIX_SUFFIX
	    b2a = dlapy2_(&b2r, &b2i);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    b2a = qlapy2(&b2r, &b2i);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    b2a = qlapy2_(&b2r, &b2i);
#endif


/*           Normalize so beta > 0, and Im( alpha1 ) > 0 */

	    BETA(ilast - 1) = b1a;
	    BETA(ilast) = b2a;
	    ALPHAR(ilast - 1) = wr * b1a * s1inv;
	    ALPHAI(ilast - 1) = wi * b1a * s1inv;
	    ALPHAR(ilast) = wr * b2a * s1inv;
	    ALPHAI(ilast) = -(wi * b2a) * s1inv;

/*           Step 3: Go to next block -- exit if finished. */

	    ilast = ifirst - 1;
	    if (ilast < *ilo) {
		goto L380;
	    }

/*           Reset counters */

	    iiter = 0;
	    eshift = 0.;
	    if (! ilschr) {
		ilastm = ilast;
		if (ifrstm > ilast) {
		    ifrstm = *ilo;
		}
	    }
	    goto L350;
	} else {

/*           Usual case: 3x3 or larger block, using Francis implic
it   
                         LONG DOUBLE-shift   

                                      2   
             Eigenvalue equation is  w  - c w + d = 0,   

                                           -1 2        -1   
             so compute 1st column of  (A B  )  - c A B   + d   
             using the formula in QZIT (from EISPACK)   

             We assume that the block is at least 3x3 */

	    ad11 = ascale * A(ilast-1,ilast-1) / (bscale * B(ilast-1,ilast-1));
	    ad21 = ascale * A(ilast,ilast-1) / (bscale * B(ilast-1,ilast-1));
	    ad12 = ascale * A(ilast-1,ilast) / (bscale * B(ilast,ilast));
	    ad22 = ascale * A(ilast,ilast) / (bscale * B(ilast,ilast));
	    u12 = B(ilast-1,ilast) / B(ilast,ilast);
	    ad11l = ascale * A(ifirst,ifirst) / (bscale * B(ifirst,ifirst));
	    ad21l = ascale * A(ifirst+1,ifirst) / (bscale * B(ifirst,ifirst));
	    ad12l = ascale * A(ifirst,ifirst+1) / (bscale * B(ifirst+1,ifirst+1));
	    ad22l = ascale * A(ifirst+1,ifirst+1) / (bscale *
		     B(ifirst+1,ifirst+1));
	    ad32l = ascale * A(ifirst+2,ifirst+1) / (bscale *
		     B(ifirst+1,ifirst+1));
	    u12l = B(ifirst,ifirst+1) / B(ifirst+1,ifirst+1);

	    V(0) = (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * u12 
		    * ad11l + (ad12l - ad11l * u12l) * ad21l;
	    V(1) = (ad22l - ad11l - ad21l * u12l - (ad11 - ad11l) - (ad22 - 
		    ad11l) + ad21 * u12) * ad21l;
	    V(2) = ad32l * ad21l;

	    istart = ifirst;


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfg_(&c__3, v, &V(1), &c__1, &tau);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfg(&c__3, v, &V(1), &c__1, &tau);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfg_(&c__3, v, &V(1), &c__1, &tau);
#endif

	    V(0) = 1.;

/*           Sweep */

	    i__2 = ilast - 2;
	    for (j = istart; j <= ilast-2; ++j) {

/*              All but last elements: use 3x3 Householder tra
nsforms.   

                Zero (j-1)st column of A */

		if (j > istart) {
		    V(0) = A(j,j-1);
		    V(1) = A(j+1,j-1);
		    V(2) = A(j+2,j-1);


#ifdef PETSC_PREFIX_SUFFIX
		    dlarfg_(&c__3, &A(j,j-1), &V(1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlarfg(&c__3, &A(j,j-1), &V(1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlarfg_(&c__3, &A(j,j-1), &V(1), &c__1, &
#endif

			    tau);
		    V(0) = 1.;
		    A(j+1,j-1) = 0.;
		    A(j+2,j-1) = 0.;
		}

		i__3 = ilastm;
		for (jc = j; jc <= ilastm; ++jc) {
		    temp = tau * (A(j,jc) + V(1) * A(j+1,jc) + V(2) * A(j+2,jc));
		    A(j,jc) -= temp;
		    A(j+1,jc) -= temp * V(1);
		    A(j+2,jc) -= temp * V(2);
		    temp2 = tau * (B(j,jc) + V(1) * B(j+1,jc) + V(2) * B(j+2,jc));
		    B(j,jc) -= temp2;
		    B(j+1,jc) -= temp2 * V(1);
		    B(j+2,jc) -= temp2 * V(2);
/* L230: */
		}
		if (ilq) {
		    i__3 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			temp = tau * (Q(jr,j) + V(1) * Q(jr,j+1) + V(2) * Q(jr,j+2)
				);
			Q(jr,j) -= temp;
			Q(jr,j+1) -= temp * V(1);
			Q(jr,j+2) -= temp * V(2);
/* L240: */
		    }
		}

/*              Zero j-th column of B (see DLAGBC for details)
   

                Swap rows to pivot */

		ilpivt = 0;
/* Computing MAX */
		d__3 = (d__1 = B(j+1,j+1), ABS(d__1)), d__4 =
			 (d__2 = B(j+1,j+2), ABS(d__2));
		temp = MAX(d__3,d__4);
/* Computing MAX */
		d__3 = (d__1 = B(j+2,j+1), ABS(d__1)), d__4 =
			 (d__2 = B(j+2,j+2), ABS(d__2));
		temp2 = MAX(d__3,d__4);
		if (MAX(temp,temp2) < safmin) {
		    scale = 0.;
		    u1 = 1.;
		    u2 = 0.;
		    goto L250;
		} else if (temp >= temp2) {
		    w11 = B(j+1,j+1);
		    w21 = B(j+2,j+1);
		    w12 = B(j+1,j+2);
		    w22 = B(j+2,j+2);
		    u1 = B(j+1,j);
		    u2 = B(j+2,j);
		} else {
		    w21 = B(j+1,j+1);
		    w11 = B(j+2,j+1);
		    w22 = B(j+1,j+2);
		    w12 = B(j+2,j+2);
		    u2 = B(j+1,j);
		    u1 = B(j+2,j);
		}

/*              Swap columns if nec. */

		if (ABS(w12) > ABS(w11)) {
		    ilpivt = 1;
		    temp = w12;
		    temp2 = w22;
		    w12 = w11;
		    w22 = w21;
		    w11 = temp;
		    w21 = temp2;
		}

/*              LU-factor */

		temp = w21 / w11;
		u2 -= temp * u1;
		w22 -= temp * w12;
		w21 = 0.;

/*              Compute SCALE */

		scale = 1.;
		if (ABS(w22) < safmin) {
		    scale = 0.;
		    u2 = 1.;
		    u1 = -w12 / w11;
		    goto L250;
		}
		if (ABS(w22) < ABS(u2)) {
		    scale = (d__1 = w22 / u2, ABS(d__1));
		}
		if (ABS(w11) < ABS(u1)) {
/* Computing MIN */
		    d__2 = scale, d__3 = (d__1 = w11 / u1, ABS(d__1));
		    scale = MIN(d__2,d__3);
		}

/*              Solve */

		u2 = scale * u2 / w22;
		u1 = (scale * u1 - w12 * u2) / w11;

L250:
		if (ilpivt) {
		    temp = u2;
		    u2 = u1;
		    u1 = temp;
		}

/*              Compute Householder Vector   

   Computing 2nd power */
		d__1 = scale;
/* Computing 2nd power */
		d__2 = u1;
/* Computing 2nd power */
		d__3 = u2;
		t = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		tau = scale / t + 1.;
		vs = -1. / (scale + t);
		V(0) = 1.;
		V(1) = vs * u1;
		V(2) = vs * u2;

/*              Apply transformations from the right.   

   Computing MIN */
		i__4 = j + 3;
		i__3 = MIN(i__4,ilast);
		for (jr = ifrstm; jr <= MIN(j+3,ilast); ++jr) {
		    temp = tau * (A(jr,j) + V(1) * A(jr,j+1) + V(2) * A(jr,j+2));
		    A(jr,j) -= temp;
		    A(jr,j+1) -= temp * V(1);
		    A(jr,j+2) -= temp * V(2);
/* L260: */
		}
		i__3 = j + 2;
		for (jr = ifrstm; jr <= j+2; ++jr) {
		    temp = tau * (B(jr,j) + V(1) * B(jr,j+1) + V(2) * B(jr,j+2));
		    B(jr,j) -= temp;
		    B(jr,j+1) -= temp * V(1);
		    B(jr,j+2) -= temp * V(2);
/* L270: */
		}
		if (ilz) {
		    i__3 = *n;
		    for (jr = 1; jr <= *n; ++jr) {
			temp = tau * (Z(jr,j) + V(1) * Z(jr,j+1) + V(2) * Z(jr,j+2)
				);
			Z(jr,j) -= temp;
			Z(jr,j+1) -= temp * V(1);
			Z(jr,j+2) -= temp * V(2);
/* L280: */
		    }
		}
		B(j+1,j) = 0.;
		B(j+2,j) = 0.;
/* L290: */
	    }

/*           Last elements: Use Givens rotations   

             Rotations from the left */

	    j = ilast - 1;
	    temp = A(j,j-1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&temp, &A(j+1,j-1), &c, &s, &A(j,j-1));
#endif

	    A(j+1,j-1) = 0.;

	    i__2 = ilastm;
	    for (jc = j; jc <= ilastm; ++jc) {
		temp = c * A(j,jc) + s * A(j+1,jc);
		A(j+1,jc) = -s * A(j,jc) + c * A(j+1,jc);
		A(j,jc) = temp;
		temp2 = c * B(j,jc) + s * B(j+1,jc);
		B(j+1,jc) = -s * B(j,jc) + c * B(j+1,jc);
		B(j,jc) = temp2;
/* L300: */
	    }
	    if (ilq) {
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    temp = c * Q(jr,j) + s * Q(jr,j+1);
		    Q(jr,j+1) = -s * Q(jr,j) + c * 
			    Q(jr,j+1);
		    Q(jr,j) = temp;
/* L310: */
		}
	    }

/*           Rotations from the right. */

	    temp = B(j+1,j+1);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&temp, &B(j+1,j), &c, &s, &B(j+1,j+1));
#endif

	    B(j+1,j) = 0.;

	    i__2 = ilast;
	    for (jr = ifrstm; jr <= ilast; ++jr) {
		temp = c * A(jr,j+1) + s * A(jr,j);
		A(jr,j) = -s * A(jr,j+1) + c * A(jr,j);
		A(jr,j+1) = temp;
/* L320: */
	    }
	    i__2 = ilast - 1;
	    for (jr = ifrstm; jr <= ilast-1; ++jr) {
		temp = c * B(jr,j+1) + s * B(jr,j);
		B(jr,j) = -s * B(jr,j+1) + c * B(jr,j);
		B(jr,j+1) = temp;
/* L330: */
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    temp = c * Z(jr,j+1) + s * Z(jr,j);
		    Z(jr,j) = -s * Z(jr,j+1) + c * 
			    Z(jr,j);
		    Z(jr,j+1) = temp;
/* L340: */
		}
	    }

/*           End of LONG DOUBLE-Shift code */

	}

	goto L350;

/*        End of iteration loop */

L350:
/* L360: */
	;
    }

/*     Drop-through = non-convergence   

   L370: */
    *info = ilast;
    goto L420;

/*     Successful completion of all QZ steps */

L380:

/*     Set Eigenvalues 1:ILO-1 */

    i__1 = *ilo - 1;
    for (j = 1; j <= *ilo-1; ++j) {
	if (B(j,j) < 0.) {
	    if (ilschr) {
		i__2 = j;
		for (jr = 1; jr <= j; ++jr) {
		    A(jr,j) = -A(jr,j);
		    B(jr,j) = -B(jr,j);
/* L390: */
		}
	    } else {
		A(j,j) = -A(j,j);
		B(j,j) = -B(j,j);
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= *n; ++jr) {
		    Z(jr,j) = -Z(jr,j);
/* L400: */
		}
	    }
	}
	ALPHAR(j) = A(j,j);
	ALPHAI(j) = 0.;
	BETA(j) = B(j,j);
/* L410: */
    }

/*     Normal Termination */

    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size 
*/

L420:
    WORK(1) = (LONG DOUBLE) (*n);
    return;

/*     End of DHGEQZ */

} /* dhgeqz_ */

