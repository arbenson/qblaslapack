#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#define SIGN(a,b)          ( ((b)>=0) ? (a) : -(a) )
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggbal_(char *job, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggbal(char *job, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggbal_(char *job, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *b, int *ldb, int *ilo, int *ihi, 
	LONG DOUBLE *lscale, LONG DOUBLE *rscale, LONG DOUBLE *work, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAL balances a pair of general real matrices (A,B).  This   
    involves, first, permuting A and B by similarity transformations to   
    isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N   
    elements on the diagonal; and second, applying a diagonal similarity 
  
    transformation to rows and columns ILO to IHI to make the rows   
    and columns as close in norm as possible. Both steps are optional.   

    Balancing may reduce the 1-norm of the matrices, and improve the   
    accuracy of the computed eigenvalues and/or eigenvectors in the   
    generalized eigenvalue problem A*x = lambda*B*x.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the operations to be performed on A and B:   
            = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0   
                    and RSCALE(I) = 1.0 for i = 1,...,N.   
            = 'P':  permute only;   
            = 'S':  scale only;   
            = 'B':  both permute and scale.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the input matrix A.   
            On exit,  A is overwritten by the balanced matrix.   
            If JOB = 'N', A is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the input matrix B.   
            On exit,  B is overwritten by the balanced matrix.   
            If JOB = 'N', B is not referenced.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,N).   

    ILO     (output) INTEGER   
    IHI     (output) INTEGER   
            ILO and IHI are set to ints such that on exit   
            A(i,j) = 0 and B(i,j) = 0 if i > j and   
            j = 1,...,ILO-1 or i = IHI+1,...,N.   
            If JOB = 'N' or 'S', ILO = 1 and IHI = N.   

    LSCALE  (output) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the left side of A and B.  If P(j) is the index of the   
            row interchanged with row j, and D(j)   
            is the scaling factor applied to row j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    RSCALE  (output) LONG DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the right side of A and B.  If P(j) is the index of the   
            column interchanged with column j, and D(j)   
            is the scaling factor applied to column j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. WARD, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b34 = 10.;
    static LONG DOUBLE c_b70 = .5;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static int lcab;
    static LONG DOUBLE beta, coef;
    static int irab, lrab;
    static LONG DOUBLE basl, cmax;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE coef2, coef5;
    static int i, j, k, l, m;
    static LONG DOUBLE P_gamma, t, alpha;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);
    static LONG DOUBLE sfmin, sfmax;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dswap_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qswap_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static int iflow;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static int kount, jc;
    static LONG DOUBLE ta, tb, tc;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static int ir, it;
    static LONG DOUBLE ew;
    static int nr;
    static LONG DOUBLE pgamma;

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
    static int lsfmin, lsfmax, ip1, jp1, lm1;
    static LONG DOUBLE cab, rab, ewc, cor, sum;
    static int nrp2, icab;



#define LSCALE(I) lscale[(I)-1]
#define RSCALE(I) rscale[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(
	    job, "S") && ! lsame_(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    } else if (*ldb < MAX(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGBAL", &i__1);
	return;
    }

    k = 1;
    l = *n;

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (lsame_(job, "N")) {
	*ilo = 1;
	*ihi = *n;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    LSCALE(i) = 1.;
	    RSCALE(i) = 1.;
/* L10: */
	}
	return;
    }

    if (k == l) {
	*ilo = 1;
	*ihi = 1;
	LSCALE(1) = 1.;
	RSCALE(1) = 1.;
	return;
    }

    if (lsame_(job, "S")) {
	goto L190;
    }

    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues.   

       Find row with one nonzero in columns 1 through L */

L20:
    l = lm1;
    if (l != 1) {
	goto L30;
    }

    RSCALE(1) = 1.;
    LSCALE(1) = 1.;
    goto L190;

L30:
    lm1 = l - 1;
    for (i = l; i >= 1; --i) {
	i__1 = lm1;
	for (j = 1; j <= lm1; ++j) {
	    jp1 = j + 1;
	    if (A(i,j) != 0. || B(i,j) != 0.) {
		goto L50;
	    }
/* L40: */
	}
	j = l;
	goto L70;

L50:
	i__1 = l;
	for (j = jp1; j <= l; ++j) {
	    if (A(i,j) != 0. || B(i,j) != 0.) {
		goto L80;
	    }
/* L60: */
	}
	j = jp1 - 1;

L70:
	m = l;
	iflow = 1;
	goto L160;
L80:
	;
    }
    goto L100;

/*     Find column with one nonzero in rows K through N */

L90:
    ++k;

L100:
    i__1 = l;
    for (j = k; j <= l; ++j) {
	i__2 = lm1;
	for (i = k; i <= lm1; ++i) {
	    ip1 = i + 1;
	    if (A(i,j) != 0. || B(i,j) != 0.) {
		goto L120;
	    }
/* L110: */
	}
	i = l;
	goto L140;
L120:
	i__2 = l;
	for (i = ip1; i <= l; ++i) {
	    if (A(i,j) != 0. || B(i,j) != 0.) {
		goto L150;
	    }
/* L130: */
	}
	i = ip1 - 1;
L140:
	m = k;
	iflow = 2;
	goto L160;
L150:
	;
    }
    goto L190;

/*     Permute rows M and I */

L160:
    LSCALE(m) = (LONG DOUBLE) i;
    if (i == m) {
	goto L170;
    }
    i__1 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&i__1, &A(i,k), lda, &A(m,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&i__1, &A(i,k), lda, &A(m,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&i__1, &A(i,k), lda, &A(m,k), lda);
#endif

    i__1 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&i__1, &B(i,k), ldb, &B(m,k), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&i__1, &B(i,k), ldb, &B(m,k), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&i__1, &B(i,k), ldb, &B(m,k), ldb);
#endif


/*     Permute columns M and J */

L170:
    RSCALE(m) = (LONG DOUBLE) j;
    if (j == m) {
	goto L180;
    }

#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&l, &A(1,j), &c__1, &A(1,m), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dswap_(&l, &B(1,j), &c__1, &B(1,m), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qswap(&l, &B(1,j), &c__1, &B(1,m), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qswap_(&l, &B(1,j), &c__1, &B(1,m), &c__1);
#endif


L180:
    switch (iflow) {
	case 1:  goto L20;
	case 2:  goto L90;
    }

L190:
    *ilo = k;
    *ihi = l;

    if (*ilo == *ihi) {
	return;
    }

    if (lsame_(job, "P")) {
	return;
    }

/*     Balance the submatrix in rows ILO to IHI. */

    nr = *ihi - *ilo + 1;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	RSCALE(i) = 0.;
	LSCALE(i) = 0.;

	WORK(i) = 0.;
	WORK(i + *n) = 0.;
	WORK(i + (*n << 1)) = 0.;
	WORK(i + *n * 3) = 0.;
	WORK(i + (*n << 2)) = 0.;
	WORK(i + *n * 5) = 0.;
/* L200: */
    }

/*     Compute right side vector in resulting linear equations */

    basl = log(c_b34);
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *ihi;
	for (j = *ilo; j <= *ihi; ++j) {
	    tb = B(i,j);
	    ta = A(i,j);
	    if (ta == 0.) {
		goto L210;
	    }
	    d__1 = ABS(ta);
	    ta = log(d__1) / basl;
L210:
	    if (tb == 0.) {
		goto L220;
	    }
	    d__1 = ABS(tb);
	    tb = log(d__1) / basl;
L220:
	    WORK(i + (*n << 2)) = WORK(i + (*n << 2)) - ta - tb;
	    WORK(j + *n * 5) = WORK(j + *n * 5) - ta - tb;
/* L230: */
	}
/* L240: */
    }

    coef = 1. / (LONG DOUBLE) (nr << 1);
    coef2 = coef * coef;
    coef5 = coef2 * .5;
    nrp2 = nr + 2;
    beta = 0.;
    it = 1;

/*     Start generalized conjugate gradient iteration */

L250:


#ifdef PETSC_PREFIX_SUFFIX
    P_gamma = ddot_(&nr, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + (*n << 2))
#endif
#ifdef Q_C_PREFIX_SUFFIX
    P_gamma = qdot(&nr, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + (*n << 2))
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    P_gamma = qdot_(&nr, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + (*n << 2))
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    , &c__1) + ddot_(&nr, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo + *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    , &c__1) + qdot(&nr, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo + *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    , &c__1) + qdot_(&nr, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo + *
#endif

	    n * 5), &c__1);

    ew = 0.;
    ewc = 0.;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	ew += WORK(i + (*n << 2));
	ewc += WORK(i + *n * 5);
/* L260: */
    }

/* Computing 2nd power */
    d__1 = ew;
/* Computing 2nd power */
    d__2 = ewc;
/* Computing 2nd power */
    d__3 = ew - ewc;
    P_gamma = coef * P_gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
    if (P_gamma == 0.) {
	goto L350;
    }
    if (it != 1) {
	beta = P_gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.);
    tc = coef5 * (ew - ewc * 3.);


#ifdef PETSC_PREFIX_SUFFIX
    dscal_(&nr, &beta, &WORK(*ilo), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qscal(&nr, &beta, &WORK(*ilo), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qscal_(&nr, &beta, &WORK(*ilo), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dscal_(&nr, &beta, &WORK(*ilo + *n), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qscal(&nr, &beta, &WORK(*ilo + *n), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qscal_(&nr, &beta, &WORK(*ilo + *n), &c__1);
#endif



#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(&nr, &coef, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + *n), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(&nr, &coef, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + *n), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(&nr, &coef, &WORK(*ilo + (*n << 2)), &c__1, &WORK(*ilo + *n), &
#endif

	    c__1);

#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(&nr, &coef, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(&nr, &coef, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(&nr, &coef, &WORK(*ilo + *n * 5), &c__1, &WORK(*ilo), &c__1);
#endif


    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	WORK(i) += tc;
	WORK(i + *n) += t;
/* L270: */
    }

/*     Apply matrix to vector */

    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (j = *ilo; j <= *ihi; ++j) {
	    if (A(i,j) == 0.) {
		goto L280;
	    }
	    ++kount;
	    sum += WORK(j);
L280:
	    if (B(i,j) == 0.) {
		goto L290;
	    }
	    ++kount;
	    sum += WORK(j);
L290:
	    ;
	}
	WORK(i + (*n << 1)) = (LONG DOUBLE) kount * WORK(i + *n) + sum;
/* L300: */
    }

    i__1 = *ihi;
    for (j = *ilo; j <= *ihi; ++j) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (i = *ilo; i <= *ihi; ++i) {
	    if (A(i,j) == 0.) {
		goto L310;
	    }
	    ++kount;
	    sum += WORK(i + *n);
L310:
	    if (B(i,j) == 0.) {
		goto L320;
	    }
	    ++kount;
	    sum += WORK(i + *n);
L320:
	    ;
	}
	WORK(j + *n * 3) = (LONG DOUBLE) kount * WORK(j) + sum;
/* L330: */
    }


#ifdef PETSC_PREFIX_SUFFIX
    sum = ddot_(&nr, &WORK(*ilo + *n), &c__1, &WORK(*ilo + (*n << 1)), &c__1) 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sum = qdot(&nr, &WORK(*ilo + *n), &c__1, &WORK(*ilo + (*n << 1)), &c__1) 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sum = qdot_(&nr, &WORK(*ilo + *n), &c__1, &WORK(*ilo + (*n << 1)), &c__1) 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    + ddot_(&nr, &WORK(*ilo), &c__1, &WORK(*ilo + *n * 3), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    + qdot(&nr, &WORK(*ilo), &c__1, &WORK(*ilo + *n * 3), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    + qdot_(&nr, &WORK(*ilo), &c__1, &WORK(*ilo + *n * 3), &c__1);
#endif

    alpha = P_gamma / sum;

/*     Determine correction to current iteration */

    cmax = 0.;
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	cor = alpha * WORK(i + *n);
	if (ABS(cor) > cmax) {
	    cmax = ABS(cor);
	}
	LSCALE(i) += cor;
	cor = alpha * WORK(i);
	if (ABS(cor) > cmax) {
	    cmax = ABS(cor);
	}
	RSCALE(i) += cor;
/* L340: */
    }
    if (cmax < .5) {
	goto L350;
    }

    d__1 = -alpha;

#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(&nr, &d__1, &WORK(*ilo + (*n << 1)), &c__1, &WORK(*ilo + (*n << 2))
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(&nr, &d__1, &WORK(*ilo + (*n << 1)), &c__1, &WORK(*ilo + (*n << 2))
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(&nr, &d__1, &WORK(*ilo + (*n << 1)), &c__1, &WORK(*ilo + (*n << 2))
#endif

	    , &c__1);
    d__1 = -alpha;

#ifdef PETSC_PREFIX_SUFFIX
    daxpy_(&nr, &d__1, &WORK(*ilo + *n * 3), &c__1, &WORK(*ilo + *n * 5), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qaxpy(&nr, &d__1, &WORK(*ilo + *n * 3), &c__1, &WORK(*ilo + *n * 5), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qaxpy_(&nr, &d__1, &WORK(*ilo + *n * 3), &c__1, &WORK(*ilo + *n * 5), &
#endif

	    c__1);

    pgamma = P_gamma;
    ++it;
    if (it <= nrp2) {
	goto L250;
    }

/*     End generalized conjugate gradient iteration */

L350:

#ifdef PETSC_PREFIX_SUFFIX
    sfmin = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    sfmin = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    sfmin = qlamch_("S");
#endif

    sfmax = 1. / sfmin;
    lsfmin = (int) (log(sfmin) / basl + 1.);
    lsfmax = (int) (log(sfmax) / basl);
    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *n - *ilo + 1;

#ifdef PETSC_PREFIX_SUFFIX
	irab = idamax_(&i__2, &A(i,*ilo), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	irab = iqamax(&i__2, &A(i,*ilo), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	irab = iqamax_(&i__2, &A(i,*ilo), lda);
#endif

	rab = (d__1 = A(i,irab+*ilo-1), ABS(d__1));
	i__2 = *n - *ilo + 1;

#ifdef PETSC_PREFIX_SUFFIX
	irab = idamax_(&i__2, &B(i,*ilo), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	irab = iqamax(&i__2, &B(i,*ilo), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	irab = iqamax_(&i__2, &B(i,*ilo), lda);
#endif

/* Computing MAX */
	d__2 = rab, d__3 = (d__1 = B(i,irab+*ilo-1), ABS(
		d__1));
	rab = MAX(d__2,d__3);
	d__1 = rab + sfmin;
	lrab = (int) (log(d__1) / basl + 1.);
	ir = (int) (LSCALE(i) + SIGN(c_b70, LSCALE(i)));
/* Computing MIN */
	i__2 = MAX(ir,lsfmin), i__2 = MIN(i__2,lsfmax), i__3 = lsfmax - lrab;
	ir = MIN(i__2,i__3);
	LSCALE(i) = pow((LONG DOUBLE)c_b34, (LONG DOUBLE)ir);

#ifdef PETSC_PREFIX_SUFFIX
	icab = idamax_(ihi, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	icab = iqamax(ihi, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	icab = iqamax_(ihi, &A(1,i), &c__1);
#endif

	cab = (d__1 = A(icab,i), ABS(d__1));

#ifdef PETSC_PREFIX_SUFFIX
	icab = idamax_(ihi, &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	icab = iqamax(ihi, &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	icab = iqamax_(ihi, &B(1,i), &c__1);
#endif

/* Computing MAX */
	d__2 = cab, d__3 = (d__1 = B(icab,i), ABS(d__1));
	cab = MAX(d__2,d__3);
	d__1 = cab + sfmin;
	lcab = (int) (log(d__1) / basl + 1.);
	jc = (int) (RSCALE(i) + SIGN(c_b70, RSCALE(i)));
/* Computing MIN */
	i__2 = MAX(jc,lsfmin), i__2 = MIN(i__2,lsfmax), i__3 = lsfmax - lcab;
	jc = MIN(i__2,i__3);
	RSCALE(i) = pow(c_b34, (LONG DOUBLE)jc);
/* L360: */
    }

/*     Row scaling of matrices A and B */

    i__1 = *ihi;
    for (i = *ilo; i <= *ihi; ++i) {
	i__2 = *n - *ilo + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__2, &LSCALE(i), &A(i,*ilo), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__2, &LSCALE(i), &A(i,*ilo), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__2, &LSCALE(i), &A(i,*ilo), lda);
#endif

	i__2 = *n - *ilo + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__2, &LSCALE(i), &B(i,*ilo), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__2, &LSCALE(i), &B(i,*ilo), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__2, &LSCALE(i), &B(i,*ilo), ldb);
#endif

/* L370: */
    }

/*     Column scaling of matrices A and B */

    i__1 = *ihi;
    for (j = *ilo; j <= *ihi; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(ihi, &RSCALE(j), &A(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(ihi, &RSCALE(j), &A(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(ihi, &RSCALE(j), &A(1,j), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dscal_(ihi, &RSCALE(j), &B(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(ihi, &RSCALE(j), &B(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(ihi, &RSCALE(j), &B(1,j), &c__1);
#endif

/* L380: */
    }

    return;

/*     End of DGGBAL */

} /* dggbal_ */

