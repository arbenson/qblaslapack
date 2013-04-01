#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaqtr_(long int *ltran, long int *lreal, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaqtr(long int *ltran, long int *lreal, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaqtr_(long int *ltran, long int *lreal, int *n, 
#endif

	LONG DOUBLE *t, int *ldt, LONG DOUBLE *b, LONG DOUBLE *w, LONG DOUBLE 
	*scale, LONG DOUBLE *x, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAQTR solves the real quasi-triangular system   

                 op(T)*p = scale*c,               if LREAL = .TRUE.   

    or the complex quasi-triangular systems   

               op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.   

    in real arithmetic, where T is upper quasi-triangular.   
    If LREAL = .FALSE., then the first diagonal block of T must be   
    1 by 1, B is the specially structured matrix   

                   B = [ b(1) b(2) ... b(n) ]   
                       [       w            ]   
                       [           w        ]   
                       [              .     ]   
                       [                 w  ]   

    op(A) = A or A', A' denotes the conjugate transpose of   
    matrix A.   

    On input, X = [ c ].  On output, X = [ p ].   
                  [ d ]                  [ q ]   

    This subroutine is designed for the condition number estimation   
    in routine DTRSNA.   

    Arguments   
    =========   

    LTRAN   (input) LOGICAL   
            On entry, LTRAN specifies the option of conjugate transpose: 
  
               = .FALSE.,    op(T+i*B) = T+i*B,   
               = .TRUE.,     op(T+i*B) = (T+i*B)'.   

    LREAL   (input) LOGICAL   
            On entry, LREAL specifies the input matrix structure:   
               = .FALSE.,    the input is complex   
               = .TRUE.,     the input is real   

    N       (input) INTEGER   
            On entry, N specifies the order of T+i*B. N >= 0.   

    T       (input) LONG DOUBLE PRECISION array, dimension (LDT,N)   
            On entry, T contains a matrix in Schur canonical form.   
            If LREAL = .FALSE., then the first diagonal block of T mu   
            be 1 by 1.   

    LDT     (input) INTEGER   
            The leading dimension of the matrix T. LDT >= MAX(1,N).   

    B       (input) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, B contains the elements to form the matrix   
            B as described above.   
            If LREAL = .TRUE., B is not referenced.   

    W       (input) LONG DOUBLE PRECISION   
            On entry, W is the diagonal element of the matrix B.   
            If LREAL = .TRUE., W is not referenced.   

    SCALE   (output) LONG DOUBLE PRECISION   
            On exit, SCALE is the scale factor.   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (2*N)   
            On entry, X contains the right hand side of the system.   
            On exit, X is overwritten by the solution.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            On exit, INFO is set to   
               0: successful exit.   
                 1: the some diagonal 1 by 1 block has been perturbed by 
  
                    a small number SMIN to keep nonsingularity.   
                 2: the some diagonal 2 by 2 block has been perturbed by 
  
                    a small number in DLALN2 to keep nonsingularity.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

    ===================================================================== 
  


       Do not test the input parameters for errors   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static long int c_false = 0;
    static int c__2 = 2;
    static LONG DOUBLE c_b21 = 1.;
    static LONG DOUBLE c_b25 = 0.;
    static long int c_true = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */

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
    static int ierr;
    static LONG DOUBLE smin, xmax, d[4]	/* was [2][2] */;
    static int i, j, k;
    static LONG DOUBLE v[4]	/* was [2][2] */, z;

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

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dasum_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qasum(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qasum_(int *, LONG DOUBLE *, int *);
#endif


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
    static int jnext, P_j1, j2;
    static LONG DOUBLE sminw;
    static int n1, n2;
    static LONG DOUBLE xnorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaln2_(long int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaln2(long int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaln2_(long int *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *
	    , LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *), dlange_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlange_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlange_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static LONG DOUBLE si, xj;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    static LONG DOUBLE scaloc, sr;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qladiv_(LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE bignum;
    static long int notran;
    static LONG DOUBLE smlnum, rec, eps, tjj, tmp;



#define V(I) v[(I)]
#define WAS(I) was[(I)]
#define B(I) b[(I)-1]
#define X(I) x[(I)-1]
#define WORK(I) work[(I)-1]

#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    notran = ! (*ltran);
    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Set constants to control overflow */


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("P");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("P");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("P");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S") / eps;
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S") / eps;
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S") / eps;
#endif

    bignum = 1. / smlnum;


#ifdef PETSC_PREFIX_SUFFIX
    xnorm = dlange_("M", n, n, &T(1,1), ldt, d);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    xnorm = qlange("M", n, n, &T(1,1), ldt, d);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    xnorm = qlange_("M", n, n, &T(1,1), ldt, d);
#endif

    if (! (*lreal)) {
/* Computing MAX */

#ifdef PETSC_PREFIX_SUFFIX
	d__1 = xnorm, d__2 = ABS(*w), d__1 = MAX(d__1,d__2), d__2 = dlange_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	d__1 = xnorm, d__2 = ABS(*w), d__1 = MAX(d__1,d__2), d__2 = qlange(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	d__1 = xnorm, d__2 = ABS(*w), d__1 = MAX(d__1,d__2), d__2 = qlange_(
#endif

		"M", n, &c__1, &B(1), n, d);
	xnorm = MAX(d__1,d__2);
    }
/* Computing MAX */
    d__1 = smlnum, d__2 = eps * xnorm;
    smin = MAX(d__1,d__2);

/*     Compute 1-norm of each column of strictly upper triangular   
       part of T to control overflow in triangular solver. */

    WORK(1) = 0.;
    i__1 = *n;
    for (j = 2; j <= *n; ++j) {
	i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	WORK(j) = dasum_(&i__2, &T(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	WORK(j) = qasum(&i__2, &T(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	WORK(j) = qasum_(&i__2, &T(1,j), &c__1);
#endif

/* L10: */
    }

    if (! (*lreal)) {
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    WORK(i) += (d__1 = B(i), ABS(d__1));
/* L20: */
	}
    }

    n2 = *n << 1;
    n1 = *n;
    if (! (*lreal)) {
	n1 = n2;
    }

#ifdef PETSC_PREFIX_SUFFIX
    k = idamax_(&n1, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    k = iqamax(&n1, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    k = iqamax_(&n1, &X(1), &c__1);
#endif

    xmax = (d__1 = X(k), ABS(d__1));
    *scale = 1.;

    if (xmax > bignum) {
	*scale = bignum / xmax;

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&n1, scale, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&n1, scale, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&n1, scale, &X(1), &c__1);
#endif

	xmax = bignum;
    }

    if (*lreal) {

	if (notran) {

/*           Solve T*p = scale*c */

	    jnext = *n;
	    for (j = *n; j >= 1; --j) {
		if (j > jnext) {
		    goto L30;
		}
		P_j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1 && T(j,j-1) != 0.) {
		    P_j1 = j - 1;
		    j2 = j;
		    jnext = j - 2;
		}

		if (P_j1 == j2) {

/*                 Meet 1 by 1 diagonal block   

                   Scale to avoid overflow when computing 
  
                       x(j) = b(j)/T(j,j) */

		    xj = (d__1 = X(P_j1), ABS(d__1));
		    tjj = (d__1 = T(P_j1,P_j1), ABS(d__1));
		    tmp = T(P_j1,P_j1);
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }

		    if (xj == 0.) {
			goto L30;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    X(P_j1) /= tmp;
		    xj = (d__1 = X(P_j1), ABS(d__1));

/*                 Scale x if necessary to avoid overflow 
when adding a   
                   multiple of column P_j1 of T. */

		    if (xj > 1.) {
			rec = 1. / xj;
			if (WORK(P_j1) > (bignum - xmax) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			}
		    }
		    if (P_j1 > 1) {
			i__1 = P_j1 - 1;
			d__1 = -X(P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif

				, &c__1);
			i__1 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
			k = idamax_(&i__1, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			k = iqamax(&i__1, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			k = iqamax_(&i__1, &X(1), &c__1);
#endif

			xmax = (d__1 = X(k), ABS(d__1));
		    }

		} else {

/*                 Meet 2 by 2 diagonal block   

                   Call 2 by 2 linear system solve, to tak
e   
                   care of possible overflow by scaling fa
ctor. */

		    d[0] = X(P_j1);
		    d[1] = X(j2);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif

			    &c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(n, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(n, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(n, &scaloc, &X(1), &c__1);
#endif

			*scale *= scaloc;
		    }
		    X(P_j1) = V(0);
		    X(j2) = V(1);

/*                 Scale V(1,1) (= X(P_J1)) and/or V(2,1) (=
X(J2))   
                   to avoid overflow in updating right-han
d side.   

   Computing MAX */
		    d__1 = ABS(V(0)), d__2 = ABS(V(1));
		    xj = MAX(d__1,d__2);
		    if (xj > 1.) {
			rec = 1. / xj;
/* Computing MAX */
			d__1 = WORK(P_j1), d__2 = WORK(j2);
			if (MAX(d__1,d__2) > (bignum - xmax) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			}
		    }

/*                 Update right-hand side */

		    if (P_j1 > 1) {
			i__1 = P_j1 - 1;
			d__1 = -X(P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif

				, &c__1);
			i__1 = P_j1 - 1;
			d__1 = -X(j2);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif

				, &c__1);
			i__1 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
			k = idamax_(&i__1, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			k = iqamax(&i__1, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			k = iqamax_(&i__1, &X(1), &c__1);
#endif

			xmax = (d__1 = X(k), ABS(d__1));
		    }

		}

L30:
		;
	    }

	} else {

/*           Solve T'*p = scale*c */

	    jnext = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (j < jnext) {
		    goto L40;
		}
		P_j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < *n && T(j+1,j) != 0.) {
		    P_j1 = j;
		    j2 = j + 1;
		    jnext = j + 2;
		}

		if (P_j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t. */

		    xj = (d__1 = X(P_j1), ABS(d__1));
		    if (xmax > 1.) {
			rec = 1. / xmax;
			if (WORK(P_j1) > (bignum - xj) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    X(P_j1) -= ddot_(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    X(P_j1) -= qdot(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    X(P_j1) -= qdot_(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif

			    c__1);

		    xj = (d__1 = X(P_j1), ABS(d__1));
		    tjj = (d__1 = T(P_j1,P_j1), ABS(d__1));
		    tmp = T(P_j1,P_j1);
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    X(P_j1) /= tmp;
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = X(P_j1), ABS(d__1));
		    xmax = MAX(d__2,d__3);

		} else {

/*                 2 by 2 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side elements by inner produ
ct.   

   Computing MAX */
		    d__3 = (d__1 = X(P_j1), ABS(d__1)), d__4 = (d__2 = X(j2), 
			    ABS(d__2));
		    xj = MAX(d__3,d__4);
		    if (xmax > 1.) {
			rec = 1. / xmax;
/* Computing MAX */
			d__1 = WORK(j2), d__2 = WORK(P_j1);
			if (MAX(d__1,d__2) > (bignum - xj) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(n, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(n, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[0] = X(P_j1) - ddot_(&i__2, &T(1,P_j1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[0] = X(P_j1) - qdot(&i__2, &T(1,P_j1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[0] = X(P_j1) - qdot_(&i__2, &T(1,P_j1), &c__1, &
#endif

			    X(1), &c__1);
		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[1] = X(j2) - ddot_(&i__2, &T(1,j2), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[1] = X(j2) - qdot(&i__2, &T(1,j2), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[1] = X(j2) - qdot_(&i__2, &T(1,j2), &c__1, &
#endif

			    X(1), &c__1);


#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__1, &smin, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, &
#endif

			    c_b25, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(n, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(n, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(n, &scaloc, &X(1), &c__1);
#endif

			*scale *= scaloc;
		    }
		    X(P_j1) = V(0);
		    X(j2) = V(1);
/* Computing MAX */
		    d__3 = (d__1 = X(P_j1), ABS(d__1)), d__4 = (d__2 = X(j2), 
			    ABS(d__2)), d__3 = MAX(d__3,d__4);
		    xmax = MAX(d__3,xmax);

		}
L40:
		;
	    }
	}

    } else {

/* Computing MAX */
	d__1 = eps * ABS(*w);
	sminw = MAX(d__1,smin);
	if (notran) {

/*           Solve (T + iB)*(p+iq) = c+id */

	    jnext = *n;
	    for (j = *n; j >= 1; --j) {
		if (j > jnext) {
		    goto L70;
		}
		P_j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1 && T(j,j-1) != 0.) {
		    P_j1 = j - 1;
		    j2 = j;
		    jnext = j - 2;
		}

		if (P_j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 division */

		    z = *w;
		    if (P_j1 == 1) {
			z = B(1);
		    }
		    xj = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(*n + P_j1), ABS(
			    d__2));
		    tjj = (d__1 = T(P_j1,P_j1), ABS(d__1)) + ABS(z);
		    tmp = T(P_j1,P_j1);
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }

		    if (xj == 0.) {
			goto L70;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }

#ifdef PETSC_PREFIX_SUFFIX
		    dladiv_(&X(P_j1), &X(*n + P_j1), &tmp, &z, &sr, &si);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qladiv(&X(P_j1), &X(*n + P_j1), &tmp, &z, &sr, &si);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qladiv_(&X(P_j1), &X(*n + P_j1), &tmp, &z, &sr, &si);
#endif

		    X(P_j1) = sr;
		    X(*n + P_j1) = si;
		    xj = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(*n + P_j1), ABS(
			    d__2));

/*                 Scale x if necessary to avoid overflow 
when adding a   
                   multiple of column P_j1 of T. */

		    if (xj > 1.) {
			rec = 1. / xj;
			if (WORK(P_j1) > (bignum - xmax) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			}
		    }

		    if (P_j1 > 1) {
			i__1 = P_j1 - 1;
			d__1 = -X(P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif

				, &c__1);
			i__1 = P_j1 - 1;
			d__1 = -X(*n + P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif

				n + 1), &c__1);

			X(1) += B(P_j1) * X(*n + P_j1);
			X(*n + 1) -= B(P_j1) * X(P_j1);

			xmax = 0.;
			i__1 = P_j1 - 1;
			for (k = 1; k <= P_j1-1; ++k) {
/* Computing MAX */
			    d__3 = xmax, d__4 = (d__1 = X(k), ABS(d__1)) + (
				    d__2 = X(k + *n), ABS(d__2));
			    xmax = MAX(d__3,d__4);
/* L50: */
			}
		    }

		} else {

/*                 Meet 2 by 2 diagonal block */

		    d[0] = X(P_j1);
		    d[1] = X(j2);
		    d[2] = X(*n + P_j1);
		    d[3] = X(*n + j2);
		    d__1 = -(*w);

#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_false, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_false, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_false, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &
#endif

			    c_b25, &d__1, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {
			i__1 = *n << 1;

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(&i__1, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(&i__1, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(&i__1, &scaloc, &X(1), &c__1);
#endif

			*scale = scaloc * *scale;
		    }
		    X(P_j1) = V(0);
		    X(j2) = V(1);
		    X(*n + P_j1) = V(2);
		    X(*n + j2) = V(3);

/*                 Scale X(P_J1), .... to avoid overflow in 
  
                   updating right hand side.   

   Computing MAX */
		    d__1 = ABS(V(0)) + ABS(V(2)), d__2 = ABS(V(1)) + ABS(V(3))
			    ;
		    xj = MAX(d__1,d__2);
		    if (xj > 1.) {
			rec = 1. / xj;
/* Computing MAX */
			d__1 = WORK(P_j1), d__2 = WORK(j2);
			if (MAX(d__1,d__2) > (bignum - xmax) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			}
		    }

/*                 Update the right-hand side. */

		    if (P_j1 > 1) {
			i__1 = P_j1 - 1;
			d__1 = -X(P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(1)
#endif

				, &c__1);
			i__1 = P_j1 - 1;
			d__1 = -X(j2);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(1)
#endif

				, &c__1);

			i__1 = P_j1 - 1;
			d__1 = -X(*n + P_j1);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,P_j1), &c__1, &X(*
#endif

				n + 1), &c__1);
			i__1 = P_j1 - 1;
			d__1 = -X(*n + j2);

#ifdef PETSC_PREFIX_SUFFIX
			daxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qaxpy(&i__1, &d__1, &T(1,j2), &c__1, &X(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qaxpy_(&i__1, &d__1, &T(1,j2), &c__1, &X(*
#endif

				n + 1), &c__1);

			X(1) = X(1) + B(P_j1) * X(*n + P_j1) + B(j2) * X(*n + j2);
			X(*n + 1) = X(*n + 1) - B(P_j1) * X(P_j1) - B(j2) * X(j2);

			xmax = 0.;
			i__1 = P_j1 - 1;
			for (k = 1; k <= P_j1-1; ++k) {
/* Computing MAX */
			    d__3 = (d__1 = X(k), ABS(d__1)) + (d__2 = X(k + *
				    n), ABS(d__2));
			    xmax = MAX(d__3,xmax);
/* L60: */
			}
		    }

		}
L70:
		;
	    }

	} else {

/*           Solve (T + iB)'*(p+iq) = c+id */

	    jnext = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (j < jnext) {
		    goto L80;
		}
		P_j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < *n && T(j+1,j) != 0.) {
		    P_j1 = j;
		    j2 = j + 1;
		    jnext = j + 2;
		}

		if (P_j1 == j2) {

/*                 1 by 1 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t. */

		    xj = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(P_j1 + *n), ABS(
			    d__2));
		    if (xmax > 1.) {
			rec = 1. / xmax;
			if (WORK(P_j1) > (bignum - xj) * rec) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    X(P_j1) -= ddot_(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    X(P_j1) -= qdot(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    X(P_j1) -= qdot_(&i__2, &T(1,P_j1), &c__1, &X(1), &
#endif

			    c__1);
		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    X(*n + P_j1) -= ddot_(&i__2, &T(1,P_j1), &c__1, &X(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    X(*n + P_j1) -= qdot(&i__2, &T(1,P_j1), &c__1, &X(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    X(*n + P_j1) -= qdot_(&i__2, &T(1,P_j1), &c__1, &X(
#endif

			    *n + 1), &c__1);
		    if (P_j1 > 1) {
			X(P_j1) -= B(P_j1) * X(*n + 1);
			X(*n + P_j1) += B(P_j1) * X(1);
		    }
		    xj = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(P_j1 + *n), ABS(
			    d__2));

		    z = *w;
		    if (P_j1 == 1) {
			z = B(1);
		    }

/*                 Scale if necessary to avoid overflow in
   
                   complex division */

		    tjj = (d__1 = T(P_j1,P_j1), ABS(d__1)) + ABS(z);
		    tmp = T(P_j1,P_j1);
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }

		    if (tjj < 1.) {
			if (xj > bignum * tjj) {
			    rec = 1. / xj;

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }
		    d__1 = -z;

#ifdef PETSC_PREFIX_SUFFIX
		    dladiv_(&X(P_j1), &X(*n + P_j1), &tmp, &d__1, &sr, &si);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qladiv(&X(P_j1), &X(*n + P_j1), &tmp, &d__1, &sr, &si);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qladiv_(&X(P_j1), &X(*n + P_j1), &tmp, &d__1, &sr, &si);
#endif

		    X(P_j1) = sr;
		    X(P_j1 + *n) = si;
/* Computing MAX */
		    d__3 = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(P_j1 + *n), 
			    ABS(d__2));
		    xmax = MAX(d__3,xmax);

		} else {

/*                 2 by 2 diagonal block   

                   Scale if necessary to avoid overflow in
 forming the   
                   right-hand side element by inner produc
t.   

   Computing MAX */
		    d__5 = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(*n + P_j1), 
			    ABS(d__2)), d__6 = (d__3 = X(j2), ABS(d__3)) + (
			    d__4 = X(*n + j2), ABS(d__4));
		    xj = MAX(d__5,d__6);
		    if (xmax > 1.) {
			rec = 1. / xmax;
/* Computing MAX */
			d__1 = WORK(P_j1), d__2 = WORK(j2);
			if (MAX(d__1,d__2) > (bignum - xj) / xmax) {

#ifdef PETSC_PREFIX_SUFFIX
			    dscal_(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			    qscal(&n2, &rec, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			    qscal_(&n2, &rec, &X(1), &c__1);
#endif

			    *scale *= rec;
			    xmax *= rec;
			}
		    }

		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[0] = X(P_j1) - ddot_(&i__2, &T(1,P_j1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[0] = X(P_j1) - qdot(&i__2, &T(1,P_j1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[0] = X(P_j1) - qdot_(&i__2, &T(1,P_j1), &c__1, &
#endif

			    X(1), &c__1);
		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[1] = X(j2) - ddot_(&i__2, &T(1,j2), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[1] = X(j2) - qdot(&i__2, &T(1,j2), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[1] = X(j2) - qdot_(&i__2, &T(1,j2), &c__1, &
#endif

			    X(1), &c__1);
		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[2] = X(*n + P_j1) - ddot_(&i__2, &T(1,P_j1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[2] = X(*n + P_j1) - qdot(&i__2, &T(1,P_j1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[2] = X(*n + P_j1) - qdot_(&i__2, &T(1,P_j1), &
#endif

			    c__1, &X(*n + 1), &c__1);
		    i__2 = P_j1 - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    d[3] = X(*n + j2) - ddot_(&i__2, &T(1,j2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    d[3] = X(*n + j2) - qdot(&i__2, &T(1,j2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    d[3] = X(*n + j2) - qdot_(&i__2, &T(1,j2), &
#endif

			    c__1, &X(*n + 1), &c__1);
		    d[0] -= B(P_j1) * X(*n + 1);
		    d[1] -= B(j2) * X(*n + 1);
		    d[2] += B(P_j1) * X(1);
		    d[3] += B(j2) * X(1);


#ifdef PETSC_PREFIX_SUFFIX
		    dlaln2_(&c_true, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlaln2(&c_true, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlaln2_(&c_true, &c__2, &c__2, &sminw, &c_b21, &T(P_j1,P_j1), ldt, &c_b21, &c_b21, d, &c__2, &c_b25, 
#endif

			    w, v, &c__2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }

		    if (scaloc != 1.) {

#ifdef PETSC_PREFIX_SUFFIX
			dscal_(&n2, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qscal(&n2, &scaloc, &X(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qscal_(&n2, &scaloc, &X(1), &c__1);
#endif

			*scale = scaloc * *scale;
		    }
		    X(P_j1) = V(0);
		    X(j2) = V(1);
		    X(*n + P_j1) = V(2);
		    X(*n + j2) = V(3);
/* Computing MAX */
		    d__5 = (d__1 = X(P_j1), ABS(d__1)) + (d__2 = X(*n + P_j1), 
			    ABS(d__2)), d__6 = (d__3 = X(j2), ABS(d__3)) + (
			    d__4 = X(*n + j2), ABS(d__4)), d__5 = MAX(d__5,
			    d__6);
		    xmax = MAX(d__5,xmax);

		}

L80:
		;
	    }

	}

    }

    return;

/*     End of DLAQTR */

} /* dlaqtr_ */

